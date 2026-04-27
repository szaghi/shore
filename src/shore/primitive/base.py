"""Primitive mesh base class and shared utilities.

A *primitive* is a parametric structured grid built from simple
geometric inputs (corners, radii, axis) — no STL, no projection, no
hyperbolic marching.  Primitives are the connective tissue of an
overset assembly: backgrounds, buffers, and the conforming alternative
to Chimera fringes.

Each concrete primitive is a subclass of :class:`PrimitiveMesh` that
populates ``self.blocks`` with one or more :class:`HexBlock` instances
at construction time.  Downstream operations (Jacobian check, ``.geo``
write, VTK export) are inherited.
"""

from __future__ import annotations

import warnings
from abc import ABC
from pathlib import Path

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.spacing import Spacing1D
from shore.volume.spacing import tanh_one_sided, tanh_two_sided, uniform

# ---------------------------------------------------------------------------
# Axis resolution: --axis "x" | "y" | "z" | "1.0 0.0 0.0"
# ---------------------------------------------------------------------------

_AXIS_LABELS = {
    "x": np.array([1.0, 0.0, 0.0]),
    "y": np.array([0.0, 1.0, 0.0]),
    "z": np.array([0.0, 0.0, 1.0]),
}


def resolve_axis(spec: str | tuple[float, float, float] | np.ndarray) -> np.ndarray:
    """Resolve an axis specifier into a unit vector in R^3.

    Parameters
    ----------
    spec :
        Either a single-character axis label (``'x'``, ``'y'``, ``'z'``,
        case-insensitive) or a 3-component sequence interpreted as a
        free axis direction (will be normalised).

    Returns
    -------
    numpy.ndarray, shape (3,)
        Unit vector in the chosen direction.

    Raises
    ------
    ValueError
        If *spec* is neither a known label nor a 3-vector, or if a
        3-vector has zero magnitude.
    """
    if isinstance(spec, str):
        key = spec.strip().lower()
        if key in _AXIS_LABELS:
            return _AXIS_LABELS[key].copy()
        # Fall through: maybe a whitespace-separated triple "1 0 0".
        parts = spec.split()
        if len(parts) == 3:
            try:
                spec = tuple(float(x) for x in parts)
            except ValueError as exc:
                raise ValueError(
                    f"axis: cannot parse {spec!r} as 'x'/'y'/'z' or '<x> <y> <z>'"
                ) from exc
        else:
            raise ValueError(f"axis: {spec!r} is not 'x', 'y', 'z', or a 3-component vector")

    arr = np.asarray(spec, dtype=float).ravel()
    if arr.shape != (3,):
        raise ValueError(f"axis: expected 3 components, got shape {arr.shape}")
    mag = float(np.linalg.norm(arr))
    if mag < 1e-15:
        raise ValueError("axis: zero-magnitude vector is not a valid direction")
    return arr / mag


def axis_frame(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an orthonormal right-handed frame ``(e1, e2, axis)`` for *axis*.

    *axis* is assumed to be a unit vector.  ``e1`` and ``e2`` span the
    plane perpendicular to *axis*; the choice of ``e1`` is deterministic
    but arbitrary (whichever world axis is "least parallel" to *axis*).
    """
    axis = np.asarray(axis, dtype=float)
    # Pick the world axis least parallel to *axis* to seed the frame.
    candidates = np.eye(3)
    dots = np.abs(candidates @ axis)
    seed = candidates[int(np.argmin(dots))]
    e1 = seed - np.dot(seed, axis) * axis
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    return e1, e2, axis


# ---------------------------------------------------------------------------
# 1D distribution dispatch — single source of truth for all primitives
# ---------------------------------------------------------------------------


def coords_1d(
    n: int,
    start: float,
    stop: float,
    law: str = "uniform",
    beta: float = 3.0,
) -> np.ndarray:
    """Build *n* sample points on ``[start, stop]`` using *law*.

    *law* is one of ``"uniform"``, ``"tanh"`` (one-sided, dense at
    *start*), or ``"tanh2"`` (two-sided, dense at both ends).

    Raises
    ------
    ValueError
        If *law* is unrecognised or *n* < 2.
    """
    if n < 2:
        raise ValueError(f"coords_1d: n must be >= 2, got {n}")
    if law == "uniform":
        return uniform(n, start, stop)
    if law == "tanh":
        return tanh_one_sided(n, start, stop, beta)
    if law == "tanh2":
        return tanh_two_sided(n, start, stop, beta)
    raise ValueError(f"coords_1d: unknown law {law!r} (expected 'uniform', 'tanh', or 'tanh2')")


# ---------------------------------------------------------------------------
# Per-axis TFI: bilinear blend between 4 parallel edges' distributions.
#
# A structured hex block has 4 edges parallel to each axis.  For axis ``i``
# they are j_lo_k_lo, j_hi_k_lo, j_lo_k_hi, j_hi_k_hi.  When all 4 carry
# the same Spacing1D the construction reduces to a tensor product and the
# axis coordinate depends only on ``i``.  When they disagree, the axis
# coordinate becomes a function of (i, j, k) via bilinear blending of the
# four parametric positions:
#
#     t(i, j, k) = (1-η)(1-ζ) * t_00(i) + η(1-ζ) * t_10(i)
#                + (1-η)ζ     * t_01(i) + η ζ    * t_11(i)
#
# with ``η = j/(nj-1)``, ``ζ = k/(nk-1)``.  Each ``t_ab(i)`` is the
# corresponding edge's normalised arc-length parameter on ``[0, 1]``.
# Linear interpolation between ``start`` and ``stop`` produces the world
# coordinate.
#
# This generalises straightforwardly to the j- and k-axes: rotate the
# argument names.
# ---------------------------------------------------------------------------


def _normalise_to_unit(coords: np.ndarray) -> np.ndarray:
    """Map ``coords`` ∈ [a, b] to ``t`` ∈ [0, 1] preserving spacing."""
    a, b = float(coords[0]), float(coords[-1])
    if b - a == 0.0:
        raise ParameterError(f"_normalise_to_unit: degenerate interval [{a}, {b}]")
    return (coords - a) / (b - a)


def axis_field_tfi(
    n_axis: int,
    n_t1: int,
    n_t2: int,
    start: float,
    stop: float,
    edge_specs: dict[tuple[int, int], Spacing1D],
    axis_label: str,
) -> np.ndarray:
    """Return ``(n_axis, n_t1, n_t2)`` axis-coordinate field via TFI.

    The 4 parallel edges of an axis are keyed by their endpoint indices in
    the two transverse directions:

    - ``(0, 0)`` — both transverse axes at lo (e.g. j_lo_k_lo for axis i)
    - ``(1, 0)`` — first transverse hi, second lo (j_hi_k_lo)
    - ``(0, 1)`` — first lo, second hi              (j_lo_k_hi)
    - ``(1, 1)`` — both hi                          (j_hi_k_hi)

    Each ``edge_specs`` value is a ``Spacing1D`` describing the node
    distribution along that edge.  All 4 keys must be present.

    Parameters
    ----------
    n_axis :
        Number of nodes along the axis being built.
    n_t1, n_t2 :
        Number of nodes along the two transverse axes.
    start, stop :
        World coordinates of the axis interval.
    edge_specs :
        Mapping ``{(0,0), (1,0), (0,1), (1,1) → Spacing1D}``.
    axis_label :
        Used in error messages (``'i'`` / ``'j'`` / ``'k'``).

    Returns
    -------
    numpy.ndarray, shape ``(n_axis, n_t1, n_t2)``
    """
    expected = {(0, 0), (1, 0), (0, 1), (1, 1)}
    missing = expected - set(edge_specs.keys())
    if missing:
        raise ParameterError(
            f"axis_field_tfi[{axis_label}]: missing edge specs for {sorted(missing)}"
        )

    # Build each edge's parametric position t in [0, 1], shape (n_axis,).
    t_edges: dict[tuple[int, int], np.ndarray] = {}
    for key, sp in edge_specs.items():
        coords = sp.build(n_axis, start, stop)
        t_edges[key] = _normalise_to_unit(coords)

    # Bilinear blend weights, shape (n_t1, n_t2).
    eta = np.linspace(0.0, 1.0, n_t1)[:, None]  # (n_t1, 1)
    zeta = np.linspace(0.0, 1.0, n_t2)[None, :]  # (1,    n_t2)
    w00 = (1.0 - eta) * (1.0 - zeta)
    w10 = eta * (1.0 - zeta)
    w01 = (1.0 - eta) * zeta
    w11 = eta * zeta

    # Combined parametric field, shape (n_axis, n_t1, n_t2).
    t_field = (
        t_edges[(0, 0)][:, None, None] * w00[None, :, :]
        + t_edges[(1, 0)][:, None, None] * w10[None, :, :]
        + t_edges[(0, 1)][:, None, None] * w01[None, :, :]
        + t_edges[(1, 1)][:, None, None] * w11[None, :, :]
    )
    return start + (stop - start) * t_field


def resolve_edge_specs(
    axis: str,
    edge_pins: dict[str, Spacing1D] | None,
    axis_default: Spacing1D,
) -> dict[tuple[int, int], Spacing1D]:
    """Resolve the 4 parallel-edge ``Spacing1D``s for one axis.

    Parameters
    ----------
    axis :
        ``'i'`` / ``'j'`` / ``'k'``.
    edge_pins :
        Optional sparse user override: ``{edge_name: Spacing1D}``.  Edges
        not in the dict (or unmentioned for this axis) fall back to
        *axis_default*.  Edge names are the canonical 12 in
        :mod:`shore.mesh.edge`.
    axis_default :
        ``Spacing1D`` used for any edge not pinned in *edge_pins*.

    Returns
    -------
    dict ``{(0,0), (1,0), (0,1), (1,1) → Spacing1D}``

    Notes
    -----
    Edge-name → axis-key mapping (from :mod:`shore.mesh.block`):

    - axis ``i`` (j and k constant): j_lo_k_lo=(0,0), j_hi_k_lo=(1,0),
      j_lo_k_hi=(0,1), j_hi_k_hi=(1,1)  — first key = j, second = k
    - axis ``j`` (i and k constant): i_lo_k_lo=(0,0), i_hi_k_lo=(1,0),
      i_lo_k_hi=(0,1), i_hi_k_hi=(1,1)  — first key = i, second = k
    - axis ``k`` (i and j constant): i_lo_j_lo=(0,0), i_hi_j_lo=(1,0),
      i_lo_j_hi=(0,1), i_hi_j_hi=(1,1)  — first key = i, second = j
    """
    edge_name_map: dict[str, tuple[tuple[int, int], ...]] = {
        "i": (
            ("j_lo_k_lo", (0, 0)),
            ("j_hi_k_lo", (1, 0)),
            ("j_lo_k_hi", (0, 1)),
            ("j_hi_k_hi", (1, 1)),
        ),
        "j": (
            ("i_lo_k_lo", (0, 0)),
            ("i_hi_k_lo", (1, 0)),
            ("i_lo_k_hi", (0, 1)),
            ("i_hi_k_hi", (1, 1)),
        ),
        "k": (
            ("i_lo_j_lo", (0, 0)),
            ("i_hi_j_lo", (1, 0)),
            ("i_lo_j_hi", (0, 1)),
            ("i_hi_j_hi", (1, 1)),
        ),
    }
    if axis not in edge_name_map:
        raise ParameterError(f"resolve_edge_specs: axis must be 'i'/'j'/'k', got {axis!r}")

    pins = edge_pins or {}
    out: dict[tuple[int, int], Spacing1D] = {}
    for name, key in edge_name_map[axis]:
        out[key] = pins.get(name, axis_default)
    return out


def stamp_edge_specs_onto_block(
    block: HexBlock,
    *,
    edge_pins: dict[str, Spacing1D] | None,
) -> None:
    """Populate ``block.edge(name).spacing`` for every name in *edge_pins*.

    After construction, downstream readers (e.g. ``Mesh._resolve_k_steps``)
    may want to learn the actual per-edge schedule used.  This stamps the
    user's pins onto the block so they survive past ``build()``.
    """
    if not edge_pins:
        return
    for name, spec in edge_pins.items():
        block.edge(name).spacing = spec


# Suppress "unused" in the interim while the helper is being threaded through
# the primitive constructors.
__all__ = [
    "PrimitiveMesh",
    "axis_field_tfi",
    "axis_frame",
    "coords_1d",
    "resolve_axis",
    "resolve_edge_specs",
    "stamp_edge_specs_onto_block",
]
_ = warnings  # touched in concrete primitives


# ---------------------------------------------------------------------------
# Primitive ABC
# ---------------------------------------------------------------------------


class PrimitiveMesh(ABC):
    """Abstract base for parametric structured grids.

    Concrete subclasses populate ``self.blocks`` (a list of
    :class:`HexBlock`) and ``self.label_default`` (the file name to
    use when ``write_geo(stem)`` is called and the writer must pick a
    suffix for multi-block primitives).
    """

    blocks: list[HexBlock]
    label_default: str = "primitive"

    # ------------------------------------------------------------------
    # I/O — primitives write a *single* .geo per call by default; multi-
    # block primitives override ``write_geo`` to disambiguate.
    # ------------------------------------------------------------------

    def write_geo(self, output: str | Path) -> None:
        """Write the primitive to one or more ``.geo`` files.

        Single-block primitives write *output* verbatim (a single file).
        Multi-block primitives append ``_<label>.geo`` to *output*'s stem.
        """
        from shore.io.geo import write_geo

        out = Path(output)
        if len(self.blocks) == 1:
            # Treat *output* as a full file path (with or without .geo
            # extension).
            target = out if out.suffix == ".geo" else out.with_suffix(".geo")
            write_geo(target, self.blocks[0].nodes)
            return

        # Multi-block: *output* is a stem; one file per block.
        stem = out.with_suffix("") if out.suffix == ".geo" else out
        for blk in self.blocks:
            write_geo(stem.parent / f"{stem.name}_{blk.label}.geo", blk.nodes)

    def write_vtk(self, output: str | Path, with_jacobian: bool = True) -> None:
        """Write the primitive as a single ``.vtm`` MultiBlock file.

        Single- and multi-block primitives both produce one ``.vtm``;
        the difference is just the number of children inside it.
        """
        from shore.io.vtk import write_caps_vtm

        path = Path(output).with_suffix(".vtm")
        write_caps_vtm(path, [blk.nodes for blk in self.blocks], with_jacobian=with_jacobian)
