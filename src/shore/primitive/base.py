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

from abc import ABC
from pathlib import Path

import numpy as np

from shore.mesh.block import HexBlock
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
