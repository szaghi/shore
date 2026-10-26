"""Meshing pipeline steps shared by the CLI and the command layer (ADR 0004).

``build_topology`` turns a topology name plus parameters into wired blocks.
It is the single code path behind ``shore mesh`` and
``Project.build_topology``, so a journal replay produces the same blocks as
the CLI run that recorded it.
"""

from __future__ import annotations

import inspect
from collections.abc import Callable
from typing import Any

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.body import Body
from shore.mesh.spacing import Spacing1D
from shore.mesh.topology import (
    CGridTopology,
    CHTopology,
    CubedSphereTopology,
    OGridTopology,
    OHTopology,
    Topology,
)

TOPOLOGIES = ("cubed_sphere", "ogrid", "cgrid", "ch", "oh")


def k_spacing(
    law: str = "geometric",
    *,
    ds: float = 1e-3,
    growth: float = 1.15,
    beta: float = 3.0,
    thickness: float | None = None,
) -> Spacing1D:
    """Wall-normal spacing law, as ``shore mesh --volume-k-*`` builds it."""
    if law in {"tanh", "tanh2"} and thickness is None:
        raise ParameterError(f"k-spacing law '{law}' needs a total thickness.")
    return Spacing1D(law=law, ds=ds, growth=growth, beta=beta, total_thickness=thickness)


def make_topology(name: str, options: dict[str, Any] | None = None) -> Topology:
    """Topology object for *name*; *options* are its constructor arguments.

    ``ch`` takes ``chord_min/max``, ``normal_min/max``, ``nk_channel``,
    ``freestream_world``, ``span_axis``, ``k_law_channel``,
    ``k_beta_channel``, ``i_split_body_upper/lower``, ``ni_outflow``,
    ``c1_seam``; ``oh`` takes ``box_vmin``, ``box_vmax``, ``nk_channel``,
    ``k_law_channel``, ``k_beta_channel``, ``c1_seam``; the others take none.
    """
    opts = dict(options or {})
    if name == "cubed_sphere":
        cls: type[Topology] = CubedSphereTopology
    elif name == "ogrid":
        cls = OGridTopology
    elif name == "cgrid":
        cls = CGridTopology
    elif name == "ch":
        cls = CHTopology
        if opts.get("freestream_world") is not None:
            opts["freestream_world"] = np.array(opts["freestream_world"], dtype=float)
    elif name == "oh":
        cls = OHTopology
    else:
        raise ParameterError(f"Unknown topology '{name}' (expected one of {TOPOLOGIES}).")
    _check_kwargs(cls, opts, f"topology '{name}' options")
    return cls(**opts)


def build_topology(
    body: Body,
    name: str,
    *,
    nk: int,
    spacing: Spacing1D,
    options: dict[str, Any] | None = None,
    **params: Any,
) -> tuple[Topology, list[HexBlock]]:
    """Build the wired blocks of topology *name* on *body*.

    *params* go to the topology's ``build``: ``ni``, ``nj``, ``theta_cap``,
    ``i_spacing``, ``i_beta``, ``j_spacing``, ``j_beta`` for the ray-cast
    topologies (``cubed_sphere``, ``ogrid``, ``oh``); ``ni_body``,
    ``ni_wake``, ``n_stations``, ``wake_length``, ``wake_growth``,
    ``span_axis``, ``te_dihedral_deg``, ``blocks`` and the ``cap_*``
    options for ``cgrid``; the C-grid subset for ``ch``.
    """
    topo = make_topology(name, options)
    _check_kwargs(topo.build, params, f"topology '{name}' build")
    if name in ("cgrid", "ch"):
        blocks = topo.build(surface=None, nk=nk, mesh=body.mesh, spacing_k=spacing, **params)
    else:
        blocks = topo.build(
            surface=None, nk=nk, spacing_k=spacing, projector=body.projector, **params
        )
    return topo, blocks


def _check_kwargs(fn: Callable[..., Any], kwargs: dict[str, Any], what: str) -> None:
    """Reject keyword names that are not named parameters of *fn*.

    The topology ``build`` methods end in a ``**kwargs`` catch-all; checking
    against the named parameters keeps a misspelt journal argument from
    being silently ignored.
    """
    named = {
        n
        for n, p in inspect.signature(fn).parameters.items()
        if p.kind not in (inspect.Parameter.VAR_KEYWORD, inspect.Parameter.VAR_POSITIONAL)
    }
    unknown = sorted(set(kwargs) - named)
    if unknown:
        raise ParameterError(f"{what}: unknown parameter(s) {', '.join(unknown)}.")


# ---------------------------------------------------------------------------
# Primitives, block files, splits, wing-tip caps
# ---------------------------------------------------------------------------

PRIMITIVES = ("box", "annulus", "flat_caps", "h_grid_rect", "h_grid_tube")


def build_primitive(kind: str, **params: Any) -> Any:
    """Build parametric primitive *kind* (``shore primitive``); returns the
    primitive object (``.blocks``, ``.write_geo``)."""
    from shore.primitive import (
        AnnularCylinderMesh,
        BoxMesh,
        FlatCapsCylinder,
        HGridRectMesh,
        HGridTubeMesh,
    )

    classes: dict[str, Any] = {
        "box": BoxMesh,
        "annulus": AnnularCylinderMesh,
        "flat_caps": FlatCapsCylinder,
        "h_grid_rect": HGridRectMesh,
        "h_grid_tube": HGridTubeMesh,
    }
    if kind not in classes:
        raise ParameterError(f"Unknown primitive '{kind}' (expected one of {PRIMITIVES}).")
    _check_kwargs(classes[kind].build, params, f"primitive '{kind}'")
    return classes[kind].build(**params)


def common_base_stem(stems: list[str]) -> str:
    """Longest common ``<base>_`` prefix of *stems*, ending at an underscore.

    A single stem gives everything up to and including its last underscore.
    """
    if len(stems) >= 2:
        common = stems[0]
        for s in stems[1:]:
            i = 0
            limit = min(len(common), len(s))
            while i < limit and common[i] == s[i]:
                i += 1
            common = common[:i]
        return common.rsplit("_", 1)[0] + "_" if "_" in common else ""
    return stems[0].rsplit("_", 1)[0] + "_" if "_" in stems[0] else ""


def load_blocks(
    geo_paths: list[str] | tuple[str, ...], adjacency_paths: list[str] | tuple[str, ...] = ()
) -> list[HexBlock]:
    """Blocks from ``.geo`` files, wired from adjacency sidecars (``shore split``).

    Labels come from the file stems with their longest common ``<base>_``
    prefix removed (``sphere_sub0.geo`` -> ``sub0``); a single file drops
    everything up to its last underscore.
    """
    from pathlib import Path

    from shore.io.adjacency import read_adjacency_json, wire_blocks_from_adjacency
    from shore.io.geo import read_geo

    paths = [Path(p) for p in geo_paths]
    stems = [p.stem for p in paths]
    base_stem = common_base_stem(stems)
    blocks = []
    for p, stem in zip(paths, stems, strict=True):
        label = stem[len(base_stem) :] if base_stem and stem.startswith(base_stem) else stem
        blocks.append(HexBlock(label=label, nodes=read_geo(p)))
    for adj in adjacency_paths:
        wire_blocks_from_adjacency(blocks, read_adjacency_json(adj), source=str(adj))
    return blocks


def apply_splits(
    blocks: list[HexBlock], splits: list[tuple[str, str, tuple[int, ...]]]
) -> list[HexBlock]:
    """Split *blocks* by ``(label, axis, at)`` directives in one pass."""
    from shore.split import BlockSplit, split_blocks

    labels = {b.label for b in blocks}
    for label, _, _ in splits:
        if label not in labels:
            raise ParameterError(
                f"split: no block named {label!r} in the input (have {sorted(labels)!r})"
            )
    return split_blocks(blocks, [BlockSplit(label=lb, axis=ax, at=at) for lb, ax, at in splits])


def build_tip_cap(
    body_geo: str,
    *,
    tip: str,
    le_frac: float,
    te_frac: float,
    span_axis: str,
    tip_inset: float | None,
    nk: int | None,
    n_cap_j: int,
    cap_j_spacing: str,
    cap_j_beta: float,
    auto_round_te: bool,
    te_thickness: float | None,
    te_thickness_frac: float | None,
) -> tuple[Any, dict[str, Any]]:
    """Wing-tip cap from a C-grid body block ``.geo`` (``shore cap-tip``).

    Returns the cap ``Mesh`` and the values resolved automatically
    (``tip_inset``: 1e-3 chord when *None*; ``nk``: the body's when *None*).
    """
    from shore.io.geo import read_geo
    from shore.mesh.mesh import Mesh
    from shore.mesh.topology import CGridTipCap
    from shore.topology.cgrid_cap import find_symmetric_le_cuts

    body = read_geo(body_geo)
    if body.ndim != 4 or body.shape[3] != 3:
        raise ParameterError(f"body geo has shape {body.shape}; expected (ni, n_stations, nk, 3).")
    if body.shape[1] < 2:
        raise ParameterError(
            f"body has {body.shape[1]} spanwise stations; "
            "need >= 2 to infer the outward span direction."
        )
    j_extreme = 0 if tip == "jlo" else body.shape[1] - 1
    j_inboard = 1 if tip == "jlo" else body.shape[1] - 2
    cross_section = body[:, j_extreme, 0, :].copy()  # (ni, 3) at the wall (k=0)
    inboard_dir = body[:, j_inboard, 0, :].mean(axis=0) - cross_section.mean(axis=0)
    inboard_norm = float(np.linalg.norm(inboard_dir))
    if inboard_norm < 1e-12:
        raise ParameterError(
            "adjacent spanwise stations have coincident centroids — cannot infer span direction."
        )
    span_dir = -inboard_dir / inboard_norm
    span_axis_idx = "xyz".index(span_axis)
    le_up, le_lo = find_symmetric_le_cuts(
        cross_section, le_frac=le_frac, span_axis_world=span_axis_idx
    )
    te_up, te_lo = find_symmetric_le_cuts(
        cross_section, le_frac=te_frac, span_axis_world=span_axis_idx
    )
    resolved: dict[str, Any] = {}
    if tip_inset is None:
        extents = cross_section.max(axis=0) - cross_section.min(axis=0)
        extents[span_axis_idx] = -np.inf
        tip_inset = 1e-3 * float(extents.max())
        resolved["tip_inset"] = tip_inset
    if nk is None:
        nk = int(body.shape[2])
        resolved["nk"] = nk
    topo = CGridTipCap(
        le_cut_upper_idx=le_up,
        le_cut_lower_idx=le_lo,
        te_cut_upper_idx=te_up,
        te_cut_lower_idx=te_lo,
        span_dir=tuple(span_dir.tolist()),
        tip_inset=tip_inset,
        n_cap_j=n_cap_j,
        cap_j_spacing=cap_j_spacing,
        cap_j_beta=cap_j_beta,
        is_blunt=False,
        auto_round_te=auto_round_te,
        te_thickness=te_thickness,
        te_thickness_frac=None if te_thickness is not None else te_thickness_frac,
    )
    return Mesh.from_array(surface=cross_section, topology=topo, nk=nk), resolved
