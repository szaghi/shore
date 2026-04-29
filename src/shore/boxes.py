"""Generate :class:`shore.io.boxes.Box` instances from an STL surface.

The Xall ``cc.par`` boxes section needs hexahedra that *cover the
body interior*: every cell of the background mesh whose centre falls
inside the solid body must lie inside at least one box, so the
chimera hole-cutter can mark it as donors-only.  See
:mod:`shore.io.boxes` for the role boxes play.

This module turns a triangulated body into a set of such hexahedra.
The default method (``"voxel-fill"``) voxelises the body's AABB at
a resolution chosen from ``max_outward_gap`` and emits each voxel
inside the body as an axis-aligned :class:`Box`.  Coverage is
guaranteed by construction; the only failure mode is excessive box
count for tiny gap budgets relative to the body size.

Methods
-------
- ``"voxel-fill"`` (default) — voxelise the body's AABB at a
  resolution such that every body-overlapping voxel has corners
  within ``max_outward_gap`` of the body surface (outward only).
  Each voxel that overlaps the body becomes one axis-aligned
  :class:`Box`; together they form a watertight cover of the
  body (interior + surface shell), suitable for chimera
  hole-cutting.
- ``"aabb"``  — single world-axis-aligned bounding box of the body.
  Best effort: emits a ``UserWarning`` when the measured outward
  gap exceeds ``max_outward_gap``.  Useful as a quick preview only.
- ``"obb"``   — single oriented bounding box (area-weighted PCA).
  Best effort with the same warning behaviour as ``aabb``.

Outward gap (best-effort methods only)
--------------------------------------
``aabb`` and ``obb`` measure the actual distance from sampled
box-surface points (8 corners + a parametric ``kxk`` grid on each
face) to the STL surface (outward only) and compare to
``max_outward_gap``.  ``voxel-fill`` enforces the gap by its voxel-
size choice and does not need to measure.

Watertightness
--------------
``voxel-fill`` calls ``trimesh.contains()`` which assumes a
watertight body STL.  Non-watertight inputs return unreliable
inside/outside classifications; :func:`shore.surface.io.load_stl`
already warns on non-watertight meshes, and :func:`boxes_from_stl`
re-checks and raises if voxel-fill is requested on one.

API
---
.. autofunction:: boxes_from_stl
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Literal

import numpy as np

from shore.exceptions import MeshInputError, ParameterError
from shore.io.boxes import Box
from shore.surface.io import load_stl

Method = Literal["aabb", "obb", "voxel-fill"]


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def boxes_from_stl(
    stl_path: str | Path,
    associated_block: str,
    *,
    max_outward_gap: float,
    method: Method = "voxel-fill",
    face_samples_per_axis: int = 4,
) -> list[Box]:
    """Build a list of :class:`Box` covering the STL body interior.

    Parameters
    ----------
    stl_path :
        Path to the body's STL file.  Loaded via
        :func:`shore.surface.io.load_stl` (with hygiene).
    associated_block :
        SHORE block label every returned box references.  Resolved to
        a 1-based index by :func:`shore.io.cc_par.write_cc_par`.
    max_outward_gap :
        Maximum allowed distance, in STL world units, from any point
        on a box's surface to the nearest STL surface point —
        measured **outward only** (points on the box surface that lie
        inside the body don't count).  Strictly positive.  Should be
        smaller than the body's near-mesh wall-normal thickness so
        the boxes fit inside the wall.
    method :
        Algorithm: ``"voxel-fill"`` (default — guaranteed body-volume
        coverage by voxelising the AABB), ``"aabb"`` / ``"obb"``
        (single-box, best-effort).
    face_samples_per_axis :
        ``"aabb"`` / ``"obb"`` post-check sampling density.  Each of
        the 6 box faces is sampled on a ``kxk`` grid (total = 8
        corners + 6*k*k samples per box) for the outward-gap
        measurement.  Default 4 (104 samples/box).  Ignored by
        ``voxel-fill``.

    Returns
    -------
    list[Box]
        One or more boxes, all sharing ``associated_block``.

    Raises
    ------
    ParameterError
        If parameters are out of range, or ``method`` is unknown.
    MeshInputError
        Forwarded from :func:`load_stl`, or raised when
        ``method='voxel-fill'`` on a non-watertight STL.
    """
    if max_outward_gap <= 0.0:
        raise ParameterError(f"boxes_from_stl: max_outward_gap must be > 0, got {max_outward_gap}")
    if face_samples_per_axis < 1:
        raise ParameterError(
            f"boxes_from_stl: face_samples_per_axis must be >= 1, got {face_samples_per_axis}"
        )
    if method not in ("aabb", "obb", "voxel-fill"):
        raise ParameterError(
            f"boxes_from_stl: method must be 'aabb' | 'obb' | 'voxel-fill', got {method!r}"
        )

    mesh = load_stl(stl_path)
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    if vertices.shape[0] < 4:
        raise MeshInputError(
            f"boxes_from_stl: STL has only {vertices.shape[0]} vertices "
            "(need >= 4 to bound a hexahedron)"
        )

    if method == "voxel-fill":
        if not mesh.is_watertight:
            raise MeshInputError(
                f"boxes_from_stl: method='voxel-fill' requires a watertight "
                f"STL (its inside/outside test relies on closed-surface "
                f"ray casting); {stl_path} is not watertight.  Repair the "
                f"STL or use method='aabb'/'obb' as a best-effort fallback."
            )
        return _voxel_fill(
            mesh,
            associated_block=associated_block,
            max_outward_gap=max_outward_gap,
        )

    # Best-effort single-box methods.
    import trimesh.proximity

    proximity = trimesh.proximity.ProximityQuery(mesh)

    if method == "aabb":
        return [
            _box_from_aabb(
                vertices,
                associated_block,
                max_outward_gap=max_outward_gap,
                proximity=proximity,
                face_samples_per_axis=face_samples_per_axis,
                source_path=str(stl_path),
            )
        ]

    # method == "obb"
    return [
        _box_from_obb(
            vertices,
            _area_weights(mesh),
            associated_block,
            max_outward_gap=max_outward_gap,
            proximity=proximity,
            face_samples_per_axis=face_samples_per_axis,
            source_path=str(stl_path),
        )
    ]


# ---------------------------------------------------------------------------
# Voxel-fill (default; guarantees body-volume coverage)
# ---------------------------------------------------------------------------


def _voxel_fill(
    mesh,
    *,
    associated_block: str,
    max_outward_gap: float,
) -> list[Box]:
    """Voxelise the body's AABB, then greedy-merge the cover.

    Stage 1 — seed grid.  Voxel side is chosen so the voxel
    half-diagonal does not exceed ``max_outward_gap``:
    ``side <= max_outward_gap / sqrt(3)``.  Inclusion test:
    ``signed_distance(centre) > -half_diag``, which admits both
    centre-inside voxels and body-shell voxels (centre outside but
    within ``half_diag`` of the surface).  Together they form a
    watertight cover of the body, and every voxel's corners sit at
    most ``max_outward_gap`` outside the body surface.

    Stage 2 — greedy AABB merge.  The seed cover is wasteful: a
    cubic grid forces all three axes to the same fine resolution
    even when the local body geometry would tolerate a coarser slab.
    :func:`_greedy_merge_grid` repeatedly merges face-adjacent grid
    blocks into larger AABBs, accepting each merge only when the
    merged block's 8 corners *all* satisfy the outward-gap budget.
    Body coverage is preserved by construction (a merge is a strict
    superset of its parts), and the gap budget is rechecked at the
    enlarged corners.  Order: merge along the longest axis first to
    favour slab-like outputs.
    """
    import trimesh.proximity

    sqrt3 = float(np.sqrt(3.0))
    side = max_outward_gap / sqrt3

    lo = np.asarray(mesh.vertices, dtype=np.float64).min(axis=0)
    hi = np.asarray(mesh.vertices, dtype=np.float64).max(axis=0)
    extent = hi - lo
    counts = np.maximum(1, np.ceil(extent / side).astype(np.int64))
    actual_side = extent / counts
    half_diag = 0.5 * float(np.linalg.norm(actual_side))

    cx = lo[0] + (np.arange(counts[0]) + 0.5) * actual_side[0]
    cy = lo[1] + (np.arange(counts[1]) + 0.5) * actual_side[1]
    cz = lo[2] + (np.arange(counts[2]) + 0.5) * actual_side[2]
    grid_x, grid_y, grid_z = np.meshgrid(cx, cy, cz, indexing="ij")
    centres = np.column_stack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()])

    proximity = trimesh.proximity.ProximityQuery(mesh)
    signed = proximity.signed_distance(centres)
    include = (signed > -half_diag).reshape(tuple(int(c) for c in counts))

    if not include.any():
        return []

    blocks = _greedy_merge_grid(
        include,
        lo=lo,
        side=actual_side,
        proximity=proximity,
        max_outward_gap=max_outward_gap,
    )

    return [Box.from_aabb(box_lo, box_hi, associated_block) for box_lo, box_hi in blocks]


def _greedy_merge_grid(
    include: np.ndarray,
    *,
    lo: np.ndarray,
    side: np.ndarray,
    proximity,
    max_outward_gap: float,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Greedy axis-aligned merge over a regular voxel cover.

    Parameters
    ----------
    include :
        Bool array of shape ``(Nx, Ny, Nz)``; True voxels form the
        seed cover and must be tiled by the returned AABBs.
    lo, side :
        World-space origin and per-axis voxel side.  World coordinates
        of the AABB ``[i0:i1, j0:j1, k0:k1]`` are
        ``lo + (i0, j0, k0) * side`` (low corner) and
        ``lo + (i1, j1, k1) * side`` (high corner).
    proximity :
        ``trimesh.proximity.ProximityQuery`` for the body.
    max_outward_gap :
        Per-vertex outward-gap budget that every merged AABB's 8
        corners must satisfy.

    Returns
    -------
    list of ``(lo3, hi3)``
        World-space low/high corners of the merged AABBs.  Their
        union equals the union of the seed voxels.

    Algorithm
    ---------
    Each block is a half-open grid range ``[i0:i1, j0:j1, k0:k1]``
    plus a flag for whether the block is fully populated (only fully
    populated blocks are merge candidates; seed voxels start as
    unit-extent fully populated blocks).

    Two blocks are merge-eligible iff:
      1. they are face-adjacent along one axis ``a`` (one's high
         coord on axis ``a`` equals the other's low coord on ``a``);
      2. their grid ranges on the other two axes match exactly;
      3. the merged AABB's 8 world-space corners *all* have outward
         gap to the body <= ``max_outward_gap``.

    Pass order: along each axis in turn, in descending order of the
    seed grid's extent on that axis.  This biases the result toward
    slab-shaped outputs when the body's elongation differs by axis.
    """
    # Active blocks as (i0, i1, j0, j1, k0, k1).  Initial seed: one
    # block per True voxel.
    blocks: list[tuple[int, int, int, int, int, int]] = []
    ii, jj, kk = np.where(include)
    for i, j, k in zip(ii.tolist(), jj.tolist(), kk.tolist(), strict=True):
        blocks.append((i, i + 1, j, j + 1, k, k + 1))

    # Merge order: longest seed-axis first.
    axes_order = sorted(range(3), key=lambda a: -include.shape[a])

    for axis in axes_order:
        blocks = _merge_pass_axis(
            blocks,
            axis=axis,
            lo=lo,
            side=side,
            proximity=proximity,
            max_outward_gap=max_outward_gap,
        )

    # Convert grid-range blocks to world-space (lo3, hi3) pairs.
    out: list[tuple[np.ndarray, np.ndarray]] = []
    for i0, i1, j0, j1, k0, k1 in blocks:
        block_lo = lo + np.array([i0, j0, k0]) * side
        block_hi = lo + np.array([i1, j1, k1]) * side
        out.append((block_lo, block_hi))
    return out


def _merge_pass_axis(
    blocks: list[tuple[int, int, int, int, int, int]],
    *,
    axis: int,
    lo: np.ndarray,
    side: np.ndarray,
    proximity,
    max_outward_gap: float,
) -> list[tuple[int, int, int, int, int, int]]:
    """One full merge sweep along ``axis``.

    Group blocks by their (other-two-axes) grid range.  Inside each
    group, sort by axis-low coord and walk a sweep, merging
    consecutive face-adjacent blocks whenever the merged corners fit
    the gap budget.  Repeat the sweep until a pass makes no merges.
    """
    other_axes = tuple(a for a in range(3) if a != axis)

    while True:
        groups: dict[tuple[int, int, int, int], list[int]] = {}
        for idx, b in enumerate(blocks):
            key_lo_a = b[2 * other_axes[0]]
            key_hi_a = b[2 * other_axes[0] + 1]
            key_lo_b = b[2 * other_axes[1]]
            key_hi_b = b[2 * other_axes[1] + 1]
            key = (key_lo_a, key_hi_a, key_lo_b, key_hi_b)
            groups.setdefault(key, []).append(idx)

        merged_any = False
        keep = np.ones(len(blocks), dtype=bool)
        new_blocks: list[tuple[int, int, int, int, int, int]] = []

        for indices in groups.values():
            if len(indices) < 2:
                continue
            # Sort group by low coord on merge axis.
            indices.sort(key=lambda i: blocks[i][2 * axis])

            # Sweep: try to extend a running block forward.
            run_idx = indices[0]
            run = blocks[run_idx]
            run_keep = [run_idx]
            for nxt_idx in indices[1:]:
                nxt = blocks[nxt_idx]
                # Face-adjacency on merge axis.
                if run[2 * axis + 1] != nxt[2 * axis]:
                    # Gap; flush run, start fresh.
                    if len(run_keep) > 1:
                        merged_any = True
                        for i in run_keep:
                            keep[i] = False
                        new_blocks.append(run)
                    run = nxt
                    run_idx = nxt_idx
                    run_keep = [nxt_idx]
                    continue

                # Merged extent: low coords from run, high coords from
                # run except along merge axis where it takes nxt's high.
                candidate = (
                    run[0],
                    nxt[1] if axis == 0 else run[1],
                    run[2],
                    nxt[3] if axis == 1 else run[3],
                    run[4],
                    nxt[5] if axis == 2 else run[5],
                )
                if _block_corners_within_gap(
                    candidate,
                    lo=lo,
                    side=side,
                    proximity=proximity,
                    max_outward_gap=max_outward_gap,
                ):
                    run = candidate
                    run_keep.append(nxt_idx)
                else:
                    if len(run_keep) > 1:
                        merged_any = True
                        for i in run_keep:
                            keep[i] = False
                        new_blocks.append(run)
                    run = nxt
                    run_idx = nxt_idx
                    run_keep = [nxt_idx]

            if len(run_keep) > 1:
                merged_any = True
                for i in run_keep:
                    keep[i] = False
                new_blocks.append(run)

        if not merged_any:
            return blocks

        blocks = [b for i, b in enumerate(blocks) if keep[i]] + new_blocks


def _block_corners_within_gap(
    block: tuple[int, int, int, int, int, int],
    *,
    lo: np.ndarray,
    side: np.ndarray,
    proximity,
    max_outward_gap: float,
) -> bool:
    """True iff every world-space corner of ``block`` has outward gap <= budget."""
    i0, i1, j0, j1, k0, k1 = block
    block_lo = lo + np.array([i0, j0, k0]) * side
    block_hi = lo + np.array([i1, j1, k1]) * side
    corners = np.array(
        [
            [block_lo[0], block_lo[1], block_lo[2]],
            [block_hi[0], block_lo[1], block_lo[2]],
            [block_lo[0], block_hi[1], block_lo[2]],
            [block_hi[0], block_hi[1], block_lo[2]],
            [block_lo[0], block_lo[1], block_hi[2]],
            [block_hi[0], block_lo[1], block_hi[2]],
            [block_lo[0], block_hi[1], block_hi[2]],
            [block_hi[0], block_hi[1], block_hi[2]],
        ]
    )
    signed = proximity.signed_distance(corners)
    outward = np.maximum(0.0, -signed)
    return bool(outward.max() <= max_outward_gap + 1e-12)


# ---------------------------------------------------------------------------
# AABB / OBB best-effort single-box methods
# ---------------------------------------------------------------------------


def _box_from_aabb(
    vertices: np.ndarray,
    associated_block: str,
    *,
    max_outward_gap: float,
    proximity,
    face_samples_per_axis: int,
    source_path: str = "",
) -> Box:
    """Single world-axis-aligned :class:`Box` covering all vertices.

    Best effort: post-checks the measured outward gap against
    ``max_outward_gap`` and emits a ``UserWarning`` when exceeded.
    """
    lo = vertices.min(axis=0)
    hi = vertices.max(axis=0)
    box = Box.from_aabb(lo, hi, associated_block)
    _post_check_gap(
        box,
        max_outward_gap=max_outward_gap,
        proximity=proximity,
        face_samples_per_axis=face_samples_per_axis,
        source_path=source_path,
        method_name="aabb",
    )
    return box


def _box_from_obb(
    vertices: np.ndarray,
    weights: np.ndarray,
    associated_block: str,
    *,
    max_outward_gap: float,
    proximity,
    face_samples_per_axis: int,
    source_path: str = "",
) -> Box:
    """Single tight OBB :class:`Box` from area-weighted PCA.

    Best effort with the same gap post-check as :func:`_box_from_aabb`.
    """
    axes, centre = _principal_axes_weighted(vertices, weights)
    centred = vertices - centre
    projected = centred @ axes
    lo = projected.min(axis=0)
    hi = projected.max(axis=0)
    half = 0.5 * (hi - lo)
    box_centre = centre + (0.5 * (lo + hi)) @ axes.T

    verts_local = np.empty((8, 3), dtype=np.float64)
    for idx in range(8):
        i = (idx & 1) * 2 - 1
        j = ((idx >> 1) & 1) * 2 - 1
        k = ((idx >> 2) & 1) * 2 - 1
        verts_local[idx] = (i * half[0], j * half[1], k * half[2])
    verts_world = verts_local @ axes.T + box_centre

    box = Box(vertices=verts_world, associated_block=associated_block)
    _post_check_gap(
        box,
        max_outward_gap=max_outward_gap,
        proximity=proximity,
        face_samples_per_axis=face_samples_per_axis,
        source_path=source_path,
        method_name="obb",
    )
    return box


def _area_weights(mesh) -> np.ndarray:
    """Per-vertex area weights: 1/3 of the sum of incident triangle areas.

    PCA on raw vertex coordinates is biased by tessellation density
    (an over-tessellated patch dominates the covariance).  Weighting
    each vertex by its incident triangle area cancels that bias.
    """
    faces = np.asarray(mesh.faces, dtype=np.int64)
    areas = np.asarray(mesh.area_faces, dtype=np.float64)
    n_vertices = mesh.vertices.shape[0]
    weights = np.zeros(n_vertices, dtype=np.float64)
    for v in range(3):
        np.add.at(weights, faces[:, v], areas / 3.0)
    total = weights.sum()
    if total <= 0.0:
        weights = np.ones(n_vertices, dtype=np.float64) / n_vertices
    else:
        weights /= total
    return weights


def _principal_axes_weighted(
    vertices: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (axes, centre) for an area-weighted PCA on *vertices*."""
    centre = (weights[:, None] * vertices).sum(axis=0)
    centred = vertices - centre
    cov = (centred * weights[:, None]).T @ centred
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = np.argsort(eigenvalues)[::-1]
    axes = eigenvectors[:, order]
    if np.linalg.det(axes) < 0:
        axes[:, 2] = -axes[:, 2]
    return axes, centre


# ---------------------------------------------------------------------------
# Outward-gap measurement (used by aabb/obb post-check)
# ---------------------------------------------------------------------------


def _max_outward_gap(
    box_vertices: np.ndarray,
    *,
    proximity,
    face_samples_per_axis: int,
) -> float:
    """Largest outward distance from any sampled box-surface point to the STL.

    Sample set: the 8 corners + a parametric ``kxk`` grid on each
    face.  ``trimesh.proximity.signed_distance`` returns positive
    values inside the body, negative outside; outward gap is
    ``max(0, -signed)``.
    """
    samples = _box_surface_samples(box_vertices, face_samples_per_axis)
    signed = proximity.signed_distance(samples)
    outward = np.maximum(0.0, -signed)
    return float(outward.max())


def _box_surface_samples(
    box_vertices: np.ndarray,
    face_samples_per_axis: int,
) -> np.ndarray:
    """8 corners + (kxk grid) per face, in world coordinates.

    Bilinear interpolation over each face's 4 corners works for any
    hexahedron orientation, including skewed boxes.
    """
    samples = [box_vertices]

    face_quads: tuple[tuple[int, int, int, int], ...] = (
        (0, 2, 4, 6),  # i_lo
        (1, 3, 5, 7),  # i_hi
        (0, 1, 4, 5),  # j_lo
        (2, 3, 6, 7),  # j_hi
        (0, 1, 2, 3),  # k_lo
        (4, 5, 6, 7),  # k_hi
    )
    k = face_samples_per_axis
    if k >= 1:
        u = np.linspace(0.0, 1.0, k + 2)[1:-1] if k >= 2 else np.array([0.5])
        for q0, q1, q2, q3 in face_quads:
            v00 = box_vertices[q0]
            v10 = box_vertices[q1]
            v01 = box_vertices[q2]
            v11 = box_vertices[q3]
            uu, vv = np.meshgrid(u, u, indexing="ij")
            uu = uu.ravel()
            vv = vv.ravel()
            face_pts = (
                (1 - uu)[:, None] * (1 - vv)[:, None] * v00
                + uu[:, None] * (1 - vv)[:, None] * v10
                + (1 - uu)[:, None] * vv[:, None] * v01
                + uu[:, None] * vv[:, None] * v11
            )
            samples.append(face_pts)
    return np.vstack(samples)


def _post_check_gap(
    box: Box,
    *,
    max_outward_gap: float,
    proximity,
    face_samples_per_axis: int,
    source_path: str,
    method_name: str,
) -> None:
    """Best-effort outward-gap check for single-box methods.

    Emits a ``UserWarning`` when the measured outward gap exceeds the
    budget; never raises.  Caller decides whether to act on the warning.
    """
    gap = _max_outward_gap(
        box.vertices,
        proximity=proximity,
        face_samples_per_axis=face_samples_per_axis,
    )
    if gap > max_outward_gap:
        prefix = f"boxes_from_stl ({source_path})" if source_path else "boxes_from_stl"
        warnings.warn(
            f"{prefix}: method={method_name!r} produced a single box whose "
            f"max outward gap is {gap:.4g} > max_outward_gap "
            f"{max_outward_gap:.4g}.  This box has corners that sit further "
            "from the body surface than the user budget; the box may not "
            "fit inside the body's near-mesh.  Switch to "
            "method='voxel-fill' for a multi-box cover that satisfies "
            "the gap budget *and* covers the body interior.",
            UserWarning,
            stacklevel=4,
        )
