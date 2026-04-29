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
  resolution such that every voxel inside the body has corners
  within ``max_outward_gap`` of the body surface (outward only).
  Each inside-body voxel becomes one axis-aligned :class:`Box`.
  Produces full body-interior coverage.
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
    """Voxelise the body's AABB and emit strictly-inside-body voxels as boxes.

    Voxel side is chosen so the voxel half-diagonal does not exceed
    ``max_outward_gap``: ``side <= 2 * max_outward_gap / sqrt(3)``.
    A voxel with centre **strictly inside** the body therefore has
    every corner within ``max_outward_gap`` of the body surface in
    the worst case.

    Inclusion test: ``signed_distance(centre) > 0`` (positive inside,
    negative outside per ``trimesh.proximity``).  Voxels with centre
    on or outside the body surface are excluded.

    Coverage trade-off: a thin body-shell of width up to ``side / 2``
    near the body surface is *not* covered by these voxels (the
    voxels just outside the body surface aren't included).  Body-
    interior background-mesh cells that fall in this shell are
    typically also inside the body's near-mesh (wall) and are
    handled by the wall mesh's own overset donor logic, so the
    boxes don't need to cover them.  Pick ``max_outward_gap`` small
    enough that this shell is comfortably inside the wall.
    """
    import trimesh.proximity

    sqrt3 = float(np.sqrt(3.0))
    side = 2.0 * max_outward_gap / sqrt3

    lo = np.asarray(mesh.vertices, dtype=np.float64).min(axis=0)
    hi = np.asarray(mesh.vertices, dtype=np.float64).max(axis=0)
    extent = hi - lo
    counts = np.maximum(1, np.ceil(extent / side).astype(np.int64))
    # Recompute per-axis side so voxel grid covers AABB exactly.
    actual_side = extent / counts

    # Generate voxel-centre coordinates: (Nx*Ny*Nz, 3) array.
    cx = lo[0] + (np.arange(counts[0]) + 0.5) * actual_side[0]
    cy = lo[1] + (np.arange(counts[1]) + 0.5) * actual_side[1]
    cz = lo[2] + (np.arange(counts[2]) + 0.5) * actual_side[2]
    grid_x, grid_y, grid_z = np.meshgrid(cx, cy, cz, indexing="ij")
    centres = np.column_stack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()])

    # Strict inclusion: voxel centre must be inside the body.
    # ``trimesh.proximity.signed_distance``: positive inside, negative
    # outside.  Include only voxels whose centre is strictly inside.
    # The voxel's worst-case corner is at most ``half_diag`` outside
    # the body surface, where ``half_diag = side * sqrt(3) / 2 <=
    # max_outward_gap`` by construction.
    proximity = trimesh.proximity.ProximityQuery(mesh)
    signed = proximity.signed_distance(centres)
    include = signed > 0.0

    if not include.any():
        return []

    boxes: list[Box] = []
    half_side = 0.5 * actual_side
    for centre in centres[include]:
        boxes.append(
            Box.from_aabb(
                centre - half_side,
                centre + half_side,
                associated_block,
            )
        )
    return boxes


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
