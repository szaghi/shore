"""Hyperbolic marching primitives.

Provides the per-layer normal estimation, Laplacian smoothing, and STL
surface-normal query used by the OOP mesh marcher (``shore.mesh.Mesh``).
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Normal estimation
# ---------------------------------------------------------------------------


def _compute_normals(
    layer: np.ndarray,
    periodic: tuple[bool, bool] = (False, True),
    ghost_rows: dict[str, np.ndarray] | None = None,
) -> np.ndarray:
    """Estimate outward unit normals for every node of a structured layer.

    Central differences in both directions; one-sided at non-periodic boundaries.
    When ``ghost_rows`` is supplied, the boundary in question uses a centred
    difference against the ghost row instead of a one-sided fallback — this
    preserves spacing continuity across multi-block seams during marching.

    Parameters
    ----------
    layer :
        Node coordinates, shape (ni, nj, 3).
    periodic :
        ``(periodic_i, periodic_j)``.  When ``True`` the corresponding
        direction uses modular (wrap-around) indexing for the end nodes;
        when ``False`` one-sided differences are used at the boundaries
        unless a ghost row is supplied for that side.
    ghost_rows :
        Optional virtual neighbour rows for any of the 4 faces.  Keys:
        ``'i_lo'`` / ``'i_hi'`` give shape ``(nj, 3)``; ``'j_lo'`` / ``'j_hi'``
        give shape ``(ni, 3)``.  When a key is present the corresponding
        boundary's tangent becomes the centred difference between the next
        interior row and the ghost row.

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
    """
    ni, nj = layer.shape[0], layer.shape[1]
    periodic_i, periodic_j = periodic
    g = ghost_rows or {}

    # --- i-direction tangent ---
    ti = np.empty_like(layer)
    if periodic_i:
        ip1 = (np.arange(ni) + 1) % ni
        im1 = (np.arange(ni) - 1) % ni
        ti = layer[ip1] - layer[im1]
    else:
        ti[1:-1] = layer[2:] - layer[:-2]
        if "i_lo" in g:
            ti[0] = layer[1] - g["i_lo"]
        else:
            ti[0] = layer[1] - layer[0]
        if "i_hi" in g:
            ti[-1] = g["i_hi"] - layer[-2]
        else:
            ti[-1] = layer[-1] - layer[-2]

    # --- j-direction tangent ---
    if periodic_j:
        jp1 = (np.arange(nj) + 1) % nj
        jm1 = (np.arange(nj) - 1) % nj
        tj = layer[:, jp1] - layer[:, jm1]
    else:
        tj = np.empty_like(layer)
        tj[:, 1:-1] = layer[:, 2:] - layer[:, :-2]
        if "j_lo" in g:
            tj[:, 0] = layer[:, 1] - g["j_lo"]
        else:
            tj[:, 0] = layer[:, 1] - layer[:, 0]
        if "j_hi" in g:
            tj[:, -1] = g["j_hi"] - layer[:, -2]
        else:
            tj[:, -1] = layer[:, -1] - layer[:, -2]

    normals = np.cross(ti, tj)
    mag = np.linalg.norm(normals, axis=-1, keepdims=True)
    mag = np.where(mag < 1e-15, 1.0, mag)
    return normals / mag


# ---------------------------------------------------------------------------
# Laplacian smoothing
# ---------------------------------------------------------------------------


def _laplacian_smooth(
    layer: np.ndarray,
    strength: float,
    periodic: tuple[bool, bool] = (False, True),
) -> np.ndarray:
    """One step of in-plane Laplacian smoothing.

    Boundary rows/columns in non-periodic directions are clamped (never moved).
    The k=0 wall row (index 0) is always clamped regardless of periodicity.

    Parameters
    ----------
    layer :
        Shape (ni, nj, 3).
    strength :
        Blend factor in [0, 1].
    periodic :
        ``(periodic_i, periodic_j)``.  Non-periodic boundaries are clamped.
        Default ``(False, True)`` preserves the original O-grid behaviour.

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
    """
    nj = layer.shape[1]
    periodic_i, periodic_j = periodic

    # j neighbours
    if periodic_j:
        jp1 = (np.arange(nj) + 1) % nj
        jm1 = (np.arange(nj) - 1) % nj
        j_fwd = layer[:, jp1]
        j_bwd = layer[:, jm1]
    else:
        j_fwd = np.empty_like(layer)
        j_bwd = np.empty_like(layer)
        j_fwd[:, :-1] = layer[:, 1:]
        j_fwd[:, -1] = layer[:, -1]  # clamp: boundary contributes itself
        j_bwd[:, 1:] = layer[:, :-1]
        j_bwd[:, 0] = layer[:, 0]

    # i neighbours
    if periodic_i:
        i_fwd = np.roll(layer, -1, axis=0)
        i_bwd = np.roll(layer, 1, axis=0)
    else:
        i_fwd = np.empty_like(layer)
        i_bwd = np.empty_like(layer)
        i_fwd[:-1] = layer[1:]
        i_fwd[-1] = layer[-1]
        i_bwd[1:] = layer[:-1]
        i_bwd[0] = layer[0]

    laplacian = (j_fwd + j_bwd + i_fwd + i_bwd) / 4.0
    smoothed = layer + strength * (laplacian - layer)

    # always clamp the wall row
    smoothed[0] = layer[0]
    # clamp non-periodic boundaries
    if not periodic_i:
        smoothed[-1] = layer[-1]
    if not periodic_j:
        smoothed[:, 0] = layer[:, 0]
        smoothed[:, -1] = layer[:, -1]

    return smoothed


SMOOTH_MODES = ("limited", "legacy")


def _limit_smoothing(
    raw: np.ndarray,
    smoothed: np.ndarray,
    normals: np.ndarray,
    step: float,
    max_pullback: float = 0.25,
    max_skew: float = 0.25,
) -> np.ndarray:
    """Limit the Laplacian smoothing motion where it would damage the layer.

    Smoothing node *positions* on a curved layer also moves them along the
    normal: the layer shrinks toward its centre of curvature by about
    ``strength * h**2 * kappa``, independently of the marching step. When
    that inward motion exceeds the step the new layer falls behind the
    previous one and the cells fold (negative Jacobian at small first
    steps); large tangential motion over a thin layer skews it.

    Per node, the smoothing motion ``smoothed - raw`` is split into normal
    and tangential parts and limited:

    - the normal part may pull the node back by at most
      ``max_pullback * step``, so the layer always advances by at least
      ``(1 - max_pullback) * step``;
    - the tangential part is capped at ``max_skew * step``.

    The defaults (0.25, 0.25) are the tightest limits that leave the
    60x60x10 cubed-sphere regression grid bit-identical to unlimited
    smoothing; on a coarse 10x12 sphere with a first step of 0.3 % of the
    cell size they keep the worst wall-cell tilt at 19-24 degrees (8 degrees
    unsmoothed, 64-68 degrees at (0.5, 1.0)).

    Nodes where neither limit binds keep ``smoothed`` exactly, so grids that
    marched cleanly before are unchanged.

    Parameters
    ----------
    raw :
        Layer after the normal step, before smoothing, (ni, nj, 3).
    smoothed :
        The same layer after the Laplacian sweeps.
    normals :
        Unit marching normals used for the step, (ni, nj, 3).
    step :
        Marching step of this layer.
    max_pullback, max_skew :
        Limits, in units of the step.
    """
    delta = smoothed - raw
    dn = np.einsum("...j,...j->...", delta, normals)
    dt = delta - dn[..., None] * normals
    mag = np.linalg.norm(dt, axis=-1)
    h = abs(step)
    bind = (dn < -max_pullback * h) | (mag > max_skew * h)
    if not bind.any():
        return smoothed
    dn_lim = np.maximum(dn, -max_pullback * h)
    scale = np.where(mag > max_skew * h, max_skew * h / np.where(mag > 0.0, mag, 1.0), 1.0)
    limited = raw + dn_lim[..., None] * normals + dt * scale[..., None]
    return np.where(bind[..., None], limited, smoothed)


# ---------------------------------------------------------------------------
# STL surface normal query
# ---------------------------------------------------------------------------


def _stl_normals(
    trimesh_obj: object,
    points: np.ndarray,
    interpolate: bool = True,
) -> np.ndarray:
    """Query surface normals at arbitrary points from a trimesh object.

    Uses closest-point projection.  When ``interpolate=True`` (default), returns
    the barycentric blend of the closest triangle's three vertex normals at the
    closest point — C0-continuous across triangle edges, which suppresses
    facet-aligned banding on coarse STLs.  When ``interpolate=False`` returns
    the flat face normal of the closest triangle (legacy behaviour).

    Vertex normals are area-weighted by trimesh and require a clean mesh
    (consistent face orientation, no degenerate faces, merged vertices).
    ``shore.surface.io.load_stl`` enforces all three by default.

    Parameters
    ----------
    trimesh_obj :
        A ``trimesh.Trimesh`` instance.
    points :
        Shape (N, 3).
    interpolate :
        When *True*, barycentric blend of vertex normals over the closest
        triangle.  When *False*, flat face normal of the closest triangle.

    Returns
    -------
    numpy.ndarray, shape (N, 3), unit normals.
    """
    import trimesh.proximity as prox
    import trimesh.triangles as tri_utils

    closest, _, tri_ids = prox.closest_point(trimesh_obj, points)

    if interpolate:
        tris = trimesh_obj.triangles[tri_ids]
        bary = tri_utils.points_to_barycentric(tris, closest)
        vn = trimesh_obj.vertex_normals[trimesh_obj.faces[tri_ids]]
        normals = (bary[..., None] * vn).sum(axis=1)
    else:
        normals = trimesh_obj.face_normals[tri_ids]

    mag = np.linalg.norm(normals, axis=-1, keepdims=True)
    mag = np.where(mag < 1e-15, 1.0, mag)
    return normals / mag
