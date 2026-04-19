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
) -> np.ndarray:
    """Estimate outward unit normals for every node of a structured layer.

    Central differences in both directions; one-sided at non-periodic boundaries.

    Parameters
    ----------
    layer :
        Node coordinates, shape (ni, nj, 3).
    periodic :
        ``(periodic_i, periodic_j)``.  When ``True`` the corresponding
        direction uses modular (wrap-around) indexing for the end nodes;
        when ``False`` one-sided differences are used at the boundaries.
        Default ``(False, True)`` preserves the original O-grid behaviour
        (i open, j periodic).

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
    """
    ni, nj = layer.shape[0], layer.shape[1]
    periodic_i, periodic_j = periodic

    # --- i-direction tangent ---
    ti = np.empty_like(layer)
    if periodic_i:
        ip1 = (np.arange(ni) + 1) % ni
        im1 = (np.arange(ni) - 1) % ni
        ti  = layer[ip1] - layer[im1]
    else:
        ti[1:-1] = layer[2:] - layer[:-2]
        ti[0]    = layer[1]  - layer[0]
        ti[-1]   = layer[-1] - layer[-2]

    # --- j-direction tangent ---
    if periodic_j:
        jp1 = (np.arange(nj) + 1) % nj
        jm1 = (np.arange(nj) - 1) % nj
        tj  = layer[:, jp1] - layer[:, jm1]
    else:
        tj = np.empty_like(layer)
        tj[:, 1:-1] = layer[:, 2:] - layer[:, :-2]
        tj[:, 0]    = layer[:, 1]  - layer[:, 0]
        tj[:, -1]   = layer[:, -1] - layer[:, -2]

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
        j_fwd[:, -1]  = layer[:, -1]   # clamp: boundary contributes itself
        j_bwd[:, 1:]  = layer[:, :-1]
        j_bwd[:, 0]   = layer[:, 0]

    # i neighbours
    if periodic_i:
        i_fwd = np.roll(layer, -1, axis=0)
        i_bwd = np.roll(layer,  1, axis=0)
    else:
        i_fwd = np.empty_like(layer)
        i_bwd = np.empty_like(layer)
        i_fwd[:-1] = layer[1:]
        i_fwd[-1]  = layer[-1]
        i_bwd[1:]  = layer[:-1]
        i_bwd[0]   = layer[0]

    laplacian = (j_fwd + j_bwd + i_fwd + i_bwd) / 4.0
    smoothed  = layer + strength * (laplacian - layer)

    # always clamp the wall row
    smoothed[0] = layer[0]
    # clamp non-periodic boundaries
    if not periodic_i:
        smoothed[-1] = layer[-1]
    if not periodic_j:
        smoothed[:, 0]  = layer[:, 0]
        smoothed[:, -1] = layer[:, -1]

    return smoothed


# ---------------------------------------------------------------------------
# STL surface normal query
# ---------------------------------------------------------------------------

def _stl_normals(trimesh_obj: object, points: np.ndarray) -> np.ndarray:
    """Query surface normals at arbitrary points from a trimesh object.

    Uses closest-point projection; returns the face normal of the nearest
    triangle (flat, not interpolated — sufficient for blending guidance).

    Parameters
    ----------
    trimesh_obj :
        A ``trimesh.Trimesh`` instance.
    points :
        Shape (N, 3).

    Returns
    -------
    numpy.ndarray, shape (N, 3), unit normals.
    """
    import trimesh.proximity as prox
    _, _, tri_ids = prox.closest_point(trimesh_obj, points)
    normals = trimesh_obj.face_normals[tri_ids]
    mag = np.linalg.norm(normals, axis=-1, keepdims=True)
    mag = np.where(mag < 1e-15, 1.0, mag)
    return normals / mag


