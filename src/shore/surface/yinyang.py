"""Yin-Yang double O-grid: two rotated standard O-grids that together
cover a full sphere without placing any grid node at a singular pole.

Mathematical foundation
-----------------------
Yin piece: standard lat-lon directions with colatitude theta in (theta_m, pi - theta_m)
    d_yin(theta, phi) = (sin theta cos phi, sin theta sin phi, cos theta)

Yang piece: same (ni, nj) grid, directions rotated by R = R_y(+pi/2)
    R = [[ 0,  0,  1],
         [ 0,  1,  0],
         [-1,  0,  0]]
    d_yang(theta, phi) = R · d_yin(theta, phi)

Yang's "poles" land at (+-1, 0, 0) — on the original equator, not at the original poles.

Full-coverage constraint: a point escapes both blocks iff |z| >= cos(theta_m) AND
|x| >= cos(theta_m). On the unit sphere x^2 + z^2 <= 1, so both constraints hold
simultaneously only if 2*cos^2(theta_m) <= 1, i.e. theta_m >= pi/4.
Full coverage is guaranteed if and only if theta_m < pi/4.

Reference: Kageyama & Sato (2004), Geochem. Geophys. Geosyst. 5(9).
"""

from __future__ import annotations

import numpy as np
import trimesh

from shore.exceptions import MeshInputError
from shore.surface.remesh import sphere_latlon_directions


YANG_ROTATION: np.ndarray = np.array(
    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]], dtype=np.float64
)

DEFAULT_THETA_MARGIN: float = np.pi / 6


def yang_latlon_directions(
    ni: int,
    nj: int,
    theta_margin: float = DEFAULT_THETA_MARGIN,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
) -> np.ndarray:
    """Build (ni, nj) unit direction vectors for the Yang O-grid block.

    The Yang block is a standard lat-lon grid rotated by 90 degrees around
    the y-axis. Its "poles" land on the original equator at (+-1, 0, 0).

    Parameters
    ----------
    ni :
        Number of points in the latitude direction (pole-to-pole, exclusive).
    nj :
        Number of points in the longitude direction (full circle, periodic).
    theta_margin :
        Polar-cap exclusion angle in radians. Must be less than pi/4 for
        full-coverage guarantee (default pi/6).

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
        Unit direction vectors for the Yang block.
    """
    d = sphere_latlon_directions(
        ni, nj, theta_margin,
        i_spacing=i_spacing, i_beta=i_beta,
        j_spacing=j_spacing, j_beta=j_beta, j_cluster_lon=j_cluster_lon,
    )  # (ni, nj, 3)
    flat = d.reshape(-1, 3)  # (ni*nj, 3)
    rotated = (YANG_ROTATION @ flat.T).T  # (ni*nj, 3)
    return rotated.reshape(ni, nj, 3)


def project_yang_to_surface(
    mesh: trimesh.Trimesh,
    ni: int,
    nj: int,
    theta_margin: float = DEFAULT_THETA_MARGIN,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
) -> np.ndarray:
    """Project a Yang O-grid sphere grid onto a triangulated surface.

    Identical ray-casting procedure as project_sphere_to_surface, but using
    yang_latlon_directions in place of sphere_latlon_directions.

    Parameters
    ----------
    mesh :
        The body surface (raw STL).
    ni :
        Latitude resolution of the structured grid.
    nj :
        Longitude resolution of the structured grid.
    theta_margin :
        Polar-cap exclusion angle in radians (default pi/6).

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
        Coordinates of the k=0 Yang surface layer.

    Raises
    ------
    ValueError
        If any ray misses the body (centroid is outside or body is not
        star-shaped).
    """
    centroid = mesh.centroid  # (3,)
    directions = yang_latlon_directions(
        ni, nj, theta_margin,
        i_spacing=i_spacing, i_beta=i_beta,
        j_spacing=j_spacing, j_beta=j_beta, j_cluster_lon=j_cluster_lon,
    )  # (ni, nj, 3)

    flat_dirs = directions.reshape(-1, 3)  # (ni*nj, 3)
    origins = np.tile(centroid, (flat_dirs.shape[0], 1))  # (ni*nj, 3)

    # trimesh ray casting: returns hit locations, ray indices, face indices
    intersector = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    locations, ray_indices, _ = intersector.intersects_location(
        ray_origins=origins,
        ray_directions=flat_dirs,
        multiple_hits=False,
    )

    # Check every ray hit
    hit_rays = set(ray_indices.tolist())
    n_rays = flat_dirs.shape[0]
    if len(hit_rays) != n_rays:
        missing = n_rays - len(hit_rays)
        raise MeshInputError(
            f"{missing} of {n_rays} rays missed the body.  "
            "The body may not be star-shaped w.r.t. its centroid, or the STL "
            "has holes.  Try repairing the mesh or shifting the centroid."
        )

    # Sort by ray index to restore (ni, nj) order
    order = np.argsort(ray_indices)
    surface_pts = locations[order].reshape(ni, nj, 3)

    return surface_pts


def coverage_mask_yin(
    directions: np.ndarray,
    theta_margin: float,
) -> np.ndarray:
    """Compute the Yin coverage mask for a set of direction vectors.

    A direction is covered by Yin if |z| < cos(theta_margin), i.e. the point
    lies within the polar-cap exclusion zone of the Yin block.

    Parameters
    ----------
    directions :
        Unit direction vectors, shape (N, 3).
    theta_margin :
        Polar-cap exclusion angle in radians.

    Returns
    -------
    numpy.ndarray, shape (N,)
        Boolean mask where True indicates the direction is covered by Yin.
    """
    return np.abs(directions[:, 2]) < np.cos(theta_margin)


def coverage_mask_yang(
    directions: np.ndarray,
    theta_margin: float,
) -> np.ndarray:
    """Compute the Yang coverage mask for a set of direction vectors.

    A direction is covered by Yang if |x| < cos(theta_margin), i.e. the point
    lies within the polar-cap exclusion zone of the Yang block (which has its
    pseudo-poles on the x-axis).

    Parameters
    ----------
    directions :
        Unit direction vectors, shape (N, 3).
    theta_margin :
        Polar-cap exclusion angle in radians.

    Returns
    -------
    numpy.ndarray, shape (N,)
        Boolean mask where True indicates the direction is covered by Yang.
    """
    return np.abs(directions[:, 0]) < np.cos(theta_margin)