"""Unstructured STL surface → structured (ni, nj) quad grid via sphere projection.

Strategy
--------
1. Place a unit sphere centered at the body centroid.
2. Build a structured (ni, nj) lat-lon grid on that sphere.
3. For each sphere point, cast a ray from the centroid outward and find
   the first intersection with the STL body.
4. The intersection points form the k=0 layer of the O-grid.

Limitation: the approach requires that every ray from the centroid hits the
body exactly once (i.e. the body is star-shaped w.r.t. its centroid).  This
covers spheres, ellipsoids, fuselages, and most convex bodies.  Concave or
multiply-connected surfaces need block decomposition — deferred to future work.
"""

from __future__ import annotations

import numpy as np
import trimesh


def sphere_latlon_directions(ni: int, nj: int) -> np.ndarray:
    """Build (ni, nj) unit direction vectors on a lat-lon sphere.

    The poles are excluded so the grid is open (no singular points).
    The grid closes in the j-direction (periodic), not in i.

    Parameters
    ----------
    ni :
        Number of points in the latitude direction (pole-to-pole, exclusive).
    nj :
        Number of points in the longitude direction (full circle, periodic).

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
        Unit direction vectors.
    """
    # latitude: (0, pi) exclusive — avoids poles
    lats = np.linspace(0.0, np.pi, ni + 2)[1:-1]
    # longitude: [0, 2*pi) — nj points, periodic
    lons = np.linspace(0.0, 2.0 * np.pi, nj, endpoint=False)

    lat_grid, lon_grid = np.meshgrid(lats, lons, indexing="ij")

    x = np.sin(lat_grid) * np.cos(lon_grid)
    y = np.sin(lat_grid) * np.sin(lon_grid)
    z = np.cos(lat_grid)

    return np.stack([x, y, z], axis=-1)  # (ni, nj, 3)


def project_sphere_to_surface(
    mesh: trimesh.Trimesh,
    ni: int,
    nj: int,
) -> np.ndarray:
    """Project a structured sphere grid onto a triangulated surface.

    Parameters
    ----------
    mesh :
        The body surface (raw STL).
    ni :
        Latitude resolution of the structured grid.
    nj :
        Longitude resolution of the structured grid.

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
        Coordinates of the k=0 structured surface layer.

    Raises
    ------
    ValueError
        If any ray misses the body (centroid is outside or body is not
        star-shaped).
    """
    centroid = mesh.centroid  # (3,)
    directions = sphere_latlon_directions(ni, nj)  # (ni, nj, 3)

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
        raise ValueError(
            f"{missing} of {n_rays} rays missed the body.  "
            "The body may not be star-shaped w.r.t. its centroid, or the STL "
            "has holes.  Try repairing the mesh or shifting the centroid."
        )

    # Sort by ray index to restore (ni, nj) order
    order = np.argsort(ray_indices)
    surface_pts = locations[order].reshape(ni, nj, 3)

    return surface_pts
