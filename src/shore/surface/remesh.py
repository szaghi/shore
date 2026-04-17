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

from shore.exceptions import MeshInputError, ParameterError
from shore.volume.spacing import build_i_coords, build_j_coords


def sphere_latlon_directions(
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    theta_margin: float | None = None,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
) -> np.ndarray:
    """Build (ni, nj) unit direction vectors on a lat-lon sphere.

    The poles are excluded so the grid is open (no singular points).
    The grid closes in the j-direction (periodic), not in i.

    Parameters
    ----------
    ni :
        Number of points in the latitude direction (pole-to-pole, exclusive).
    nj :
        Number of points in the longitude direction (full circle, periodic).
    theta_cap :
        Polar exclusion angle in radians: the equatorial band covers
        ``[theta_cap, pi - theta_cap]``.  The cap boundary parallel has
        radius ``sin(theta_cap)`` on the unit sphere.  Must be in
        ``(0, pi/2)``.  Default ``0.05*pi`` (~9°).
    theta_cap_deg :
        Convenience alias for ``theta_cap`` in degrees.  Mutually exclusive
        with ``theta_cap`` and ``theta_margin``.
    theta_margin :
        Deprecated alias for ``theta_cap`` (radians).

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
        Unit direction vectors.

    Raises
    ------
    ValueError
        If more than one of ``theta_cap``, ``theta_cap_deg``, ``theta_margin``
        are supplied, or if the resolved angle is outside ``(0, pi/2)``.
    """
    # --- resolve theta_cap ---
    n_supplied = sum(x is not None for x in (theta_cap, theta_cap_deg, theta_margin))
    if n_supplied > 1:
        raise ParameterError(
            "Supply at most one of theta_cap, theta_cap_deg, theta_margin."
        )
    if theta_cap_deg is not None:
        theta_cap = np.deg2rad(theta_cap_deg)
    elif theta_margin is not None:
        theta_cap = theta_margin

    if theta_cap is None:
        # Default: 0.05*pi (~9°)
        theta_cap = 0.05 * np.pi

    if not (0.0 < theta_cap < np.pi / 2):
        raise ParameterError(
            f"theta_cap={theta_cap:.6f} must be in (0, pi/2). "
            "0 would collapse the cap to a point; pi/2 would put the "
            "cap edge at the equator."
        )

    lats = build_i_coords(ni, theta_cap, np.pi - theta_cap, spacing=i_spacing, beta=i_beta)

    # longitude: [0, 2*pi) — nj points, periodic
    lons = build_j_coords(nj, spacing=j_spacing, beta=j_beta, cluster_lon=j_cluster_lon)

    lat_grid, lon_grid = np.meshgrid(lats, lons, indexing="ij")

    x = np.sin(lat_grid) * np.cos(lon_grid)
    y = np.sin(lat_grid) * np.sin(lon_grid)
    z = np.cos(lat_grid)

    return np.stack([x, y, z], axis=-1)  # (ni, nj, 3)


def project_sphere_to_surface(
    mesh: trimesh.Trimesh,
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    theta_margin: float | None = None,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
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
    theta_cap :
        Polar exclusion angle in radians: the equatorial band covers
        ``[theta_cap, pi - theta_cap]``.  The cap boundary parallel has
        radius ``sin(theta_cap)`` on a unit sphere.  Must be in
        ``(0, pi/2)``.  Default ``0.05*pi`` (~9°).
    theta_cap_deg :
        Convenience alias for ``theta_cap`` in degrees.  Exactly one of
        ``theta_cap`` and ``theta_cap_deg`` may be supplied.
    theta_margin :
        Deprecated alias for ``theta_cap`` (radians).  Will be removed in a
        future release.

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
        Coordinates of the k=0 structured surface layer.

    Raises
    ------
    ValueError
        If more than one of ``theta_cap``, ``theta_cap_deg``, and
        ``theta_margin`` are supplied, if ``theta_cap`` is outside
        ``(0, pi/2)``, or if any ray misses the body.
    """
    # --- resolve theta_cap ---
    n_supplied = sum(x is not None for x in (theta_cap, theta_cap_deg, theta_margin))
    if n_supplied > 1:
        raise ParameterError(
            "Supply at most one of theta_cap, theta_cap_deg, theta_margin."
        )
    if theta_cap_deg is not None:
        theta_cap = np.deg2rad(theta_cap_deg)
    elif theta_margin is not None:
        theta_cap = theta_margin

    centroid = mesh.centroid  # (3,)
    directions = sphere_latlon_directions(
        ni, nj, theta_cap=theta_cap,
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
