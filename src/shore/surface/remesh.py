"""Unstructured STL surface → structured (ni, nj) quad grid via sphere projection.

Strategy
--------
1. Place a unit sphere centered at the body centroid.
2. Build a structured (ni, nj) lat-lon grid on that sphere.
3. For each sphere point, cast a ray from the centroid outward and find
   the first intersection with the STL body.
4. The intersection points form the k=0 layer of the O-grid.

For the 6-block caps topology a unified projection covers the full sphere in
one ray-cast (``project_caps_to_surface``).  The polar caps and equatorial band
share the same ``di`` by construction — no Laplacian seed needed.

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


def project_caps_to_surface(
    mesh: trimesh.Trimesh,
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
    pole_eps_deg: float = 0.5,
) -> tuple[np.ndarray, int]:
    """Project a unified lat-lon grid (caps + equatorial) onto a surface.

    Unlike ``project_sphere_to_surface``, this function covers the full sphere
    in a single ray-cast.  The equatorial band and both polar caps share the
    same angular step ``di`` by construction — the seam between equatorial and
    cap blocks is spacing-continuous.

    The grid is parameterised over ``[pole_eps, pi - pole_eps]`` with
    ``ni_total`` rows, where ``pole_eps`` is a small exclusion near the
    singular poles.  The split index ``nc_cap`` is chosen so that the cap
    boundary row lands as close as possible to the requested ``theta_cap``.

    Parameters
    ----------
    mesh :
        Triangulated body surface.
    ni :
        Number of rows in the equatorial band (unchanged from the single-block
        interface).
    nj :
        Number of longitude points (multiple of 4 for caps topology).
    theta_cap :
        Requested polar exclusion angle in radians.  The actual split is
        snapped to the nearest integer multiple of ``di_eq`` within 5°.
        Default ``0.05 * pi`` (~9°).
    theta_cap_deg :
        Convenience alias in degrees.  Mutually exclusive with ``theta_cap``.
    i_spacing :
        Latitude spacing law applied uniformly across the full colatitude range.
        ``'uniform'`` (default), ``'tanh'``, or ``'tanh2'``.
    i_beta :
        Clustering strength for tanh/tanh2 i-spacing.
    j_spacing :
        Longitude spacing law.  ``'uniform'`` or ``'tanh2'``.
    j_beta :
        Clustering strength for tanh2 j-spacing.
    j_cluster_lon :
        Target longitude for tanh2 j-clustering (radians).
    pole_eps_deg :
        Angular exclusion at each pole in degrees (default 0.5°).

    Returns
    -------
    surface : numpy.ndarray, shape (ni_total, nj, 3)
        Full projected surface.  Rows ``0 .. nc_cap-1`` are the south cap
        (colatitude near 0), rows ``nc_cap .. ni_total-nc_cap`` are the
        equatorial band, rows ``ni_total-nc_cap .. ni_total-1`` are the north
        cap (colatitude near pi).
    nc_cap : int
        Number of rows in each polar cap block (including the shared boundary
        row with the equatorial band).

    Raises
    ------
    ParameterError
        If theta_cap and theta_cap_deg are both supplied.
    MeshInputError
        If any ray misses the body.
    """
    if theta_cap is not None and theta_cap_deg is not None:
        raise ParameterError(
            "Supply at most one of theta_cap and theta_cap_deg."
        )
    if theta_cap_deg is not None:
        theta_cap = np.deg2rad(theta_cap_deg)
    if theta_cap is None:
        theta_cap = 0.05 * np.pi

    pole_eps = np.deg2rad(pole_eps_deg)

    # Target equatorial step (used only to estimate nc_cap).
    di_eq = (np.pi - 2.0 * theta_cap) / max(ni - 1, 1)

    # Number of cap rows: how many steps fit between pole_eps and theta_cap.
    nc_cap = max(2, round((theta_cap - pole_eps) / di_eq) + 1)

    # Total rows: equatorial ni + 2*(nc_cap-1) cap rows (boundary shared).
    ni_total = ni + 2 * (nc_cap - 1)

    # Build a uniform latitude array over the full [pole_eps, pi-pole_eps] range.
    # All ni_total rows share the same step; the equatorial band is the middle ni
    # rows.  The actual step di_actual may differ slightly from di_eq due to the
    # integer rounding of nc_cap — that is acceptable and unavoidable.
    theta_start = pole_eps
    theta_stop  = np.pi - pole_eps
    lats = build_i_coords(
        ni_total, theta_start, theta_stop,
        spacing=i_spacing, beta=i_beta,
    )

    lons = build_j_coords(nj, spacing=j_spacing, beta=j_beta, cluster_lon=j_cluster_lon)

    lat_grid, lon_grid = np.meshgrid(lats, lons, indexing="ij")
    x = np.sin(lat_grid) * np.cos(lon_grid)
    y = np.sin(lat_grid) * np.sin(lon_grid)
    z = np.cos(lat_grid)
    directions = np.stack([x, y, z], axis=-1)  # (ni_total, nj, 3)

    centroid  = mesh.centroid
    flat_dirs = directions.reshape(-1, 3)
    origins   = np.tile(centroid, (flat_dirs.shape[0], 1))

    intersector = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
    locations, ray_indices, _ = intersector.intersects_location(
        ray_origins=origins,
        ray_directions=flat_dirs,
        multiple_hits=False,
    )

    hit_rays = set(ray_indices.tolist())
    n_rays   = flat_dirs.shape[0]
    if len(hit_rays) != n_rays:
        missing = n_rays - len(hit_rays)
        raise MeshInputError(
            f"{missing} of {n_rays} rays missed the body.  "
            "The body may not be star-shaped w.r.t. its centroid, or the STL "
            "has holes.  Try repairing the mesh or shifting the centroid."
        )

    order       = np.argsort(ray_indices)
    surface_pts = locations[order].reshape(ni_total, nj, 3)

    return surface_pts, nc_cap


def project_gnomonic_cap_to_surface(
    mesh: trimesh.Trimesh,
    eq_rows: list[np.ndarray],
    nc: int,
    pole: str = "north",
) -> np.ndarray:
    """Ray-cast a (nc, nc) gnomonic cap patch onto the STL surface.

    Builds ``nc*nc`` ray directions by transfinite interpolation (TFI) of the
    4 boundary edges (the equatorial sub-block pole rows), then casts all rays
    from the body centroid.  Boundary nodes are copied from the input rows to
    guarantee exact conformity with adjacent equatorial sub-blocks.

    The boundary layout follows the same convention as
    ``_overwrite_cap_boundaries``:

    North pole::

        cap[0,  :] ← eq_rows[0]            (sub0 pole row)
        cap[:, -1] ← eq_rows[1]            (sub1 pole row)
        cap[-1, :] ← eq_rows[2][::-1]      (sub2 pole row reversed)
        cap[:,  0] ← eq_rows[3][::-1]      (sub3 pole row reversed)

    South pole::

        cap[:,  0] ← eq_rows[0]            (sub0 pole row)
        cap[-1, :] ← eq_rows[1]            (sub1 pole row)
        cap[:, -1] ← eq_rows[2][::-1]      (sub2 pole row reversed)
        cap[0,  :] ← eq_rows[3][::-1]      (sub3 pole row reversed)

    Parameters
    ----------
    mesh :
        Triangulated body surface.
    eq_rows :
        ``[sub0_row, sub1_row, sub2_row, sub3_row]``, each shape ``(nc, 3)``.
        Boundary edge positions on the STL surface.
    nc :
        Edge size of the cap block (= nj//4 + 1).
    pole :
        ``'north'`` or ``'south'``.

    Returns
    -------
    numpy.ndarray, shape (nc, nc, 3)

    Raises
    ------
    MeshInputError
        If any interior ray misses the body.
    """
    # Build the 4 boundary edge arrays in cap (i, j) order:
    # edge_i0[j]  = cap[0,  j]
    # edge_i1[j]  = cap[-1, j]
    # edge_j0[i]  = cap[i,  0]
    # edge_j1[i]  = cap[i, -1]
    if pole == "north":
        edge_i0 = eq_rows[0]            # cap[0, :]
        edge_j1 = eq_rows[1]            # cap[:, -1]
        edge_i1 = eq_rows[2][::-1]      # cap[-1, :] ← sub2 reversed
        edge_j0 = eq_rows[3][::-1]      # cap[:, 0]  ← sub3 reversed
    else:
        edge_j0 = eq_rows[0]            # cap[:, 0]
        edge_i1 = eq_rows[1]            # cap[-1, :]
        edge_j1 = eq_rows[2][::-1]      # cap[:, -1] ← sub2 reversed
        edge_i0 = eq_rows[3][::-1]      # cap[0, :]  ← sub3 reversed

    u = np.linspace(0.0, 1.0, nc)   # u = i/(nc-1)
    v = np.linspace(0.0, 1.0, nc)   # v = j/(nc-1)
    ug, vg = np.meshgrid(u, v, indexing="ij")  # (nc, nc)

    # Transfinite interpolation (Gordon-Hall):
    # P(u,v) = (1-u)*e_i0(v) + u*e_i1(v) + (1-v)*e_j0(u) + v*e_j1(u)
    #        - corner blend (bilinear)
    pts = (
        (1.0 - ug[..., None]) * edge_i0[np.newaxis, :, :]          # (1,nc,3) broadcast
        +       ug[..., None]  * edge_i1[np.newaxis, :, :]
        + (1.0 - vg[..., None]) * edge_j0[:, np.newaxis, :]         # (nc,1,3) broadcast
        +       vg[..., None]  * edge_j1[:, np.newaxis, :]
        - (1.0 - ug[..., None]) * (1.0 - vg[..., None]) * edge_i0[0]
        - (1.0 - ug[..., None]) *        vg[..., None]  * edge_i0[-1]
        -        ug[..., None]  * (1.0 - vg[..., None]) * edge_i1[0]
        -        ug[..., None]  *        vg[..., None]  * edge_i1[-1]
    )  # (nc, nc, 3)

    # Normalise to unit directions for ray casting
    mag = np.linalg.norm(pts, axis=-1, keepdims=True)
    mag = np.where(mag < 1e-15, 1.0, mag)
    dirs = pts / mag  # (nc, nc, 3)

    # Only ray-cast interior nodes; boundary nodes use input positions exactly
    # Build a mask for interior points (not on any boundary edge)
    interior_mask = np.ones((nc, nc), dtype=bool)
    interior_mask[0,  :] = False
    interior_mask[-1, :] = False
    interior_mask[:,  0] = False
    interior_mask[:, -1] = False

    interior_dirs = dirs[interior_mask]   # (n_interior, 3)
    n_int = interior_dirs.shape[0]

    cap = np.empty((nc, nc, 3), dtype=np.float64)

    if n_int > 0:
        centroid  = mesh.centroid
        origins   = np.tile(centroid, (n_int, 1))

        intersector = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
        locations, ray_indices, _ = intersector.intersects_location(
            ray_origins=origins,
            ray_directions=interior_dirs,
            multiple_hits=False,
        )

        hit_rays = set(ray_indices.tolist())
        if len(hit_rays) != n_int:
            missing = n_int - len(hit_rays)
            raise MeshInputError(
                f"{missing} of {n_int} gnomonic cap interior rays missed the body. "
                "The STL may have holes near the pole or the centroid is off-centre."
            )

        order = np.argsort(ray_indices)
        cap[interior_mask] = locations[order]

    # Fill boundary edges from input rows (exact conformity, no ray-cast error)
    cap[0,  :, :] = edge_i0
    cap[-1, :, :] = edge_i1
    cap[:,  0, :] = edge_j0
    cap[:, -1, :] = edge_j1

    return cap
