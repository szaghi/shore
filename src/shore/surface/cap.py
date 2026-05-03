"""Polar cap k=0 seed face via gnomonic projection of a flat tangent square.

The cap is constructed by placing a flat planar square (a diamond, with its 4
corners on the tangent-plane cardinal axes at distance ``tan(theta_cap)``) at
the cap pole and then projecting every node radially onto the body surface
through the cached ``Projector``.  The 4 corners land at colatitude
``theta_cap`` exactly (on a sphere); cap edges are great-circle arcs between
adjacent corners; the interior is filled by bilinear interpolation in the
tangent plane (which produces equispaced corner-to-corner straight lines and,
after radial projection, equispaced great-circle arcs on the sphere).

Why a diamond and not an axis-aligned square: the 4 cap corners must coincide
with the 4 sub-block-junction longitudes ``phi ∈ {0, π/2, π, 3π/2}`` to
satisfy the existing ``_overwrite_cap_boundaries`` contract.  A diamond with
corners on the cardinal axes places the cap-corners exactly at those
longitudes.  An axis-aligned square would put corners at ``phi ∈ {±π/4,
±3π/4}`` — a 45° offset from the sub-block junctions.

Spacing inheritance
-------------------
The tangent-plane parametric grid uses ``build_j_coords``-style 1D laws on
``[-tan(theta_cap), +tan(theta_cap)]``, so the cap edges' node distribution
matches the equator's ``j_spacing`` law along the seam.  Cap edges are then
the equator pole rows (bit-exact) and the equator inherits this distribution
on its first/last latitude row.

Pole convention
---------------
- ``pole='south'`` (shore: ``cap_south``, sphere top, z > 0): no rotation.
- ``pole='north'`` (shore: ``cap_north``, sphere bottom, z < 0): cap-local
  ``+z'`` axis maps to world ``-z``; the (u, v) tangent-plane coords are
  reflected so the cap's index orientation matches the existing override
  convention (``cap[0, :]`` → sub0 pole row, etc.).
"""

from __future__ import annotations

import numpy as np

from shore.volume.spacing import tanh_two_sided


def _tangent_coords(
    nc: int,
    half_side: float,
    spacing: str,
    beta: float,
) -> np.ndarray:
    """Return (nc,) coordinates in [-half_side, +half_side] under the given law.

    ``spacing='uniform'`` is a plain linspace.  ``spacing='tanh2'`` uses
    ``tanh_two_sided`` to cluster nodes near the corners (and thus near the
    sub-block junctions on the seam).
    """
    if spacing == "uniform":
        return np.linspace(-half_side, half_side, nc)
    if spacing == "tanh2":
        return tanh_two_sided(nc, -half_side, half_side, beta)
    # 'tanh' is not a meaningful periodic-axis law and is rejected upstream
    # for j_spacing; treat any other value as uniform fallback.
    return np.linspace(-half_side, half_side, nc)


def build_cap_k0(
    theta_cap: float,
    nc: int,
    pole: str,
    projector: object,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    phi_offset: float = 0.0,
) -> np.ndarray:
    """Build the ``(nc, nc, 3)`` cap k=0 seed face via gnomonic projection.

    Parameters
    ----------
    theta_cap :
        Polar exclusion angle in radians.  The 4 cap corners sit at colatitude
        ``theta_cap`` on a sphere; cap-edge midpoints sit at colatitude
        ``atan(tan(theta_cap)/sqrt(2))`` (smaller than ``theta_cap`` — edges
        dip toward the pole, which is the natural shape of a great-circle arc
        between two parallel-circle corners).
    nc :
        Cap resolution.  Must equal ``nj // 4 + 1`` so the cap edges have the
        same node count as the sub-block pole rows.
    pole :
        ``'north'`` (cap_north, ``z<0`` of the body) or ``'south'`` (cap_south,
        ``z>0``).
    projector :
        ``Projector`` with cached BVH.
    j_spacing :
        ``'uniform'`` (default) or ``'tanh2'``.  The tangent-plane node
        distribution along the cap edges; inherited by the equator seam.
    j_beta :
        Clustering strength for ``'tanh2'`` (ignored for ``'uniform'``).
    phi_offset :
        Azimuthal rotation (radians) applied to the gnomonic diamond about
        the polar axis.  Default ``0.0`` places cap corners at world azimuths
        ``{0, pi/2, pi, 3*pi/2}`` (the CubedSphereTopology convention).
        :class:`shore.mesh.OHTopology` passes ``pi/4`` so the four sub-block
        *centres* (rather than corners) face the world cardinal axes — which
        is what makes a 1:1 sub-block-to-box-face mapping work for the
        wrap-channel hybrid.  The cap array's parametric (a, b) layout is
        unchanged; only the world-space orientation rotates.

    Returns
    -------
    numpy.ndarray, shape ``(nc, nc, 3)``
        Cap k=0 layer with nodes lying on the body surface.
    """
    if pole not in ("north", "south"):
        raise ValueError(f"pole must be 'north' or 'south', got {pole!r}")
    if not (0.0 < theta_cap < np.pi / 2):
        raise ValueError(f"theta_cap must be in (0, pi/2), got {theta_cap}")

    half_side = float(np.tan(theta_cap))

    # 1D tangent-plane coordinates in (u, v) ∈ [-half_side, half_side].
    u_lin = _tangent_coords(nc, half_side, j_spacing, j_beta)
    v_lin = _tangent_coords(nc, half_side, j_spacing, j_beta)

    # The cap is a diamond in the tangent plane (z = 1):
    #   C0 = (+h, 0, 1)   at cap[0, 0]   (azimuth   0 deg)
    #   C1 = (0, +h, 1)   at cap[-1, 0]  (azimuth  90 deg)
    #   C2 = (-h, 0, 1)   at cap[-1, -1] (azimuth 180 deg)
    #   C3 = (0, -h, 1)   at cap[0, -1]  (azimuth 270 deg)
    # with h = half_side = tan(theta_cap).  Parametrise by (a, b) in [0, 1]^2
    # where a = (u_lin + h) / (2h), b similar.  Bilinear interpolation gives:
    #   x_local = h * (1 - a - b)
    #   y_local = h * (a - b)
    a = (u_lin + half_side) / (2.0 * half_side)  # (nc,)
    b = (v_lin + half_side) / (2.0 * half_side)  # (nc,)
    a_grid, b_grid = np.meshgrid(a, b, indexing="ij")  # (nc, nc) each

    x_local = half_side * (1.0 - a_grid - b_grid)
    y_local = half_side * (a_grid - b_grid)
    z_local = np.ones_like(x_local)

    # World rotation: south pole keeps z = +1; north pole flips z = -1.  For
    # the north pole we additionally negate y so the cap's index winding
    # matches ``_overwrite_cap_boundaries`` (cap[0, :] traverses sub0's pole
    # row in the same direction the override expects).
    if pole == "south":
        x_world = x_local
        y_world = y_local
        z_world = z_local
    else:  # north
        x_world = x_local
        y_world = -y_local
        z_world = -z_local

    # Stack into (nc, nc, 3) directions and normalise (gnomonic projection
    # onto the unit sphere of directions).
    dirs = np.stack([x_world, y_world, z_world], axis=-1)
    dirs = dirs / np.linalg.norm(dirs, axis=-1, keepdims=True)

    # Optional phi_offset: rotate the (x, y) tangent-plane components about
    # the polar axis.  z stays fixed.  OHTopology uses phi_offset=pi/4 so
    # sub-block centres face cardinal axes.
    if phi_offset != 0.0:
        c, s = float(np.cos(phi_offset)), float(np.sin(phi_offset))
        x_rot = c * dirs[..., 0] - s * dirs[..., 1]
        y_rot = s * dirs[..., 0] + c * dirs[..., 1]
        dirs = np.stack([x_rot, y_rot, dirs[..., 2]], axis=-1)

    # Ray-cast through the cached Projector.
    result = projector.project(dirs.reshape(-1, 3), on_failure="raise")  # type: ignore[attr-defined]
    return result.locations.reshape(nc, nc, 3)
