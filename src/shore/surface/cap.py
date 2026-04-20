"""Polar cap k=0 seed face via Shirley concentric map + stereographic projection.

Phase 4 rewrite: replaces ``project_gnomonic_cap_to_surface`` (Cartesian TFI in
world coords) and the Laplacian-based ``_build_cap_k0`` / ``_pull_cap_interior``
helpers with a single mathematically clean construction that achieves C1 spacing
continuity at the cap/equator seam (exact at edge midpoints, bounded at corners).

Pipeline
--------
1. Compute the equator's first latitude step ``dtheta_seam = lats_eq[1] - lats_eq[0]``.
2. Build a radial profile of rings in colatitude: ``ring_theta[0] = theta_cap``,
   ``ring_theta[1] = theta_cap - dtheta_seam`` (C1 by construction),
   ``ring_theta[2:]`` from ``build_i_coords`` matching the equator's spacing law.
3. Map the logical ``(nc, nc)`` grid to rings via Shirley's concentric
   square-to-disk map, rotated so the 4 cap corners align with the 4 sub-block
   seam corners.
4. Invert the stereographic projection ``rho = 2 tan(theta/2)`` to get ray
   directions on the sphere.
5. Ray-cast through the ``Projector``; override the boundary nodes with the
   exact equator pole-row coordinates.

Shore naming convention
-----------------------
This module accepts ``pole='north'`` / ``pole='south'`` consistent with shore's
existing ``_overwrite_cap_boundaries`` convention:

- ``pole='north'`` → the block labelled ``cap_north`` at ``z < 0`` of the sphere.
- ``pole='south'`` → the block labelled ``cap_south`` at ``z > 0`` of the sphere.

The z-sign in the ray direction is selected accordingly.
"""

from __future__ import annotations

import warnings

import numpy as np

from shore.volume.spacing import build_i_coords

# ---------------------------------------------------------------------------
# Shirley concentric square-to-disk
# ---------------------------------------------------------------------------


def _shirley_rotated(
    u: np.ndarray,
    v: np.ndarray,
    phi_0: float,
    pole: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Shirley concentric map ``[-1,1]^2 -> unit disk``, returning ``(r_disk, phi_disk)``.

    The signed-r form of Shirley & Chiu (1997) is used so the quadrant is
    preserved correctly.  For ``pole='north'`` the ``v`` axis is flipped so the
    cap's boundary traversal direction matches the physical longitude
    direction of the equator (the ``_overwrite_cap_boundaries`` mapping for
    the north cap has reversed orientation relative to the south cap).

    Parameters
    ----------
    u, v :
        Grid arrays in ``[-1, 1]``.  Same shape.
    phi_0 :
        Rotation applied to ``phi_disk`` so the cap's (0, 0) corner lands at
        the physical longitude of the sub-block junction.
    pole :
        ``'north'`` or ``'south'``.

    Returns
    -------
    r_disk, phi_disk :
        Disk radius ``in [0, 1]`` and angle (radians).
    """
    if pole == "north":
        v = -v

    abs_u = np.abs(u)
    abs_v = np.abs(v)
    use_u = abs_u >= abs_v

    safe_u = np.where(u == 0.0, 1.0, u)
    safe_v = np.where(v == 0.0, 1.0, v)
    r_signed = np.where(use_u, u, v)
    phi_sgn = np.where(
        use_u,
        (np.pi / 4.0) * (v / safe_u),
        (np.pi / 2.0) - (np.pi / 4.0) * (u / safe_v),
    )
    is_origin = (abs_u == 0.0) & (abs_v == 0.0)
    r_signed = np.where(is_origin, 0.0, r_signed)
    phi_sgn = np.where(is_origin, 0.0, phi_sgn)

    x = r_signed * np.cos(phi_sgn)
    y = r_signed * np.sin(phi_sgn)

    r_disk = np.hypot(x, y)
    phi_disk = np.arctan2(y, x) + phi_0
    return r_disk, phi_disk


# ---------------------------------------------------------------------------
# Boundary override (same mapping as sphere_cap._overwrite_cap_boundaries)
# ---------------------------------------------------------------------------


def _override_boundary(
    cap: np.ndarray,
    eq_pole_rows: list[np.ndarray],
    pole: str,
) -> np.ndarray:
    """Replace the 4 boundary edges of the cap with the equator pole rows."""
    if pole == "north":
        cap[0, :, :] = eq_pole_rows[0]
        cap[:, -1, :] = eq_pole_rows[1]
        cap[-1, :, :] = eq_pole_rows[2][::-1]
        cap[:, 0, :] = eq_pole_rows[3][::-1]
    else:
        cap[:, 0, :] = eq_pole_rows[0]
        cap[-1, :, :] = eq_pole_rows[1]
        cap[:, -1, :] = eq_pole_rows[2][::-1]
        cap[0, :, :] = eq_pole_rows[3][::-1]
    return cap


# ---------------------------------------------------------------------------
# Corner longitude
# ---------------------------------------------------------------------------


def _longitude_of_cap_00(eq_pole_rows: list[np.ndarray]) -> float:
    """Physical longitude of the shared corner at cap ``(0, 0)``.

    For both poles, ``cap[0, 0]`` is always ``eq_pole_rows[0][0]`` (the first
    node of the first sub-block's pole row).
    """
    p = eq_pole_rows[0][0]
    return float(np.arctan2(p[1], p[0]))


# ---------------------------------------------------------------------------
# Ring-theta profile
# ---------------------------------------------------------------------------


def _ring_theta_profile(
    n_rings: int,
    theta_cap: float,
    dtheta_seam: float,
    pole_eps_rad: float,
    i_spacing: str,
    i_beta: float,
) -> np.ndarray:
    """Return (n_rings,) colatitudes for the cap rings, from seam inward.

    ``ring_theta[0] = theta_cap`` (boundary) and ``ring_theta[-1]`` sits near
    ``pole_eps_rad``.  The seam step is exact (C1 at edge midpoints) when
    ``theta_cap - dtheta_seam > pole_eps_rad``; otherwise a uniform fallback
    is used and a warning is emitted.
    """
    ring_theta = np.empty(n_rings)
    ring_theta[0] = theta_cap

    if n_rings == 1:
        return ring_theta

    if theta_cap - dtheta_seam > pole_eps_rad:
        ring_theta[1] = theta_cap - dtheta_seam
        if n_rings > 2:
            # Distribute the remaining rings smoothly from ring_theta[1] down to
            # pole_eps_rad using the inherited spacing law, then drop the first
            # point (which duplicates ring_theta[1]).  This preserves C1 at the
            # seam and gives monotone interior rings even when
            # theta_cap - 2 dtheta_seam <= pole_eps.
            tail = build_i_coords(
                n_rings - 1,
                theta_cap - dtheta_seam,
                pole_eps_rad,
                spacing=i_spacing,
                beta=i_beta,
            )
            ring_theta[2:] = tail[1:]
        return ring_theta

    warnings.warn(
        f"Cap C1 seam match not achievable: dtheta_eq_seam={dtheta_seam:.4f} rad "
        f"exceeds theta_cap - pole_eps={theta_cap - pole_eps_rad:.4f} rad. "
        "Falling back to uniform cap spacing (spacing continuity only approximate).",
        stacklevel=3,
    )
    return build_i_coords(n_rings, theta_cap, pole_eps_rad, spacing="uniform")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def build_cap_k0(
    eq_pole_rows: list[np.ndarray],
    theta_cap: float,
    ni: int,
    pole: str,
    projector: object | None = None,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
    pole_eps_rad: float | None = None,
) -> np.ndarray:
    """Build the ``(nc, nc, 3)`` cap k=0 seed face via Shirley + stereographic.

    Parameters
    ----------
    eq_pole_rows :
        ``[sub0_row, sub1_row, sub2_row, sub3_row]`` — each shape ``(nc, 3)``.
        Equatorial sub-block pole rows at k=0.
    theta_cap :
        Polar exclusion angle in radians used by the equator projection.
    ni :
        Number of rows in the equatorial band (used only to compute the
        equator's seam step ``lats_eq[1] - lats_eq[0]``).
    pole :
        ``'north'`` (shore cap_north, ``z < 0``) or ``'south'`` (``z > 0``).
    projector :
        ``Projector`` with cached BVH.  When *None*, the cap is constructed
        analytically on a sphere of radius ``mean(||eq_pole_rows||)`` — used
        only by unit tests that build a cap from a synthetic sphere layer.
    i_spacing :
        Latitude spacing law used by the equator (``'uniform'``, ``'tanh'``,
        ``'tanh2'``).  Inherited so the cap inherits the equator's density.
    i_beta :
        Clustering strength for tanh/tanh2 (same semantics as the equator).
    pole_eps_rad :
        Angular exclusion at the geographic pole (default 0.5 deg).

    Returns
    -------
    numpy.ndarray, shape ``(nc, nc, 3)``
    """
    if pole not in ("north", "south"):
        raise ValueError(f"pole must be 'north' or 'south', got {pole!r}")
    if pole_eps_rad is None:
        pole_eps_rad = float(np.deg2rad(0.5))

    nc = eq_pole_rows[0].shape[0]

    # --- equator's seam step ---
    lats_eq = build_i_coords(
        ni,
        theta_cap,
        np.pi - theta_cap,
        spacing=i_spacing,
        beta=i_beta,
    )
    dtheta_seam = float(lats_eq[1] - lats_eq[0])

    # --- ring radial profile ---
    n_rings = (nc + 1) // 2
    ring_theta = _ring_theta_profile(
        n_rings,
        theta_cap,
        dtheta_seam,
        pole_eps_rad,
        i_spacing,
        i_beta,
    )
    ring_rho = 2.0 * np.tan(ring_theta / 2.0)

    # --- Shirley mapping ---
    u_lin = np.linspace(-1.0, 1.0, nc)
    v_lin = np.linspace(-1.0, 1.0, nc)
    u, v = np.meshgrid(u_lin, v_lin, indexing="ij")

    phi_00 = _longitude_of_cap_00(eq_pole_rows)
    # For north: v-flipped Shirley puts cap(0,0) at disk angle 3*pi/4.
    # For south: no flip, disk angle at cap(0,0) = -3*pi/4 (= 5*pi/4 mod 2*pi).
    phi_0 = phi_00 - 3.0 * np.pi / 4.0 if pole == "north" else phi_00 + 3.0 * np.pi / 4.0
    r_disk, phi_disk = _shirley_rotated(u, v, phi_0, pole)

    # --- ring index → rho ---
    ring_idx = np.rint((1.0 - r_disk) * (nc - 1) / 2.0).astype(int)
    ring_idx = np.clip(ring_idx, 0, n_rings - 1)
    rho = ring_rho[ring_idx]

    # --- stereographic inverse → sphere directions ---
    theta_node = 2.0 * np.arctan(rho / 2.0)
    # Shore inversion: 'north' label = z<0, 'south' label = z>0.
    z_sign = -1.0 if pole == "north" else 1.0
    dirs = np.stack(
        [
            np.sin(theta_node) * np.cos(phi_disk),
            np.sin(theta_node) * np.sin(phi_disk),
            z_sign * np.cos(theta_node),
        ],
        axis=-1,
    )  # (nc, nc, 3)

    # --- project or analytical fallback ---
    if projector is not None:
        result = projector.project(dirs.reshape(-1, 3), on_failure="raise")  # type: ignore[attr-defined]
        cap = result.locations.reshape(nc, nc, 3)
    else:
        all_boundary = np.concatenate(list(eq_pole_rows), axis=0)
        sphere_r = float(np.mean(np.linalg.norm(all_boundary, axis=-1)))
        cap = dirs * sphere_r

    # --- bit-exact boundary override ---
    _override_boundary(cap, eq_pole_rows, pole)
    return cap
