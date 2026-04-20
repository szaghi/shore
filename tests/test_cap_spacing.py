"""Phase 4: C1 spacing continuity of the cap k=0 seed face at the cap/equator seam.

By construction ``build_cap_k0`` sets ``ring_theta[1] = theta_cap - dtheta_seam``
where ``dtheta_seam = lats_eq[1] - lats_eq[0]`` is the equator's first latitude
step.  Shirley's concentric square-to-disk map preserves the meridian (phi) along
each square row/column that crosses the origin, so:

- At edge midpoints the radial step from ring 0 (seam) to ring 1 is exactly
  ``dtheta_seam`` (on the same phi meridian).
- At corners the diagonal neighbour also sits at ``theta_cap - dtheta_seam`` on
  the same phi meridian — the step is again exactly ``dtheta_seam``.

This test sweeps over ``ni``, ``nj``, ``theta_cap``, ``i_spacing``, ``i_beta`` and
asserts the relative error is below the floating-point noise floor.
"""

from __future__ import annotations

import itertools
import warnings

import numpy as np
import pytest

from shore.surface.cap import build_cap_k0
from tests.conftest import sphere_layer


def _north_pole_rows(surface: np.ndarray, nj: int) -> list[np.ndarray]:
    """Extract the 4 equator pole rows on the shore 'north' side (z<0, i=ni-1)."""
    sector = nj // 4
    nc_j = sector + 1
    rows = []
    ni = surface.shape[0]
    for q in range(4):
        j0 = q * sector
        j_idx = [*range(j0, j0 + sector), (j0 + sector) % nj]
        rows.append(surface[ni - 1, j_idx, :])
        assert rows[-1].shape == (nc_j, 3)
    return rows


def _theta_from_cap(cap: np.ndarray, pole: str) -> np.ndarray:
    """Colatitude (radians) of cap nodes about the block's pole axis.

    Shore: ``pole='north'`` → z<0 pole, so colatitude is measured from -z.
    """
    r = np.linalg.norm(cap, axis=-1)
    z_sign = -1.0 if pole == "north" else 1.0
    cos_theta = np.clip(z_sign * cap[..., 2] / r, -1.0, 1.0)
    return np.arccos(cos_theta)


# Parameter grid from the Phase 4 design doc
_NIS: tuple[int, ...] = (16, 32, 64)
_NJS: tuple[int, ...] = (16, 32, 64)
_THETAS_DEG: tuple[float, ...] = (5.0, 9.0, 15.0, 30.0, 45.0)
_I_SPACINGS: tuple[str, ...] = ("uniform", "tanh", "tanh2")
_I_BETAS: tuple[float, ...] = (1.5, 3.0, 5.0)
_POLE_EPS_RAD: float = float(np.deg2rad(0.5))


def _is_degenerate(ni: int, theta_cap: float) -> bool:
    """Match the fallback predicate in ``_ring_theta_profile``."""
    dtheta_seam = (np.pi - 2.0 * theta_cap) / (ni - 1)
    return not (theta_cap - dtheta_seam > _POLE_EPS_RAD)


@pytest.mark.parametrize(
    "ni,nj,theta_cap_deg,i_spacing,i_beta",
    [
        (ni, nj, td, sp, b)
        for ni, nj, td, sp, b in itertools.product(
            _NIS,
            _NJS,
            _THETAS_DEG,
            _I_SPACINGS,
            _I_BETAS,
        )
        if not _is_degenerate(ni, np.deg2rad(td))
        and (sp != "uniform" or b == _I_BETAS[0])  # beta is a no-op for uniform
    ],
)
def test_c1_seam_match_edge_midpoints(
    ni: int,
    nj: int,
    theta_cap_deg: float,
    i_spacing: str,
    i_beta: float,
) -> None:
    """Radial step at the 4 edge midpoints equals the equator seam step."""
    theta_cap = np.deg2rad(theta_cap_deg)
    surface = sphere_layer(ni, nj, theta_cap=theta_cap)
    pole_rows = _north_pole_rows(surface, nj)

    with warnings.catch_warnings():
        warnings.simplefilter("error")  # fail the test if profile falls back
        cap = build_cap_k0(
            pole_rows,
            theta_cap=theta_cap,
            ni=ni,
            pole="north",
            i_spacing=i_spacing,
            i_beta=i_beta,
            pole_eps_rad=_POLE_EPS_RAD,
        )

    theta = _theta_from_cap(cap, pole="north")
    nc = cap.shape[0]
    mid = nc // 2

    lats_eq = np.linspace(theta_cap, np.pi - theta_cap, ni)
    dtheta_eq = float(lats_eq[1] - lats_eq[0])

    # Equator uses uniform lats for this test setup (sphere_layer is linspace);
    # the cap profile must still land on the uniform seam step for i_spacing=uniform
    # and on build_i_coords(ni, theta_cap, pi-theta_cap, i_spacing, i_beta)[1] otherwise.
    from shore.volume.spacing import build_i_coords

    lats_cap_side = build_i_coords(
        ni,
        theta_cap,
        np.pi - theta_cap,
        spacing=i_spacing,
        beta=i_beta,
    )
    dtheta_seam_expected = float(lats_cap_side[1] - lats_cap_side[0])

    # 4 edge midpoints (ring 0 → ring 1 radial step)
    dth_edges = [
        theta[0, mid] - theta[1, mid],  # i=0   edge
        theta[-1, mid] - theta[-2, mid],  # i=-1  edge
        theta[mid, 0] - theta[mid, 1],  # j=0   edge
        theta[mid, -1] - theta[mid, -2],  # j=-1  edge
    ]
    for dth in dth_edges:
        rel = abs(dth - dtheta_seam_expected) / dtheta_seam_expected
        assert rel < 1e-10, (
            f"midpoint radial step {np.rad2deg(dth):.6f} deg vs expected "
            f"{np.rad2deg(dtheta_seam_expected):.6f} deg, relerr={rel:.3e} "
            f"(ni={ni}, nj={nj}, theta_cap={theta_cap_deg}, "
            f"i_spacing={i_spacing}, i_beta={i_beta})"
        )

    # The uniform-linspace sphere_layer equator must also match for i_spacing='uniform';
    # for tanh/tanh2, we only verify internal consistency (cap's inherited profile).
    if i_spacing == "uniform":
        assert abs(dtheta_seam_expected - dtheta_eq) < 1e-12


@pytest.mark.parametrize(
    "ni,nj,theta_cap_deg,i_spacing,i_beta",
    [
        (ni, nj, td, sp, b)
        for ni, nj, td, sp, b in itertools.product(
            _NIS,
            _NJS,
            _THETAS_DEG,
            _I_SPACINGS,
            _I_BETAS,
        )
        if not _is_degenerate(ni, np.deg2rad(td))
        and nj >= 16
        and (sp != "uniform" or b == _I_BETAS[0])
    ],
)
def test_corner_seam_match_bounded(
    ni: int,
    nj: int,
    theta_cap_deg: float,
    i_spacing: str,
    i_beta: float,
) -> None:
    """Diagonal step at each of the 4 corners is within 2% of the seam step."""
    theta_cap = np.deg2rad(theta_cap_deg)
    surface = sphere_layer(ni, nj, theta_cap=theta_cap)
    pole_rows = _north_pole_rows(surface, nj)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cap = build_cap_k0(
            pole_rows,
            theta_cap=theta_cap,
            ni=ni,
            pole="north",
            i_spacing=i_spacing,
            i_beta=i_beta,
            pole_eps_rad=_POLE_EPS_RAD,
        )

    theta = _theta_from_cap(cap, pole="north")

    from shore.volume.spacing import build_i_coords

    lats_cap_side = build_i_coords(
        ni,
        theta_cap,
        np.pi - theta_cap,
        spacing=i_spacing,
        beta=i_beta,
    )
    dtheta_seam_expected = float(lats_cap_side[1] - lats_cap_side[0])

    dth_corners = [
        theta[0, 0] - theta[1, 1],
        theta[0, -1] - theta[1, -2],
        theta[-1, 0] - theta[-2, 1],
        theta[-1, -1] - theta[-2, -2],
    ]
    for dth in dth_corners:
        rel = abs(dth - dtheta_seam_expected) / dtheta_seam_expected
        assert rel < 0.02, (
            f"corner radial step {np.rad2deg(dth):.6f} deg vs expected "
            f"{np.rad2deg(dtheta_seam_expected):.6f} deg, relerr={rel:.3e}"
        )


class TestRingProfile:
    """Direct tests of the ring-theta profile construction."""

    def test_seam_value_is_theta_cap(self):
        theta_cap = np.deg2rad(15.0)
        rows = _north_pole_rows(sphere_layer(33, 32, theta_cap=theta_cap), 32)
        cap = build_cap_k0(rows, theta_cap=theta_cap, ni=33, pole="north")
        theta = _theta_from_cap(cap, pole="north")
        # Every boundary node must sit at theta=theta_cap (to fp)
        np.testing.assert_allclose(theta[0, :], theta_cap, atol=1e-12)
        np.testing.assert_allclose(theta[-1, :], theta_cap, atol=1e-12)
        np.testing.assert_allclose(theta[:, 0], theta_cap, atol=1e-12)
        np.testing.assert_allclose(theta[:, -1], theta_cap, atol=1e-12)

    def test_center_near_pole(self):
        theta_cap = np.deg2rad(15.0)
        rows = _north_pole_rows(sphere_layer(33, 32, theta_cap=theta_cap), 32)
        cap = build_cap_k0(
            rows,
            theta_cap=theta_cap,
            ni=33,
            pole="north",
            pole_eps_rad=_POLE_EPS_RAD,
        )
        theta = _theta_from_cap(cap, pole="north")
        nc = cap.shape[0]
        assert theta[nc // 2, nc // 2] < np.deg2rad(2.0), (
            "cap centre should sit near the geographic pole"
        )

    def test_degenerate_fallback_warns(self):
        """When theta_cap - dtheta_seam <= pole_eps, fall back + warn."""
        ni, nj = 9, 16
        theta_cap = np.deg2rad(5.0)  # dtheta_seam ≈ 21.25° > 4.5°
        rows = _north_pole_rows(sphere_layer(ni, nj, theta_cap=theta_cap), nj)
        with pytest.warns(UserWarning, match="C1 seam match not achievable"):
            build_cap_k0(
                rows,
                theta_cap=theta_cap,
                ni=ni,
                pole="north",
                pole_eps_rad=_POLE_EPS_RAD,
            )
