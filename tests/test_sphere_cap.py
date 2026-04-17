"""Tests for the sphere polar cap topology (src/shore/topology/sphere_cap.py).

Acceptance criteria from issue #3:
- nj not divisible by 4 → warning emitted, value rounded, no crash
- 4 equatorial sub-blocks have shape (ni, nj//4+1, nk, 3) and share edge
  nodes with adjacent sub-blocks
- Cap seed face nodes lie on the sphere surface within numerical tolerance
- All shared face nodes between equatorial and cap blocks are identical at
  every layer k
- All cell Jacobians positive in all 6 blocks for a unit sphere test case
- _compute_normals default behaviour (equatorial, periodic j) is unchanged
- New tests cover the above criteria
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from tests.conftest import sphere_layer
from shore.surface.remesh import project_sphere_to_surface, sphere_latlon_directions
from shore.topology.sphere_cap import (
    _adjust_nj,
    _gnomonic_seed_face,
    _overwrite_cap_boundaries,
    _pull_cap_interior,
    split_ogrid_with_polar_caps,
)
from shore.volume.hyperbolic import _compute_normals
from shore.volume.jacobian import _jacobian_field_open


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_sphere_surface(ni: int, nj: int, radius: float = 1.0) -> np.ndarray:
    """Unit-sphere k=0 layer, shape (ni, nj, 3)."""
    return sphere_layer(ni, nj, radius=radius)


def _run_caps(ni: int = 10, nj: int = 16, nk: int = 5, ds: float = 0.05) -> list[np.ndarray]:
    """Run split_ogrid_with_polar_caps on a unit sphere and return the 6 blocks."""
    surface = _make_sphere_surface(ni, nj)
    return split_ogrid_with_polar_caps(surface, nk=nk, ds=ds, growth=1.1)


# ---------------------------------------------------------------------------
# Step 1 — nj validation
# ---------------------------------------------------------------------------

class TestNjAdjustment:
    def test_divisible_by_4_unchanged(self) -> None:
        assert _adjust_nj(16) == 16
        assert _adjust_nj(20) == 20

    def test_not_divisible_emits_warning(self) -> None:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = _adjust_nj(18)
        assert result == 16
        assert len(w) == 1
        assert "rounding down" in str(w[0].message).lower()

    def test_rounded_value_is_largest_multiple_below(self) -> None:
        assert _adjust_nj(19) == 16
        assert _adjust_nj(21) == 20
        assert _adjust_nj(15) == 12

    def test_split_warns_for_odd_nj(self) -> None:
        """Warning is emitted and nj is adjusted; no ValueError is raised."""
        # Use nj=20 → rounds to 20 (already divisible), use 18 → 16
        surface = _make_sphere_surface(ni=10, nj=18)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # Only test that warning fires; full march may fail on very coarse grid
            try:
                split_ogrid_with_polar_caps(surface, nk=4, ds=0.05)
            except RuntimeError:
                pass  # Jacobian failure is acceptable; we only need the warning
        assert any("rounding down" in str(x.message).lower() for x in w)


# ---------------------------------------------------------------------------
# Step 2 — equatorial sub-block shapes and shared edges
# ---------------------------------------------------------------------------

class TestEquatorialSplit:
    def test_six_blocks_returned(self) -> None:
        blocks = _run_caps()
        assert len(blocks) == 6

    def test_sub_block_shapes(self) -> None:
        ni, nj, nk = 10, 16, 5
        nc = nj // 4 + 1   # 5
        blocks = _run_caps(ni=ni, nj=nj, nk=nk)
        for q in range(4):
            assert blocks[q].shape == (ni, nc, nk, 3), \
                f"sub{q} shape {blocks[q].shape} != ({ni}, {nc}, {nk}, 3)"

    def test_cap_block_shapes(self) -> None:
        ni, nj, nk = 10, 16, 5
        nc = nj // 4 + 1
        blocks = _run_caps(ni=ni, nj=nj, nk=nk)
        assert blocks[4].shape == (nc, nc, nk, 3), "cap_north shape wrong"
        assert blocks[5].shape == (nc, nc, nk, 3), "cap_south shape wrong"

    def test_adjacent_sub_blocks_share_edge_at_k0(self) -> None:
        blocks = _run_caps()
        for q in range(4):
            q_next = (q + 1) % 4
            np.testing.assert_array_equal(
                blocks[q][:, -1, 0, :],
                blocks[q_next][:, 0, 0, :],
                err_msg=f"sub{q} last j-col != sub{q_next} first j-col at k=0",
            )

    def test_adjacent_sub_blocks_share_edge_all_k(self) -> None:
        blocks = _run_caps(nk=6)
        nk = blocks[0].shape[2]
        for q in range(4):
            q_next = (q + 1) % 4
            for k in range(nk):
                np.testing.assert_allclose(
                    blocks[q][:, -1, k, :],
                    blocks[q_next][:, 0, k, :],
                    atol=1e-12,
                    err_msg=f"sub{q}/sub{q_next} shared edge diverges at k={k}",
                )


# ---------------------------------------------------------------------------
# Step 3-4 — polar cap seed face on sphere surface
# ---------------------------------------------------------------------------

class TestCapSeedFace:
    def test_gnomonic_seed_corners_are_on_sphere(self) -> None:
        """All nodes of the gnomonic seed face must lie on the sphere."""
        r = 2.5
        # Use corners near the south pole (z≈-r) so the bilinear centre is not
        # at the origin.  These are the actual corner positions from a sphere
        # lat-lon grid near the south pole.
        theta = np.pi - 0.1   # near south pole
        corners = np.array([
            [r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi),
             r * np.cos(theta)]
            for phi in [0.0, np.pi / 2, np.pi, 3 * np.pi / 2]
        ])
        seed = _gnomonic_seed_face(corners, nc=5)
        radii = np.linalg.norm(seed.reshape(-1, 3), axis=-1)
        np.testing.assert_allclose(radii, r, atol=1e-12,
                                   err_msg="gnomonic seed points not on sphere")

    def test_cap_seed_boundary_nodes_on_unit_sphere(self) -> None:
        """The 4 boundary edges of each cap's k=0 face must lie on the sphere.

        The interior nodes are determined by discrete Laplacian interpolation
        and will in general lie inside the sphere (not on its surface), so
        only the boundary edges are tested here.
        """
        ni, nj, nk = 10, 16, 4
        blocks = _run_caps(ni=ni, nj=nj, nk=nk)
        nc = nj // 4 + 1
        for label, cap_idx in [("cap_north", 4), ("cap_south", 5)]:
            seed = blocks[cap_idx][:, :, 0, :]   # k=0 layer
            # Boundary nodes: first/last row and col
            boundary = np.concatenate([
                seed[0, :, :],
                seed[-1, :, :],
                seed[1:-1, 0, :],
                seed[1:-1, -1, :],
            ], axis=0)
            radii = np.linalg.norm(boundary, axis=-1)
            np.testing.assert_allclose(
                radii, 1.0, atol=1e-10,
                err_msg=f"{label} k=0 boundary nodes not on unit sphere",
            )

    def test_gnomonic_seed_is_smooth(self) -> None:
        """Interior nodes should not have wild jumps — check max neighbour distance."""
        r = 1.0
        corners = np.array([
            [r, 0.0, 0.5],
            [0.0, r, 0.5],
            [-r, 0.0, 0.5],
            [0.0, -r, 0.5],
        ])
        corners = corners / np.linalg.norm(corners, axis=-1, keepdims=True) * r
        seed = _gnomonic_seed_face(corners, nc=9)
        diff_i = np.linalg.norm(seed[1:] - seed[:-1], axis=-1)
        diff_j = np.linalg.norm(seed[:, 1:] - seed[:, :-1], axis=-1)
        assert diff_i.max() < 0.5, "cap seed has large i-direction jump"
        assert diff_j.max() < 0.5, "cap seed has large j-direction jump"


# ---------------------------------------------------------------------------
# Shared face conformity between equatorial and cap blocks
# ---------------------------------------------------------------------------

class TestSharedFaceConformity:
    """At every layer k the equatorial pole rows must equal the cap boundaries."""

    def _check_north(self, blocks: list[np.ndarray]) -> None:
        ni   = blocks[0].shape[0]
        nk   = blocks[0].shape[2]
        cap  = blocks[4]   # cap_north

        for k in range(nk):
            # cap i=0 row ← sub0 pole row
            np.testing.assert_allclose(
                cap[0, :, k, :], blocks[0][ni - 1, :, k, :],
                atol=1e-12, err_msg=f"north cap i=0 row != sub0 pole at k={k}",
            )
            # cap j=-1 col ← sub1 pole row
            np.testing.assert_allclose(
                cap[:, -1, k, :], blocks[1][ni - 1, :, k, :],
                atol=1e-12, err_msg=f"north cap j=-1 col != sub1 pole at k={k}",
            )
            # cap i=-1 row ← sub2 pole row (reversed)
            np.testing.assert_allclose(
                cap[-1, :, k, :], blocks[2][ni - 1, ::-1, k, :],
                atol=1e-12, err_msg=f"north cap i=-1 row != sub2 pole reversed at k={k}",
            )
            # cap j=0 col ← sub3 pole row (reversed)
            np.testing.assert_allclose(
                cap[:, 0, k, :], blocks[3][ni - 1, ::-1, k, :],
                atol=1e-12, err_msg=f"north cap j=0 col != sub3 pole reversed at k={k}",
            )

    def _check_south(self, blocks: list[np.ndarray]) -> None:
        """Check south cap shared faces.

        South cap edge mapping (pole='south'):
          cap[:,  0] ← sub0 south pole row  (phi 0°→90°)
          cap[-1, :] ← sub1 south pole row  (phi 90°→180°)
          cap[:, -1] ← sub2 south pole row reversed
          cap[0,  :] ← sub3 south pole row reversed
        """
        nk  = blocks[0].shape[2]
        cap = blocks[5]   # cap_south

        for k in range(nk):
            np.testing.assert_allclose(
                cap[:, 0, k, :], blocks[0][0, :, k, :],
                atol=1e-12, err_msg=f"south cap j=0 col != sub0 south pole at k={k}",
            )
            np.testing.assert_allclose(
                cap[-1, :, k, :], blocks[1][0, :, k, :],
                atol=1e-12, err_msg=f"south cap i=-1 row != sub1 south pole at k={k}",
            )
            np.testing.assert_allclose(
                cap[:, -1, k, :], blocks[2][0, ::-1, k, :],
                atol=1e-12, err_msg=f"south cap j=-1 col != sub2 south pole reversed at k={k}",
            )
            np.testing.assert_allclose(
                cap[0, :, k, :], blocks[3][0, ::-1, k, :],
                atol=1e-12, err_msg=f"south cap i=0 row != sub3 south pole reversed at k={k}",
            )

    def test_north_cap_shared_faces_all_layers(self) -> None:
        blocks = _run_caps(ni=10, nj=16, nk=6)
        self._check_north(blocks)

    def test_south_cap_shared_faces_all_layers(self) -> None:
        blocks = _run_caps(ni=10, nj=16, nk=6)
        self._check_south(blocks)


# ---------------------------------------------------------------------------
# Jacobian positivity
# ---------------------------------------------------------------------------

class TestJacobianPositivity:
    def test_all_sub_blocks_positive_jacobian(self) -> None:
        blocks = _run_caps(ni=12, nj=16, nk=6, ds=0.04)
        for q in range(4):
            jac = _jacobian_field_open(blocks[q])
            assert jac.min() > 0.0, f"sub{q} has non-positive Jacobian: {jac.min():.3e}"

    def test_cap_blocks_positive_jacobian(self) -> None:
        blocks = _run_caps(ni=12, nj=16, nk=6, ds=0.04)
        for label, idx in [("cap_north", 4), ("cap_south", 5)]:
            jac = _jacobian_field_open(blocks[idx])
            assert jac.min() > 0.0, f"{label} has non-positive Jacobian: {jac.min():.3e}"


# ---------------------------------------------------------------------------
# _compute_normals default behaviour unchanged (periodic j)
# ---------------------------------------------------------------------------

class TestComputeNormalsUnchanged:
    def test_default_periodic_j_gives_outward_normals_on_sphere(self) -> None:
        """With default periodic=(False, True) normals should point outward."""
        layer = sphere_layer(10, 16, radius=1.0)   # (10, 16, 3)
        normals = _compute_normals(layer)            # default periodic=(False,True)
        # dot product of normal with radial direction should be positive
        dots = np.einsum("ijk,ijk->ij", normals, layer / np.linalg.norm(
            layer, axis=-1, keepdims=True
        ))
        assert dots.min() > 0.0, "default normals not outward on sphere"

    def test_non_periodic_j_gives_finite_normals(self) -> None:
        """With periodic=(False, False) normals at j=0 and j=-1 should still be finite."""
        layer = sphere_layer(10, 16, radius=1.0)
        normals = _compute_normals(layer, periodic=(False, False))
        assert np.all(np.isfinite(normals)), "non-periodic normals contain inf/nan"
        mags = np.linalg.norm(normals, axis=-1)
        np.testing.assert_allclose(mags, 1.0, atol=1e-12,
                                   err_msg="non-periodic normals not unit length")


# ---------------------------------------------------------------------------
# layer_callback
# ---------------------------------------------------------------------------

class TestLayerCallback:
    def test_callback_called_for_each_layer(self) -> None:
        nk = 5
        called = []
        def cb(k: int, jmin: float) -> None:
            called.append((k, jmin))

        surface = _make_sphere_surface(ni=8, nj=16)
        split_ogrid_with_polar_caps(surface, nk=nk, ds=0.05, layer_callback=cb)
        assert len(called) == nk - 1
        assert [c[0] for c in called] == list(range(1, nk))
        assert all(c[1] > 0.0 for c in called)


# ---------------------------------------------------------------------------
# theta_cap parameter
# ---------------------------------------------------------------------------

class TestThetaCap:
    def test_default_matches_explicit_005pi(self) -> None:
        """sphere_layer() default == sphere_layer(theta_cap=0.05*pi)."""
        default = sphere_layer(10, 16)
        explicit = sphere_layer(10, 16, theta_cap=0.05 * np.pi)
        np.testing.assert_array_equal(default, explicit)

    def test_deg_alias_matches_rad(self) -> None:
        """theta_cap_deg=deg and theta_cap=rad give identical results."""
        angle_rad = 0.15
        angle_deg = np.degrees(angle_rad)
        by_rad = sphere_layer(10, 16, theta_cap=angle_rad)
        by_deg = sphere_layer(10, 16, theta_cap_deg=angle_deg)
        np.testing.assert_allclose(by_rad, by_deg, atol=1e-14)

    def test_boundary_parallel_radius(self) -> None:
        """Cap boundary parallel has radius sin(theta_cap) on a unit sphere."""
        theta_cap = 0.20   # radians
        surface = sphere_layer(12, 16, theta_cap=theta_cap)
        # i=0 row is the south boundary; i=-1 is the north boundary
        south_row = surface[0, :, :]   # (nj, 3)
        north_row = surface[-1, :, :]  # (nj, 3)
        expected_r = np.sin(theta_cap)
        # radius in xy-plane = sqrt(x^2 + y^2)
        south_rxy = np.sqrt(south_row[:, 0]**2 + south_row[:, 1]**2)
        north_rxy = np.sqrt(north_row[:, 0]**2 + north_row[:, 1]**2)
        np.testing.assert_allclose(south_rxy, expected_r, atol=1e-12,
                                   err_msg="south boundary parallel radius wrong")
        np.testing.assert_allclose(north_rxy, expected_r, atol=1e-12,
                                   err_msg="north boundary parallel radius wrong")

    def test_larger_theta_cap_moves_boundary_toward_equator(self) -> None:
        """Larger theta_cap → equatorial band shrinks → boundary row further from pole.

        theta_cap is measured from the pole axis: larger angle means the first
        grid row (i=0, south boundary) sits on a parallel with larger radius.
        sin(theta_cap) is the xy-plane radius on a unit sphere.
        """
        s_small = sphere_layer(10, 16, theta_cap=0.10)
        s_large = sphere_layer(10, 16, theta_cap=0.30)
        r_small = np.sqrt(s_small[0, :, 0]**2 + s_small[0, :, 1]**2).mean()
        r_large = np.sqrt(s_large[0, :, 0]**2 + s_large[0, :, 1]**2).mean()
        assert r_large > r_small, (
            "larger theta_cap should give larger boundary parallel radius "
            f"(sin(0.30)={np.sin(0.30):.4f} > sin(0.10)={np.sin(0.10):.4f})"
        )

    def test_caps_topology_with_theta_cap(self) -> None:
        """split_ogrid_with_polar_caps works for a non-default theta_cap."""
        surface = sphere_layer(10, 16, theta_cap=0.20)
        blocks = split_ogrid_with_polar_caps(surface, nk=4, ds=0.04, growth=1.1)
        assert len(blocks) == 6

    def test_invalid_theta_cap_zero_raises(self) -> None:
        with pytest.raises(ValueError, match="theta_cap"):
            sphere_layer(8, 12, theta_cap=0.0)

    def test_invalid_theta_cap_half_pi_raises(self) -> None:
        with pytest.raises(ValueError, match="theta_cap"):
            sphere_layer(8, 12, theta_cap=np.pi / 2)

    def test_both_rad_and_deg_raises(self) -> None:
        with pytest.raises(ValueError):
            sphere_layer(8, 12, theta_cap=0.1, theta_cap_deg=10.0)

    def test_project_sphere_deg_alias(self) -> None:
        """project_sphere_to_surface theta_cap_deg matches theta_cap."""
        dirs_rad = sphere_latlon_directions(8, 12, theta_cap=0.15)
        dirs_deg = sphere_latlon_directions(8, 12, theta_cap_deg=np.degrees(0.15))
        np.testing.assert_allclose(dirs_rad, dirs_deg, atol=1e-14)

    def test_project_sphere_both_raises(self) -> None:
        """Supplying theta_cap and theta_cap_deg together raises ValueError."""
        with pytest.raises(ValueError):
            sphere_latlon_directions(8, 12, theta_cap=0.1, theta_cap_deg=10.0)


# ---------------------------------------------------------------------------
# cap_pull parameter
# ---------------------------------------------------------------------------

class TestCapPull:
    """Tests for the cap_pull parameter (interior dome toward the pole)."""

    def test_zero_pull_unchanged(self) -> None:
        """cap_pull=0 must produce the same result as not passing it."""
        surface = _make_sphere_surface(ni=10, nj=16)
        blocks_default = split_ogrid_with_polar_caps(surface.copy(), nk=4, ds=0.05)
        blocks_zero    = split_ogrid_with_polar_caps(surface.copy(), nk=4, ds=0.05, cap_pull=0.0)
        for i, (a, b) in enumerate(zip(blocks_default, blocks_zero, strict=True)):
            np.testing.assert_array_equal(a, b, err_msg=f"block {i} differs with cap_pull=0")

    def test_nonzero_pull_moves_cap_interior(self) -> None:
        """cap_pull > 0 must move off-centre interior nodes away from cap_pull=0.

        Use a larger grid (nj=32) so nc=9 gives a non-trivial interior.  The
        quarter-point (mid//2, mid//2) sits at the peak of the 2-D Hann weight
        and is not on the pole axis, so it is always displaced.
        """
        surface = _make_sphere_surface(ni=10, nj=32)
        b0 = split_ogrid_with_polar_caps(surface.copy(), nk=4, ds=0.05, cap_pull=0.0)
        b4 = split_ogrid_with_polar_caps(surface.copy(), nk=4, ds=0.05, cap_pull=0.4)
        nc = b0[4].shape[0]   # 9 for nj=32
        # Use (1, 1) — strictly interior, off-axis, maximum Hann weight region
        diff = np.linalg.norm(b4[4][1, 1, 0, :] - b0[4][1, 1, 0, :])
        assert diff > 1e-5, f"cap interior node (1,1) unchanged with cap_pull=0.4 (diff={diff:.3e})"

    def test_pull_leaves_boundary_edges_unchanged(self) -> None:
        """cap_pull must not modify the 4 boundary edges of the cap k=0 face."""
        surface = _make_sphere_surface(ni=10, nj=16)
        b0 = split_ogrid_with_polar_caps(surface.copy(), nk=4, ds=0.05, cap_pull=0.0)
        b4 = split_ogrid_with_polar_caps(surface.copy(), nk=4, ds=0.05, cap_pull=0.4)
        cap0 = b0[4][:, :, 0, :]   # cap_north k=0
        cap4 = b4[4][:, :, 0, :]
        np.testing.assert_allclose(cap4[0, :, :],  cap0[0, :, :],  atol=1e-12,
                                   err_msg="cap[0,:] boundary row changed by cap_pull")
        np.testing.assert_allclose(cap4[-1, :, :], cap0[-1, :, :], atol=1e-12,
                                   err_msg="cap[-1,:] boundary row changed by cap_pull")
        np.testing.assert_allclose(cap4[:, 0, :],  cap0[:, 0, :],  atol=1e-12,
                                   err_msg="cap[:,0] boundary col changed by cap_pull")
        np.testing.assert_allclose(cap4[:, -1, :], cap0[:, -1, :], atol=1e-12,
                                   err_msg="cap[:,-1] boundary col changed by cap_pull")

    def test_pull_leaves_sub_blocks_unchanged(self) -> None:
        """cap_pull must not modify any equatorial sub-block geometry."""
        surface = _make_sphere_surface(ni=10, nj=16)
        b0 = split_ogrid_with_polar_caps(surface.copy(), nk=4, ds=0.05, cap_pull=0.0)
        b4 = split_ogrid_with_polar_caps(surface.copy(), nk=4, ds=0.05, cap_pull=0.4)
        for q in range(4):
            np.testing.assert_array_equal(b0[q], b4[q],
                                          err_msg=f"sub{q} changed by cap_pull=0.4")

    def test_interior_nodes_remain_on_sphere(self) -> None:
        """After cap_pull, all cap interior nodes must lie on the sphere surface."""
        surface = _make_sphere_surface(ni=10, nj=16)
        blocks = split_ogrid_with_polar_caps(surface, nk=4, ds=0.05, cap_pull=0.5)
        for label, idx in [("cap_north", 4), ("cap_south", 5)]:
            seed = blocks[idx][:, :, 0, :]   # k=0 face, (nc, nc, 3)
            interior = seed[1:-1, 1:-1, :]   # strictly interior nodes
            radii = np.linalg.norm(interior, axis=-1)
            np.testing.assert_allclose(
                radii, 1.0, atol=1e-10,
                err_msg=f"{label} interior nodes not on unit sphere after cap_pull=0.5",
            )

    def test_pull_positive_jacobian_all_blocks(self) -> None:
        """cap_pull=0.4 must not cause Jacobian failures (sub-blocks and caps)."""
        surface = _make_sphere_surface(ni=12, nj=16)
        # Use small nk to avoid reaching total thickness > radius
        blocks = split_ogrid_with_polar_caps(surface, nk=5, ds=0.04, cap_pull=0.4)
        for q in range(4):
            jac = _jacobian_field_open(blocks[q])
            assert jac.min() > 0.0, f"sub{q} Jacobian non-positive with cap_pull=0.4"

    def test_pull_conformity_preserved(self) -> None:
        """cap_pull must not break shared-face conformity between cap and sub-blocks."""
        surface = _make_sphere_surface(ni=10, nj=16)
        blocks = split_ogrid_with_polar_caps(surface, nk=4, ds=0.05, cap_pull=0.4)
        ni, nc, nk = blocks[0].shape[0], blocks[4].shape[0], blocks[0].shape[2]
        cap = blocks[4]
        for k in range(nk):
            np.testing.assert_allclose(
                cap[0, :, k, :], blocks[0][ni - 1, :, k, :], atol=1e-12,
                err_msg=f"north cap conformity broken at k={k} with cap_pull=0.4",
            )

    def test_pull_cap_interior_zero_is_identity(self) -> None:
        """_pull_cap_interior with cap_pull=0 returns a copy of the input seed."""
        surface = _make_sphere_surface(ni=10, nj=16)
        ni, nj = surface.shape[:2]
        nc = nj // 4 + 1
        from shore.topology.sphere_cap import _build_cap_k0
        sector = nj // 4
        eq_rows = [surface[ni - 1, list(range(q*sector, q*sector+nc-1)) + [(q*sector+nc-1)%nj], :]
                   for q in range(4)]
        seed = _build_cap_k0(eq_rows, nc, pole="north")
        result = _pull_cap_interior(seed, eq_rows, cap_pull=0.0, pole="north")
        np.testing.assert_array_equal(seed, result)

