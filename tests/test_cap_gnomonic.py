"""Gnomonic flat-square cap projection tests.

The cap k=0 seed is built by placing a flat square in the tangent plane at
the cap pole (corners on the cardinal axes at distance ``tan(theta_cap)``)
and ray-casting onto the body surface.  On a unit sphere, the 4 cap corners
land at colatitude ``theta_cap`` exactly; cap edges are great-circle arcs
between adjacent corners; cells are roughly square everywhere by
construction.
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.mesh import CubedSphereTopology, Mesh, Spacing1D
from shore.surface.cap import build_cap_k0


def _theta(p: np.ndarray, pole: str) -> float:
    """Colatitude (radians) of a point about the body's pole axis.

    Shore convention: ``pole='north'`` → cap_north → ``z<0`` of body.
    ``pole='south'`` → cap_south → ``z>0`` of body.
    """
    r = float(np.linalg.norm(p))
    z_sign = -1.0 if pole == "north" else 1.0
    return float(np.arccos(np.clip(z_sign * p[2] / r, -1.0, 1.0)))


def _cap_quad_aspect(layer: np.ndarray) -> np.ndarray:
    """Per-quad max-edge / min-edge aspect ratio."""
    p00 = layer[:-1, :-1]
    p10 = layer[1:, :-1]
    p01 = layer[:-1, 1:]
    p11 = layer[1:, 1:]
    e_top = np.linalg.norm(p01 - p00, axis=-1)
    e_bot = np.linalg.norm(p11 - p10, axis=-1)
    e_left = np.linalg.norm(p10 - p00, axis=-1)
    e_right = np.linalg.norm(p11 - p01, axis=-1)
    edges = np.stack([e_top, e_bot, e_left, e_right])
    return edges.max(axis=0) / np.maximum(edges.min(axis=0), 1e-30)


# ---------------------------------------------------------------------------
# Direct tests on build_cap_k0
# ---------------------------------------------------------------------------


class TestCapCornersOnSphere:
    """The 4 cap corners must land at colatitude theta_cap on a unit sphere."""

    @pytest.mark.parametrize("theta_deg", [3.0, 9.0, 15.0, 30.0, 45.0])
    @pytest.mark.parametrize("pole", ["north", "south"])
    def test_four_corners_at_theta_cap(self, sphere_stl, theta_deg, pole):
        from shore.mesh import Body

        body = Body.from_stl(sphere_stl)
        theta_cap = float(np.deg2rad(theta_deg))
        nc = 9  # nj // 4 + 1 with nj=32
        cap = build_cap_k0(
            theta_cap=theta_cap,
            nc=nc,
            pole=pole,
            projector=body.projector,
        )
        corners = [cap[0, 0], cap[-1, 0], cap[-1, -1], cap[0, -1]]
        for c in corners:
            theta_actual = _theta(c, pole)
            # On a sphere of radius 1, corner colatitude must match theta_cap
            # to within ray-cast / icosphere tessellation noise (~1e-3 rad
            # for an icosphere with ~5000 triangles).
            assert abs(theta_actual - theta_cap) < 5e-3, (
                f"corner {c} at colatitude {np.rad2deg(theta_actual):.3f} deg, "
                f"expected {theta_deg:.3f} deg (pole={pole})"
            )

    def test_cap_corners_on_body_surface(self, sphere_stl):
        """Every cap node must lie on the unit sphere (radius ≈ 1)."""
        from shore.mesh import Body

        body = Body.from_stl(sphere_stl)
        cap = build_cap_k0(
            theta_cap=float(np.deg2rad(9.0)),
            nc=9,
            pole="south",
            projector=body.projector,
        )
        radii = np.linalg.norm(cap, axis=-1)
        # icosphere is not a perfect sphere; allow ~3% radial deviation
        assert radii.min() > 0.95
        assert radii.max() < 1.05


# ---------------------------------------------------------------------------
# Cell quality on the resulting cap
# ---------------------------------------------------------------------------


class TestCapCellQuality:
    """Gnomonic projection of a uniform square produces near-square cells."""

    def test_cap_aspect_ratio_low(self, sphere_stl):
        """Worst cap-quad aspect ratio should be modest (< 2.5).

        The gnomonic-square projection is only mildly anisotropic — corner
        quads stretch ~sqrt(2) along the diagonal — so AR_max should sit
        well below the topology-bad.png Shirley AR (which exceeded 6).
        """
        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=20,
            nj=32,
            nk=3,
            spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
        )
        for cap in (m.blocks[4], m.blocks[5]):
            ar = _cap_quad_aspect(cap.nodes[:, :, 0, :])
            assert ar.max() < 2.5, f"{cap.label}: AR_max={ar.max():.2f}, expected < 2.5"

    def test_cap_no_concentric_pole_pattern(self, sphere_stl):
        """The cap centre must not be a single repeated point.

        Shirley/concentric maps collapsed the cap center to a single ring
        (the ``rho=0`` shell); the gnomonic-square has a well-defined
        2D-grid centre with finite cell size.
        """
        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=20,
            nj=32,
            nk=3,
        )
        cap = m.blocks[4].nodes[:, :, 0, :]
        nc = cap.shape[0]
        c = nc // 2
        # The 4 cells around the cap centre must each have positive area.
        for di, dj in [(-1, -1), (-1, 0), (0, -1), (0, 0)]:
            i, j = c + di, c + dj
            d1 = cap[i + 1, j + 1] - cap[i, j]
            d2 = cap[i, j + 1] - cap[i + 1, j]
            area = 0.5 * np.linalg.norm(np.cross(d1, d2))
            assert area > 1e-6, f"degenerate cap-centre quad at ({i}, {j}), area={area:.2e}"


# ---------------------------------------------------------------------------
# Seam C0 conformity
# ---------------------------------------------------------------------------


class TestSeamConformity:
    """Cap edges must equal the equator's pole rows bit-exactly."""

    def test_north_seam_c0(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=12,
            nj=20,
            nk=3,
        )
        sub_blocks = m.blocks[:4]
        cap_north = m.blocks[4]
        ni = sub_blocks[0].nodes.shape[0]

        np.testing.assert_array_equal(
            cap_north.nodes[0, :, 0, :],
            sub_blocks[0].nodes[ni - 1, :, 0, :],
        )
        np.testing.assert_array_equal(
            cap_north.nodes[:, -1, 0, :],
            sub_blocks[1].nodes[ni - 1, :, 0, :],
        )
        np.testing.assert_array_equal(
            cap_north.nodes[-1, :, 0, :],
            sub_blocks[2].nodes[ni - 1, ::-1, 0, :],
        )
        np.testing.assert_array_equal(
            cap_north.nodes[:, 0, 0, :],
            sub_blocks[3].nodes[ni - 1, ::-1, 0, :],
        )

    def test_south_seam_c0(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=12,
            nj=20,
            nk=3,
        )
        sub_blocks = m.blocks[:4]
        cap_south = m.blocks[5]

        np.testing.assert_array_equal(
            cap_south.nodes[:, 0, 0, :],
            sub_blocks[0].nodes[0, :, 0, :],
        )
        np.testing.assert_array_equal(
            cap_south.nodes[-1, :, 0, :],
            sub_blocks[1].nodes[0, :, 0, :],
        )
        np.testing.assert_array_equal(
            cap_south.nodes[:, -1, 0, :],
            sub_blocks[2].nodes[0, ::-1, 0, :],
        )
        np.testing.assert_array_equal(
            cap_south.nodes[0, :, 0, :],
            sub_blocks[3].nodes[0, ::-1, 0, :],
        )


# ---------------------------------------------------------------------------
# Equator latitude monotonicity
# ---------------------------------------------------------------------------


class TestEquatorMonotone:
    """The equator's colatitude along i must be monotone for each j."""

    @pytest.mark.parametrize("i_spacing", ["uniform", "tanh", "tanh2"])
    def test_meridian_colatitude_monotone(self, sphere_stl, i_spacing):
        """Colatitude (measured from cap_south pole axis) must increase
        monotonically as i runs 0 → ni-1 along each meridian."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=20,
            nj=32,
            nk=3,
            i_spacing=i_spacing,
        )
        sub0 = m.blocks[0].nodes[:, :, 0, :]
        for j in range(sub0.shape[1]):
            thetas = np.array([_theta(sub0[i, j], "south") for i in range(sub0.shape[0])])
            # i=0 is the south seam (smaller theta_south); i=ni-1 is the north
            # seam (larger theta_south, near pi).  Diffs must all be positive
            # (small negative tolerance for floating-point noise).
            assert np.all(np.diff(thetas) > -1e-3), f"j={j} non-monotone: diffs={np.diff(thetas)}"


# ---------------------------------------------------------------------------
# Seam C1 step-matching
# ---------------------------------------------------------------------------


class TestSeamStepMatching:
    """All multi-block edges meeting at a cap-seam vertex must agree on
    first-cell length at that vertex.  Three edges meet at each cap corner:
    two cap-seam edges (cap[:, 0, :] and cap[0, :, :] for sub0/south) and
    one equator meridional edge (sub_q[:, j_local, 0, :]).  The cap-seam
    first cell length at the vertex is determined by the gnomonic-square
    geometry; the equator's meridional first cell must equal it.
    """

    def _arc_anchor(self, p, q, anchor=None):
        """Angle (rad) subtended at *anchor* between p and q."""
        if anchor is None:
            anchor = np.zeros(3)
        a = (p - anchor) / np.linalg.norm(p - anchor)
        b = (q - anchor) / np.linalg.norm(q - anchor)
        return float(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0)))

    def _meridian_first_arc_south(self, sub_blocks, j_global, sector):
        """Equator's south-seam first cell length (great-circle arc on body)."""
        q = j_global // sector
        j_local = j_global % sector
        return self._arc_anchor(
            sub_blocks[q].nodes[0, j_local, 0, :],
            sub_blocks[q].nodes[1, j_local, 0, :],
        )

    def _cap_seam_first_arc_south(self, cap_south, j_global, sector):
        """Cap-seam first cell length at the south seam at global j-index.

        The seam ring concatenates the 4 cap edges (per sub-block layout in
        equator._seam_ring); this is the arc from seam[j] to seam[j+1] (mod
        nj), i.e. the cap-edge first cell *along the seam* at that vertex.
        """
        nc = cap_south.shape[0]
        # Build the seam ring inline (mirrors equator._seam_ring south-pole layout).
        sub_seams = [
            cap_south[:, 0, :],
            cap_south[-1, :, :],
            cap_south[:, -1, :][::-1],
            cap_south[0, :, :][::-1],
        ]
        nj = 4 * (nc - 1)
        ring = np.empty((nj, 3))
        for q, s in enumerate(sub_seams):
            ring[q * (nc - 1) : (q + 1) * (nc - 1)] = s[:-1]
        return self._arc_anchor(ring[j_global], ring[(j_global + 1) % nj])

    def test_first_step_matches_cap_uniform(self, sphere_stl):
        """Equator's first meridional cell length matches the cap-seam first
        cell length at every seam vertex, for i_spacing='uniform' at compatible
        sizes.  At ni=60, nj=60, theta_cap=30° the cap and equator radial steps
        are within the pinning compatibility band, so the C1 match is exact.
        """
        import warnings as _w

        with _w.catch_warnings():
            _w.simplefilter("error")  # any fallback would be a UserWarning
            m = Mesh.from_stl(
                sphere_stl,
                topology=CubedSphereTopology(),
                ni=60,
                nj=60,
                nk=3,
                theta_cap=float(np.deg2rad(30.0)),
                i_spacing="uniform",
            )
        sub_blocks = m.blocks[:4]
        cap_south = m.blocks[5].nodes[:, :, 0, :]
        nc = cap_south.shape[0]
        # Each sub-block "owns" nc-1 of the global ring (the corner with the
        # next sub-block is shared) → sector = nc - 1.
        sector = nc - 1
        nj_total = 4 * sector

        for jg in range(0, nj_total, max(1, nj_total // 16)):
            eq_arc = self._meridian_first_arc_south(sub_blocks, jg, sector)
            cap_arc = self._cap_seam_first_arc_south(cap_south, jg, sector)
            assert abs(eq_arc - cap_arc) < 5e-3, (
                f"j={jg}: eq_first_arc={np.rad2deg(eq_arc):.4f} deg, "
                f"cap_seam_first_arc={np.rad2deg(cap_arc):.4f} deg"
            )

    def test_fallback_warns_when_incompatible(self, sphere_stl):
        """When cap and equator radial steps are incompatible, the equator
        falls back to unpinned spacing and emits a single UserWarning."""
        with pytest.warns(UserWarning, match="incompatible"):
            Mesh.from_stl(
                sphere_stl,
                topology=CubedSphereTopology(),
                # Cap step ~ tan(9 deg)/3 ~ 3 deg; equator natural step
                # (uniform) ~ 162/9 ~ 18 deg.  Ratio ~6x = incompatible.
                # (9 deg was the default cap before ADR 0054: now explicit.)
                ni=10,
                nj=12,
                nk=3,
                theta_cap=float(0.05 * np.pi),
            )


# ---------------------------------------------------------------------------
# Outer-layer seam continuity (volume marching with seam-aware normals)
# ---------------------------------------------------------------------------


class TestVolumeSeamContinuity:
    """The j-step across each sub-sub seam should stay close to 1:1 at every
    marching layer.  Seam-aware normals (centred differences across block
    seams during marching) ensure both sides advance with the same outward
    velocity, so the j-spacing on the two sides agrees to within a few percent
    even at the outermost layer.
    """

    def test_sub_sub_j_seam_ratio(self, sphere_stl):
        """j-step ratio across sub_q/sub_{q+1} seam stays within 5% of 1.0
        at every k layer."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=60,
            nj=60,
            nk=20,
            theta_cap=float(np.pi / 6),
            spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
        )
        m.march(ds=0.02, growth=1.1)
        for q in range(4):
            q_next = (q + 1) % 4
            for k in range(20):
                sub_q = m.blocks[q].nodes[:, :, k, :]
                sub_n = m.blocks[q_next].nodes[:, :, k, :]
                s_left = np.linalg.norm(sub_q[:, -1, :] - sub_q[:, -2, :], axis=-1)
                s_right = np.linalg.norm(sub_n[:, 1, :] - sub_n[:, 0, :], axis=-1)
                ratio = s_right / np.maximum(s_left, 1e-30)
                max_dev = float(np.max(np.abs(ratio - 1.0)))
                assert max_dev < 0.05, (
                    f"sub{q}/sub{q_next} k={k}: max j-step ratio deviation {max_dev:.4f} > 0.05"
                )


# ---------------------------------------------------------------------------
# Marching robustness
# ---------------------------------------------------------------------------


class TestMarchingNoNegativeJacobian:
    """Multi-layer march on a unit sphere produces no negative Jacobians."""

    @pytest.mark.parametrize("nk", [5, 10, 20])
    def test_march(self, sphere_stl, nk):
        # a 30 deg cap: with the default 9 deg cap the nk=20 march folds its
        # cap-corner cells, which ADR 0052 refuses (tests/spec/test_cap_folds.py)
        from shore.quality import scaled_jacobian

        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=20,
            nj=32,
            nk=nk,
            theta_cap=float(np.deg2rad(30.0)),
            spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
        )
        m.march(ds=0.02, growth=1.1)
        assert min(float(scaled_jacobian(b.nodes).min()) for b in m.blocks) > 0.0
