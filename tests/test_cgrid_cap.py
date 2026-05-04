"""Tests for shore.topology.cgrid_cap.build_cap_blocks — pure-numpy cap fill.

The cap-fill builder takes a body's j=0 cross-section polyline (the
airfoil profile after the TE slivers have been absorbed into the wake
blocks) plus geometric parameters, and returns the cap sub-block node
arrays at one wing tip.

These tests exercise the builder in isolation — no STL, no SHORE-mesh
integration, no marching.  They pin down:

- Block shapes and dtypes.
- Bit-exact node coincidence at SHARED seams between cap sub-blocks.
- The cap_jlo plane's j_lo face = wing-tip wall (offset by tip_inset
  along span_dir from the body cross-section).
- The cap_jlo plane's k_hi face = body cross-section (= wrap envelope's
  j_lo, before the body is split).
- Sharp vs blunt switch — sharp omits cap_te; blunt includes it.
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.topology.cgrid_cap import (
    auto_round_sharp_te,
    build_cap_blocks,
    detect_te_thickness,
    find_symmetric_le_cuts,
)


def _naca0012_cross_section(ni: int = 41) -> tuple[np.ndarray, int, int]:
    """Synthetic NACA 0012-style cross section, post-extension.

    Sweeps upper-TE-cut → suction-side → LE → pressure-side →
    lower-TE-cut.  Returns (cross_section, le_cut_upper_idx,
    le_cut_lower_idx).  The cuts sit at the i indices closest to
    chord = 0.05 from each side of the LE point.
    """
    # Use a cosine distribution so the LE point is sharply resolved.
    n_half = ni // 2 + 1
    chord = np.linspace(0.0, 1.0, n_half)
    # NACA 0012 thickness distribution (open trailing edge — nonzero at x=1).
    yt = 0.6 * (
        0.2969 * np.sqrt(chord)
        - 0.1260 * chord
        - 0.3516 * chord**2
        + 0.2843 * chord**3
        - 0.1015 * chord**4
    )
    # The "post-extension" body sweeps from upper-TE-cut to lower-TE-cut,
    # so we need to clip both ends at, say, 95% chord.
    te_cut_frac = 0.95
    le_cut_frac = 0.05
    upper_x = chord[chord <= te_cut_frac]
    upper_y = yt[chord <= te_cut_frac]
    # Build a closed-style sweep: upper-TE-cut (rightmost upper) →
    # decreasing x along suction → LE (x=0) → increasing x along
    # pressure → lower-TE-cut.
    # Reverse the upper to start at TE-cut.
    upper = np.column_stack([upper_x[::-1], +upper_y[::-1]])
    lower = np.column_stack([upper_x[1:], -upper_y[1:]])  # skip duplicate LE
    pts2d = np.vstack([upper, lower])
    # Lift to 3D (z = 0).
    pts3d = np.column_stack([pts2d, np.zeros(len(pts2d))])

    # Find the LE-cut indices (chord = le_cut_frac on each side).
    # The LE point is where x is minimum.
    le_idx = int(np.argmin(pts3d[:, 0]))
    # Search backward (toward i=0) for the first node with x >= le_cut_frac.
    le_cut_upper = next(i for i in range(le_idx, -1, -1) if pts3d[i, 0] >= le_cut_frac)
    le_cut_lower = next(i for i in range(le_idx, len(pts3d)) if pts3d[i, 0] >= le_cut_frac)
    return pts3d, le_cut_upper, le_cut_lower


def _box_cross_section(ni: int = 41) -> tuple[np.ndarray, int, int]:
    """Trivial closed-box cross section for sanity checks.

    Sweeps a flat-walled prism in 3D so the cap geometry is
    deterministic and checkable by hand.  The "airfoil" here is just
    a rectangle from x = 0.05 to x = 0.95 with y = ±0.05.  i goes
    upper-TE-cut → upper-edge → LE-cut → lower-edge → lower-TE-cut.
    """
    n_half = ni // 2 + 1
    xs = np.linspace(0.05, 0.95, n_half)
    upper = np.column_stack([xs[::-1], +0.05 * np.ones(n_half)])
    lower = np.column_stack([xs[1:], -0.05 * np.ones(n_half - 1)])
    pts2d = np.vstack([upper, lower])
    pts3d = np.column_stack([pts2d, np.zeros(len(pts2d))])
    # LE point is at x=0.05 (i = n_half - 1 in the upper half).
    le_idx = int(np.argmin(pts3d[:, 0]))
    # The "LE-cut" for this rectangle: pick a chord-position cut at, say,
    # x = 0.2 (= 5% from x=0.05 along chord-of-0.9).
    le_cut_x = 0.20
    le_cut_upper = next(i for i in range(le_idx, -1, -1) if pts3d[i, 0] >= le_cut_x)
    le_cut_lower = next(i for i in range(le_idx, len(pts3d)) if pts3d[i, 0] >= le_cut_x)
    return pts3d, le_cut_upper, le_cut_lower


# ---------------------------------------------------------------------------
# Shape & dtype
# ---------------------------------------------------------------------------


class TestCapShapes:
    def test_default_returns_five_blocks(self):
        """5-block H-fill: le + upper + lower + te_upper + te_lower."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        assert set(caps.keys()) == {"le", "upper", "lower", "te_upper", "te_lower"}

    def test_blunt_returns_six_blocks(self):
        """When the body has wake_base too: 5 cap + 1 cap_wake_base."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        wake_base_cs = np.array(
            [
                [[0.95, +0.01, 0.0], [0.95, 0.0, 0.0], [0.95, -0.01, 0.0]],
                [[1.20, +0.01, 0.0], [1.20, 0.0, 0.0], [1.20, -0.01, 0.0]],
            ]
        )  # (ni_wake=2, nk_w=3, 3)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=True,
            wake_base_cross_section=wake_base_cs,
        )
        assert set(caps.keys()) == {"le", "upper", "lower", "te_upper", "te_lower", "wake_base"}

    def test_block_rank_4_with_k_eq_2(self):
        """Cap blocks have shape (ni, nj, 2, 3) — 1 cell thick in k."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        for blk in caps.values():
            assert blk.ndim == 4
            assert blk.shape[2] == 2  # 1 cell thick in k (spanwise)
            assert blk.shape[3] == 3
            assert blk.dtype == np.float64

    def test_upper_lower_have_same_ni_and_nj(self):
        """SHARED chord-cut between upper/lower needs matching i and j."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        u = caps["upper"]
        lo = caps["lower"]
        assert u.shape[0] == lo.shape[0], "i-extent must match for chord-cut SHARED"
        assert u.shape[1] == lo.shape[1], "j-extent must match"

    def test_nk_cap_greater_than_two(self):
        """Cap blocks support multi-cell spanwise resolution."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            nk_cap=8,
            is_blunt=False,
        )
        for blk in caps.values():
            assert blk.shape[2] == 8

    def test_nk_cap_layers_uniform_spaced_along_span(self):
        """Interior k-layers are uniformly spaced between outboard and inboard planes.

        k=0 is the **outboard** wing-tip wall (offset by +tip_inset *
        span_dir from the cross-section); k=nk_cap-1 is the cross-section
        plane.  So each k-step moves *inboard* (= -span_dir * inset/(nk-1)).
        """
        cs, le_up, le_lo = _naca0012_cross_section(41)
        nk_cap = 6
        inset = 0.1
        span = np.array([0.0, 0.0, 1.0])
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=span,
            tip_inset=inset,
            n_cap_j=10,
            nk_cap=nk_cap,
            is_blunt=False,
        )
        u = caps["upper"]
        for k in range(nk_cap - 1):
            step = u[:, :, k + 1, :] - u[:, :, k, :]
            expected = -(inset / (nk_cap - 1)) * span
            np.testing.assert_allclose(step, np.broadcast_to(expected, step.shape), atol=1e-12)

    def test_le_block_wraps_around_le(self):
        """cap_le's i-axis spans from upper-LE-cut to lower-LE-cut."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        ni_le = le_lo - le_up + 1
        assert caps["le"].shape[0] == ni_le


# ---------------------------------------------------------------------------
# Wing-tip wall (k_lo) lies on the cross section, offset by tip_inset
# ---------------------------------------------------------------------------


class TestCapKAxis:
    def test_klo_is_outboard_wall(self):
        """cap.k=0 (outboard wall) = cross-section + tip_inset * span_dir.

        With the chimera-overset cap convention, k=0 is the **outboard
        wing-tip wall** offset *outward* from the body's cross-section
        plane, and k=nk-1 sits on the cross-section plane (the chimera-
        overlap face).
        """
        cs, le_up, le_lo = _naca0012_cross_section(41)
        span = np.array([0.0, 0.0, 1.0])
        inset = 0.012345
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=span,
            tip_inset=inset,
            n_cap_j=8,
            is_blunt=False,
        )
        u = caps["upper"]
        # u[0, -1, -1, :] is on the cross-section plane (chimera-overlap).
        # u[0, -1, 0,  :] is the same node offset outward by +inset * span.
        offset = u[0, -1, 0, :] - u[0, -1, -1, :]
        np.testing.assert_allclose(
            offset,
            inset * span,
            atol=1e-12,
            err_msg="k_lo - k_hi should equal +tip_inset * span_dir",
        )

    def test_klo_to_khi_uses_span_direction(self):
        """For every node, k_lo - k_hi equals +tip_inset * span_dir exactly."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        span = np.array([0.6, 0.8, 0.0])
        span /= np.linalg.norm(span)
        inset = 0.05
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=span,
            tip_inset=inset,
            n_cap_j=8,
            is_blunt=False,
        )
        for blk in caps.values():
            offset = blk[..., 0, :] - blk[..., -1, :]
            expected = np.broadcast_to(inset * span, offset.shape)
            np.testing.assert_allclose(
                offset,
                expected,
                atol=1e-12,
                err_msg="k_lo - k_hi isn't +inset * span_dir",
            )


# ---------------------------------------------------------------------------
# k_hi face equals the input cross-section / wrap envelope j_lo
# ---------------------------------------------------------------------------


class TestKHiOnCrossSection:
    def test_upper_khi_jhi_traces_suction_side(self):
        """cap_upper.k_hi[i, -1, :] == cross_section[upper-LE-cut → upper-TE-cut]."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        u = caps["upper"]
        # u's i goes from upper-LE-cut (i=0) to upper-TE-cut (i=ni-1).
        # j_hi (j=nj-1) is the airfoil suction wall.
        # k_hi (k=1) is the body plane.
        # So u[i, -1, 1, :] should equal cs[le_up - i] for i = 0..le_up
        # (body i goes upper-TE-cut to lower-TE-cut; u i is reversed).
        for i in range(u.shape[0]):
            body_i = le_up - i
            np.testing.assert_allclose(u[i, -1, 1, :], cs[body_i, :], atol=1e-12)

    def test_lower_khi_jhi_traces_pressure_side(self):
        """cap_lower.k_hi[i, -1, :] == cross_section[le_lo → ni-1] forward."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        lo = caps["lower"]
        # lo's i goes from lower-LE-cut (i=0) to lower-TE-cut (i=ni-1).
        for i in range(lo.shape[0]):
            body_i = le_lo + i
            np.testing.assert_allclose(lo[i, -1, 1, :], cs[body_i, :], atol=1e-12)

    def test_le_khi_jhi_wraps_le(self):
        """cap_le.k_hi[i, -1, :] == cross_section[le_up + i] for i = 0..ni_le-1."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        le = caps["le"]
        for i in range(le.shape[0]):
            body_i = le_up + i
            np.testing.assert_allclose(le[i, -1, 1, :], cs[body_i, :], atol=1e-12)


# ---------------------------------------------------------------------------
# SHARED-seam coincidence between cap sub-blocks (bit-exact node match)
# ---------------------------------------------------------------------------


class TestCapSharedSeams:
    def test_le_ilo_matches_upper_ilo(self):
        """cap_le.i_lo (column at upper-LE-cut) == cap_upper.i_lo."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        le = caps["le"]
        u = caps["upper"]
        # le.i_lo = le[0, :, :, :], shape (nj, 2, 3).
        # u.i_lo = u[0, :, :, :], shape (nj, 2, 3).
        # Both columns should be at the upper-LE-cut chord position with
        # j running from chord-cut/chord-perp inner (j=0) to airfoil wall (j=-1).
        np.testing.assert_allclose(le[0, :, :, :], u[0, :, :, :], atol=1e-12)

    def test_le_ihi_matches_lower_ilo(self):
        """cap_le.i_hi (column at lower-LE-cut) == cap_lower.i_lo."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        le = caps["le"]
        lo = caps["lower"]
        np.testing.assert_allclose(le[-1, :, :, :], lo[0, :, :, :], atol=1e-12)

    def test_upper_jlo_matches_lower_jlo(self):
        """Chord-cut SHARED: cap_upper.j_lo == cap_lower.j_lo (both at chord cut)."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        u = caps["upper"]
        lo = caps["lower"]
        # Both blocks have j=0 = chord cut, j=nj-1 = airfoil wall (each on
        # its own side: upper.j_hi = suction wall, lower.j_hi = pressure
        # wall).  So the chord-cut SHARED face is upper.j_lo ↔ lower.j_lo,
        # but with i-reversal: upper.i goes upper-LE-cut → upper-TE-cut,
        # lower.i goes lower-LE-cut → lower-TE-cut.  At the chord cut,
        # the i-coordinates are the same (chord positions match), so
        # the j_lo columns match index-for-index.
        np.testing.assert_allclose(u[:, 0, :, :], lo[:, 0, :, :], atol=1e-12)


# ---------------------------------------------------------------------------
# j-axis distribution under the cap_j_spacing law
# ---------------------------------------------------------------------------


class TestCapJSpacing:
    def test_uniform_distribution(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=11,
            j_spacing="uniform",
            is_blunt=False,
        )
        # Pick a mid-i column of cap_upper; the j coords should be linear.
        u = caps["upper"]
        ic = u.shape[0] // 2
        col = u[ic, :, 1, :]  # k_hi (body plane), j running cut to wall
        # Distance along the line from j=0 to j=k should be linear in k/n.
        ds = np.linalg.norm(np.diff(col, axis=0), axis=-1)
        # All step sizes should be equal under uniform.
        np.testing.assert_allclose(ds, ds[0], rtol=1e-9)

    def test_tanh_clusters_near_wall(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=20,
            j_spacing="tanh",
            j_beta=4.0,
            is_blunt=False,
        )
        u = caps["upper"]
        ic = u.shape[0] // 2
        col = u[ic, :, 1, :]
        ds = np.linalg.norm(np.diff(col, axis=0), axis=-1)
        # Cluster near the wall (j_hi).  ds[-1] (near wall) << ds[0] (near cut).
        assert ds[-1] < 0.5 * ds[0]


# ---------------------------------------------------------------------------
# Parameter validation
# ---------------------------------------------------------------------------


class TestCapValidation:
    def test_le_cut_ordering(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="le_cut_upper_idx"):
            build_cap_blocks(
                cross_section=cs,
                le_cut_upper_idx=le_lo,  # swapped
                le_cut_lower_idx=le_up,
                span_dir=np.array([0.0, 0.0, 1.0]),
                tip_inset=1e-3,
                n_cap_j=10,
                is_blunt=False,
            )

    def test_negative_tip_inset(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="tip_inset"):
            build_cap_blocks(
                cross_section=cs,
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                span_dir=np.array([0.0, 0.0, 1.0]),
                tip_inset=-0.01,
                n_cap_j=10,
                is_blunt=False,
            )

    def test_n_cap_j_too_small(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="n_cap_j"):
            build_cap_blocks(
                cross_section=cs,
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                span_dir=np.array([0.0, 0.0, 1.0]),
                tip_inset=1e-3,
                n_cap_j=1,
                is_blunt=False,
            )

    def test_blunt_requires_wake_base(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="wake_base"):
            build_cap_blocks(
                cross_section=cs,
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                span_dir=np.array([0.0, 0.0, 1.0]),
                tip_inset=1e-3,
                n_cap_j=10,
                is_blunt=True,
                wake_base_cross_section=None,
            )

    def test_unknown_j_spacing(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="j_spacing"):
            build_cap_blocks(
                cross_section=cs,
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                span_dir=np.array([0.0, 0.0, 1.0]),
                tip_inset=1e-3,
                n_cap_j=10,
                j_spacing="bogus",
                is_blunt=False,
            )

    def test_span_dir_not_unit(self):
        """Non-unit span_dir is rejected — caller must normalise."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="span_dir"):
            build_cap_blocks(
                cross_section=cs,
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                span_dir=np.array([0.0, 0.0, 2.0]),  # not unit
                tip_inset=1e-3,
                n_cap_j=10,
                is_blunt=False,
            )


# ---------------------------------------------------------------------------
# find_symmetric_le_cuts
# ---------------------------------------------------------------------------


class TestFindSymmetricLeCuts:
    def test_finds_indices_symmetric_around_le(self):
        cs, _le_up, _le_lo = _naca0012_cross_section(41)
        # Reference: the helper used by the test fixture itself.
        le_up_helper, le_lo_helper = find_symmetric_le_cuts(cs, le_frac=0.05, span_axis_world=2)
        # Helper must be symmetric (LE point at midpoint).
        le_idx = int(np.argmin(cs[:, 0]))
        assert le_up_helper + le_lo_helper == 2 * le_idx
        # And both indices map to chord ~ 0.05.
        for idx in (le_up_helper, le_lo_helper):
            assert abs(cs[idx, 0] - 0.05) < 0.05

    def test_invalid_frac_rejected(self):
        cs, _, _ = _naca0012_cross_section(41)
        with pytest.raises(ParameterError):
            find_symmetric_le_cuts(cs, le_frac=0.0, span_axis_world=2)
        with pytest.raises(ParameterError):
            find_symmetric_le_cuts(cs, le_frac=1.0, span_axis_world=2)

    def test_invalid_span_axis_rejected(self):
        cs, _, _ = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="span_axis_world"):
            find_symmetric_le_cuts(cs, le_frac=0.05, span_axis_world=3)


# ---------------------------------------------------------------------------
# TE-cut: avoids the singular axis at sharp TE
# ---------------------------------------------------------------------------


def _sharp_te_cross_section(ni: int = 41) -> tuple[np.ndarray, int, int]:
    """Sharp-TE NACA 0012-style cross-section running TE point → suction → LE → pressure → TE point.

    Unlike ``_naca0012_cross_section`` (which truncates at 95% chord),
    this fixture's i=0 and i=ni-1 *both* sit at the sharp TE point —
    realistic input for the C-grid body block.  Used to verify the
    cap-builder's TE-cut singularity-avoidance.
    """
    n_half = ni // 2 + 1
    chord = np.linspace(0.0, 1.0, n_half)
    yt = 0.6 * (
        0.2969 * np.sqrt(chord)
        - 0.1260 * chord
        - 0.3516 * chord**2
        + 0.2843 * chord**3
        - 0.1015 * chord**4
    )
    # Force a sharp TE: y(x=1) = 0 by setting the last-point thickness to zero.
    yt[-1] = 0.0
    upper = np.column_stack([chord[::-1], +yt[::-1]])
    lower = np.column_stack([chord[1:], -yt[1:]])
    pts2d = np.vstack([upper, lower])
    pts3d = np.column_stack([pts2d, np.zeros(len(pts2d))])
    le_idx = int(np.argmin(pts3d[:, 0]))
    le_cut_upper = next(i for i in range(le_idx, -1, -1) if pts3d[i, 0] >= 0.05)
    le_cut_lower = next(i for i in range(le_idx, len(pts3d)) if pts3d[i, 0] >= 0.05)
    return pts3d, le_cut_upper, le_cut_lower


class TestCapTECut:
    def test_sharp_te_te_columns_collapse(self):
        """On a truly sharp TE input, the cap_te_*.i_hi columns collapse.

        This is the structural reality of the H-fill on a sharp TE,
        and is the reason ``CGridTipCap`` rejects sharp inputs unless
        auto-rounding is enabled.
        """
        cs, le_up, le_lo = _sharp_te_cross_section(41)
        assert np.allclose(cs[0, :], cs[-1, :])
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        # cap_te_upper.i_hi is the TE-base column (= body's i=0 node).
        # On sharp TE, body's i=0 == ni-1 == sharp point → degenerate.
        col = caps["te_upper"][-1, :, 0, :]
        span = float(np.linalg.norm(col.max(axis=0) - col.min(axis=0)))
        assert span < 1e-6

    def test_blunt_te_columns_non_degenerate(self):
        """On a blunt TE input, the cap_te_*.i_hi columns are non-degenerate."""
        cs, le_up, le_lo = _naca0012_cross_section(41)  # already blunt-truncated at 95% chord
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        col = caps["te_upper"][-1, :, 0, :]
        span = float(np.linalg.norm(col.max(axis=0) - col.min(axis=0)))
        assert span > 1e-3

    def test_te_cut_inboard_avoids_degeneracy(self):
        """TE-cut indices inboard of the body endpoints → te_upper.i_lo non-degenerate."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        te_up = 3
        te_lo = cs.shape[0] - 1 - te_up
        caps = build_cap_blocks(
            cross_section=cs,
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            te_cut_upper_idx=te_up,
            te_cut_lower_idx=te_lo,
            span_dir=np.array([0.0, 0.0, 1.0]),
            tip_inset=1e-3,
            n_cap_j=10,
            is_blunt=False,
        )
        # cap_upper.i_hi sits at body i = te_up — must be non-degenerate.
        col = caps["upper"][-1, :, 0, :]
        span = float(np.linalg.norm(col.max(axis=0) - col.min(axis=0)))
        assert span > 1e-3

    def test_te_cut_asymmetric_rejected(self):
        """LE/TE-cuts must be symmetric for the chord-cut SHARED face to be 1:1."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="symmetric"):
            build_cap_blocks(
                cross_section=cs,
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                te_cut_upper_idx=3,
                te_cut_lower_idx=cs.shape[0] - 5,  # asymmetric
                span_dir=np.array([0.0, 0.0, 1.0]),
                tip_inset=1e-3,
                n_cap_j=10,
                is_blunt=False,
            )


# ---------------------------------------------------------------------------
# detect_te_thickness + auto_round_sharp_te
# ---------------------------------------------------------------------------


class TestDetectTEThickness:
    def test_sharp_returns_zero(self):
        cs, _, _ = _sharp_te_cross_section(41)
        assert detect_te_thickness(cs, span_axis_world=2) < 1e-9

    def test_blunt_returns_thickness(self):
        cs, _, _ = _naca0012_cross_section(41)
        # _naca0012_cross_section is already truncated at 95% chord,
        # so i=0 (top) and i=ni-1 (bottom) sit at y = ±yt(0.95).
        te = detect_te_thickness(cs, span_axis_world=2)
        # Expected: 2 * 0.6 * (NACA polynomial at chord=0.95).
        assert te > 1e-3


class TestAutoRoundSharpTE:
    def test_round_a_sharp_cross_section(self):
        cs, _, _ = _sharp_te_cross_section(41)
        rounded = auto_round_sharp_te(cs, te_thickness=0.002, span_axis_world=2)
        # i=0 and i=ni-1 should now differ by exactly 0.002 in y.
        delta = rounded[-1, :] - rounded[0, :]
        assert abs(delta[1] + 0.002) < 1e-12  # i=ni-1 moved -te/2 = -0.001
        assert abs(delta[0]) < 1e-12  # x unchanged
        assert abs(delta[2]) < 1e-12  # z unchanged

    def test_round_preserves_midpoint(self):
        cs, _, _ = _sharp_te_cross_section(41)
        rounded = auto_round_sharp_te(cs, te_thickness=0.005, span_axis_world=2)
        midpoint = 0.5 * (rounded[0, :] + rounded[-1, :])
        np.testing.assert_allclose(midpoint, cs[0, :], atol=1e-12)

    def test_round_preserves_interior_nodes(self):
        cs, _, _ = _sharp_te_cross_section(41)
        rounded = auto_round_sharp_te(cs, te_thickness=0.003, span_axis_world=2)
        np.testing.assert_array_equal(rounded[1:-1], cs[1:-1])

    def test_blunt_input_rejected(self):
        cs, _, _ = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="not sharp"):
            auto_round_sharp_te(cs, te_thickness=0.002, span_axis_world=2)

    def test_negative_thickness_rejected(self):
        cs, _, _ = _sharp_te_cross_section(41)
        with pytest.raises(ParameterError, match="te_thickness"):
            auto_round_sharp_te(cs, te_thickness=-0.001, span_axis_world=2)
