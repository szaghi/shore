"""Tests for shore.topology.cgrid_cap: the O-H tip cap (ADR 0059).

The builder takes a wing's tip-station cross-section (upper TE -> LE ->
lower TE, blunt or rounded TE) and returns five right-handed blocks: a
core and the ring ``upper`` / ``le`` / ``lower`` / ``te`` joining it to
the section's wall, extruded outwards. These tests exercise it in
isolation: shapes, the wall, the layers, the flare, the SHARED wiring,
the validation and ``check_cap_blocks`` (ADR 0058).
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.face import FaceBC
from shore.quality import scaled_jacobian
from shore.topology.cgrid_cap import (
    CAP_BLOCKS,
    auto_round_sharp_te,
    build_cap_blocks,
    check_cap_blocks,
    detect_te_thickness,
    find_symmetric_le_cuts,
    section_normal,
    wire_cap_blocks,
)

DIST = np.array([0.0, 0.01, 0.025, 0.05, 0.1])
UP = np.array([0.0, 0.0, 1.0])


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


def _caps(cs=None, le_up=None, le_lo=None, **kw):
    if cs is None:
        cs, le_up, le_lo = _naca0012_cross_section(41)
    kw.setdefault("outward", UP)
    kw.setdefault("k_distances", DIST)
    kw.setdefault("n_cap_j", 8)
    return build_cap_blocks(cs, le_cut_upper_idx=le_up, le_cut_lower_idx=le_lo, **kw)


def _hex(arrays):
    return [
        HexBlock(label=n, nodes=arrays[n], periodic_i=False, periodic_j=False) for n in CAP_BLOCKS
    ]


class TestShapes:
    def test_five_blocks(self):
        assert tuple(_caps()) == CAP_BLOCKS

    def test_counts(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = _caps(cs, le_up, le_lo)
        n_side, n_nose = le_up + 1, le_lo - le_up + 1
        assert sorted(caps["core"].shape[:2]) == sorted((n_side, n_nose))
        for name in ("upper", "lower"):
            assert n_side in caps[name].shape[:2] and 8 in caps[name].shape[:2]
        for name in ("le", "te"):
            assert n_nose in caps[name].shape[:2] and 8 in caps[name].shape[:2]
        assert {c.shape[2] for c in caps.values()} == {DIST.size}

    @pytest.mark.parametrize(
        "kw",
        [{}, {"flare_deg": 0.0}, {"flare_deg": 60.0}, {"core_frac": 0.2}, {"core_frac": 0.8},
         {"j_spacing": "uniform"}, {"j_spacing": "tanh2"}, {"outward": -UP}],
    )  # fmt: skip
    def test_right_handed_and_valid(self, kw):
        for name, x in _caps(**kw).items():
            assert float(scaled_jacobian(x).min()) > 0.2, name

    def test_box_section(self):
        """A rectangle: its square nose is the worst case, still valid."""
        cs, le_up, le_lo = _box_cross_section(41)
        for name, x in _caps(cs, le_up, le_lo).items():
            assert float(scaled_jacobian(x).min()) > 0.1, name


class TestWall:
    def test_edge_is_the_section(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        caps = _caps(cs, le_up, le_lo)
        wall = {p.tobytes() for x in caps.values() for p in x[:, :, 0].reshape(-1, 3)}
        assert all(p.tobytes() in wall for p in cs)

    def test_wall_in_the_section_plane(self):
        for x in _caps().values():
            np.testing.assert_array_equal(x[:, :, 0, 2], 0.0)

    def test_tip_inset_moves_the_wall_outwards(self):
        for sign in (1.0, -1.0):
            for x in _caps(outward=sign * UP, tip_inset=0.1).values():
                np.testing.assert_allclose(x[:, :, 0, 2], sign * 0.1, atol=1e-15)

    def test_wall_area_is_the_section_area(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        x, y = cs[:, 0], cs[:, 1]
        section = 0.5 * abs(float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))
        area = 0.0
        for c in _caps(cs, le_up, le_lo).values():
            f = c[:, :, 0]
            a, b, cc, d = f[:-1, :-1], f[1:, :-1], f[1:, 1:], f[:-1, 1:]
            area += 0.5 * np.linalg.norm(np.cross(cc - a, d - b), axis=-1).sum()
        assert abs(area - section) < 1e-12 * section


class TestLayers:
    def test_core_lines_along_the_normal(self):
        core = _caps()["core"]
        np.testing.assert_allclose(
            core[:, :, :, 2] - core[:, :, :1, 2], np.broadcast_to(DIST, core.shape[:3])
        )
        np.testing.assert_array_equal(
            core[:, :, :, :2], np.broadcast_to(core[:, :, :1, :2], core[..., :2].shape)
        )

    @pytest.mark.parametrize("flare", [0.0, 30.0, 45.0])
    def test_flare_at_the_edge(self, flare):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        edge = {p.tobytes() for p in cs}
        for x in _caps(cs, le_up, le_lo, flare_deg=flare).values():
            w, nxt = x[:, :, 0].reshape(-1, 3), x[:, :, 1].reshape(-1, 3)
            for a, b in zip(w, nxt, strict=True):
                if a.tobytes() in edge:
                    d = (b - a) / np.linalg.norm(b - a)
                    assert abs(np.degrees(np.arccos(d[2])) - flare) < 1e-9

    def test_section_normal_side(self):
        cs, _, _ = _naca0012_cross_section(41)
        np.testing.assert_allclose(section_normal(cs, UP), UP, atol=1e-12)
        np.testing.assert_allclose(section_normal(cs, -UP), -UP, atol=1e-12)


class TestWiring:
    def _seams(self, blocks):
        return [(b, f) for b in blocks for f in b.faces() if f.bc is FaceBC.SHARED]

    def test_sixteen_shared_faces(self):
        blocks = _hex(_caps())
        wire_cap_blocks(blocks)
        assert len(self._seams(blocks)) == 16
        for b in blocks:
            assert b.face("k_lo").bc is FaceBC.WALL
            assert b.face("k_hi").bc is FaceBC.FREE

    def test_seams_match_through_axis_maps(self):
        blocks = _hex(_caps())
        wire_cap_blocks(blocks)
        owner = {id(f): b for b in blocks for f in b.faces()}
        along = {"i": ("j", "k"), "j": ("i", "k"), "k": ("i", "j")}
        for b, f in self._seams(blocks):
            partner = owner[id(f.partner)]
            assert f.partner.partner is f
            theirs = f.partner.nodes(partner)
            (a0, flip0), (_, flip1) = f.axis_map.along
            t = theirs if along[f.partner.name[0]].index(a0) == 0 else theirs.transpose(1, 0, 2)
            t = t[::-1] if flip0 else t
            t = t[:, ::-1] if flip1 else t
            np.testing.assert_array_equal(f.nodes(b), t)

    def test_ring_edges_stay_free(self):
        blocks = _hex(_caps())
        wire_cap_blocks(blocks)
        free = sum(f.bc is FaceBC.FREE for b in blocks for f in b.faces())
        assert free == 5 + 4  # every k_hi, and the ring's outer side


class TestValidation:
    def test_sharp_te_rejected(self):
        cs, le_up, le_lo = _sharp_te_cross_section(41)
        with pytest.raises(ParameterError, match="sharp"):
            _caps(cs, le_up, le_lo)

    def test_asymmetric_cuts_rejected(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="symmetric"):
            _caps(cs, le_up, le_lo + 1)

    def test_cuts_must_straddle_the_le(self):
        cs, _, _ = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="LE cuts"):
            _caps(cs, 0, cs.shape[0] - 1)

    @pytest.mark.parametrize(
        "kw", [{"core_frac": 0.0}, {"core_frac": 1.0}, {"flare_deg": 90.0}, {"flare_deg": -1.0}]
    )
    def test_bad_parameters(self, kw):
        with pytest.raises(ParameterError):
            _caps(**kw)

    @pytest.mark.parametrize("dist", [[0.0], [0.1, 0.2], [0.0, 0.2, 0.1]])
    def test_bad_distances(self, dist):
        with pytest.raises(ParameterError, match="k_distances"):
            _caps(k_distances=np.array(dist))


class TestCheckCapBlocks:
    def test_valid_cap_passes(self):
        check_cap_blocks(_hex(_caps()))

    def test_folded_block_refused(self):
        arrays = _caps()
        arrays["le"] = arrays["le"][::-1].copy()  # mirrored: every cell negative
        with pytest.raises(ParameterError, match="tip cap block le"):
            check_cap_blocks(_hex(arrays))

    def test_block_without_cells_refused(self):
        arrays = _caps()
        arrays["te"] = arrays["te"][:1].copy()
        with pytest.raises(ParameterError, match="no cells"):
            check_cap_blocks(_hex(arrays))


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
