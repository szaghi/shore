"""Tests for CGridTipCap topology + Mesh.from_array constructor.

CGridTipCap is a Topology subclass that builds the tip-cap component
of a C-grid wing — a 3 (sharp) or 4 (blunt) block H-style fill of the
wing-tip airfoil cross-section, designed to compose against the C-grid
body wrap as a chimera component.

These tests pin down:
- Block count and labels for sharp / blunt.
- Internal SHARED seams between cap sub-blocks (LE-cut, chord-cut, +
  TE-cut for blunt).
- Default BCs (k_lo = WALL outboard wing-tip, k_hi = FREE chimera fringe,
  exterior i / j faces FREE).
- Mesh.from_array constructor wiring + no-op march on CGridTipCap.
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh.face import FaceBC
from shore.mesh.mesh import Mesh
from shore.mesh.topology import CGridTipCap


def _naca0012_cross_section(ni: int = 41) -> tuple[np.ndarray, int, int]:
    """Synthetic NACA 0012-style cross section, post-extension.

    Sweeps upper-TE-cut → suction → LE → pressure → lower-TE-cut.
    Returns ``(cross_section, le_cut_upper_idx, le_cut_lower_idx)``.
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
    upper_x = chord[chord <= 0.95]
    upper_y = yt[chord <= 0.95]
    upper = np.column_stack([upper_x[::-1], +upper_y[::-1]])
    lower = np.column_stack([upper_x[1:], -upper_y[1:]])
    pts2d = np.vstack([upper, lower])
    pts3d = np.column_stack([pts2d, np.zeros(len(pts2d))])
    le_idx = int(np.argmin(pts3d[:, 0]))
    le_cut_upper = next(i for i in range(le_idx, -1, -1) if pts3d[i, 0] >= 0.05)
    le_cut_lower = next(i for i in range(le_idx, len(pts3d)) if pts3d[i, 0] >= 0.05)
    return pts3d, le_cut_upper, le_cut_lower


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


class TestCGridTipCapConstruction:
    def test_returns_five_blocks_default(self):
        """Default 5-block H-fill: le, upper, lower, te_upper, te_lower."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        m = Mesh.from_array(surface=cs, topology=topo, nk=8)
        assert len(m.blocks) == 5
        labels = sorted(b.label for b in m.blocks)
        assert labels == [
            "cap_le",
            "cap_lower",
            "cap_te_lower",
            "cap_te_upper",
            "cap_upper",
        ]

    def test_returns_six_blocks_with_body_wake_base(self):
        """5 cap blocks + 1 cap_wake_base when the body wrap is blunt."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        wb = np.array(
            [
                [[0.95, +0.005, 0.0], [0.95, 0.0, 0.0], [0.95, -0.005, 0.0]],
                [[1.20, +0.005, 0.0], [1.20, 0.0, 0.0], [1.20, -0.005, 0.0]],
            ]
        )  # (ni_wake=2, nk_w=3, 3)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=True,
            wake_base_cross_section=wb,
        )
        m = Mesh.from_array(surface=cs, topology=topo, nk=8)
        assert len(m.blocks) == 6
        labels = sorted(b.label for b in m.blocks)
        assert labels == [
            "cap_le",
            "cap_lower",
            "cap_te_lower",
            "cap_te_upper",
            "cap_upper",
            "cap_wake_base",
        ]

    def test_block_shapes_sharp(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=12,
            is_blunt=False,
        )
        m = Mesh.from_array(surface=cs, topology=topo, nk=8)
        for blk in m.blocks:
            assert blk.nodes.ndim == 4
            assert blk.nodes.shape[2] == 8  # nk_cap = nk
            assert blk.nodes.shape[1] == 12  # nj_cap = n_cap_j
            assert blk.nodes.shape[3] == 3

    def test_nk_cap_equals_nk(self):
        """When using Mesh.from_array, nk parameter sets the cap's k extent."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        m = Mesh.from_array(surface=cs, topology=topo, nk=15)
        for blk in m.blocks:
            assert blk.nodes.shape[2] == 15


# ---------------------------------------------------------------------------
# Default BCs
# ---------------------------------------------------------------------------


class TestCGridTipCapBCs:
    def _build(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        return Mesh.from_array(surface=cs, topology=topo, nk=6)

    def test_k_lo_is_wall_on_every_block(self):
        """k_lo is the outboard wing-tip wall — WALL on every cap block."""
        m = self._build()
        for blk in m.blocks:
            assert blk.face("k_lo").bc == FaceBC.WALL

    def test_k_hi_is_free_on_every_block(self):
        """k_hi is the inboard chimera-overlap plane — FREE (chimera fringe receiver)."""
        m = self._build()
        for blk in m.blocks:
            assert blk.face("k_hi").bc == FaceBC.FREE

    def test_j_hi_is_wall_on_chord_band_blocks(self):
        """cap_upper.j_hi and cap_lower.j_hi are airfoil walls (suction / pressure surfaces)."""
        m = self._build()
        labels = {b.label: b for b in m.blocks}
        assert labels["cap_upper"].face("j_hi").bc == FaceBC.WALL
        assert labels["cap_lower"].face("j_hi").bc == FaceBC.WALL

    def test_j_hi_is_wall_on_le_block(self):
        """cap_le.j_hi wraps the LE airfoil profile — WALL."""
        m = self._build()
        labels = {b.label: b for b in m.blocks}
        assert labels["cap_le"].face("j_hi").bc == FaceBC.WALL


# ---------------------------------------------------------------------------
# Internal SHARED seams between cap sub-blocks (sharp)
# ---------------------------------------------------------------------------


class TestCGridTipCapSharedSeams:
    def _build(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        return Mesh.from_array(surface=cs, topology=topo, nk=6)

    def test_le_cut_upper_seam(self):
        """cap_le.i_lo ↔ cap_upper.i_lo — LE-cut on the suction side."""
        m = self._build()
        labels = {b.label: b for b in m.blocks}
        face_le = labels["cap_le"].face("i_lo")
        face_up = labels["cap_upper"].face("i_lo")
        assert face_le.bc == FaceBC.SHARED
        assert face_up.bc == FaceBC.SHARED
        assert face_le.partner is face_up
        assert face_up.partner is face_le

    def test_le_cut_lower_seam(self):
        """cap_le.i_hi ↔ cap_lower.i_lo — LE-cut on the pressure side."""
        m = self._build()
        labels = {b.label: b for b in m.blocks}
        face_le = labels["cap_le"].face("i_hi")
        face_lo = labels["cap_lower"].face("i_lo")
        assert face_le.bc == FaceBC.SHARED
        assert face_lo.bc == FaceBC.SHARED
        assert face_le.partner is face_lo
        assert face_lo.partner is face_le

    def test_chord_cut_seam(self):
        """cap_upper.j_lo ↔ cap_lower.j_lo — chord cut at y=0."""
        m = self._build()
        labels = {b.label: b for b in m.blocks}
        face_up = labels["cap_upper"].face("j_lo")
        face_lo = labels["cap_lower"].face("j_lo")
        assert face_up.bc == FaceBC.SHARED
        assert face_lo.bc == FaceBC.SHARED
        assert face_up.partner is face_lo
        assert face_lo.partner is face_up

    def test_seam_nodes_coincide_bit_exact(self):
        """SHARED seam nodes coincide bit-exactly between partner blocks.

        The cap-fill builder allocates each sub-block in its own
        contiguous array (no shared backing buffer), so SHARED seams
        rely on **bit-exact node coincidence at construction time**
        rather than numpy view aliasing.  Mutations don't propagate;
        equality holds initially.
        """
        m = self._build()
        labels = {b.label: b for b in m.blocks}
        # Chord cut: cap_upper.j_lo == cap_lower.j_lo.
        np.testing.assert_array_equal(
            labels["cap_upper"].nodes[:, 0, :, :],
            labels["cap_lower"].nodes[:, 0, :, :],
        )
        # LE-cut suction: cap_le.i_lo == cap_upper.i_lo.
        np.testing.assert_array_equal(
            labels["cap_le"].nodes[0, :, :, :],
            labels["cap_upper"].nodes[0, :, :, :],
        )
        # LE-cut pressure: cap_le.i_hi == cap_lower.i_lo.
        np.testing.assert_array_equal(
            labels["cap_le"].nodes[-1, :, :, :],
            labels["cap_lower"].nodes[0, :, :, :],
        )


# ---------------------------------------------------------------------------
# No-op march
# ---------------------------------------------------------------------------


class TestCGridTipCapNoOpMarch:
    def test_march_does_not_modify_blocks(self):
        """Mesh.march() on a CGridTipCap mesh is a no-op (cap is volumetric at build)."""
        cs, le_up, le_lo = _naca0012_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        m = Mesh.from_array(surface=cs, topology=topo, nk=6)
        # Snapshot the pre-march block nodes.
        before = [blk.nodes.copy() for blk in m.blocks]
        m.march()  # no-op for CGridTipCap
        for blk, snap in zip(m.blocks, before, strict=True):
            np.testing.assert_array_equal(blk.nodes, snap)


# ---------------------------------------------------------------------------
# Mesh.from_array
# ---------------------------------------------------------------------------


class TestMeshFromArray:
    def test_body_is_none(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        m = Mesh.from_array(surface=cs, topology=topo, nk=6)
        assert m.body is None

    def test_topology_attached(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        m = Mesh.from_array(surface=cs, topology=topo, nk=6)
        assert m.topology is topo


# ---------------------------------------------------------------------------
# Parameter validation
# ---------------------------------------------------------------------------


class TestCGridTipCapValidation:
    def test_nk_too_small(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        with pytest.raises(ParameterError):
            Mesh.from_array(surface=cs, topology=topo, nk=1)

    def test_blunt_requires_wake_base(self):
        _cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="wake_base"):
            CGridTipCap(
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                span_dir=(0.0, 0.0, -1.0),
                tip_inset=0.05,
                n_cap_j=10,
                is_blunt=True,
                wake_base_cross_section=None,
            )

    def test_invalid_span_dir(self):
        _cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match="span_dir"):
            CGridTipCap(
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                span_dir=(0.0, 0.0, 2.0),  # not unit
                tip_inset=0.05,
                n_cap_j=10,
                is_blunt=False,
            )


# ---------------------------------------------------------------------------
# Auto-round on sharp TE
# ---------------------------------------------------------------------------


def _sharp_te_naca_cross_section(ni: int = 41) -> tuple[np.ndarray, int, int]:
    """A *sharp* NACA 0012 cross-section (i=0 and i=ni-1 coincide at TE).

    Distinct from the default _naca0012_cross_section fixture which is
    pre-truncated at 95% chord (and therefore already 'blunt').
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
    yt[-1] = 0.0  # force sharp TE
    upper = np.column_stack([chord[::-1], +yt[::-1]])
    lower = np.column_stack([chord[1:], -yt[1:]])
    pts2d = np.vstack([upper, lower])
    pts3d = np.column_stack([pts2d, np.zeros(len(pts2d))])
    le_idx = int(np.argmin(pts3d[:, 0]))
    le_cut_upper = next(i for i in range(le_idx, -1, -1) if pts3d[i, 0] >= 0.05)
    le_cut_lower = next(i for i in range(le_idx, len(pts3d)) if pts3d[i, 0] >= 0.05)
    return pts3d, le_cut_upper, le_cut_lower


class TestCGridTipCapAutoRound:
    def test_auto_round_default_thickness_frac(self, capsys):
        """Sharp input + default auto_round_te=True → cap_te_upper.i_hi non-degenerate."""
        cs, le_up, le_lo = _sharp_te_naca_cross_section(41)
        assert np.allclose(cs[0, :], cs[-1, :])  # sanity: sharp
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        m = Mesh.from_array(surface=cs, topology=topo, nk=6)
        # Verify that auto-round produced a non-degenerate cap_te_upper.i_hi.
        labels = {b.label: b for b in m.blocks}
        col = labels["cap_te_upper"].nodes[-1, :, 0, :]
        span = float(np.linalg.norm(col.max(axis=0) - col.min(axis=0)))
        assert span > 1e-4
        # The verbose stderr message was emitted.
        assert "auto-rounding" in capsys.readouterr().err

    def test_auto_round_explicit_thickness(self):
        """Explicit te_thickness parameter overrides the chord-frac default."""
        cs, le_up, le_lo = _sharp_te_naca_cross_section(41)
        te_user = 0.005
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
            te_thickness=te_user,
            te_thickness_frac=None,
        )
        m = Mesh.from_array(surface=cs, topology=topo, nk=6)
        labels = {b.label: b for b in m.blocks}
        # cap_te_upper at i=ni-1, j=nj-1 is the airfoil-suction TE-base
        # corner at the upper-side; cap_te_lower's mirror sits at -te/2.
        upper_corner = labels["cap_te_upper"].nodes[-1, -1, -1, :]
        lower_corner = labels["cap_te_lower"].nodes[-1, -1, -1, :]
        # The thickness gap between them = te_user (anti-symmetric split).
        gap = float(np.linalg.norm(upper_corner - lower_corner))
        np.testing.assert_allclose(gap, te_user, rtol=1e-9)

    def test_auto_round_disabled_rejects_sharp(self):
        cs, le_up, le_lo = _sharp_te_naca_cross_section(41)
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
            auto_round_te=False,
        )
        with pytest.raises(ParameterError, match="sharp TE"):
            Mesh.from_array(surface=cs, topology=topo, nk=6)

    def test_blunt_input_not_rounded(self, capsys):
        """A blunt input is used as-is — no auto-round message."""
        cs, le_up, le_lo = _naca0012_cross_section(41)  # already blunt-truncated
        topo = CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
        )
        Mesh.from_array(surface=cs, topology=topo, nk=6)
        # No auto-round message on blunt input.
        assert "auto-rounding" not in capsys.readouterr().err

    def test_te_thickness_and_frac_mutually_exclusive(self):
        _cs, le_up, le_lo = _sharp_te_naca_cross_section(41)
        with pytest.raises(ParameterError, match="not both"):
            CGridTipCap(
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                span_dir=(0.0, 0.0, -1.0),
                tip_inset=0.05,
                n_cap_j=10,
                is_blunt=False,
                te_thickness=0.001,
                te_thickness_frac=0.002,
            )
