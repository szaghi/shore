"""Tests for CGridTipCap: parameter validation and the build refusal.

The caps of this construction fold or have no cells (ADR 0057), so every
build is refused by ``check_cap_blocks`` (ADR 0058) until the redesign.
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
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


class TestCGridTipCapRefused:
    """The caps of this construction fold or have no cells (ADR 0057), so
    every build is refused (ADR 0058) and no mesh is made."""

    def _topo(self, cs, le_up, le_lo, **kw):
        return CGridTipCap(
            le_cut_upper_idx=le_up,
            le_cut_lower_idx=le_lo,
            span_dir=(0.0, 0.0, -1.0),
            tip_inset=0.05,
            n_cap_j=10,
            is_blunt=False,
            **kw,
        )

    def test_blunt_section_refused(self):
        cs, le_up, le_lo = _naca0012_cross_section(41)
        with pytest.raises(ParameterError, match=r"tip cap block cap_"):
            Mesh.from_array(surface=cs, topology=self._topo(cs, le_up, le_lo), nk=6)

    def test_sharp_section_refused_after_rounding(self, capsys):
        cs, le_up, le_lo = _sharp_te_naca_cross_section(41)
        with pytest.raises(ParameterError, match=r"tip cap block cap_"):
            Mesh.from_array(surface=cs, topology=self._topo(cs, le_up, le_lo), nk=6)
        assert "auto-rounding" in capsys.readouterr().err


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
