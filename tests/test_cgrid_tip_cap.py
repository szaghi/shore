"""Tests for CGridTipCap (ADR 0059): the tip cap as its own component.

CGridTipCap builds the O-H cap of ``shore.topology.cgrid_cap`` from a
cross-section polyline: five blocks ``cap_core`` / ``cap_upper`` /
``cap_le`` / ``cap_lower`` / ``cap_te``, wall at ``k_lo``, SHARED between
cap blocks, FREE elsewhere; the march is a no-op.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh.face import FaceBC
from shore.mesh.mesh import Mesh
from shore.mesh.topology import CGridTipCap
from shore.quality import scaled_jacobian

LABELS = ["cap_core", "cap_upper", "cap_le", "cap_lower", "cap_te"]


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


def _topo(cs_fn=_naca0012_cross_section, **kw):
    cs, le_up, le_lo = cs_fn(41)
    kw.setdefault("span_dir", (0.0, 0.0, -1.0))
    kw.setdefault("tip_inset", 0.05)
    kw.setdefault("n_cap_j", 10)
    return cs, CGridTipCap(le_cut_upper_idx=le_up, le_cut_lower_idx=le_lo, **kw)


def _mesh(nk=6, **kw):
    cs, topo = _topo(**kw)
    return Mesh.from_array(surface=cs, topology=topo, nk=nk)


class TestConstruction:
    def test_five_labelled_blocks(self):
        assert [b.label for b in _mesh().blocks] == LABELS

    def test_valid_cells(self):
        for b in _mesh().blocks:
            assert float(scaled_jacobian(b.nodes).min()) > 0.2, b.label

    def test_nk(self):
        assert {b.nodes.shape[2] for b in _mesh(nk=9).blocks} == {9}

    def test_wall_offset_by_tip_inset(self):
        for b in _mesh().blocks:  # span_dir -z, inset 0.05: the wall at z = -0.05
            np.testing.assert_allclose(b.nodes[:, :, 0, 2], -0.05, atol=1e-15)
            assert np.all(b.nodes[:, :, 1:, 2] < -0.05)

    def test_k_distances(self):
        dist = np.array([0.0, 0.1, 0.3])
        core = _mesh(nk=3, k_distances=dist).blocks[0].nodes
        np.testing.assert_allclose(
            -0.05 - core[:, :, :, 2], np.broadcast_to(dist, core.shape[:3]), atol=1e-15
        )

    def test_k_distances_must_match_nk(self):
        cs, topo = _topo(k_distances=np.array([0.0, 0.1, 0.3]))
        with pytest.raises(ParameterError, match="k_distances"):
            Mesh.from_array(surface=cs, topology=topo, nk=4)


class TestBCs:
    def test_wall_shared_free(self):
        blocks = _mesh().blocks
        for b in blocks:
            assert b.face("k_lo").bc is FaceBC.WALL
            assert b.face("k_hi").bc is FaceBC.FREE
        shared = sum(f.bc is FaceBC.SHARED for b in blocks for f in b.faces())
        assert shared == 16
        walls = sum(f.bc is FaceBC.WALL for b in blocks for f in b.faces())
        assert walls == 5


class TestMarch:
    def test_march_does_not_modify_blocks(self):
        m = _mesh()
        before = [blk.nodes.copy() for blk in m.blocks]
        m.march()
        for blk, snap in zip(m.blocks, before, strict=True):
            np.testing.assert_array_equal(blk.nodes, snap)

    def test_mesh_from_array(self):
        cs, topo = _topo()
        m = Mesh.from_array(surface=cs, topology=topo, nk=6)
        assert m.body is None
        assert m.topology is topo


class TestValidation:
    def test_nk_too_small(self):
        with pytest.raises(ParameterError):
            _mesh(nk=1)

    def test_invalid_span_dir(self):
        with pytest.raises(ParameterError, match="span_dir"):
            _topo(span_dir=(0.0, 0.0, 2.0))

    def test_negative_tip_inset(self):
        with pytest.raises(ParameterError, match="tip_inset"):
            _topo(tip_inset=-0.1)

    def test_old_te_arguments_warn(self):
        with pytest.warns(DeprecationWarning, match="ADR 0059"):
            _topo(is_blunt=True, wake_base_cross_section=np.zeros((3, 3, 3)))

    def test_te_thickness_and_frac_mutually_exclusive(self):
        with pytest.raises(ParameterError, match="not both"):
            _topo(te_thickness=0.01, te_thickness_frac=0.002)


class TestAutoRound:
    def test_sharp_section_rounded(self, capsys):
        m = _mesh(cs_fn=_sharp_te_naca_cross_section)
        assert "auto-rounding" in capsys.readouterr().err
        assert [b.label for b in m.blocks] == LABELS
        for b in m.blocks:
            assert float(scaled_jacobian(b.nodes).min()) > 0.2, b.label

    def test_explicit_thickness_is_the_base(self):
        m = _mesh(cs_fn=_sharp_te_naca_cross_section, te_thickness=0.005, te_thickness_frac=None)
        te = m.blocks[-1].nodes[:, :, 0]
        wall = te[:, -1] if te.shape[1] == 10 else te[-1]  # the base row, j = n_cap_j - 1
        np.testing.assert_allclose(np.ptp(wall[:, 1]), 0.005, rtol=1e-9)

    def test_disabled_rejects_sharp(self):
        with pytest.raises(ParameterError, match="sharp TE"):
            _mesh(cs_fn=_sharp_te_naca_cross_section, auto_round_te=False)

    def test_blunt_input_not_rounded(self, capsys):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            _mesh()
        assert "auto-rounding" not in capsys.readouterr().err
