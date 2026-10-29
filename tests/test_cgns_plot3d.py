"""Unit tests of the CGNS and Plot3D writers (ADR 0035)."""

from __future__ import annotations

import sys
import warnings

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.io.cgns import _inverse, transform_vector, write_cgns
from shore.io.plot3d import write_plot3d


def _flat_caps():
    from shore.project import Project

    p = Project()
    p.build_primitive("flat_caps", r_out=3.0, z0=-1.0, z1=4.0, ni=5, nc=6, nk=5, inner_frac=0.4)
    return p.mesh.blocks


class TestTransform:
    def test_same_axes_across_lo_hi(self):
        assert transform_vector("i_hi", "i_lo", (("j", False), ("k", False))) == (1, 2, 3)

    def test_wake_cut_both_lo_reversed(self):
        # C-grid wake cut: upper.k_lo <-> lower.k_lo, i reversed
        assert transform_vector("k_lo", "k_lo", (("i", True), ("j", False))) == (-1, 2, -3)

    def test_cross_axis(self):
        # an i_lo face whose j runs along the partner's i, onto the partner's j_lo
        assert transform_vector("i_lo", "j_lo", (("i", False), ("k", False))) == (-2, 1, 3)

    @pytest.mark.parametrize("t", [(1, 2, 3), (-1, 2, -3), (-2, 1, 3), (3, -1, 2)])
    def test_inverse(self, t):
        inv = _inverse(t)
        assert _inverse(inv) == t
        for a, v in enumerate(t):  # a -> |v|-1 -> a
            assert abs(inv[abs(v) - 1]) == a + 1


class TestCgnsErrors:
    def test_needs_h5py(self, tmp_path, monkeypatch):
        monkeypatch.setitem(sys.modules, "h5py", None)
        with pytest.raises(ParameterError, match=r"\[cgns\]"):
            write_cgns({"b": _flat_caps()}, tmp_path / "x.cgns")

    def test_connection_outside_base(self, tmp_path):
        blocks = _flat_caps()
        with pytest.raises(ParameterError, match="outside"):
            write_cgns({"a": blocks[:1], "b": blocks[1:]}, tmp_path / "x.cgns")

    def test_repeated_base_or_zone_names(self, tmp_path):
        (block,) = [b for b in _flat_caps() if b.label == "center"]
        with pytest.raises(ParameterError):
            write_cgns({"b": [block, block]}, tmp_path / "x.cgns")


class TestPlot3dLayout:
    def test_records(self, tmp_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks = _flat_caps()[:2]
        write_plot3d(blocks, tmp_path / "g.xyz")
        raw = (tmp_path / "g.xyz").read_bytes()

        def record(offset):
            n = int(np.frombuffer(raw, "<i4", 1, offset)[0])
            assert int(np.frombuffer(raw, "<i4", 1, offset + 4 + n)[0]) == n
            return raw[offset + 4 : offset + 4 + n], offset + 8 + n

        head, at = record(0)
        assert np.frombuffer(head, "<i4").tolist() == [2]
        dims, at = record(at)
        assert np.frombuffer(dims, "<i4").reshape(2, 3).tolist() == [
            list(b.nodes.shape[:3]) for b in blocks
        ]
        for b in blocks:
            data, at = record(at)
            xyz = np.frombuffer(data, "<f8").reshape(3, -1)
            for c in range(3):
                assert np.array_equal(xyz[c], b.nodes[..., c].ravel(order="F"))
        assert at == len(raw)

    def test_no_blocks(self, tmp_path):
        with pytest.raises(ParameterError):
            write_plot3d([], tmp_path / "g.xyz")
