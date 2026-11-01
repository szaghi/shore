"""Tests for shore.primitive.HGridTubeMesh — the O-H circular-bore H-grid (ADR 0060)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh.face import FaceBC
from shore.primitive.h_grid_tube import TUBE_BLOCKS, HGridTubeMesh
from shore.quality import scaled_jacobian


def _tube(**kw):
    params = {"length": 1.0, "radius": 0.05, "n_axial": 10, "n_circ": 16, "n_radial": 6}
    return HGridTubeMesh.build(**{**params, **kw})


class TestHGridTubeShape:
    def test_labels(self):
        assert tuple(b.label for b in _tube().blocks) == TUBE_BLOCKS

    def test_shapes(self):
        blocks = _tube(n_axial=20, n_circ=16, n_radial=6).blocks
        for b in blocks[:4]:
            assert b.nodes.shape == (20, 5, 6, 3)
        assert blocks[4].nodes.shape == (20, 5, 5, 3)

    @pytest.mark.parametrize("axis, idx", [("x", 0), ("y", 1), ("z", 2)])
    def test_streamwise_axis(self, axis, idx):
        for b in _tube(length=2.0, streamwise_axis=axis).blocks:
            assert abs(b.nodes[..., idx].max() - 2.0) < 1e-12
            assert abs(b.nodes[..., idx].min()) < 1e-12

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_right_handed_valid(self, axis):
        for b in _tube(streamwise_axis=axis).blocks:
            assert float(scaled_jacobian(b.nodes).min()) > 0.7, b.label

    def test_outer_radius(self):
        for b in _tube(radius=0.07).blocks[:4]:
            outer = b.nodes[2, :, -1, :]
            np.testing.assert_allclose(np.hypot(outer[:, 1], outer[:, 2]), 0.07, rtol=1e-12)

    def test_core_size(self):
        core = _tube(radius=1.0, core_frac=0.4).blocks[4].nodes[0]
        corner = np.abs(core[[0, 0, -1, -1], [0, -1, 0, -1]][:, 1:])
        np.testing.assert_allclose(np.hypot(corner[:, 0], corner[:, 1]), 0.4, rtol=1e-12)

    def test_no_singular_axis(self):
        for b in _tube().blocks:
            edges = np.linalg.norm(np.diff(b.nodes, axis=2), axis=-1)
            assert edges.min() > 0.0, b.label


class TestHGridTubeBCs:
    def test_default_bcs(self):
        blocks = _tube().blocks
        for b in blocks:
            assert b.face("i_lo").bc == FaceBC.FREE
            assert b.face("i_hi").bc == FaceBC.FREE
            assert not b.periodic_i and not b.periodic_j
        for b in blocks[:4]:
            assert b.face("k_hi").bc == FaceBC.WALL
        walls = sum(f.bc == FaceBC.WALL for b in blocks for f in b.faces())
        assert walls == 4

    def test_sixteen_shared_faces(self):
        blocks = _tube().blocks
        shared = [f for b in blocks for f in b.faces() if f.bc == FaceBC.SHARED]
        assert len(shared) == 16
        assert all(f.partner.partner is f for f in shared)

    def test_label_deprecated(self):
        with pytest.warns(DeprecationWarning, match="ADR 0060"):
            _tube(label="duct")


class TestHGridTubeErrors:
    def test_invalid_streamwise_axis_rejected(self):
        with pytest.raises(ParameterError, match="streamwise_axis"):
            HGridTubeMesh.build(
                length=1.0,
                radius=0.05,
                n_axial=10,
                n_circ=8,
                n_radial=4,
                streamwise_axis="w",
            )

    def test_negative_length_rejected(self):
        with pytest.raises(ParameterError, match="length"):
            HGridTubeMesh.build(length=-1.0, radius=0.05, n_axial=10, n_circ=8, n_radial=4)

    def test_negative_radius_rejected(self):
        with pytest.raises(ParameterError, match="radius"):
            HGridTubeMesh.build(length=1.0, radius=0.0, n_axial=10, n_circ=8, n_radial=4)

    def test_too_few_axial_rejected(self):
        with pytest.raises(ParameterError, match="n_axial"):
            HGridTubeMesh.build(length=1.0, radius=0.05, n_axial=1, n_circ=8, n_radial=4)

    def test_too_few_circ_rejected(self):
        with pytest.raises(ParameterError, match="n_circ"):
            HGridTubeMesh.build(length=1.0, radius=0.05, n_axial=10, n_circ=10, n_radial=4)

    def test_too_few_radial_rejected(self):
        with pytest.raises(ParameterError, match="n_radial"):
            HGridTubeMesh.build(length=1.0, radius=0.05, n_axial=10, n_circ=8, n_radial=1)

    def test_core_frac_rejected(self):
        with pytest.raises(ParameterError, match="core_frac"):
            _tube(core_frac=1.0)


class TestHGridTubeSpacing:
    def test_tanh_clusters_radial_at_the_wall(self):
        b = _tube(radius=0.1, n_circ=8, n_radial=20, k_law="tanh", k_beta=4.0).blocks[0]
        line = b.nodes[0, 1]  # the ring's middle line, core -> wall
        dr = np.linalg.norm(np.diff(line, axis=0), axis=-1)
        assert dr[-1] < dr[0] * 0.5
