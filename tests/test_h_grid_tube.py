"""Tests for shore.primitive.HGridTubeMesh — circular-bore H-grid."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh.face import FaceBC
from shore.primitive.h_grid_tube import HGridTubeMesh


class TestHGridTubeShape:
    def test_default_shape(self):
        ht = HGridTubeMesh.build(length=1.0, radius=0.05, n_axial=20, n_circ=8, n_radial=6)
        assert ht.blocks[0].nodes.shape == (20, 8, 6, 3)

    def test_streamwise_z(self):
        ht = HGridTubeMesh.build(
            length=2.0,
            radius=0.05,
            n_axial=10,
            n_circ=8,
            n_radial=4,
            streamwise_axis="z",
        )
        b = ht.blocks[0]
        assert abs(b.nodes[..., 2].max() - 2.0) < 1e-12
        assert abs(b.nodes[..., 2].min()) < 1e-12

    def test_streamwise_y(self):
        ht = HGridTubeMesh.build(
            length=3.0,
            radius=0.05,
            n_axial=10,
            n_circ=8,
            n_radial=4,
            streamwise_axis="y",
        )
        b = ht.blocks[0]
        assert abs(b.nodes[..., 1].max() - 3.0) < 1e-12

    def test_outer_radius_correct(self):
        ht = HGridTubeMesh.build(
            length=1.0,
            radius=0.07,
            n_axial=4,
            n_circ=12,
            n_radial=5,
            streamwise_axis="x",
        )
        b = ht.blocks[0]
        # Cross-section is in (y, z); outer wall is at k=n_radial-1.
        outer = b.nodes[2, :, -1, :]
        r = np.sqrt(outer[:, 1] ** 2 + outer[:, 2] ** 2)
        assert np.allclose(r, 0.07, atol=1e-12)


class TestHGridTubeCenterline:
    def test_k_lo_collapses_to_point(self):
        # The k=0 row is the polar singular axis: all n_circ angular
        # nodes collapse to the same point at the duct centerline.
        ht = HGridTubeMesh.build(
            length=1.0,
            radius=0.05,
            n_axial=4,
            n_circ=12,
            n_radial=5,
            streamwise_axis="x",
        )
        b = ht.blocks[0]
        # At any axial station, k=0 nodes are bit-coincident.
        for i in range(b.nodes.shape[0]):
            row = b.nodes[i, :, 0, :]
            np.testing.assert_array_equal(row, np.broadcast_to(row[0], row.shape))


class TestHGridTubeBCs:
    def test_default_bcs(self):
        ht = HGridTubeMesh.build(length=1.0, radius=0.05, n_axial=10, n_circ=8, n_radial=4)
        b = ht.blocks[0]
        # Streamwise: FREE.
        assert b.face("i_lo").bc == FaceBC.FREE
        assert b.face("i_hi").bc == FaceBC.FREE
        # Walls: centerline + outer.
        assert b.face("k_lo").bc == FaceBC.WALL
        assert b.face("k_hi").bc == FaceBC.WALL

    def test_periodic_j(self):
        ht = HGridTubeMesh.build(length=1.0, radius=0.05, n_axial=10, n_circ=8, n_radial=4)
        b = ht.blocks[0]
        assert b.periodic_j is True
        assert b.periodic_i is False


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
            HGridTubeMesh.build(length=1.0, radius=0.05, n_axial=10, n_circ=3, n_radial=4)

    def test_too_few_radial_rejected(self):
        with pytest.raises(ParameterError, match="n_radial"):
            HGridTubeMesh.build(length=1.0, radius=0.05, n_axial=10, n_circ=8, n_radial=1)


class TestHGridTubeSpacing:
    def test_tanh_clusters_radial_near_wall(self):
        # tanh on k clusters near the centerline (k=0); we want
        # clustering near the outer wall (k=n_radial-1) instead, so
        # users would typically use tanh2 + appropriate beta, but
        # here we just sanity-check that tanh produces non-uniform
        # radial spacing.
        ht = HGridTubeMesh.build(
            length=1.0,
            radius=0.1,
            n_axial=4,
            n_circ=8,
            n_radial=20,
            k_law="tanh",
            k_beta=4.0,
        )
        b = ht.blocks[0]
        # Radial coords at one (i, j).
        r = np.sqrt(b.nodes[0, 0, :, 1] ** 2 + b.nodes[0, 0, :, 2] ** 2)
        # Spacings.
        dr = np.diff(r)
        # tanh-1-sided clusters near r=0; first dr should be smaller
        # than the last dr.
        assert dr[0] < dr[-1] * 0.5
