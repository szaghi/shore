"""Tests for shore.primitive.HGridRectMesh — rectangular-duct H-grid."""

from __future__ import annotations

import pytest

from shore.exceptions import ParameterError
from shore.mesh.face import FaceBC
from shore.primitive.h_grid_rect import HGridRectMesh


class TestHGridRectShape:
    def test_default_shape(self):
        hg = HGridRectMesh.build(length=1.0, width=0.1, height=0.2, ni=20, nj=8, nk=6)
        assert hg.blocks[0].nodes.shape == (20, 8, 6, 3)

    def test_streamwise_x_extent(self):
        hg = HGridRectMesh.build(
            length=1.5, width=0.2, height=0.3, ni=10, nj=5, nk=4, streamwise_axis="x"
        )
        b = hg.blocks[0]
        assert abs(b.nodes[..., 0].max() - 1.5) < 1e-12
        assert abs(b.nodes[..., 0].min()) < 1e-12

    def test_streamwise_y_extent(self):
        hg = HGridRectMesh.build(
            length=2.0, width=0.1, height=0.2, ni=10, nj=5, nk=4, streamwise_axis="y"
        )
        b = hg.blocks[0]
        assert abs(b.nodes[..., 1].max() - 2.0) < 1e-12

    def test_streamwise_z_extent(self):
        hg = HGridRectMesh.build(
            length=3.0, width=0.1, height=0.2, ni=10, nj=5, nk=4, streamwise_axis="z"
        )
        b = hg.blocks[0]
        assert abs(b.nodes[..., 2].max() - 3.0) < 1e-12

    def test_cross_section_extent_streamwise_x(self):
        # streamwise=x → j is along y (width), k is along z (height).
        hg = HGridRectMesh.build(
            length=1.0, width=0.4, height=0.7, ni=10, nj=5, nk=4, streamwise_axis="x"
        )
        b = hg.blocks[0]
        assert abs(b.nodes[..., 1].max() - 0.4) < 1e-12  # width
        assert abs(b.nodes[..., 2].max() - 0.7) < 1e-12  # height


class TestHGridRectBCs:
    def test_default_bcs(self):
        hg = HGridRectMesh.build(length=1.0, width=0.1, height=0.2, ni=10, nj=5, nk=4)
        b = hg.blocks[0]
        # Streamwise faces: FREE.
        assert b.face("i_lo").bc == FaceBC.FREE
        assert b.face("i_hi").bc == FaceBC.FREE
        # Duct walls: WALL.
        assert b.face("j_lo").bc == FaceBC.WALL
        assert b.face("j_hi").bc == FaceBC.WALL
        assert b.face("k_lo").bc == FaceBC.WALL
        assert b.face("k_hi").bc == FaceBC.WALL


class TestHGridRectErrors:
    def test_invalid_streamwise_axis_rejected(self):
        with pytest.raises(ParameterError, match="streamwise_axis"):
            HGridRectMesh.build(
                length=1.0, width=0.1, height=0.2, ni=10, nj=5, nk=4, streamwise_axis="w"
            )

    def test_negative_dimension_rejected(self):
        with pytest.raises(ParameterError, match="length"):
            HGridRectMesh.build(length=-1.0, width=0.1, height=0.2, ni=10, nj=5, nk=4)
        with pytest.raises(ParameterError, match="width"):
            HGridRectMesh.build(length=1.0, width=0.0, height=0.2, ni=10, nj=5, nk=4)
        with pytest.raises(ParameterError, match="height"):
            HGridRectMesh.build(length=1.0, width=0.1, height=-0.5, ni=10, nj=5, nk=4)

    def test_too_few_nodes_rejected(self):
        # BoxMesh requires ni / nj / nk >= 2.
        with pytest.raises(ParameterError):
            HGridRectMesh.build(length=1.0, width=0.1, height=0.2, ni=1, nj=5, nk=4)


class TestHGridRectSpacing:
    def test_tanh2_clusters_near_walls(self):
        # tanh2 on j with strong beta should put the j-spacing tightly
        # at both j=0 and j=nj-1 (the duct walls).
        hg = HGridRectMesh.build(
            length=1.0,
            width=1.0,
            height=0.1,
            ni=4,
            nj=20,
            nk=4,
            j_law="tanh2",
            j_beta=4.0,
            streamwise_axis="x",
        )
        b = hg.blocks[0]
        # j-coords along the duct (any i, any k).
        y = b.nodes[0, :, 0, 1]
        first_dy = y[1] - y[0]
        last_dy = y[-1] - y[-2]
        mid_dy = y[len(y) // 2] - y[len(y) // 2 - 1]
        # First / last spacings should be much smaller than the middle.
        assert first_dy < 0.5 * mid_dy
        assert last_dy < 0.5 * mid_dy
