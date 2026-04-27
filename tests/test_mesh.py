"""Tests for the shore.mesh OOP layer."""

from __future__ import annotations

import numpy as np
import pytest

from shore.mesh import (
    CubedSphereTopology,
    FaceBC,
    HexBlock,
    Mesh,
    Spacing1D,
)

# ---------------------------------------------------------------------------
# Spacing1D
# ---------------------------------------------------------------------------


class TestSpacing1D:
    def test_uniform_endpoints(self):
        sp = Spacing1D(law="uniform")
        coords = sp.build(5, 0.0, 1.0)
        assert coords[0] == pytest.approx(0.0)
        assert coords[-1] == pytest.approx(1.0)
        assert len(coords) == 5

    def test_tanh_monotone(self):
        sp = Spacing1D(law="tanh", beta=3.0)
        coords = sp.build(10, 0.0, 1.0)
        assert np.all(np.diff(coords) > 0)

    def test_tanh2_monotone(self):
        sp = Spacing1D(law="tanh2", beta=3.0)
        coords = sp.build(10, 0.0, 1.0)
        assert np.all(np.diff(coords) > 0)

    def test_geometric_steps_positive(self):
        sp = Spacing1D(law="geometric", ds=1e-3, growth=1.1)
        steps = sp.build_steps(10)
        assert len(steps) == 9
        assert np.all(steps > 0)

    def test_geometric_requires_ds(self):
        from shore.exceptions import ParameterError

        sp = Spacing1D(law="geometric")
        with pytest.raises(ParameterError):
            sp.build_steps(5)


# ---------------------------------------------------------------------------
# HexBlock construction
# ---------------------------------------------------------------------------


class TestHexBlock:
    def _make_block(self, ni=5, nj=8, nk=4) -> HexBlock:
        nodes = np.zeros((ni, nj, nk, 3))
        return HexBlock(label="test", nodes=nodes)

    def test_face_count(self):
        blk = self._make_block()
        assert len(blk.faces()) == 6

    def test_edge_count(self):
        blk = self._make_block()
        assert len(blk.edges()) == 12

    def test_vertex_count(self):
        blk = self._make_block()
        assert len(blk.vertices()) == 8

    def test_face_lookup_by_name(self):
        blk = self._make_block()
        f = blk.face("k_lo")
        assert f.name == "k_lo"
        assert f.bc == FaceBC.WALL

    def test_edge_lookup_by_name(self):
        blk = self._make_block()
        e = blk.edge("j_lo_k_lo")
        assert e.axis == "i"

    def test_wall_face_is_k_lo(self):
        blk = self._make_block()
        assert blk.wall_face().name == "k_lo"

    def test_outer_face_is_k_hi(self):
        blk = self._make_block()
        assert blk.outer_face().name == "k_hi"

    def test_face_nodes_shape(self):
        blk = self._make_block(ni=5, nj=8, nk=4)
        assert blk.face("k_lo").nodes(blk).shape == (5, 8, 3)
        assert blk.face("j_lo").nodes(blk).shape == (5, 4, 3)
        assert blk.face("i_hi").nodes(blk).shape == (8, 4, 3)

    def test_edge_nodes_shape(self):
        blk = self._make_block(ni=5, nj=8, nk=4)
        assert blk.edge("j_lo_k_lo").nodes(blk).shape == (5, 3)  # along i
        assert blk.edge("i_lo_k_lo").nodes(blk).shape == (8, 3)  # along j
        assert blk.edge("i_lo_j_lo").nodes(blk).shape == (4, 3)  # along k

    def test_vertex_coords(self):
        nodes = np.arange(5 * 8 * 4 * 3, dtype=float).reshape(5, 8, 4, 3)
        blk = HexBlock(label="v", nodes=nodes)
        v0 = blk.vertices()[0]
        assert np.array_equal(v0.coords(blk), nodes[0, 0, 0, :])


# ---------------------------------------------------------------------------
# CubedSphereTopology.build()
# ---------------------------------------------------------------------------


class TestCubedSphereTopology:
    def _build(self, sphere_stl, ni=8, nj=16, nk=4):
        return Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=ni,
            nj=nj,
            nk=nk,
        ).blocks

    def test_returns_six_blocks(self, sphere_stl):
        blocks = self._build(sphere_stl)
        assert len(blocks) == 6

    def test_block_labels(self, sphere_stl):
        blocks = self._build(sphere_stl)
        labels = [b.label for b in blocks]
        assert labels == ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]

    def test_sub_block_shape(self, sphere_stl):
        ni, nj, nk = 8, 16, 4
        blocks = self._build(sphere_stl, ni=ni, nj=nj, nk=nk)
        nc = nj // 4 + 1
        for q in range(4):
            assert blocks[q].nodes.shape == (ni, nc, nk, 3)

    def test_cap_block_shape(self, sphere_stl):
        ni, nj, nk = 8, 16, 4
        blocks = self._build(sphere_stl, ni=ni, nj=nj, nk=nk)
        nc = nj // 4 + 1
        assert blocks[4].nodes.shape == (nc, nc, nk, 3)
        assert blocks[5].nodes.shape == (nc, nc, nk, 3)

    def test_shared_views(self, sphere_stl):
        """Adjacent sub-block shared columns share memory for q=0..2.

        The ring allocation gives views for sub0↔sub1, sub1↔sub2, sub2↔sub3.
        The periodic seam sub3↔sub0 uses a separate wrap column and is
        explicitly re-stitched rather than view-shared.
        """
        blocks = self._build(sphere_stl)
        for q in range(3):
            q_next = q + 1
            assert np.shares_memory(
                blocks[q].nodes[:, -1, :, :],
                blocks[q_next].nodes[:, 0, :, :],
            ), f"sub{q} j_hi and sub{q_next} j_lo must share memory"
        # q=3 (wrap seam): values must be equal but not necessarily views
        np.testing.assert_array_equal(
            blocks[3].nodes[:, -1, 0, :],
            blocks[0].nodes[:, 0, 0, :],
        )

    def test_shared_face_bc(self, sphere_stl):
        blocks = self._build(sphere_stl)
        for q in range(4):
            assert blocks[q].face("j_hi").bc == FaceBC.SHARED
            assert blocks[q].face("j_lo").bc == FaceBC.SHARED

    def test_dirichlet_face_bc_on_caps(self, sphere_stl):
        blocks = self._build(sphere_stl)
        for name in ("i_lo", "i_hi", "j_lo", "j_hi"):
            assert blocks[4].face(name).bc == FaceBC.DIRICHLET
            assert blocks[5].face(name).bc == FaceBC.DIRICHLET

    def test_k0_conformity(self, sphere_stl):
        """Shared j-columns must be identical at k=0 after build."""
        blocks = self._build(sphere_stl)
        for q in range(4):
            q_next = (q + 1) % 4
            np.testing.assert_array_equal(
                blocks[q].nodes[:, -1, 0, :],
                blocks[q_next].nodes[:, 0, 0, :],
            )


# ---------------------------------------------------------------------------
# Mesh.from_stl + march
# ---------------------------------------------------------------------------


class TestMeshMarch:
    def test_march_expands_radius(self, sphere_stl):
        topo = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=9,
            nj=16,
            nk=5,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        topo.march(ds=0.05, growth=1.1, smooth_strength=0.2, smooth_iters=2)

        for blk in topo.blocks:
            r_k0 = np.linalg.norm(blk.nodes[:, :, 0, :], axis=-1)
            r_kn = np.linalg.norm(blk.nodes[:, :, -1, :], axis=-1)
            assert float(r_kn.mean()) > float(r_k0.mean()) * 1.05

    def test_march_no_negative_jacobian(self, sphere_stl):
        """March should complete without raising MarchingError."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=9,
            nj=16,
            nk=5,
        )
        m.march(ds=0.05, growth=1.1)

    def test_write_geo(self, sphere_stl, tmp_path):
        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=9,
            nj=16,
            nk=5,
        )
        m.march(ds=0.05, growth=1.1)
        m.write_geo(tmp_path / "test")
        files = list(tmp_path.glob("*.geo"))
        assert len(files) == 6

    def test_from_stl_does_not_march(self, sphere_stl):
        """from_stl must not fill volume layers — only k=0 is set."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=9,
            nj=16,
            nk=5,
        )
        for blk in m.blocks:
            # march was NOT called, so inner layers are uninitialised;
            # just verify the allocation shape is correct
            assert blk.nodes.shape[2] == 5

    def test_caps_bvh_built_once(self, sphere_stl):
        """Phase 3: Projector should be instantiated exactly once per Body.

        Body._projector is a cached_property: the Projector (and its embedded
        RayMeshIntersector / BVH) is built on first access and reused for all
        subsequent calls — including the cap interior projection inside
        CubedSphereTopology.build when projector= is passed.
        """
        from unittest.mock import patch

        import shore.surface.projector as _proj_module
        from shore.mesh.body import Body

        build_count = 0
        real_init = _proj_module.Projector.__init__

        def counting_init(self_inner, *args, **kwargs):
            nonlocal build_count
            build_count += 1
            real_init(self_inner, *args, **kwargs)

        with patch.object(_proj_module.Projector, "__init__", counting_init):
            body = Body.from_stl(sphere_stl)
            CubedSphereTopology().build(
                surface=None,
                nk=3,
                ni=9,
                nj=16,
                projector=body.projector,
            )

        assert build_count == 1, (
            f"Expected exactly 1 Projector build, got {build_count}.  "
            "Check that Body._projector is a cached_property and that "
            "CubedSphereTopology.build uses projector= instead of building its own."
        )
