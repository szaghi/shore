"""Tests for block-to-geometry association (ADR 0003, roadmap phase 2)."""

from __future__ import annotations

import json
import warnings

import numpy as np
import pytest
import trimesh

from shore.exceptions import MeshInputError, ParameterError
from shore.geometry import BODY, EntityKind, EntityRef, StlGeometry
from shore.geometry.keys import EntityKey, match_keys
from shore.io.adjacency import (
    read_adjacency_json,
    wire_blocks_from_adjacency,
    write_adjacency_json,
)
from shore.mesh import Body, CGridTopology, CubedSphereTopology, Mesh
from shore.mesh.association import (
    Association,
    AssociationResolver,
    associate_wall,
    snap_block,
)
from shore.mesh.block import HexBlock
from shore.primitive.box import BoxMesh
from shore.split import BlockSplit, split_blocks

WALL_EDGES = {"j_lo_k_lo", "j_hi_k_lo", "i_lo_k_lo", "i_hi_k_lo"}


def _assert_wall_associated(block: HexBlock, assoc: Association) -> None:
    for f in block.faces():
        assert f.association == (assoc if f.name == "k_lo" else None), f.name
    for e in block.edges():
        assert e.association == (assoc if e.name in WALL_EDGES else None), e.name
    for v in block.vertices():
        assert v.association == (assoc if v.ijk[2] == 0 else None), v.ijk


@pytest.fixture(scope="module")
def sphere_mesh(sphere_stl) -> Mesh:
    return Mesh.from_file(sphere_stl, CubedSphereTopology(), ni=9, nj=16, nk=4)


# ---------------------------------------------------------------------------
# Value types
# ---------------------------------------------------------------------------


class TestAssociationType:
    def test_body(self):
        a = Association.body()
        assert a.geometry == "body"
        assert a.kind is EntityKind.BODY
        assert a.keys == ()

    def test_body_with_keys_rejected(self):
        key = StlGeometry(trimesh.creation.icosphere()).fingerprint(BODY)
        with pytest.raises(ParameterError):
            Association("body", EntityKind.BODY, (key,))

    def test_surface_without_keys_rejected(self):
        with pytest.raises(ParameterError):
            Association("body", EntityKind.SURFACE)

    def test_mixed_key_kinds_rejected(self):
        key = StlGeometry(trimesh.creation.icosphere()).fingerprint(BODY)
        with pytest.raises(ParameterError):
            Association("body", EntityKind.SURFACE, (key,))

    def test_json_round_trip(self):
        key = EntityKey(EntityKind.SURFACE, "Plane", (0.1, 0.2, 0.3), 0.5, (1.0, 0.0, 1.0))
        a = Association("wing", EntityKind.SURFACE, (key,))
        assert Association.from_json(json.loads(json.dumps(a.to_json()))) == a

    def test_malformed_json(self):
        with pytest.raises(MeshInputError):
            Association.from_json({"geometry": "body", "kind": "nonsense"})


class TestStlFingerprint:
    def test_stable_across_reload(self, sphere_stl):
        a = StlGeometry.from_file(sphere_stl).fingerprint(BODY)
        b = StlGeometry.from_file(sphere_stl).fingerprint(BODY)
        assert a == b
        assert match_keys([a], StlGeometry.from_file(sphere_stl)).count("same") == 1


# ---------------------------------------------------------------------------
# Templates associate their body wall
# ---------------------------------------------------------------------------


class TestWallAssociation:
    def test_cubed_sphere_walls(self, sphere_mesh):
        for block in sphere_mesh.blocks:
            _assert_wall_associated(block, Association.body())
            assert block.face("k_lo").body_projected
            assert not block.face("k_hi").body_projected

    def test_cgrid_wake_cut_not_associated(self, naca0012_stl):
        body = Body.from_stl(naca0012_stl, skip_hygiene=True)
        topo = CGridTopology()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks = topo.build(
                surface=None, nk=6, mesh=body.mesh, ni_body=40, ni_wake=10,
                n_stations=3, wake_length=5.0, blocks=3,
            )  # fmt: skip
        mesh = Mesh(body=body, blocks=blocks, topology=topo)
        by = {b.label: b for b in mesh.blocks}
        _assert_wall_associated(by["body"], Association.body())
        for wake in ("upper", "lower"):  # k_lo is the SHARED wake cut, not the body
            assert by[wake].face("k_lo").association is None

    def test_primitive_has_no_association(self):
        box = BoxMesh.build((0, 0, 0), (1, 1, 1), 3, 3, 3)
        for block in box.blocks:
            assert all(f.association is None for f in block.faces())
            assert not any(f.body_projected for f in block.faces())


# ---------------------------------------------------------------------------
# Split propagation
# ---------------------------------------------------------------------------


class TestSplitPropagation:
    def test_split_keeps_wall_and_leaves_cuts_free(self):
        block = HexBlock("b", np.random.default_rng(0).random((9, 9, 9, 3)), periodic_j=False)
        assoc = Association.body()
        associate_wall(block, assoc)
        chunks = split_blocks([block], [BlockSplit("b", "i", (4,)), BlockSplit("b", "k", (4,))])
        by = {c.label: c for c in chunks}
        assert set(by) == {"b__i0__k0", "b__i1__k0", "b__i0__k1", "b__i1__k1"}
        lo, hi = by["b__i0__k0"], by["b__i1__k0"]
        for c in (lo, hi):
            assert c.face("k_lo").association == assoc  # piece of the source wall
            assert c.face("k_hi").association is None  # internal cut
        # the i-cut runs across the wall: its k_lo edge lies on the source face
        assert lo.edge("i_hi_k_lo").association == assoc
        assert lo.vertices()[4].association == assoc  # (-1, 0, 0): on the source wall edge
        assert lo.face("i_hi").association is None  # internal cut face
        for c in (by["b__i0__k1"], by["b__i1__k1"]):  # away from the wall
            assert all(f.association is None for f in c.faces())
            assert all(v.association is None for v in c.vertices())

    def test_split_of_curve_associated_edge(self):
        block = HexBlock("b", np.random.default_rng(1).random((9, 5, 5, 3)), periodic_j=False)
        key = EntityKey(EntityKind.CURVE, "Line", (0.5, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0))
        curve = Association("cad", EntityKind.CURVE, (key,))
        block.edge("j_lo_k_lo").association = curve
        a, b = split_blocks([block], [BlockSplit("b", "i", (4,))])
        assert a.edge("j_lo_k_lo").association == curve
        assert b.edge("j_lo_k_lo").association == curve


# ---------------------------------------------------------------------------
# Adjacency sidecar (version 3)
# ---------------------------------------------------------------------------


class TestAdjacencyRoundTrip:
    def test_round_trip(self, sphere_mesh, tmp_path):
        path = write_adjacency_json(sphere_mesh.blocks, tmp_path / "adj")
        payload = read_adjacency_json(path)
        assert payload["version"] == 3
        assert payload["associations"] == [Association.body()]  # deduplicated
        fresh = [HexBlock(b.label, b.nodes.copy(), periodic_j=False) for b in sphere_mesh.blocks]
        wire_blocks_from_adjacency(fresh, payload)
        for old, new in zip(sphere_mesh.blocks, fresh, strict=True):
            for a, b in zip(old.faces(), new.faces(), strict=True):
                assert a.association == b.association
            for a, b in zip(old.edges(), new.edges(), strict=True):
                assert a.association == b.association
            for a, b in zip(old.vertices(), new.vertices(), strict=True):
                assert a.association == b.association

    def test_version_2_still_readable(self, tmp_path):
        box = BoxMesh.build((0, 0, 0), (1, 1, 1), 3, 3, 3)
        path = write_adjacency_json(box.blocks, tmp_path / "adj")
        data = json.loads(path.read_text())
        data["version"] = 2
        path.write_text(json.dumps(data))
        payload = read_adjacency_json(path)
        assert payload["associations"] == []

    def test_bad_association_index(self, sphere_mesh, tmp_path):
        path = write_adjacency_json(sphere_mesh.blocks, tmp_path / "adj")
        data = json.loads(path.read_text())
        data["blocks"][0]["associated"]["faces"]["k_lo"] = 7
        path.write_text(json.dumps(data))
        with pytest.raises(MeshInputError, match="out of range"):
            read_adjacency_json(path)


# ---------------------------------------------------------------------------
# Projection follows association
# ---------------------------------------------------------------------------


class TestSnapStl:
    def test_snap_returns_wall_to_body(self, sphere_stl):
        mesh = Mesh.from_file(sphere_stl, CubedSphereTopology(), ni=9, nj=16, nk=4)
        k0 = [b.nodes[:, :, 0].copy() for b in mesh.blocks]
        for b in mesh.blocks:
            b.nodes[:, :, 0] *= 1.05  # push the wall off the body
        mesh.snap()
        for b, ref in zip(mesh.blocks, k0, strict=True):
            np.testing.assert_allclose(b.nodes[:, :, 0], ref, atol=0.02)
            # closest point on the facets: on the body surface
            d = mesh.body.geometry.project(b.nodes[:, :, 0]).distances
            assert d.max() < 1e-12

    def test_unknown_geometry(self):
        block = HexBlock("b", np.zeros((3, 3, 3, 3)), periodic_j=False)
        block.face("k_lo").association = Association.body("nope")
        with pytest.raises(MeshInputError, match="nope"):
            snap_block(block, AssociationResolver({}))


@pytest.fixture(scope="module")
def bd():
    return pytest.importorskip("build123d")


def _closest_entity(g, kind: EntityKind, point) -> EntityRef:
    p = np.asarray(point, dtype=float)[None]
    return min(g.entities(kind), key=lambda r: g.project(p, r, exact=True).distances[0])


class TestSnapCad:
    def test_vertex_edge_face_priority_on_box(self, bd):
        from shore.geometry.occt import OcctGeometry

        g = OcctGeometry(bd.Box(1, 1, 1).wrapped)
        # a block whose k_lo face lies (roughly) on the box's +x face
        y, z = np.meshgrid(np.linspace(-0.5, 0.5, 5), np.linspace(-0.5, 0.5, 5), indexing="ij")
        nodes = np.empty((5, 5, 3, 3))
        for k in range(3):
            nodes[:, :, k] = np.stack([np.full_like(y, 0.5 + 0.2 * k), y, z], axis=-1)
        nodes += 0.01 * np.random.default_rng(2).standard_normal(nodes.shape)
        block = HexBlock("b", nodes, periodic_j=False)
        face = _closest_entity(g, EntityKind.SURFACE, (2.0, 0.0, 0.0))
        block.face("k_lo").association = Association.of("cad", g, [face])
        for name, mid in {
            "j_lo_k_lo": (0.5, 0.0, -0.5),
            "j_hi_k_lo": (0.5, 0.0, 0.5),
            "i_lo_k_lo": (0.5, -0.5, 0.0),
            "i_hi_k_lo": (0.5, 0.5, 0.0),
        }.items():
            edge = _closest_entity(g, EntityKind.CURVE, mid)
            block.edge(name).association = Association.of("cad", g, [edge])
        for v in block.vertices():
            if v.ijk[2] == 0:
                corner = (0.5, -0.5 if v.ijk[0] == 0 else 0.5, -0.5 if v.ijk[1] == 0 else 0.5)
                ref = _closest_entity(g, EntityKind.POINT, corner)
                v.association = Association.of("cad", g, [ref])
        snap_block(block, AssociationResolver({"cad": g}, exact=True))
        wall = block.nodes[:, :, 0]
        np.testing.assert_allclose(wall[..., 0], 0.5, atol=1e-12)  # on the +x face
        np.testing.assert_allclose(wall[:, 0, 2], -0.5, atol=1e-12)  # edges on box edges
        np.testing.assert_allclose(wall[:, -1, 2], 0.5, atol=1e-12)
        np.testing.assert_allclose(wall[0, :, 1], -0.5, atol=1e-12)
        np.testing.assert_allclose(wall[-1, :, 1], 0.5, atol=1e-12)
        np.testing.assert_allclose(wall[0, 0], [0.5, -0.5, -0.5], atol=1e-12)  # corners
        np.testing.assert_allclose(wall[-1, -1], [0.5, 0.5, 0.5], atol=1e-12)
        assert np.all(block.nodes[:, :, 1:, 0] > 0.6)  # unassociated layers untouched

    def test_associations_survive_span_edit(self, bd, tmp_path):
        from shore.geometry import load_geometry
        from tests.conftest import _naca0012_section

        def wing(span: float, name: str):
            face = bd.make_face(bd.Polyline(*[tuple(p) for p in _naca0012_section(40)], close=True))
            path = tmp_path / name
            bd.export_step(bd.extrude(face, amount=span), str(path))
            return load_geometry(path)

        g0 = wing(0.5, "w0.step")
        refs = [r for k in (EntityKind.SURFACE, EntityKind.CURVE) for r in g0.entities(k)]
        assocs = [Association.of("wing", g0, [r]) for r in refs]
        # unchanged re-read: every association resolves to the same entity, "same"
        same = AssociationResolver({"wing": wing(0.5, "w1.step")})
        for a, r in zip(assocs, refs, strict=True):
            assert same.refs(a) == [r]
            assert same.reports[a].count("same") == 1
        # +5 % span: every association still resolves, reported "changed"
        edited = AssociationResolver({"wing": wing(0.525, "w2.step")})
        for a in assocs:
            assert len(edited.refs(a)) == 1
            assert edited.reports[a].count("changed") == 1
        # +30 % span is beyond the matching tolerance for some entities
        far = AssociationResolver({"wing": wing(0.65, "w3.step")})
        with pytest.raises(MeshInputError, match="no longer resolves"):
            for a in assocs:
                far.refs(a)
