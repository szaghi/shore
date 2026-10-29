"""FROZEN acceptance spec: associating block entities with geometry (ADR 0039).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (CAD CI job, needs ``shore[cad]``; without build123d
the module is skipped, as the other CAD specs are):
``pytest tests/spec/test_associate.py`` reports 17 passed, 0 skipped, 0
xfailed. At freeze: 1 passed (A7, a guard), 16 xfailed.

Interface under test
--------------------
- ``Project.associate(block, entity, geometry="body", kind="body",
  indices=())``, the command ``Associate``. *entity* is a face (``"k_lo"``),
  an edge (``"i_lo_k_lo"``) or a vertex (``"v0"`` .. ``"v7"``, ``HexBlock.
  vertices()`` order) of a block of the current component. *kind* is
  ``"body"``, ``"surface"``, ``"curve"``, ``"point"``, or ``"none"`` to
  clear.
  - The association is ``Association.of(geometry, ..., [EntityRef(kind,
    i)])``: fingerprints, so it survives a re-import (ADR 0003).
  - An entity takes geometry of its own dimension or higher: a face a
    surface or the body; an edge also a curve; a vertex also a point.
  - Journal line ``p.associate('box', 'k_lo', 'part', 'surface', (2,))``;
    undo restores the previous association; nodes move only on ``snap``.
  - An unknown block, entity, geometry or kind, an index out of range, or
    an entity of higher dimension than the geometry raises
    ``ParameterError``, leaving project and journal unchanged.
- ``Project.snap()`` resolves associations against every project geometry
  by name. ``"body"`` is the mesh's own body, or, for a mesh loaded without
  replay (ADR 0033), the geometry its component was built on.

Cases
-----
A STEP part (a 2 x 1 x 1 box with a 0.3-radius hole) and a box primitive
above its top face (z 0.55 to 0.9).
- A1 (x 3) a face to the top surface, an edge to a top edge, a vertex to a
  corner: journal line, the association, and after ``snap`` the nodes on
  that entity.
- A2 the vertex -> edge -> face priority at shared nodes.
- A3 undo and redo; ``kind="none"`` clears.
- A4 (x 8) refusals: bad names, kinds, indices (x 7), dimensions; project
  and journal unchanged.
- A5 (x 2) saved and loaded, with and without replay: the associations,
  and the snap's result.
- A6 a cubed sphere loaded without replay snaps (its walls' ``"body"``
  association), within 1e-12 of the replayed project's snap (its seam nodes
  are no longer shared views, so they are projected once, not twice).
- A7 a snap with only wall associations equals the snap before ADR 0039
  (a guard).

Out of scope of this spec: the GUI (``tests_gui/spec/test_gui_associate.
py``), re-import of a modified CAD file (ADR 0003's key matching), exact
(tier-2) projection, association of split blocks.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from shore.exceptions import ParameterError

pytest.importorskip("build123d")  # the CAD job runs this module; the main job skips it

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

BOX = {"vmin": (-0.8, -0.4, 0.55), "vmax": (0.8, 0.4, 0.9), "ni": 5, "nj": 4, "nk": 3}


@pytest.fixture(scope="module")
def part(tmp_path_factory) -> Path:
    import build123d as bd

    path = tmp_path_factory.mktemp("cad") / "part.step"
    bd.export_step(bd.Box(2.0, 1.0, 1.0) - bd.Cylinder(0.3, 2.0), str(path))
    return path


def setup(part):
    from shore.project import Project

    p = Project()
    p.load_geometry(part, name="part")
    p.build_primitive("box", **BOX)
    return p


def entities(p):
    """(top surface, a top edge along x at y=-0.5, the corner (-1, -0.5, 0.5))."""
    from shore.geometry.base import EntityKind, EntityRef

    g = p.geometries["part"].geometry
    surfaces = [i for i in range(len(g.entities(EntityKind.SURFACE))) if _is_top(g, i)]
    curves = [i for i in range(len(g.entities(EntityKind.CURVE))) if _on_line(g, i, y=-0.5, z=0.5)]
    corner = np.array([[-1.0, -0.5, 0.5]])
    points = [
        i
        for i in range(len(g.entities(EntityKind.POINT)))
        if np.allclose(g.project(corner, EntityRef(EntityKind.POINT, i)).points[0], corner[0])
    ]
    assert len(surfaces) == 1 and len(curves) == 1 and len(points) == 1
    return surfaces[0], curves[0], points[0]


def _is_top(g, i) -> bool:
    from shore.geometry.base import EntityKind, EntityRef

    probe = g.project(np.array([[0.7, 0.3, 0.8]]), EntityRef(EntityKind.SURFACE, i))
    return bool(np.allclose(probe.points[0], (0.7, 0.3, 0.5)))


def _on_line(g, i, *, y: float, z: float) -> bool:
    from shore.geometry.base import EntityKind, EntityRef

    probe = g.project(np.array([[0.2, y - 0.1, z + 0.1], [-0.6, y - 0.1, z + 0.1]]),
                      EntityRef(EntityKind.CURVE, i)).points  # fmt: skip
    return bool(np.allclose(probe[:, 1:], (y, z)) and abs(probe[0, 0] - probe[1, 0]) > 0.5)


def ref_distance(p, nodes, kind, index) -> float:
    from shore.geometry.base import EntityKind, EntityRef

    g = p.geometries["part"].geometry
    proj = g.project(np.asarray(nodes).reshape(-1, 3), EntityRef(EntityKind(kind), index))
    return float(proj.distances.max())


# ---------------------------------------------------------------------------
# A1-A3: the command
# ---------------------------------------------------------------------------


class TestAssociate:
    @PENDING
    @pytest.mark.parametrize("entity", ["face", "edge", "vertex"])
    def test_a1_associate_and_snap(self, part, entity):
        from shore.geometry.base import EntityKind, EntityRef
        from shore.mesh.association import Association

        p = setup(part)
        top, edge, corner = entities(p)
        block = p.mesh.blocks[0]
        name, kind, index, target, view = {
            "face": ("k_lo", "surface", top, block.face("k_lo"), lambda n: n[:, :, 0]),
            "edge": ("j_lo_k_lo", "curve", edge, block.edge("j_lo_k_lo"), lambda n: n[:, 0, 0]),
            "vertex": ("v0", "point", corner, block.vertices()[0], lambda n: n[0, 0, 0]),
        }[entity]
        p.associate("box", name, "part", kind, [index])
        assert (
            p.journal().splitlines()[-1]
            == f"p.associate('box', {name!r}, 'part', {kind!r}, ({index},))"
        )
        g = p.geometries["part"].geometry
        assert target.association == Association.of("part", g, [EntityRef(EntityKind(kind), index)])
        before = block.nodes.copy()
        assert ref_distance(p, view(before), kind, index) > 1e-3  # nothing moved yet
        p.snap()
        assert ref_distance(p, view(block.nodes), kind, index) <= 1e-9

    @PENDING
    def test_a2_priority(self, part):
        p = setup(part)
        top, edge, corner = entities(p)
        p.associate("box", "k_lo", "part", "surface", [top])
        p.associate("box", "j_lo_k_lo", "part", "curve", [edge])
        p.associate("box", "v0", "part", "point", [corner])
        p.snap()
        n = p.mesh.blocks[0].nodes
        assert np.allclose(n[0, 0, 0], (-1.0, -0.5, 0.5), atol=1e-12)  # the vertex wins
        assert np.allclose(n[:, 0, 0, 1:], (-0.5, 0.5), atol=1e-12)  # the edge over the face
        assert np.allclose(n[:, 1:, 0, 2], 0.5, atol=1e-12)  # the face elsewhere

    @PENDING
    def test_a3_undo_redo_clear(self, part):
        p = setup(part)
        top, _, _ = entities(p)
        face = p.mesh.blocks[0].face("k_lo")
        before = face.association
        p.associate("box", "k_lo", "part", "surface", [top])
        after = face.association
        p.undo()
        assert face.association == before
        p.redo()
        assert face.association == after
        p.associate("box", "k_lo", "part", "none")
        assert face.association is None
        assert p.journal().splitlines()[-1] == "p.associate('box', 'k_lo', 'part', 'none')"


class TestRefusals:
    @PENDING
    @pytest.mark.parametrize(
        "args",
        [
            ("nope", "k_lo", "part", "surface", [0]),
            ("box", "nope", "part", "surface", [0]),
            ("box", "k_lo", "nope", "surface", [0]),
            ("box", "k_lo", "part", "volume", [0]),
            ("box", "k_lo", "part", "surface", [99]),
            ("box", "k_lo", "part", "surface", []),
            ("box", "k_lo", "part", "body", [1]),
        ],
    )
    def test_a4_bad_arguments(self, part, args):
        p = setup(part)
        journal = p.journal()
        before = [f.association for f in p.mesh.blocks[0].faces()]
        with pytest.raises(ParameterError):
            p.associate(*args)
        assert p.journal() == journal
        assert [f.association for f in p.mesh.blocks[0].faces()] == before

    @PENDING
    def test_a4_dimensions(self, part):
        p = setup(part)
        top, edge, corner = entities(p)
        with pytest.raises(ParameterError):
            p.associate("box", "k_lo", "part", "curve", [edge])
        with pytest.raises(ParameterError):
            p.associate("box", "k_lo", "part", "point", [corner])
        with pytest.raises(ParameterError):
            p.associate("box", "j_lo_k_lo", "part", "point", [corner])
        p.associate("box", "j_lo_k_lo", "part", "surface", [top])  # higher dimension: fine
        p.associate("box", "v3", "part", "curve", [edge])
        p.associate("box", "v4", "part", "body")


# ---------------------------------------------------------------------------
# A5-A7: persistence and snap
# ---------------------------------------------------------------------------


class TestPersistence:
    @PENDING
    @pytest.mark.parametrize("replay", [True, False])
    def test_a5_save_load(self, part, tmp_path, replay):
        from shore.project import Project

        p = setup(part)
        top, edge, corner = entities(p)
        p.associate("box", "k_lo", "part", "surface", [top])
        p.associate("box", "j_lo_k_lo", "part", "curve", [edge])
        p.associate("box", "v0", "part", "point", [corner])
        p.save(tmp_path / "saved")
        q = Project.load(tmp_path / "saved", replay=replay)
        a, b = p.mesh.blocks[0], q.mesh.blocks[0]
        assert [f.association for f in a.faces()] == [f.association for f in b.faces()]
        assert [e.association for e in a.edges()] == [e.association for e in b.edges()]
        assert [v.association for v in a.vertices()] == [v.association for v in b.vertices()]
        p.snap()
        q.snap()
        assert np.array_equal(p.mesh.blocks[0].nodes, q.mesh.blocks[0].nodes)

    @PENDING
    def test_a6_snap_after_load_without_replay(self, tmp_path):
        import trimesh

        from shore.project import Project

        trimesh.creation.icosphere(subdivisions=3).export(str(tmp_path / "s.stl"))
        p = Project()
        p.load_geometry(tmp_path / "s.stl")
        p.build_topology("cubed_sphere", nk=5, ds=0.05, growth=1.1, ni=9, nj=16)
        p.march()
        p.save(tmp_path / "saved")
        q = Project.load(tmp_path / "saved", replay=False)
        p.snap()
        q.snap()  # raised before ADR 0039: the loaded mesh has no Body
        # blocks read back from grids no longer share seam nodes as views, so a
        # seam node is projected once per block instead of twice: round-off apart
        for a, b in zip(p.mesh.blocks, q.mesh.blocks, strict=True):
            np.testing.assert_allclose(b.nodes, a.nodes, rtol=0.0, atol=1e-12)

    def test_a7_wall_snap_unchanged(self, tmp_path):  # already holds: a guard
        import trimesh

        from shore.mesh.association import AssociationResolver, snap_block
        from shore.project import Project

        trimesh.creation.icosphere(subdivisions=3).export(str(tmp_path / "s.stl"))
        p = Project()
        p.load_geometry(tmp_path / "s.stl")
        p.build_topology("ogrid", nk=5, ds=0.05, growth=1.1, ni=9, nj=16)
        p.march()
        expected = [b.nodes.copy() for b in p.mesh.blocks]
        resolver = AssociationResolver({"body": p.mesh.body.geometry})
        for b, nodes in zip(p.mesh.blocks, expected, strict=True):
            twin = type(b)(label=b.label, nodes=nodes)
            twin.face("k_lo").association = b.face("k_lo").association
            for e, f in zip(twin.edges(), b.edges(), strict=True):
                e.association = f.association
            for v, w in zip(twin.vertices(), b.vertices(), strict=True):
                v.association = w.association
            snap_block(twin, resolver)
        p.snap()
        for b, nodes in zip(p.mesh.blocks, expected, strict=True):
            assert np.array_equal(b.nodes, nodes)
