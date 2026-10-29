"""FROZEN acceptance spec: blocking, skeleton and premesh (ADR 0042).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_blocking.py``
reports 17 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 17 xfailed.

Interface under test
--------------------
A component holds a blocking: skeleton vertices, hexahedral blocks over them
(8 vertex ids, ``HexBlock.vertices()`` order, index ``4 i + 2 j + k``), a
node count per edge, associations of its faces, edges and vertices. Commands
(journaled, undoable; ``Project.blocking`` is the current component's):
- ``p.create_blocking(ni, nj, nk, *, geometry="body", margin=0.0,
  vmin=None, vmax=None)``: one block ``b0``, the bounding box of geometry
  *geometry* with each side moved out by *margin* times the box's extent
  along that axis, or the box between *vmin* and *vmax*.
- ``p.split_blocking(block, axis, fraction)``: cut across *axis* at node
  ``round(fraction * (n - 1))`` of the block's *n* nodes along it, and every
  block sharing a cut edge, at the same node (the sheet). The low part keeps
  its label; the high part of each cut block is ``b<n>``, *n* counting up
  from 1 over the project's blocking, cut blocks taken in blocking order.
- ``p.delete_blocks(labels)``: the faces left alone become boundary faces.
- ``p.associate_blocking(block, entity, geometry, kind, indices)``: as
  ``associate`` (ADR 0039), on a blocking entity.
- ``p.premesh(exact=False)``: the component's mesh, one ``HexBlock`` per
  blocking block, same labels:
  - vertices projected onto their associations;
  - edges: equally spaced nodes on the line between their vertices,
    projected onto their association;
  - faces and blocks: TFI, faces projected onto their association;
  - shared faces SHARED with axis maps, seams bit for bit; boundary faces
    WALL when associated, else FREE;
  - the blocks carry the associations, so ``snap`` keeps the nodes.
- An entity without an association of its own inherits one: an edge that of
  its associated faces, a vertex that of its associated edges, else of its
  associated faces, when they all agree.
- Queries: ``Blocking.labels()``, ``.block(label)`` (vertex ids),
  ``.shape(label)`` (node counts), ``.points`` (vertex positions).
- ``ParameterError`` for a bad block, axis, fraction (no interior node),
  entity or association, deleting every block, and blocking commands without
  a blocking. Project and journal are then unchanged.

Cases
-----
- B1 create around a geometry; B2 between corners, the journal line.
- B3 (x 2) split: node counts, cut position, labels; rounding to a node.
- B4 a cut carried through a row of three blocks (the sheet).
- B5 premesh without associations: a Cartesian grid, seams, BCs.
- B6 the 26-block blocking around a sphere (a 3 x 3 x 3 lattice, centre
  deleted, its faces associated): valid, wall on the sphere, seams closed.
- B7 snap after premesh moves nothing; B8 replay and save/load bit for bit;
  B9 the edited body (an ellipsoid file in place of the sphere) replayed.
- B10 inheritance: agreeing and disagreeing faces, an own association.
- B11 (x 5) refusals; B12 undo.

Out of scope of this spec: the GUI (``tests_gui/spec/test_gui_blocking.py``),
merge, collapse, O-grid insertion and vertex moves (later slices), spacing
laws on blocking edges (``redistribute`` after premesh), smoothing, CAD
curve and surface associations beyond what ``associate`` already covers.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")
ALONG = {"i": ("j", "k"), "j": ("i", "k"), "k": ("i", "j")}


def face_nodes(nodes, face):
    axis = "ijk".index(face[0])
    side = 0 if face.endswith("lo") else -1
    return np.take(nodes, side, axis=axis)  # over the face's along-axes, in order


def seam_gap(blocks):
    """Largest node distance over SHARED faces, mapped through their axis maps."""
    owner = {id(f): b for b in blocks for f in b.faces()}
    worst = 0.0
    for b in blocks:
        for f in b.faces():
            if f.bc.name != "SHARED":
                continue
            partner = owner[id(f.partner)]
            theirs = face_nodes(partner.nodes, f.partner.name)
            (a0, flip0), (_, flip1) = f.axis_map.along
            pax = ALONG[f.partner.name[0]]
            t = theirs if pax.index(a0) == 0 else theirs.transpose(1, 0, 2)
            t = t[::-1] if flip0 else t
            t = t[:, ::-1] if flip1 else t
            worst = max(worst, float(np.abs(face_nodes(b.nodes, f.name) - t).max()))
    return worst


def min_jacobian(blocks):
    from shore.quality import scaled_jacobian

    return min(float(scaled_jacobian(b.nodes).min()) for b in blocks)


def surface_distance(mesh, points):
    _, dist, _ = trimesh.proximity.closest_point(mesh, np.asarray(points).reshape(-1, 3))
    return float(dist.max())


def block_at(p, point):
    """The blocking block whose skeleton box holds *point* (axis-aligned lattices)."""
    b = p.blocking
    for label in b.labels():
        pts = b.points[list(b.block(label))]
        if np.all(pts.min(axis=0) < point) and np.all(point < pts.max(axis=0)):
            return label
    raise AssertionError(f"no block holds {point}")


def lattice(p, cuts=(-1.0, 1.0)):
    """Cut every axis of a cube blocking at the coordinates *cuts*."""
    for axis in "ijk":
        d = "ijk".index(axis)
        for c in cuts:
            point = np.zeros(3)
            point[d] = c
            point[[x for x in range(3) if x != d]] = -2.5  # a corner column of blocks
            label = block_at(p, point)
            pts = p.blocking.points[list(p.blocking.block(label))]
            lo, hi = pts.min(axis=0)[d], pts.max(axis=0)[d]
            p.split_blocking(label, axis, (c - lo) / (hi - lo))


def sphere_project(sphere_stl):
    """The 26-block blocking around the unit sphere, premeshed."""
    from shore.project import Project

    p = Project()
    p.load_geometry(sphere_stl)
    p.create_blocking(13, 13, 13, vmin=(-3, -3, -3), vmax=(3, 3, 3))
    lattice(p)
    centre = block_at(p, np.zeros(3))
    for face in FACES:
        p.associate_blocking(centre, face)
    p.delete_blocks([centre])
    p.premesh()
    return p, centre


def walls(blocks):
    return [face_nodes(b.nodes, f.name) for b in blocks for f in b.faces() if f.bc.name == "WALL"]


class TestCreateSplit:
    def test_b1_around_geometry(self, sphere_stl):
        from shore.project import Project

        p = Project()
        p.load_geometry(sphere_stl)
        p.create_blocking(9, 7, 5, margin=0.5)
        lo, hi = trimesh.load(sphere_stl).bounds
        extent = hi - lo
        assert p.blocking.labels() == ["b0"] and p.blocking.shape("b0") == (9, 7, 5)
        corners = p.blocking.points[list(p.blocking.block("b0"))]
        expected = [[(hi + 0.5 * extent if v >> s & 1 else lo - 0.5 * extent)[d]
                     for d, s in enumerate((2, 1, 0))] for v in range(8)]  # fmt: skip
        np.testing.assert_allclose(corners, expected, atol=1e-12)
        assert p.mesh is None  # the mesh is made by premesh only

    def test_b2_between_corners(self):
        from shore.project import Project

        p = Project()
        p.create_blocking(9, 5, 5, vmin=(0, 0, 0), vmax=(2, 1, 1))
        assert p.journal().splitlines()[-1] == (
            "p.create_blocking(9, 5, 5, vmin=(0.0, 0.0, 0.0), vmax=(2.0, 1.0, 1.0))"
        )
        corners = p.blocking.points[list(p.blocking.block("b0"))]
        np.testing.assert_array_equal(corners[7], [2.0, 1.0, 1.0])
        np.testing.assert_array_equal(corners[4], [2.0, 0.0, 0.0])

    def test_b3_split(self):
        from shore.project import Project

        p = Project()
        p.create_blocking(9, 5, 5, vmin=(0, 0, 0), vmax=(2, 1, 1))
        p.split_blocking("b0", "i", 0.5)
        assert p.journal().splitlines()[-1] == "p.split_blocking('b0', 'i', 0.5)"
        assert sorted(p.blocking.labels()) == ["b0", "b1"]
        assert p.blocking.shape("b0") == (5, 5, 5) and p.blocking.shape("b1") == (5, 5, 5)
        lo = p.blocking.points[list(p.blocking.block("b0"))]
        hi = p.blocking.points[list(p.blocking.block("b1"))]
        assert np.all(lo[4:, 0] == 1.0) and np.all(hi[:4, 0] == 1.0) and np.all(hi[4:, 0] == 2.0)

    def test_b3_rounds_to_a_node(self):
        from shore.project import Project

        p = Project()
        p.create_blocking(9, 5, 5, vmin=(0, 0, 0), vmax=(2, 1, 1))
        p.split_blocking("b0", "i", 0.3)  # node round(2.4) = 2 of 0..8
        assert p.blocking.shape("b0") == (3, 5, 5) and p.blocking.shape("b1") == (7, 5, 5)
        hi = p.blocking.points[list(p.blocking.block("b1"))]
        np.testing.assert_allclose(hi[:4, 0], 0.5, atol=1e-15)

    def test_b4_sheet(self):
        from shore.project import Project

        p = Project()
        p.create_blocking(13, 5, 5, vmin=(0, 0, 0), vmax=(3, 1, 1))
        p.split_blocking("b0", "i", 1 / 3)
        p.split_blocking("b1", "i", 0.5)
        p.split_blocking("b1", "j", 0.5)  # the middle block: the cut runs through the row
        b = p.blocking
        assert sorted(b.labels()) == ["b0", "b1", "b2", "b3", "b4", "b5"]
        for low, high in (("b0", "b3"), ("b1", "b4"), ("b2", "b5")):  # in blocking order
            assert b.shape(low) == (5, 3, 5) and b.shape(high) == (5, 3, 5)
            assert np.all(b.points[list(b.block(low))][2::4, 1] == 0.5)
        p.premesh()
        assert len(p.mesh.blocks) == 6 and seam_gap(p.mesh.blocks) == 0.0
        assert min_jacobian(p.mesh.blocks) > 0.999


class TestPremesh:
    def test_b5_cartesian(self):
        from shore.project import Project

        p = Project()
        p.create_blocking(9, 5, 5, vmin=(0, 0, 0), vmax=(2, 1, 1))
        p.split_blocking("b0", "i", 0.25)
        p.split_blocking("b0", "k", 0.5)
        p.premesh()
        blocks = {b.label: b for b in p.mesh.blocks}
        assert sorted(blocks) == sorted(p.blocking.labels())
        x, y, z = np.linspace(0, 2, 9), np.linspace(0, 1, 5), np.linspace(0, 1, 5)
        grid = np.stack(np.meshgrid(x, y, z, indexing="ij"), axis=-1)
        for label, b in blocks.items():
            corner = p.blocking.points[p.blocking.block(label)[0]]
            i0, j0, k0 = (round(c / h) for c, h in zip(corner, (0.25, 0.25, 0.25), strict=True))
            ni, nj, nk = b.nodes.shape[:3]
            np.testing.assert_allclose(b.nodes, grid[i0 : i0 + ni, j0 : j0 + nj, k0 : k0 + nk],
                                       atol=1e-12)  # fmt: skip
        assert seam_gap(p.mesh.blocks) == 0.0
        bcs = sorted(f.bc.name for b in p.mesh.blocks for f in b.faces())
        assert bcs.count("SHARED") == 8 and bcs.count("FREE") == 16 and "WALL" not in bcs

    def test_b6_sphere(self, sphere_stl):
        p, centre = sphere_project(sphere_stl)
        blocks = p.mesh.blocks
        assert len(blocks) == 26 and centre not in {b.label for b in blocks}
        bcs = [f.bc.name for b in blocks for f in b.faces()]
        assert bcs.count("WALL") == 6 and bcs.count("FREE") == 54 and bcs.count("SHARED") == 96
        assert surface_distance(trimesh.load(sphere_stl), np.concatenate(
            [w.reshape(-1, 3) for w in walls(blocks)])) < 1e-9  # fmt: skip
        assert min_jacobian(blocks) > 0.5
        assert seam_gap(blocks) == 0.0

    def test_b7_snap_keeps_nodes(self, sphere_stl):
        p, _ = sphere_project(sphere_stl)
        before = [b.nodes.copy() for b in p.mesh.blocks]
        p.snap()
        moved = max(float(np.abs(b.nodes - o).max()) for b, o in zip(p.mesh.blocks, before,
                                                                     strict=True))  # fmt: skip
        assert moved < 1e-12

    def test_b8_replay_and_load(self, sphere_stl, tmp_path):
        from shore.io.grd import read_grd_metadata
        from shore.project import Project
        from shore.project.project import run_journal

        p, _ = sphere_project(sphere_stl)
        p.save_journal(tmp_path / "journal.py")
        q = run_journal(tmp_path / "journal.py")
        for a, b in zip(p.mesh.blocks, q.mesh.blocks, strict=True):
            assert a.label == b.label and np.array_equal(a.nodes, b.nodes)
        p.save(tmp_path / "proj")
        r = Project.load(tmp_path / "proj")
        assert r.blocking.labels() == p.blocking.labels()
        p.write_xall(tmp_path / "xall")
        assert len(read_grd_metadata(tmp_path / "xall" / "mesh.grd")) == 26

    def test_b9_edited_body(self, sphere_stl, tmp_path):
        from shore.project.project import run_journal

        p, _ = sphere_project(sphere_stl)
        edited = trimesh.load(sphere_stl)
        edited.apply_scale((1.15, 0.85, 0.95))
        (tmp_path / "in").mkdir()
        edited.export(tmp_path / "in" / "sphere.stl")
        journal = p.journal().replace(repr(str(sphere_stl)), repr("sphere.stl"))
        (tmp_path / "journal.py").write_text(journal)
        q = run_journal(tmp_path / "journal.py", input_dir=tmp_path / "in")
        blocks = q.mesh.blocks
        edited = trimesh.load(tmp_path / "in" / "sphere.stl")  # as stored: float32
        assert surface_distance(edited, np.concatenate(
            [w.reshape(-1, 3) for w in walls(blocks)])) < 1e-9  # fmt: skip
        assert min_jacobian(blocks) > 0.3
        assert seam_gap(blocks) == 0.0

    def test_b10_inheritance(self, sphere_stl, tmp_path):
        from shore.project import Project

        big = trimesh.load(sphere_stl)
        big.apply_scale(2.0)
        big.export(tmp_path / "big.stl")
        p = Project()
        p.load_geometry(sphere_stl)
        p.load_geometry(tmp_path / "big.stl", name="big")
        p.create_blocking(5, 5, 5, vmin=(-0.5, -0.5, -0.5), vmax=(0.5, 0.5, 0.5))
        p.associate_blocking("b0", "i_lo")  # the unit sphere
        p.associate_blocking("b0", "j_lo", "big", "body")
        p.premesh()
        nodes = p.mesh.blocks[0].nodes
        unit = trimesh.load(sphere_stl)
        # i_lo alone: its face, its other edges and corners on the unit sphere
        assert surface_distance(unit, nodes[0, 1:, :]) < 1e-9
        # the edge both faces hold disagrees: a straight line between unmoved corners
        edge = nodes[0, 0, :]
        np.testing.assert_allclose(edge[[0, -1]], [[-0.5, -0.5, -0.5], [-0.5, -0.5, 0.5]],
                                   atol=1e-15)  # fmt: skip
        np.testing.assert_allclose(edge, np.linspace(edge[0], edge[-1], 5), atol=1e-12)
        assert surface_distance(big, nodes[1:, 0, :]) < 1e-9
        # an association of its own wins: the edge and its corners on the big sphere
        p.associate_blocking("b0", "i_lo_j_lo", "big", "body")
        p.premesh()
        assert surface_distance(big, p.mesh.blocks[0].nodes[0, 0, :]) < 1e-9
        b = p.mesh.blocks[0]
        assert b.edge("i_lo_j_lo").association.geometry == "big"
        assert b.edge("i_lo_k_lo").association.geometry == "body"
        assert b.face("i_lo").bc.name == "WALL" and b.face("i_hi").bc.name == "FREE"


class TestRefusals:
    @pytest.mark.parametrize(
        "call",
        [
            ("split_blocking", ("b0", "i", 0.01)),  # node 0: no interior node
            ("split_blocking", ("b0", "x", 0.5)),
            ("split_blocking", ("b9", "i", 0.5)),
            ("delete_blocks", (["b0"],)),
            ("associate_blocking", ("b0", "i_lo_j_lo", "body", "point", (0,))),
        ],
    )
    def test_b11_refused(self, sphere_stl, call):
        from shore.project import Project

        p = Project()
        p.load_geometry(sphere_stl)
        p.create_blocking(9, 5, 5, vmin=(0, 0, 0), vmax=(2, 1, 1))
        journal, labels = p.journal(), p.blocking.labels()
        method, args = call
        with pytest.raises(ParameterError):
            getattr(p, method)(*args)
        assert p.journal() == journal and p.blocking.labels() == labels
        with pytest.raises(ParameterError):
            Project().premesh()  # no blocking

    def test_b12_undo(self):
        from shore.project import Project

        p = Project()
        p.create_blocking(9, 5, 5, vmin=(0, 0, 0), vmax=(2, 1, 1))
        p.split_blocking("b0", "i", 0.5)
        p.premesh()
        nodes = [b.nodes.copy() for b in p.mesh.blocks]
        p.delete_blocks(["b1"])
        p.premesh()
        assert len(p.mesh.blocks) == 1
        p.undo()
        p.undo()
        assert sorted(p.blocking.labels()) == ["b0", "b1"]
        assert all(np.array_equal(b.nodes, n) for b, n in zip(p.mesh.blocks, nodes, strict=True))
        p.undo()
        p.undo()
        assert p.blocking.labels() == ["b0"] and p.mesh is None
