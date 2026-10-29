"""FROZEN acceptance spec: O-grid insertion in a blocking (ADR 0043).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_ogrid.py``
reports 13 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 13 xfailed.

Interface under test
--------------------
``p.ogrid(blocks=(), *, faces=(), exclude=(), factor=0.5, nodes=5)``, the
command ``Ogrid`` (journaled, undoable), on the current component's blocking
(ADR 0042):
- the block set is *blocks*, or every block when empty;
- the layered faces are *faces* (``(block, face)`` pairs), or every face on
  the set's boundary (a face of a set block not shared with another set
  block), less *exclude*; all must lie on that boundary;
- the set's other boundary faces stay open: the O-grid passes through them,
  and they must not be shared with a block outside the set;
- every vertex ``v`` of a layered face gets an inner copy
  ``v + (1 - factor) * d / 2``, where ``d`` sums, over the distinct edges of
  set blocks that leave a layered face at ``v`` into its block, the edge
  vector from ``v`` (inside one block: the block scaled by *factor* about its
  centre);
- the set's blocks take the inner copies and keep their labels; each layered
  face gets a block between it and its copy, labelled ``b<n>`` in block and
  face order (``i_lo, i_hi, j_lo, j_hi, k_lo, k_hi``), with *nodes* nodes
  across and the face's counts along it;
- associations: a layered face keeps its own (the layer's outer face); a
  face or edge of the set that the copies replace passes its association to
  its replacement; an open face also to the layer faces lying on it;
- ``ParameterError`` for a factor outside (0, 1), fewer than 2 nodes, a bad
  block, a face off the set's boundary, an open face shared with a block
  outside the set, no face left to layer. Project and journal are then
  unchanged.

Cases
-----
- O1 one block: core, layers, labels, counts, premesh.
- O2 an open face (half-O): the core's face in its plane, associations
  passed to the faces on it.
- O3 two blocks, the whole blocking: the shared face stays inside.
- O4 the faces of a hole: ADR 0042's 26 blocks around a sphere get a layer
  on the hole. Grid lines leave the wall within 10 degrees of its normal
  (without the O-grid, over 30).
- O5 the centre block O-gridded, its core associated and deleted.
- O6 a layered face keeps its association on the layer only.
- O7 journal, replay, undo; O8 (x 6) refusals.

Out of scope of this spec: picking blocks and faces in the GUI (forms only,
``tests_gui/spec/test_gui_ogrid.py``), O-grids through faces shared with
blocks outside the set, smoothing of the layer, collapse, merge and vertex
moves.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")
ALONG = {"i": ("j", "k"), "j": ("i", "k"), "k": ("i", "j")}
BIT = {"i": 4, "j": 2, "k": 1}


def corners(face):
    """Local vertex indices of *face*, over its along-axes."""
    axis, hi = face[0], face.endswith("hi")
    p, q = ALONG[axis]
    base = BIT[axis] if hi else 0
    return [[base + a * BIT[p] + b * BIT[q] for b in (0, 1)] for a in (0, 1)]


def face_key(p, label, face):
    verts = p.blocking.block(label)
    return ("f", *sorted(verts[c] for row in corners(face) for c in row))


def face_nodes(nodes, face):
    return np.take(nodes, 0 if face.endswith("lo") else -1, axis="ijk".index(face[0]))


def seam_gap(blocks):
    owner = {id(f): b for b in blocks for f in b.faces()}
    worst = 0.0
    for b in blocks:
        for f in b.faces():
            if f.bc.name != "SHARED":
                continue
            theirs = face_nodes(owner[id(f.partner)].nodes, f.partner.name)
            (a0, flip0), (_, flip1) = f.axis_map.along
            t = theirs if ALONG[f.partner.name[0]].index(a0) == 0 else theirs.transpose(1, 0, 2)
            t = t[::-1] if flip0 else t
            t = t[:, ::-1] if flip1 else t
            worst = max(worst, float(np.abs(face_nodes(b.nodes, f.name) - t).max()))
    return worst


def min_jacobian(blocks):
    from shore.quality import scaled_jacobian

    return min(float(scaled_jacobian(b.nodes).min()) for b in blocks)


def wall_angle(blocks):
    """Largest angle (deg) between a wall node's first grid segment and the radius."""
    worst = 0.0
    for b in blocks:
        for f in b.faces():
            if f.bc.name != "WALL":
                continue
            a, lo = "ijk".index(f.name[0]), f.name.endswith("lo")
            wall = np.take(b.nodes, 0 if lo else -1, axis=a)
            seg = np.take(b.nodes, 1 if lo else -2, axis=a) - wall
            cos = np.einsum("...i,...i", seg, wall) / (
                np.linalg.norm(seg, axis=-1) * np.linalg.norm(wall, axis=-1)
            )
            worst = max(worst, float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))).max()))
    return worst


def box(ni=9, lo=(-1, -1, -1), hi=(1, 1, 1)):
    from shore.project import Project

    p = Project()
    p.create_blocking(ni, ni, ni, vmin=lo, vmax=hi)
    return p


def corner_points(p, label):
    return p.blocking.points[list(p.blocking.block(label))]


def block_at(p, point):
    b = p.blocking
    for label in b.labels():
        pts = b.points[list(b.block(label))]
        if np.all(pts.min(axis=0) < point) and np.all(point < pts.max(axis=0)):
            return label
    raise AssertionError(f"no block holds {point}")


def lattice(p, lo, cut):
    for axis in "ijk":
        d = "ijk".index(axis)
        for c in (-cut, cut):
            point = np.zeros(3)
            point[d] = c
            point[[x for x in range(3) if x != d]] = -0.5 * (lo + cut)
            label = block_at(p, point)
            pts = corner_points(p, label)
            a, z = pts.min(axis=0)[d], pts.max(axis=0)[d]
            p.split_blocking(label, axis, (c - a) / (z - a))


def sphere_lattice(sphere_stl, lo, cut):
    """ADR 0042's blocking: a lattice of [-lo, lo] cut at +-cut, the centre's
    faces associated with the sphere; the centre's label."""
    from shore.project import Project

    p = Project()
    p.load_geometry(sphere_stl)
    p.create_blocking(13, 13, 13, vmin=(-lo,) * 3, vmax=(lo,) * 3)
    lattice(p, lo, cut)
    centre = block_at(p, np.zeros(3))
    return p, centre


def hole_faces(p, centre):
    """The (block, face) pairs of the set's blocks on the deleted centre's faces."""
    keys = {face_key(p, centre, f) for f in FACES}
    return [(label, f) for label in p.blocking.labels() if label != centre for f in FACES
            if face_key(p, label, f) in keys]  # fmt: skip


class TestOgrid:
    @PENDING
    def test_o1_one_block(self):
        p = box()
        p.ogrid(["b0"], factor=0.5, nodes=4)
        b = p.blocking
        assert b.labels() == ["b0", "b1", "b2", "b3", "b4", "b5", "b6"]
        np.testing.assert_allclose(corner_points(p, "b0"), 0.5 * corner_points(box(), "b0"),
                                   atol=1e-15)  # fmt: skip
        assert b.shape("b0") == (9, 9, 9)
        assert b.shape("b1") == (4, 9, 9) and b.shape("b2") == (4, 9, 9)  # i_lo, i_hi layers
        assert b.shape("b5") == (9, 9, 4) and b.shape("b6") == (9, 9, 4)  # k_lo, k_hi
        # each layer between its outer face and the core's
        outer = corner_points(p, "b2")
        assert np.all(outer[4:, 0] == 1.0) and np.all(outer[:4, 0] == 0.5)
        p.premesh()
        blocks = p.mesh.blocks
        assert len(blocks) == 7 and seam_gap(blocks) == 0.0 and min_jacobian(blocks) > 0.5
        bcs = [f.bc.name for x in blocks for f in x.faces()]
        assert bcs.count("FREE") == 6 and bcs.count("SHARED") == 36

    @PENDING
    def test_o2_open_face(self, tmp_path):
        floor = trimesh.creation.box(extents=(10.0, 10.0, 1.0))
        floor.apply_translation((0.0, 0.0, -1.5))  # its top face: the plane z = -1
        floor.export(tmp_path / "floor.stl")
        p = box()
        p.load_geometry(tmp_path / "floor.stl", name="floor")
        p.associate_blocking("b0", "k_lo", "floor", "body")
        p.ogrid(["b0"], exclude=[("b0", "k_lo")], factor=0.5, nodes=4)
        b = p.blocking
        assert b.labels() == ["b0", "b1", "b2", "b3", "b4", "b5"]
        core = corner_points(p, "b0")
        np.testing.assert_allclose(core[0::2, 2], -1.0, atol=1e-15)  # k_lo corners: in plane
        np.testing.assert_allclose(core[1::2, 2], 0.5, atol=1e-15)
        np.testing.assert_allclose(np.abs(core[:, :2]), 0.5, atol=1e-15)
        # the open face's association: on the core's face and the four layer faces on it
        for label in b.labels()[:5]:
            assert b.associations[face_key(p, label, "k_lo")].geometry == "floor"
        p.premesh()
        assert seam_gap(p.mesh.blocks) == 0.0 and min_jacobian(p.mesh.blocks) > 0.5
        walls = [(x.label, f.name) for x in p.mesh.blocks for f in x.faces() if f.bc.name == "WALL"]
        assert sorted(walls) == [(label, "k_lo") for label in ("b0", "b1", "b2", "b3", "b4")]

    @PENDING
    def test_o3_two_blocks(self):
        p = box(17, lo=(-2, -1, -1), hi=(2, 1, 1))
        p.split_blocking("b0", "i", 0.5)
        p.ogrid(factor=0.5, nodes=4)
        b = p.blocking
        assert len(b.labels()) == 12
        # the middle top vertex: moved by (1 - 0.5) * ((0, -2, 0) + (0, 0, -2)) / 2
        top = corner_points(p, "b0")[7]
        np.testing.assert_allclose(top, [0.0, 0.5, 0.5], atol=1e-15)
        np.testing.assert_allclose(corner_points(p, "b0")[0], [-1.5, -0.5, -0.5], atol=1e-15)
        assert face_key(p, "b0", "i_hi") == face_key(p, "b1", "i_lo")  # still shared
        p.premesh()
        assert seam_gap(p.mesh.blocks) == 0.0 and min_jacobian(p.mesh.blocks) > 0.5

    @PENDING
    def test_o4_hole(self, sphere_stl):
        p, centre = sphere_lattice(sphere_stl, 3.0, 1.0)
        for face in FACES:
            p.associate_blocking(centre, face)
        hole = hole_faces(p, centre)
        assert len(hole) == 6
        p.delete_blocks([centre])
        p.premesh()
        without = wall_angle(p.mesh.blocks)
        p.ogrid(faces=hole, factor=0.5, nodes=5)
        p.premesh()
        blocks = p.mesh.blocks
        assert len(blocks) == 32
        assert sum(f.bc.name == "WALL" for x in blocks for f in x.faces()) == 6
        mesh = trimesh.load(sphere_stl)
        wall = np.concatenate([face_nodes(x.nodes, f.name).reshape(-1, 3) for x in blocks
                               for f in x.faces() if f.bc.name == "WALL"])  # fmt: skip
        assert float(trimesh.proximity.closest_point(mesh, wall)[1].max()) < 1e-9
        assert seam_gap(blocks) == 0.0 and min_jacobian(blocks) > 0.3
        assert without > 30.0 and wall_angle(blocks) < 10.0

    @PENDING
    def test_o5_centre_block(self, sphere_stl):
        p, centre = sphere_lattice(sphere_stl, 4.0, 2.0)
        p.ogrid([centre], factor=0.5, nodes=5)
        for face in FACES:
            p.associate_blocking(centre, face)  # the core's faces
        p.delete_blocks([centre])
        p.premesh()
        blocks = p.mesh.blocks
        assert len(blocks) == 32 and seam_gap(blocks) == 0.0 and min_jacobian(blocks) > 0.3
        assert wall_angle(blocks) < 10.0

    @PENDING
    def test_o6_layered_face_association(self, sphere_stl):
        p = box()
        p.load_geometry(sphere_stl)
        p.associate_blocking("b0", "i_hi")
        outer = face_key(p, "b0", "i_hi")
        p.ogrid(["b0"], factor=0.5, nodes=4)
        b = p.blocking
        assert face_key(p, "b2", "i_hi") == outer and outer in b.associations
        assert face_key(p, "b0", "i_hi") not in b.associations  # the core's face
        assert face_key(p, "b2", "i_lo") == face_key(p, "b0", "i_hi")

    @PENDING
    def test_o7_journal_replay_undo(self, tmp_path):
        from shore.project.project import run_journal

        p = box()
        p.ogrid(["b0"], exclude=[("b0", "k_lo")], factor=0.6, nodes=3)
        assert p.journal().splitlines()[-1] == (
            "p.ogrid(('b0',), exclude=(('b0', 'k_lo'),), factor=0.6, nodes=3)"
        )
        p.premesh()
        p.save_journal(tmp_path / "journal.py")
        q = run_journal(tmp_path / "journal.py")
        for a, c in zip(p.mesh.blocks, q.mesh.blocks, strict=True):
            assert a.label == c.label and np.array_equal(a.nodes, c.nodes)
        p.undo()
        p.undo()
        assert p.blocking.labels() == ["b0"]
        np.testing.assert_array_equal(corner_points(p, "b0"), corner_points(box(), "b0"))


class TestRefusals:
    @PENDING
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"blocks": ["b0"], "factor": 1.0},
            {"blocks": ["b0"], "nodes": 1},
            {"blocks": ["b9"]},
            {"blocks": ["b0", "b1"], "faces": [("b0", "i_hi")]},  # inside the set
            {"blocks": ["b0"], "exclude": [("b0", "i_hi")]},  # open, shared with b1
            {"blocks": ["b0"], "exclude": [("b0", f) for f in FACES]},  # nothing left
        ],
    )
    def test_o8_refused(self, kwargs):
        p = box(17, lo=(-2, -1, -1), hi=(2, 1, 1))
        p.split_blocking("b0", "i", 0.5)
        journal, labels = p.journal(), p.blocking.labels()
        with pytest.raises(ParameterError):
            p.ogrid(**kwargs)
        assert p.journal() == journal and p.blocking.labels() == labels
