"""FROZEN acceptance spec: grid lines leaving the walls along their normal (ADR 0048).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_orthogonalize.py``
reports 7 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 7 xfailed.

Interface under test
--------------------
``p.orthogonalize(strength=1.0)``, the command ``Orthogonalize`` (journaled,
undoable; journal ``p.orthogonalize()``), on the current component's mesh
(after ``premesh``, ADR 0042):
- every block with exactly one WALL face has its grid lines across the wall
  reshaped to leave the wall along its normal: the wall nodes and the
  opposite face do not move, and each line's nodes keep their order;
- the normal of a wall node is common to every block holding it, so seams
  stay closed; blocks whose faces the seams change are re-filled;
- ``ParameterError`` for a strength <= 0, or one that would invert a block
  that was valid; the mesh and the journal are then unchanged.

Wall angle: the largest (or median) angle between a wall node's first grid
segment and the wall normal (the sphere's radius; on the wing, the STL
normal at the closest point, nodes within 0.05 of the TE, the tip or the
root left out: sharp edges have no single normal).

Cases
-----
- R1 the sphere with ADR 0043's hole layer: wall angle 5.7 -> below 3 deg,
  minimum scaled Jacobian kept, walls and outer faces fixed.
- R2 ADR 0042's 26-block lattice: wall angle 40.8 -> below 25 deg, valid,
  seams closed (re-filled neighbours included).
- R3 (x 2) refusals: an inverting strength, a strength of 0.
- R4 the wing on a plane (ADR 0046 route B): median wall angle below 10 deg
  and p95 below 20 deg where the plain premesh gives a median over 25.
- R5 a spacing law applied after (``redistribute``, ADR 0041) moves the
  nodes along the reshaped lines: every node stays on its line.
- R6 undo, replay bit for bit.

Out of scope of this spec: elliptic smoothing of block interiors, sharp
edges of the wall (no single normal), wall-normal marching.
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
    return np.take(nodes, 0 if face.endswith("lo") else -1, axis="ijk".index(face[0]))


def next_nodes(nodes, face):
    return np.take(nodes, 1 if face.endswith("lo") else -2, axis="ijk".index(face[0]))


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


def wall_angles(blocks, normal):
    """Per wall node: (its position, the angle deg of its first segment to normal(points))."""
    pts, angles = [], []
    for b in blocks:
        for f in b.faces():
            if f.bc.name != "WALL":
                continue
            wall = face_nodes(b.nodes, f.name).reshape(-1, 3)
            seg = next_nodes(b.nodes, f.name).reshape(-1, 3) - wall
            n = normal(wall)
            cos = np.abs(np.einsum("ij,ij->i", seg, n)) / np.linalg.norm(seg, axis=1)
            pts.append(wall)
            angles.append(np.degrees(np.arccos(np.clip(cos, 0.0, 1.0))))
    return np.concatenate(pts), np.concatenate(angles)


def radial(points):
    return points / np.linalg.norm(points, axis=1, keepdims=True)


def block_at(p, point):
    b = p.blocking
    for label in b.labels():
        pts = b.points[list(b.block(label))]
        if np.all(pts.min(axis=0) < point) and np.all(point < pts.max(axis=0)):
            return label
    raise AssertionError(f"no block holds {point}")


def cut(p, axis, c, other):
    d = "ijk".index(axis)
    point = np.array(other, float)
    point[d] = c
    label = block_at(p, point)
    pts = p.blocking.points[list(p.blocking.block(label))]
    a, z = pts.min(axis=0)[d], pts.max(axis=0)[d]
    p.split_blocking(label, axis, (c - a) / (z - a))


def hole_route(p, n, box, xs, ys, zs, centre_point, faces=FACES, geometry="body", nodes=5):
    """ADR 0046 route B: a lattice, the centre block's faces associated and the
    block deleted, the hole layered (factor 0.5, *nodes* across); premeshed."""
    p.create_blocking(*n, vmin=box[0], vmax=box[1])
    corner = np.array(box[0], float) + 0.01
    for axis, cuts in zip("ijk", (xs, ys, zs), strict=True):
        for c in cuts:
            cut(p, axis, c, corner)
    centre = block_at(p, np.asarray(centre_point, float))
    for f in faces:
        p.associate_blocking(centre, f, geometry, "body")
    keys = {p.blocking.entity(centre, f)[0] for f in faces}
    p.delete_blocks([centre])
    hole = [(lab, f) for lab in p.blocking.labels() for f in FACES
            if p.blocking.entity(lab, f)[0] in keys]  # fmt: skip
    p.ogrid(faces=hole, factor=0.5, nodes=nodes)
    p.premesh()


def sphere(sphere_stl, layered):
    from shore.project import Project

    p = Project()
    p.load_geometry(sphere_stl)
    if layered:
        hole_route(p, (13, 13, 13), ((-3,) * 3, (3,) * 3), (-1, 1), (-1, 1), (-1, 1), (0, 0, 0))
        return p
    p.create_blocking(13, 13, 13, vmin=(-3, -3, -3), vmax=(3, 3, 3))
    corner = np.full(3, -2.99)
    for axis in "ijk":
        for c in (-1.0, 1.0):
            cut(p, axis, c, corner)
    centre = block_at(p, np.zeros(3))
    for f in FACES:
        p.associate_blocking(centre, f)
    p.delete_blocks([centre])
    p.premesh()
    return p


def wing(tmp_path):
    """A blunt-TE NACA 0012 wing, chord 1, span 2 along z, root on z = 0 (no root cap)."""
    from shore.profiles import naca
    from shore.shapes import extrude_section

    w = extrude_section(naca("0012", 60, closed_te=False), span=2.0, n_stations=21,
                        span_axis="z", end_caps=True)  # fmt: skip
    root = np.all(np.abs(w.vertices[w.faces][:, :, 2]) < 1e-12, axis=1)
    w.update_faces(~root)
    w.remove_unreferenced_vertices()
    path = tmp_path / "wing.stl"
    w.export(path)
    return path


class TestOrthogonalize:
    @PENDING
    def test_r1_hole_layer(self, sphere_stl):
        p = sphere(sphere_stl, layered=True)
        _, before = wall_angles(p.mesh.blocks, radial)
        jac = min_jacobian(p.mesh.blocks)
        fixed = {}  # each wall block's wall face and opposite face
        for b in p.mesh.blocks:
            for f in b.faces():
                if f.bc.name == "WALL":
                    opposite = f.name[:2] + ("hi" if f.name.endswith("lo") else "lo")
                    fixed[(b.label, f.name)] = face_nodes(b.nodes, f.name).copy()
                    fixed[(b.label, opposite)] = face_nodes(b.nodes, opposite).copy()
        assert len(fixed) == 12
        p.orthogonalize()
        assert p.journal().splitlines()[-1] == "p.orthogonalize()"
        _, after = wall_angles(p.mesh.blocks, radial)
        assert before.max() > 5.0 and after.max() < 3.0
        assert min_jacobian(p.mesh.blocks) >= jac - 1e-9
        assert seam_gap(p.mesh.blocks) == 0.0
        blocks = {b.label: b for b in p.mesh.blocks}
        for (label, face), old in fixed.items():
            np.testing.assert_array_equal(face_nodes(blocks[label].nodes, face), old)

    @PENDING
    def test_r2_lattice(self, sphere_stl):
        p = sphere(sphere_stl, layered=False)
        _, before = wall_angles(p.mesh.blocks, radial)
        p.orthogonalize()
        _, after = wall_angles(p.mesh.blocks, radial)
        assert before.max() > 40.0 and after.max() < 25.0
        assert min_jacobian(p.mesh.blocks) > 0.2
        assert seam_gap(p.mesh.blocks) == 0.0

    @PENDING
    @pytest.mark.parametrize("strength", [2.0, 0.0])
    def test_r3_refused(self, sphere_stl, strength):
        p = sphere(sphere_stl, layered=False)
        journal = p.journal()
        nodes = [b.nodes.copy() for b in p.mesh.blocks]
        with pytest.raises(ParameterError):
            p.orthogonalize(strength)
        assert p.journal() == journal
        for b, old in zip(p.mesh.blocks, nodes, strict=True):
            np.testing.assert_array_equal(b.nodes, old)

    @PENDING
    def test_r4_wing(self, tmp_path):
        from shore.project import Project

        path = wing(tmp_path)
        mesh = trimesh.load(path)

        def normal(points):
            _, _, tri = trimesh.proximity.closest_point(mesh, points)
            return mesh.face_normals[tri]

        p = Project()
        p.load_geometry(path, name="wing", skip_hygiene=True)
        hole_route(p, (41, 21, 21), ((-3, -3, 0), (4, 3, 5)), (-0.2, 1.2), (-0.3, 0.3), (2.3,),
                   (0.5, 0.0, 1.0), faces=("i_lo", "i_hi", "j_lo", "j_hi", "k_hi"),
                   geometry="wing", nodes=9)  # fmt: skip

        def smooth_region(blocks):
            pts, ang = wall_angles(blocks, normal)
            keep = (pts[:, 2] > 0.05) & (pts[:, 2] < 1.95) & (pts[:, 0] < 0.95)
            return ang[keep]

        before = smooth_region(p.mesh.blocks)
        p.orthogonalize()
        after = smooth_region(p.mesh.blocks)
        assert np.median(before) > 25.0
        assert np.median(after) < 10.0 and np.percentile(after, 95) < 20.0
        assert min_jacobian(p.mesh.blocks) > 0.3 and seam_gap(p.mesh.blocks) == 0.0

    @PENDING
    def test_r5_then_redistribute(self, sphere_stl):
        p = sphere(sphere_stl, layered=True)
        p.orthogonalize()
        layer = next(b for b in p.mesh.blocks if any(f.bc.name == "WALL" for f in b.faces()))
        wall = next(f.name for f in layer.faces() if f.bc.name == "WALL")
        axis = "ijk".index(wall[0])
        lines = np.moveaxis(layer.nodes, axis, -2).reshape(-1, layer.nodes.shape[axis], 3).copy()
        across = {"i": "j_lo_k_lo", "j": "i_lo_k_lo", "k": "i_lo_j_lo"}[wall[0]]
        p.redistribute(layer.label, across, "tanh", beta=3.0)
        layer = next(b for b in p.mesh.blocks if b.label == layer.label)
        new = np.moveaxis(layer.nodes, axis, -2).reshape(lines.shape)
        assert not np.array_equal(new, lines)
        for old_line, new_line in zip(lines, new, strict=True):  # each node on its old line
            a, b = old_line[:-1], old_line[1:]
            ab = b - a
            t = np.clip(np.einsum("sij,ij->si", new_line[:, None] - a, ab) / (ab * ab).sum(1), 0, 1)
            d = np.linalg.norm(a + t[..., None] * ab - new_line[:, None], axis=-1).min(axis=1)
            assert d.max() < 1e-12
        assert min_jacobian(p.mesh.blocks) > 0.0 and seam_gap(p.mesh.blocks) == 0.0

    @PENDING
    def test_r6_undo_replay(self, sphere_stl, tmp_path):
        from shore.project.project import run_journal

        p = sphere(sphere_stl, layered=True)
        nodes = [b.nodes.copy() for b in p.mesh.blocks]
        p.orthogonalize(0.8)
        p.save_journal(tmp_path / "journal.py")
        q = run_journal(tmp_path / "journal.py")
        for a, b in zip(p.mesh.blocks, q.mesh.blocks, strict=True):
            assert np.array_equal(a.nodes, b.nodes)
        p.undo()
        for b, old in zip(p.mesh.blocks, nodes, strict=True):
            np.testing.assert_array_equal(b.nodes, old)
