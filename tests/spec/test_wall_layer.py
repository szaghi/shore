"""FROZEN acceptance spec: wall layers with an absolute first cell (ADR 0050).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_wall_layer.py``
reports 6 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 6 xfailed.

Interface under test
--------------------
``p.wall_layer(first_step)``, the command ``WallLayer`` (journaled,
undoable; journal ``p.wall_layer(0.001)``), on the current mesh: on every
line across the wall of each block with exactly one WALL face, geometric
cells from the wall, the first exactly *first_step* long (absolute), the
ratio solved per line from its length; nodes stay on the line's polyline
(after ``orthogonalize``, ADR 0048: along the wall normal), the wall node
and the far node do not move. Seams stay closed; blocks whose faces the
seams change are re-filled; ``ParameterError`` for a first step <= 0, one
not below a line's uniform spacing, or one that would invert a valid block
(project and journal unchanged).

Cases
-----
- L1 the sphere's hole layer (17 nodes across), orthogonalized: first cells
  1e-3, growing, ends fixed, seams closed, Jacobian and wall angle kept.
- L2 the wing on a plane (ADR 0046 route B, 25 nodes across),
  orthogonalized: first cells 1e-4, valid, median wall angle below 5 deg.
- L3 (x 3) refusals: 0, a step above a line's uniform spacing, a lattice
  without a layer (ADR 0042) whose corner blocks would invert.
- L4 journal, replay, undo.

Out of scope of this spec: hyperbolic marching of the wall layer (the
algebraic layer here is its substitute, ADR 0050), growth-ratio limits,
blocks with several walls.
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


def layered_sphere(sphere_stl):
    from shore.project import Project

    p = Project()
    p.load_geometry(sphere_stl)
    hole_route(p, (13, 13, 13), ((-3,) * 3, (3,) * 3), (-1, 1), (-1, 1), (-1, 1), (0, 0, 0),
               nodes=17)  # fmt: skip
    p.orthogonalize()
    return p


def wall_lines(blocks):
    """Every line across the wall of a block with one WALL face, wall node first."""
    out = []
    for b in blocks:
        walls = [f.name for f in b.faces() if f.bc.name == "WALL"]
        if len(walls) != 1:
            continue
        axis = "ijk".index(walls[0][0])
        x = np.moveaxis(b.nodes, axis, -2)
        if walls[0].endswith("hi"):
            x = x[..., ::-1, :]
        out.append(x.reshape(-1, x.shape[-2], 3))
    return np.concatenate(out)


class TestWallLayer:
    def test_l1_sphere(self, sphere_stl):
        p = layered_sphere(sphere_stl)
        before = wall_lines(p.mesh.blocks)
        jac = min_jacobian(p.mesh.blocks)
        _, angle = wall_angles(p.mesh.blocks, radial)
        p.wall_layer(1e-3)
        assert p.journal().splitlines()[-1] == "p.wall_layer(0.001)"
        after = wall_lines(p.mesh.blocks)
        cells = np.linalg.norm(np.diff(after, axis=1), axis=-1)
        np.testing.assert_allclose(cells[:, 0], 1e-3, rtol=1e-9)
        assert np.all(np.diff(cells, axis=1) > 0.0)
        np.testing.assert_array_equal(after[:, [0, -1]], before[:, [0, -1]])
        assert seam_gap(p.mesh.blocks) == 0.0
        assert min_jacobian(p.mesh.blocks) >= jac - 1e-9
        _, new_angle = wall_angles(p.mesh.blocks, radial)
        assert new_angle.max() < 3.0 and angle.max() < 3.0

    def test_l2_wing(self, tmp_path):
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
                   geometry="wing", nodes=25)  # fmt: skip
        p.orthogonalize()
        p.wall_layer(1e-4)
        cells = np.linalg.norm(np.diff(wall_lines(p.mesh.blocks), axis=1), axis=-1)
        np.testing.assert_allclose(cells[:, 0], 1e-4, rtol=1e-9)
        assert min_jacobian(p.mesh.blocks) > 0.3 and seam_gap(p.mesh.blocks) == 0.0
        pts, ang = wall_angles(p.mesh.blocks, normal)
        keep = (pts[:, 2] > 0.05) & (pts[:, 2] < 1.95) & (pts[:, 0] < 0.95)
        assert np.median(ang[keep]) < 5.0 and np.percentile(ang[keep], 95) < 12.0

    @pytest.mark.parametrize("case", ["zero", "too_large", "lattice"])
    def test_l3_refused(self, sphere_stl, case):
        p = sphere(sphere_stl, layered=False) if case == "lattice" else layered_sphere(sphere_stl)
        journal = p.journal()
        nodes = [b.nodes.copy() for b in p.mesh.blocks]
        with pytest.raises(ParameterError):
            p.wall_layer({"zero": 0.0, "too_large": 0.2, "lattice": 1e-3}[case])
        assert p.journal() == journal
        for b, old in zip(p.mesh.blocks, nodes, strict=True):
            np.testing.assert_array_equal(b.nodes, old)

    def test_l4_journal_undo(self, sphere_stl, tmp_path):
        from shore.project.project import run_journal

        p = layered_sphere(sphere_stl)
        nodes = [b.nodes.copy() for b in p.mesh.blocks]
        p.wall_layer(5e-3)
        p.save_journal(tmp_path / "journal.py")
        q = run_journal(tmp_path / "journal.py")
        for a, b in zip(p.mesh.blocks, q.mesh.blocks, strict=True):
            assert np.array_equal(a.nodes, b.nodes)
        p.undo()
        for b, old in zip(p.mesh.blocks, nodes, strict=True):
            np.testing.assert_array_equal(b.nodes, old)
