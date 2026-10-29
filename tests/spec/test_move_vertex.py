"""FROZEN acceptance spec: moving blocking vertices (ADR 0045).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_move_vertex.py``
reports 9 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 9 xfailed.

Interface under test
--------------------
``p.move_vertex(block, vertex, xyz)``, the command ``MoveVertex`` (journaled,
undoable; journal ``p.move_vertex('b0', 'v7', (1.2, 1.3, 1.1))``): vertex
``v0``..``v7`` of blocking block *block* (ADR 0042) moves to *xyz*. When the
vertex has an association, its own or inherited (ADR 0042), it moves to the
projection of *xyz* on it instead, so an associated vertex slides on its
geometry. Every block sharing the vertex follows. ``ParameterError`` for a
bad block or vertex name or an *xyz* that is not three finite numbers;
project and journal are then unchanged.

Cases
-----
- M1 an unassociated vertex shared by two blocks: the point, the journal,
  premesh corners, seams closed.
- M2 an associated vertex (inherited from its face): onto the sphere's STL;
  premesh keeps it there.
- M3 the edited body replayed: the stored point is re-projected onto the new
  body by premesh.
- M4 undo; replay bit for bit.
- M5 (x 5) refusals.

Out of scope of this spec: dragging with the mouse (the GUI moves a selected
vertex to typed coordinates or a geometry pick,
``tests_gui/spec/test_gui_move_vertex.py``), moving edges or faces, keeping
blocks valid (a move may invert one; premesh reports nothing).
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")


def row():
    from shore.project import Project

    p = Project()
    p.create_blocking(9, 5, 5, vmin=(0, 0, 0), vmax=(2, 1, 1))
    p.split_blocking("b0", "i", 0.5)
    return p


def point(p, block, n):
    return p.blocking.points[p.blocking.block(block)[n]]


def on_sphere(sphere_stl):
    from shore.project import Project

    p = Project()
    p.load_geometry(sphere_stl)
    p.create_blocking(5, 5, 5, vmin=(-0.5, -0.5, -0.5), vmax=(0.5, 0.5, 0.5))
    p.associate_blocking("b0", "k_hi")  # v7 inherits it
    return p


def distance(path, pts):
    return float(
        trimesh.proximity.closest_point(trimesh.load(path), np.reshape(pts, (-1, 3)))[1].max()
    )


class TestMove:
    def test_m1_free_vertex(self):
        p = row()
        p.move_vertex("b0", "v7", (1.2, 1.3, 1.1))
        assert p.journal().splitlines()[-1] == "p.move_vertex('b0', 'v7', (1.2, 1.3, 1.1))"
        np.testing.assert_array_equal(point(p, "b0", 7), [1.2, 1.3, 1.1])
        np.testing.assert_array_equal(point(p, "b1", 3), [1.2, 1.3, 1.1])  # shared
        p.premesh()
        a, b = p.mesh.blocks
        np.testing.assert_array_equal(a.nodes[-1, -1, -1], [1.2, 1.3, 1.1])
        np.testing.assert_array_equal(b.nodes[0, -1, -1], [1.2, 1.3, 1.1])
        np.testing.assert_array_equal(a.nodes[-1], b.nodes[0])  # the shared face
        from shore.quality import scaled_jacobian

        assert min(float(scaled_jacobian(x.nodes).min()) for x in (a, b)) > 0.5

    def test_m2_associated_vertex(self, sphere_stl):
        p = on_sphere(sphere_stl)
        p.move_vertex("b0", "v7", (0.3, 0.9, 0.9))
        moved = point(p, "b0", 7)
        assert distance(sphere_stl, moved) < 1e-9
        assert np.linalg.norm(moved - np.array([0.3, 0.9, 0.9])) < 0.5
        p.premesh()
        np.testing.assert_allclose(p.mesh.blocks[0].nodes[-1, -1, -1], moved, atol=1e-12)

    def test_m3_edited_body(self, sphere_stl, tmp_path):
        from shore.project.project import run_journal

        p = on_sphere(sphere_stl)
        p.move_vertex("b0", "v7", (0.3, 0.9, 0.9))
        p.premesh()
        edited = trimesh.load(sphere_stl)
        edited.apply_scale(1.2)
        (tmp_path / "in").mkdir()
        edited.export(tmp_path / "in" / "sphere.stl")
        journal = p.journal().replace(repr(str(sphere_stl)), repr("sphere.stl"))
        (tmp_path / "journal.py").write_text(journal)
        q = run_journal(tmp_path / "journal.py", input_dir=tmp_path / "in")
        assert distance(tmp_path / "in" / "sphere.stl", q.mesh.blocks[0].nodes[-1, -1, -1]) < 1e-9

    def test_m4_undo_replay(self, tmp_path):
        from shore.project.project import run_journal

        p = row()
        before = point(p, "b0", 7).copy()
        p.move_vertex("b0", "v7", (1.2, 1.3, 1.1))
        p.premesh()
        p.save_journal(tmp_path / "journal.py")
        q = run_journal(tmp_path / "journal.py")
        for a, b in zip(p.mesh.blocks, q.mesh.blocks, strict=True):
            assert np.array_equal(a.nodes, b.nodes)
        p.undo()
        p.undo()
        np.testing.assert_array_equal(point(p, "b0", 7), before)


class TestRefusals:
    @pytest.mark.parametrize(
        "args",
        [
            ("b9", "v0", (0, 0, 0)),
            ("b0", "v8", (0, 0, 0)),
            ("b0", "i_lo", (0, 0, 0)),
            ("b0", "v0", (0, 0)),
            ("b0", "v0", (0, float("nan"), 0)),
        ],
    )
    def test_m5_refused(self, args):
        p = row()
        journal, points = p.journal(), p.blocking.points.copy()
        with pytest.raises(ParameterError):
            p.move_vertex(*args)
        assert p.journal() == journal
        np.testing.assert_array_equal(p.blocking.points, points)
