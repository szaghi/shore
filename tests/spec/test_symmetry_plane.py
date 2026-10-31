"""FROZEN acceptance spec: symmetry planes of a blocking (ADR 0053).

Frozen 2026-10-30, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (CAD job, needs ``shore[cad]``):
``pytest tests/spec/test_symmetry_plane.py`` reports 7 passed, 0 skipped, 0
xfailed. At freeze: 1 passed, 6 xfailed.

Interface under test
--------------------
``p.symmetry_plane(point, normal)``, the command ``SymmetryPlane``
(journaled, undoable; journal ``p.symmetry_plane((0.0, 0.0, 0.0), (0.0,
0.0, 1.0))``), on the current component's blocking: declares the plane
through *point* normal to *normal*. ``premesh`` then makes every boundary
face without an association whose corners lie on the plane a ``SYMMETRY``
face; ``orthogonalize`` keeps the lines leaving the wall from a node on the
plane in the plane, so the faces on it stay on it. ``ParameterError`` for a
zero normal, a point or normal that is not three numbers, or no blocking
(project and journal unchanged). Without a declared plane nothing changes.

The case
--------
ADR 0051's half wing-body (``tests/spec/test_wing_body.py``: the model, the
cuts, the hole, the associations, O-grid, premesh exact, orthogonalize 0.5,
wall layer 2e-3), with the plane z = 0 declared right after
``create_blocking``. ADR 0051's blocking left those faces up to 0.026 off
the plane after ``orthogonalize``.

Cases
-----
- Y1 after premesh the faces on z = 0 are SYMMETRY, no other face is, and
  no FREE face lies on the plane.
- Y2 after orthogonalize and the wall layer every node of a SYMMETRY face
  is on the plane (|z| <= 1e-12).
- Y3 the rest holds as in ADR 0051: valid (min scaled Jacobian > 0.1),
  seams closed bit for bit, walls on the body (1e-8), first cells 2e-3.
- Y4 the outputs: CGNS ``BCSymmetryPlane`` on every SYMMETRY face; the
  face graph (``adjacency.json``) records them (``cc.par`` maps them to
  natural code 2, ADR 0036).
- Y5 journal line, replay bit for bit, undo (the blocking has no plane).
- Y6 refusals: zero normal, two numbers, no blocking.
- Y7 without the declaration: no SYMMETRY face (the faces on z = 0 FREE),
  as before ADR 0053 (a guard; holds at freeze).

Out of scope of this spec: planes through associated faces, a body that
crosses the plane, quarter models (two planes; allowed, not checked), the
GUI form beyond its existence.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("build123d")

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-30; implementation pending")

_spec = importlib.util.spec_from_file_location(
    "wing_body_spec", Path(__file__).parent / "test_wing_body.py"
)
wb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(wb)

PLANE = ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))


def build(step, plane=True):
    """ADR 0051's journal with the symmetry plane declared after create_blocking."""
    from shore.project import Project

    original = Project.create_blocking

    def create_then_plane(self, *args, **kwargs):
        original(self, *args, **kwargs)
        if plane:
            self.symmetry_plane(*PLANE)

    Project.create_blocking = create_then_plane
    try:
        return wb.build(step)
    finally:
        Project.create_blocking = original


@pytest.fixture(scope="module")
def step(tmp_path_factory):
    return wb.wing_body(tmp_path_factory.mktemp("sym") / "wb.step")


@pytest.fixture(scope="module")
def built(step):
    return build(step)


def on_plane_faces(blocks):
    """(block, face name) of faces whose nodes all have z == 0 within 1e-9."""
    return [(b, f) for b in blocks for f in wb.FACES
            if np.abs(wb.face_nodes(b.nodes, f)[..., 2]).max() <= 1e-9]  # fmt: skip


def bc(block, name):
    return block.face(name).bc.name


class TestSymmetryPlane:
    def test_y1_symmetry_faces(self, step):
        p = build(step)
        p.undo()  # wall_layer
        p.undo()  # orthogonalize: the premeshed blocking
        blocks = p.mesh.blocks
        flat = on_plane_faces(blocks)
        assert len(flat) >= 9
        sym = [(b, f.name) for b in blocks for f in b.faces() if f.bc.name == "SYMMETRY"]
        assert {(id(b), f) for b, f in sym} == {(id(b), f) for b, f in flat}
        assert all(bc(b, f) != "FREE" for b, f in flat)

    def test_y2_stays_on_plane(self, built):
        worst = 0.0
        for b in built.mesh.blocks:
            for f in b.faces():
                if f.bc.name == "SYMMETRY":
                    worst = max(worst, float(np.abs(wb.face_nodes(b.nodes, f.name)[..., 2]).max()))
        assert worst <= 1e-12

    def test_y3_quality_as_adr_0051(self, built):
        from shore.geometry.base import BODY

        blocks = built.mesh.blocks
        assert len(blocks) == 90
        assert wb.min_jacobian(blocks) > 0.1
        assert wb.seam_gap(blocks) == 0.0
        g = built.geometries["wb"].geometry
        for b in blocks:
            for f in b.faces():
                if f.bc.name == "WALL":
                    w = wb.face_nodes(b.nodes, f.name).reshape(-1, 3)
                    assert float(g.project(w, BODY, exact=True).distances.max()) < 1e-8
                    nxt = wb.next_nodes(b.nodes, f.name).reshape(-1, 3)
                    np.testing.assert_allclose(np.linalg.norm(nxt - w, axis=1), 2e-3, rtol=1e-9)

    def test_y4_outputs(self, built, tmp_path):
        import CGNS.MAP

        n_sym = sum(f.bc.name == "SYMMETRY" for b in built.mesh.blocks for f in b.faces())
        assert n_sym >= 9
        built.write_cgns(tmp_path / "wb.cgns")
        found = []

        def walk(node):
            if node[3] == "BC_t":
                found.append(node[1].tobytes().decode())
            for child in node[2]:
                walk(child)

        walk(CGNS.MAP.load(str(tmp_path / "wb.cgns"))[0])
        assert found.count("BCSymmetryPlane") == n_sym
        built.write_adjacency(tmp_path / "wb.adjacency.json")
        from shore.io.adjacency import read_adjacency_json

        faces = read_adjacency_json(tmp_path / "wb.adjacency.json")
        assert str(faces).count("SYMMETRY") >= n_sym

    def test_y5_journal_replay_undo(self, built, tmp_path):
        from shore.project.project import run_journal

        journal = built.journal()
        assert "p.symmetry_plane((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))" in journal
        built.save_journal(tmp_path / "j.py")
        q = run_journal(tmp_path / "j.py")
        for a, b in zip(built.mesh.blocks, q.mesh.blocks, strict=True):
            assert a.label == b.label and np.array_equal(a.nodes, b.nodes)

        from shore.project import Project

        r = Project()
        r.create_blocking(5, 5, 5, vmin=(0, 0, 0), vmax=(1, 1, 1))
        r.symmetry_plane(*PLANE)
        assert len(r.blocking.planes) == 1
        r.undo()
        assert r.blocking.planes == ()

    def test_y6_refusals(self):
        from shore.project import Project

        p = Project()
        with pytest.raises(ParameterError):
            p.symmetry_plane((0, 0, 0), (0, 0, 1))  # no blocking
        p.create_blocking(5, 5, 5, vmin=(0, 0, 0), vmax=(1, 1, 1))
        journal = p.journal()
        with pytest.raises(ParameterError):
            p.symmetry_plane((0, 0, 0), (0, 0, 0))
        with pytest.raises(ParameterError):
            p.symmetry_plane((0, 0), (0, 0, 1))
        assert p.journal() == journal

    def test_y7_unchanged_without_plane(self, step):
        p = wb.build(step)
        assert not [f for b in p.mesh.blocks for f in b.faces() if f.bc.name == "SYMMETRY"]
        p.undo()  # wall_layer
        p.undo()  # orthogonalize: the premeshed blocking, its faces on the plane
        flat = on_plane_faces(p.mesh.blocks)
        assert len(flat) >= 9
        assert all(bc(b, f) == "FREE" for b, f in flat)
