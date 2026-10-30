"""FROZEN acceptance spec: the wing-body blocking built in the GUI (ADR 0051).

Roadmap phase 7's acceptance, interactive half. Frozen 2026-10-30; no
implementation is pending. It may not be edited (``tests/spec/SHA256SUMS``,
``scripts/check_spec.py``).

Verdict for acceptance (GUI job, needs ``shore[gui,cad]``):
``pytest tests_gui/spec/test_gui_wing_body.py`` reports 2 passed, 0
skipped, 0 xfailed. At freeze: 2 passed.

The case
--------
``tests/spec/test_wing_body.py``'s half wing-body model and blocking, built
through the window: the geometry and the box from their forms, the cuts from
the split form, the hole's faces and the junction edges selected in the
Blocking dock (ADR 0044) and associated with the geometry selected in the
Association dock (the body, then the junction curves), the hole's blocks
selected and deleted, the hole's faces selected and O-gridded from the dock
(factor 0.7, 17 nodes), then the premesh (exact), orthogonal-lines (0.5)
and wall-layer (2e-3) forms.

Cases
-----
- I1 every action succeeds, and the session's journal, replayed headless,
  gives the session's grid bit for bit; the grid is valid, its seams
  closed, its walls on the body with first cells 2e-3.
- I2 the same journal replayed on the edited model (span 2.7, fuselage
  radius 0.55): valid, on the new body.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("build123d")

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")
ALONG = {"i": ("j", "k"), "j": ("i", "k"), "k": ("i", "j")}
HOLE = [(x, y, 0.4) for x in (-1.5, 0.5, 2.0) for y in (-0.55, 0.0, 0.55)] + [(0.5, 0.0, 2.0)]
BASE_EDGES = ("j_lo_k_lo", "j_hi_k_lo", "i_lo_k_lo", "i_hi_k_lo")


def wing_body(path, span=3.0, radius=0.5, length=5.0):
    from build123d import (
        Box,
        Cylinder,
        Line,
        Pos,
        Sphere,
        Spline,
        Wire,
        export_step,
        extrude,
        make_face,
    )

    from shore.profiles import naca

    body = Cylinder(radius, length, rotation=(0, 90, 0))
    body = body + Pos(-length / 2, 0, 0) * Sphere(radius) + Pos(length / 2, 0, 0) * Sphere(radius)
    sec = naca("0012", 40, closed_te=False)
    le = int(np.argmin(sec[:, 0]))
    upper = [(float(x), float(y), 0.0) for x, y in sec[: le + 1]]
    lower = [(float(x), float(y), 0.0) for x, y in sec[le:]]
    wire = Wire([Spline(*upper), Spline(*lower), Line(lower[-1], upper[0])])
    wing = extrude(make_face(wire), amount=span)
    export_step((body + wing) & (Pos(0, 0, 10) * Box(40, 40, 20)), str(path))
    return path


def block_at(blocking, point):
    for label in blocking.labels():
        pts = blocking.points[list(blocking.block(label))]
        if np.all(pts.min(axis=0) < point) and np.all(point < pts.max(axis=0)):
            return label
    raise AssertionError(f"no block holds {point}")


def face_nodes(nodes, face):
    return np.take(nodes, 0 if face.endswith("lo") else -1, axis="ijk".index(face[0]))


def checks(project):
    from shore.geometry.base import BODY
    from shore.quality import scaled_jacobian

    blocks = project.mesh.blocks
    assert min(float(scaled_jacobian(b.nodes).min()) for b in blocks) > 0.1
    owner = {id(f): b for b in blocks for f in b.faces()}
    g = project.geometries["wb"].geometry
    for b in blocks:
        for f in b.faces():
            if f.bc.name == "SHARED":
                theirs = face_nodes(owner[id(f.partner)].nodes, f.partner.name)
                (a0, flip0), (_, flip1) = f.axis_map.along
                t = theirs if ALONG[f.partner.name[0]].index(a0) == 0 else theirs.transpose(1, 0, 2)
                t = t[::-1] if flip0 else t
                t = t[:, ::-1] if flip1 else t
                assert np.array_equal(face_nodes(b.nodes, f.name), t)
            if f.bc.name == "WALL":
                w = face_nodes(b.nodes, f.name).reshape(-1, 3)
                assert float(g.project(w, BODY, exact=True).distances.max()) < 1e-8
                nxt = np.take(b.nodes, 1 if f.name.endswith("lo") else -2,
                              axis="ijk".index(f.name[0])).reshape(-1, 3)  # fmt: skip
                np.testing.assert_allclose(np.linalg.norm(nxt - w, axis=1), 2e-3, rtol=1e-9)


@pytest.fixture(scope="module")
def session_build(qapp, tmp_path_factory):
    """The blocking built through the window; (window, model path)."""
    from shore.geometry.base import EntityKind
    from shore.gui.session import SkeletonKey
    from shore.gui.window import MainWindow

    d = tmp_path_factory.mktemp("gui_wing_body")
    model = wing_body(d / "wb.step")
    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)

    def ok():
        w.wait_idle()
        assert w.last_error is None, w.log.toPlainText()

    w.run_form("load_geometry", {"path": str(model), "name": "wb"})
    ok()
    w.run_form("create_blocking", {"ni": "53", "nj": "33", "nk": "25",
                                   "vmin": "(-6, -4, 0)", "vmax": "(7, 4, 6)"})  # fmt: skip
    ok()
    corner = np.array([-5.99, -3.99, 0.01])
    for axis, cuts in (("i", (-3.0, -0.3, 1.3, 3.0)), ("j", (-0.8, -0.3, 0.3, 0.8)),
                       ("k", (0.8, 3.3))):  # fmt: skip
        dim = "ijk".index(axis)
        for c in cuts:
            point = corner.copy()
            point[dim] = c
            blocking = w.session.project.blocking
            label = block_at(blocking, point)
            pts = blocking.points[list(blocking.block(label))]
            a, z = pts.min(axis=0)[dim], pts.max(axis=0)[dim]
            w.run_form("split_blocking", {"block": label, "axis": axis,
                                          "fraction": repr(float((c - a) / (z - a)))})  # fmt: skip
            ok()
    blocking = w.session.project.blocking
    hole = [block_at(blocking, np.array(c)) for c in HOLE]
    owners: dict = {}
    for label in hole:
        for f in FACES:
            owners.setdefault(blocking.entity(label, f)[0], []).append((label, f))
    on_plane = {k for k in owners if np.allclose(blocking.points[list(k[1:])][:, 2], 0.0)}
    walls = {k: v[0] for k, v in owners.items() if len(v) == 1 and k not in on_plane}
    # the hole's faces, associated with the body
    w.blocking_mode.setCurrentText("face")
    w.clear_blocking_selection()
    for label, f in walls.values():
        w.select_skeleton(SkeletonKey("main", label, f), add=True)
    w.selected_geometry = ("wb", "body", [])
    w.associate_blocking_selection()
    ok()
    # the base of the wing's block, associated with the junction curves
    g = w.session.project.geometries["wb"].geometry
    curves = g.display().curves
    junction = [k for k in range(len(g.entities(EntityKind.CURVE)))
                if np.allclose(np.hypot(curves[k][:, 1], curves[k][:, 2]), 0.5, atol=1e-3)
                and np.all(np.abs(curves[k][:, 1]) < 0.07)]  # fmt: skip
    w.blocking_mode.setCurrentText("edge")
    w.clear_blocking_selection()
    for e in BASE_EDGES:
        w.select_skeleton(SkeletonKey("main", hole[-1], e), add=True)
    w.selected_geometry = ("wb", "curve", junction)
    w.associate_blocking_selection()
    ok()
    # delete the hole's blocks
    w.blocking_mode.setCurrentText("block")
    w.clear_blocking_selection()
    for label in hole:
        w.select_skeleton(SkeletonKey("main", label, ""), add=True)
    w.delete_selection()
    ok()
    # O-grid the hole's faces
    blocking = w.session.project.blocking
    w.blocking_mode.setCurrentText("face")
    w.clear_blocking_selection()
    for k, _ in w.session.skeleton_markers("face"):
        if blocking.entity(k.block, k.edge)[0] in walls:
            w.select_skeleton(k, add=True)
    w.ogrid_factor.setText("0.7")
    w.ogrid_nodes.setText("17")
    w._ogrid_from_panel()
    ok()
    w.run_form("premesh", {"exact": True})
    ok()
    w.run_form("orthogonalize", {"strength": "0.5"})
    ok()
    w.run_form("wall_layer", {"first_step": "2e-3"})
    ok()
    return w, model


class TestGuiWingBody:
    def test_i1_session_and_replay(self, session_build, tmp_path):
        from shore.project.project import run_journal

        w, _ = session_build
        project = w.session.project
        assert len(project.mesh.blocks) == 90
        checks(project)
        project.save_journal(tmp_path / "journal.py")
        q = run_journal(tmp_path / "journal.py")
        for a, b in zip(project.mesh.blocks, q.mesh.blocks, strict=True):
            assert a.label == b.label and np.array_equal(a.nodes, b.nodes)

    def test_i2_edited_model(self, session_build, tmp_path):
        from shore.project.project import run_journal

        w, model = session_build
        (tmp_path / "in").mkdir()
        wing_body(tmp_path / "in" / "wb.step", span=2.7, radius=0.55)
        journal = w.session.project.journal().replace(repr(str(model)), repr("wb.step"))
        (tmp_path / "journal.py").write_text(journal)
        q = run_journal(tmp_path / "journal.py", input_dir=tmp_path / "in")
        assert len(q.mesh.blocks) == 90
        checks(q)
