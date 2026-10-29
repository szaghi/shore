"""FROZEN acceptance spec: GUI slice 3, association editing (ADR 0039).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job, ``.[dev,gui,cad]``):
``pytest tests_gui/spec/test_gui_associate.py`` reports 9 passed, 0 skipped,
0 xfailed. At freeze: 0 passed, 9 xfailed.

Interface under test
--------------------
- ``MainWindow.selection_mode`` (a combo box: ``"face"``, ``"edge"``,
  ``"vertex"``) and ``pick_at(x, y, add=False)``, what a double click (shift
  for *add*) does:
  - on a block it selects the picked face, the face's nearest edge or its
    nearest corner (``selected_block = (BlockKey, "k_hi" | "j_lo_k_hi" |
    "v5")``);
  - on a geometry it selects the entity (``selected_geometry = (name,
    kind, [indices])``); *add* appends an entity of the same geometry and
    kind.
- ``associate_selection()`` runs ``Project.associate`` on the selection;
  ``clear_association()`` runs it with ``kind="none"``. If the block is in
  another component, ``use_component`` runs first. Errors go to the log;
  nothing changes.
- The model tree's third column shows each face's association (``"body"``,
  ``"part: surface 2, 5"``), and the associated edges and vertices follow
  the six faces as extra rows.
- ``shore.gui.viewport``: ``nearest_edge(nodes, face, xyz)``,
  ``nearest_vertex(nodes, face, xyz)``, ``edge_nodes``, ``VERTEX_IJK``.

Cases
-----
A STEP part (a 2 x 1 x 1 box with a hole) and a box primitive above it.
- U1 (x 3) selection modes on a box face: the face, the nearest edge, the
  nearest corner (computed here independently).
- U2 geometry picks select, and shift adds.
- U3 associate through the window: journal line, tree column, undo.
- U4 then snap through the form: the face on the picked surfaces.
- U5 refusals logged (a face with a curve, nothing selected), journal
  unchanged.
- U6 wall associations of a built topology show as ``"body"``.
- U7 a block of another component: ``use_component`` first.

Out of scope of this spec: mouse events, the look of highlights,
association of slices, multi-block selection.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

BOX = {"kind": "box",
       "params": "{'vmin': (-0.8, -0.4, 0.55), 'vmax': (0.8, 0.4, 0.9), 'ni': 5, 'nj': 4, "
                 "'nk': 3}"}  # fmt: skip
CUBED = {"topology": "cubed_sphere", "nk": "5", "ds": "0.05", "growth": "1.1",
         "params": "{'ni': 9, 'nj': 16}"}  # fmt: skip


@pytest.fixture(scope="module")
def part(tmp_path_factory) -> Path:
    import build123d as bd

    path = tmp_path_factory.mktemp("cad") / "part.step"
    bd.export_step(bd.Box(2.0, 1.0, 1.0) - bd.Cylinder(0.3, 2.0), str(path))
    return path


def window(qapp, part):
    from shore.gui.window import MainWindow

    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    w.run_form("load_geometry", {"path": str(part), "name": "part"})
    w.run_form("build_primitive", BOX)
    w.wait_idle()
    assert w.last_error is None
    return w


def only(w, kind):
    for k in ("block", "geometry", "slice"):
        w.set_visible(k, k == kind)
    w.viewport.show_all()
    w.viewport.draw_now()


def top_quad_pixel(w):
    """A k_hi quad near the j_lo edge and the v5 corner: its centre and pixel."""
    from shore.gui.viewport import quad_centre

    key, block = w.session.blocks()[0]
    centre = quad_centre(block.nodes, "k_hi", (3, 0, 1))
    return key, block, centre, w.viewport.project(centre)


def geometry_pixel(w, face_index):
    """Pixel of a camera-facing triangle of CAD face *face_index*."""
    vp = w.viewport
    d = w.session.project.geometries["part"].geometry.display()
    eye = vp.camera_position()
    for t in np.flatnonzero(d.triangle_entity == face_index):
        corners = d.points[d.triangles[t]]
        n = np.cross(corners[1] - corners[0], corners[2] - corners[0])
        c = corners.mean(axis=0)
        if n @ (eye - c) > 0.3 * np.linalg.norm(n) * np.linalg.norm(eye - c):
            return vp.project(c)
    raise AssertionError(f"face {face_index} not facing the camera")


def top_and_side(w):
    """(top surface z = 0.5, the side surface y = +0.5) of the part."""
    d = w.session.project.geometries["part"].geometry.display()
    found = {}
    for f in np.unique(d.triangle_entity):
        pts = d.points[d.triangles[d.triangle_entity == f]].reshape(-1, 3)
        if np.allclose(pts[:, 2], 0.5):
            found["top"] = int(f)
        if np.allclose(pts[:, 1], 0.5):
            found["side"] = int(f)
    return found["top"], found["side"]


# ---------------------------------------------------------------------------
# U1-U2: selection
# ---------------------------------------------------------------------------


class TestSelection:
    @pytest.mark.parametrize("mode", ["face", "edge", "vertex"])
    def test_u1_block_selection(self, qapp, part, mode):
        from shore.gui.viewport import VERTEX_IJK, edge_nodes

        w = window(qapp, part)
        only(w, "block")
        key, block, centre, (x, y) = top_quad_pixel(w)
        w.selection_mode.setCurrentText(mode)
        w.pick_at(x, y)
        if mode == "face":
            expected = "k_hi"
        elif mode == "edge":

            def dist(e):
                pts = edge_nodes(block.nodes, e)
                return min(np.linalg.norm(pts - centre, axis=1))  # coarse: node distance

            candidates = ("j_lo_k_hi", "j_hi_k_hi", "i_lo_k_hi", "i_hi_k_hi")
            expected = min(candidates, key=dist)
        else:
            corners = [v for v, ijk in enumerate(VERTEX_IJK) if ijk[2] == -1]
            nearest = min(
                corners, key=lambda v: np.linalg.norm(block.nodes[VERTEX_IJK[v]] - centre)
            )
            expected = f"v{nearest}"
        assert w.selected_block == (key, expected)

    def test_u2_geometry_selection(self, qapp, part):
        w = window(qapp, part)
        only(w, "geometry")
        top, side = top_and_side(w)
        w.pick_at(*geometry_pixel(w, top))
        assert w.selected_geometry == ("part", "surface", [top])
        w.pick_at(*geometry_pixel(w, side), add=True)
        assert w.selected_geometry == ("part", "surface", [top, side])
        w.pick_at(*geometry_pixel(w, side))
        assert w.selected_geometry == ("part", "surface", [side])


# ---------------------------------------------------------------------------
# U3-U7: associating
# ---------------------------------------------------------------------------


def face_row(w, component_index, block_index, face):
    block = w.tree.topLevelItem(component_index).child(block_index)
    rows = {block.child(i).text(0): block.child(i) for i in range(block.childCount())}
    return rows[face]


def select_top_face_and_surfaces(w):
    only(w, "block")
    _, _, _, (x, y) = top_quad_pixel(w)
    w.selection_mode.setCurrentText("face")
    w.pick_at(x, y)
    only(w, "geometry")
    top, side = top_and_side(w)
    w.pick_at(*geometry_pixel(w, top))
    w.pick_at(*geometry_pixel(w, side), add=True)
    return top, side


class TestAssociating:
    def test_u3_associate_and_undo(self, qapp, part):
        w = window(qapp, part)
        top, side = select_top_face_and_surfaces(w)
        lines = w.session.journal_lines()
        w.associate_selection()
        assert w.last_error is None
        assert w.session.journal_lines() == [
            *lines,
            f"p.associate('box', 'k_hi', 'part', 'surface', ({top}, {side}))",
        ]
        assert face_row(w, 0, 0, "k_hi").text(2) == f"part: surface {top}, {side}"
        w.undo()
        assert w.session.journal_lines() == lines
        assert face_row(w, 0, 0, "k_hi").text(2) == ""

    def test_u4_snap_after_associating(self, qapp, part):
        w = window(qapp, part)
        select_top_face_and_surfaces(w)
        w.associate_selection()
        w.run_form("snap", {})
        w.wait_idle()
        assert w.last_error is None and w.session.journal_lines()[-1] == "p.snap()"
        nodes = w.session.blocks()[0][1].nodes[:, :, -1]
        on_top = np.isclose(nodes[..., 2], 0.5, atol=1e-9)
        on_side = np.isclose(nodes[..., 1], 0.5, atol=1e-9)
        assert np.all(on_top | on_side)

    def test_u5_refusals(self, qapp, part):
        from shore.gui.window import MainWindow

        empty = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
        empty.associate_selection()
        assert empty.last_error is not None and "Error" in empty.log.toPlainText()
        w = window(qapp, part)
        only(w, "block")
        _, _, _, (x, y) = top_quad_pixel(w)
        w.pick_at(x, y)
        w.selected_geometry = ("part", "curve", [0])  # a face cannot take a curve
        lines = w.session.journal_lines()
        w.associate_selection()
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert w.session.journal_lines() == lines

    def test_u6_wall_associations_shown(self, qapp, sphere_stl):
        from shore.gui.window import MainWindow

        w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.run_form("build_topology", CUBED)
        w.wait_idle()
        for b in range(6):
            assert face_row(w, 0, b, "k_lo").text(2) == "body"
            assert face_row(w, 0, b, "k_hi").text(2) == ""

    def test_u7_other_component(self, qapp, part):
        w = window(qapp, part)
        top, side = select_top_face_and_surfaces(w)
        w.run_form("use_component", {"name": "other"})
        lines = w.session.journal_lines()
        w.associate_selection()
        assert w.last_error is None
        assert w.session.journal_lines()[len(lines) :] == [
            "p.use_component('main')",
            f"p.associate('box', 'k_hi', 'part', 'surface', ({top}, {side}))",
        ]
