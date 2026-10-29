"""FROZEN acceptance spec: moving blocking vertices from the GUI (ADR 0045).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job):
``pytest tests_gui/spec/test_gui_move_vertex.py`` reports 5 passed, 0
skipped, 0 xfailed. At freeze: 0 passed, 5 xfailed.

Interface under test
--------------------
- The Blocking dock's selection mode ``"vertex"``: a marker per blocking
  vertex (once, named ``v0``..``v7`` by its first block,
  ``Session.skeleton_markers("vertex")``); a pick selects it (ADR 0044).
- ``move_selection(xyz)`` (the panel: ``move_xyz`` as ``"x, y, z"``,
  ``_move_from_panel``) runs ``Project.move_vertex`` on the one selected
  vertex.
- ``move_selection_to_pick()`` moves it to the last picked geometry point
  (``last_geometry_point``, set by a geometry pick).
- Errors go to the log and nothing changes.

Cases
-----
A blocking row (two blocks, 2 x 1 x 1).
- X1 vertex markers.
- X2 a picked vertex moved from the panel: the journal line, the skeleton
  redrawn.
- X3 a vertex moved to a picked geometry point.
- X4 (x 2) refusals: two entities selected; a bad coordinate text.

Out of scope of this spec: mouse dragging, the move itself
(``tests/spec/test_move_vertex.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

ROW = {"ni": "9", "nj": "5", "nk": "5", "vmin": "(3, 3, 3)", "vmax": "(5, 4, 4)"}


def window(qapp):
    from shore.gui.window import MainWindow

    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    w.run_form("create_blocking", ROW)
    w.wait_idle()
    w.session.run("split_blocking", "b0", "i", 0.5)
    w.viewport.show_all()
    w.blocking_mode.setCurrentText("vertex")
    return w


def key(block, entity):
    from shore.gui.session import SkeletonKey

    return SkeletonKey("main", block, entity)


def pick_vertex(w, k):
    w.viewport.draw_now()
    x, y = w.viewport.project(dict(w.session.skeleton_markers("vertex"))[k])
    return w.pick_at(x, y)


class TestMoveVertex:
    @PENDING
    def test_x1_markers(self, qapp):
        w = window(qapp)
        keys = w.viewport.marker_keys
        assert len(keys) == 12 and keys == [k for k, _ in w.session.skeleton_markers("vertex")]
        assert key("b0", "v7") in keys and key("b1", "v3") not in keys  # once, first block
        np.testing.assert_array_equal(dict(w.session.skeleton_markers("vertex"))[key("b0", "v7")],
                                      [4.0, 4.0, 4.0])  # fmt: skip

    @PENDING
    def test_x2_move_from_panel(self, qapp):
        w = window(qapp)
        pick = pick_vertex(w, key("b0", "v7"))
        assert pick.key == key("b0", "v7") and w.blocking_selection == [key("b0", "v7")]
        w.move_xyz.setText("4.2, 4.3, 4.1")
        w._move_from_panel()
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == "p.move_vertex('b0', 'v7', (4.2, 4.3, 4.1))"
        markers = dict(w.session.skeleton_markers("vertex"))
        assert markers[key("b0", "v7")].tolist() == [4.2, 4.3, 4.1]
        uploaded = w.viewport.uploaded_positions("skeleton")
        assert np.any(np.all(np.abs(uploaded - [4.2, 4.3, 4.1]) < 1e-5, axis=1))

    @PENDING
    def test_x3_move_to_geometry_pick(self, qapp, sphere_stl):
        from shore.gui.viewport import GeometryPick

        w = window(qapp)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.set_visible("skeleton", False)
        w.viewport.look_at_point(np.zeros(3), 6.0, view_dir=np.array([0.0, 0.0, -1.0]))
        w.viewport.draw_now()
        x, y = w.viewport.project(np.array([0.0, 0.0, 1.0]))
        pick = w.pick_at(x, y)
        assert isinstance(pick, GeometryPick)
        w.select_skeleton(key("b1", "v6"))
        w.move_selection_to_pick()
        assert w.last_error is None
        moved = w.session.project.blocking.points[w.session.project.blocking.block("b1")[6]]
        np.testing.assert_array_equal(moved, pick.xyz)

    @PENDING
    def test_x4_two_selected(self, qapp):
        w = window(qapp)
        lines = w.session.journal_lines()
        w.select_skeleton(key("b0", "v7"))
        w.select_skeleton(key("b0", "v0"), add=True)
        w.move_selection((0.0, 0.0, 0.0))
        assert w.last_error is not None and w.session.journal_lines() == lines

    @PENDING
    def test_x4_bad_text(self, qapp):
        w = window(qapp)
        lines = w.session.journal_lines()
        w.select_skeleton(key("b0", "v7"))
        w.move_xyz.setText("1, 2")
        w._move_from_panel()
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert w.session.journal_lines() == lines
