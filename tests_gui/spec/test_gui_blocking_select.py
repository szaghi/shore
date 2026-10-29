"""FROZEN acceptance spec: selecting blocking faces and blocks in the GUI (ADR 0044).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job):
``pytest tests_gui/spec/test_gui_blocking_select.py`` reports 9 passed, 0
skipped, 0 xfailed. At freeze: 0 passed, 9 xfailed.

Interface under test
--------------------
- The Blocking dock's selection mode ``blocking_mode``: ``"edge"`` (the
  skeleton, ADR 0042), ``"face"`` or ``"block"``. In face and block mode the
  viewport draws a marker per blocking face (named by its first block) or
  block, at its skeleton centre (``Session.skeleton_markers(kind)``,
  ``viewport.marker_keys``), over everything; the Display dock's blocking
  toggle hides them with the skeleton.
- A pick on a marker returns ``MarkerPick(key, xyz)``, with ``key`` a
  ``SkeletonKey(component, block, face)`` (``""`` for a block). A pick on a
  marker or a skeleton edge selects it: ``blocking_selection`` becomes
  ``[key]``; with shift (``add=True``) the key is added, or removed if
  already selected. ``select_skeleton(key, add=)`` does the same without a
  pick; ``clear_blocking_selection()`` empties it.
- ``ogrid_selection(factor, nodes)`` (the panel: ``ogrid_factor``,
  ``ogrid_nodes``, ``_ogrid_from_panel``) runs ``Project.ogrid``:
  - the selected blocks are the set (every block when none is selected);
  - the selected faces are the layered faces, each named by its block in the
    set.
- ``delete_selection()`` runs ``Project.delete_blocks`` on the selected
  blocks.
- ``associate_blocking_selection()`` runs ``Project.associate_blocking`` for
  each selected face or edge, with the Association dock's selected geometry
  (``selected_geometry``, ADR 0039).
- Errors go to the log and nothing changes.

Cases
-----
- K1 markers per mode; hidden with the blocking.
- K2 (x 2) picking markers: select, add, toggle, replace; a skeleton edge
  pick selects too.
- K3 an O-grid of a picked block from the panel.
- K4 the hole of a 3 x 3 x 3 lattice: its six faces selected, one O-grid.
- K5 delete and associate the selection.
- K6 (x 3) refusals: nothing selected; faces given to delete; no geometry.

Out of scope of this spec: rubber-band selection, selecting through hidden
blocks, the O-grid itself (``tests/spec/test_ogrid.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

ROW = {"ni": "9", "nj": "5", "nk": "5", "vmin": "(0, 0, 0)", "vmax": "(2, 1, 1)"}
FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")


def window(qapp):
    from shore.gui.window import MainWindow

    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    w.run_form("create_blocking", ROW)
    w.wait_idle()
    w.session.run("split_blocking", "b0", "i", 0.5)
    w.viewport.show_all()
    return w


def pick_marker(w, key, add=False):
    w.viewport.draw_now()
    markers = dict(w.session.skeleton_markers("block" if key.edge == "" else "face"))
    x, y = w.viewport.project(markers[key])
    return w.pick_at(x, y, add=add)


def key(block, entity=""):
    from shore.gui.session import SkeletonKey

    return SkeletonKey("main", block, entity)


class TestSelect:
    def test_k1_markers(self, qapp):
        w = window(qapp)
        assert w.viewport.marker_keys == []  # edge mode
        w.blocking_mode.setCurrentText("face")
        assert len(w.viewport.marker_keys) == 11
        assert w.viewport.marker_keys == [k for k, _ in w.session.skeleton_markers("face")]
        w.blocking_mode.setCurrentText("block")
        assert w.viewport.marker_keys == [key("b0"), key("b1")]
        np.testing.assert_allclose(dict(w.session.skeleton_markers("block"))[key("b1")],
                                   [1.5, 0.5, 0.5], atol=1e-15)  # fmt: skip
        w.set_visible("skeleton", False)
        w.viewport.draw_now()
        markers = dict(w.session.skeleton_markers("block"))
        x, y = w.viewport.project(markers[key("b1")])
        pick = w.pick_at(x, y)
        assert pick is None or type(pick).__name__ == "PickResult"

    def test_k2_pick_blocks(self, qapp):
        from shore.gui.viewport import MarkerPick

        w = window(qapp)
        w.blocking_mode.setCurrentText("block")
        pick = pick_marker(w, key("b1"))
        assert isinstance(pick, MarkerPick) and pick.key == key("b1")
        assert w.blocking_selection == [key("b1")]
        pick_marker(w, key("b0"), add=True)
        assert w.blocking_selection == [key("b1"), key("b0")]
        pick_marker(w, key("b1"), add=True)  # toggles off
        assert w.blocking_selection == [key("b0")]
        pick_marker(w, key("b1"))  # replaces
        assert w.blocking_selection == [key("b1")]
        w.clear_blocking_selection()
        assert w.blocking_selection == []

    def test_k2_pick_faces_and_edges(self, qapp):
        w = window(qapp)
        w.blocking_mode.setCurrentText("face")
        pick_marker(w, key("b0", "i_lo"))
        assert w.blocking_selection == [key("b0", "i_lo")]
        w.blocking_mode.setCurrentText("edge")
        w.viewport.draw_now()
        edges = dict(w.session.skeleton())
        edge = next(k for k in edges if (k.block, k.edge) == ("b1", "j_lo_k_hi"))
        x, y = w.viewport.project(edges[edge][2])
        w.pick_at(x, y, add=True)
        assert w.blocking_selection == [key("b0", "i_lo"), edge]

    def test_k3_ogrid_block(self, qapp):
        w = window(qapp)
        w.blocking_mode.setCurrentText("block")
        pick_marker(w, key("b1"))
        w.ogrid_factor.setText("0.6")
        w.ogrid_nodes.setText("4")
        w._ogrid_from_panel()
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == "p.ogrid(('b1',), factor=0.6, nodes=4)"
        assert len(w.session.project.blocking.labels()) == 8

    def test_k4_hole(self, qapp):
        from shore.gui.window import MainWindow

        w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
        run = w.session.run
        run("create_blocking", 13, 13, 13, vmin=(-3, -3, -3), vmax=(3, 3, 3))
        run("split_blocking", "b0", "i", 1 / 3)
        run("split_blocking", "b1", "i", 0.5)
        run("split_blocking", "b0", "j", 1 / 3)
        run("split_blocking", "b3", "j", 0.5)
        run("split_blocking", "b0", "k", 1 / 3)
        run("split_blocking", "b9", "k", 0.5)
        blocking = w.session.project.blocking
        centres = {lab: blocking.points[list(blocking.block(lab))].mean(axis=0)
                   for lab in blocking.labels()}  # fmt: skip
        centre = next(lab for lab, c in centres.items() if np.allclose(c, 0))
        hole = {blocking.entity(centre, f)[0] for f in FACES}
        run("delete_blocks", [centre])
        blocking = w.session.project.blocking
        w.blocking_mode.setCurrentText("face")
        for k, _ in w.session.skeleton_markers("face"):
            if blocking.entity(k.block, k.edge)[0] in hole:
                w.select_skeleton(k, add=True)
        assert len(w.blocking_selection) == 6
        w.ogrid_selection(0.5, 5)
        assert w.last_error is None
        line = w.session.journal_lines()[-1]
        assert line.startswith("p.ogrid((), faces=((") and line.count("_") == 6
        assert len(w.session.project.blocking.labels()) == 32

    def test_k5_delete_and_associate(self, qapp, sphere_stl):
        w = window(qapp)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.blocking_mode.setCurrentText("face")
        w.select_skeleton(key("b0", "i_lo"))
        w.select_skeleton(key("b1", "i_hi"), add=True)
        w.selected_geometry = ("body", "body", [])
        w.associate_blocking_selection()
        assert w.last_error is None
        assert w.session.journal_lines()[-2:] == [
            "p.associate_blocking('b0', 'i_lo', 'body', 'body')",
            "p.associate_blocking('b1', 'i_hi', 'body', 'body')",
        ]
        w.blocking_mode.setCurrentText("block")
        pick_marker(w, key("b0"))
        w.delete_selection()
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == "p.delete_blocks(('b0',))"
        assert w.blocking_selection == []
        assert w.viewport.marker_keys == [key("b1")]

    @pytest.mark.parametrize("case", ["empty", "faces_to_delete", "no_geometry"])
    def test_k6_refused(self, qapp, case):
        w = window(qapp)
        lines = w.session.journal_lines()
        if case == "empty":
            w.ogrid_selection()
            assert w.last_error is not None
            w.delete_selection()
        elif case == "faces_to_delete":
            w.select_skeleton(key("b0", "i_lo"))
            w.delete_selection()
        else:
            w.select_skeleton(key("b0", "i_lo"))
            w.selected_geometry = None
            w.associate_blocking_selection()
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert w.session.journal_lines() == lines
