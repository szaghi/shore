"""FROZEN acceptance spec: GUI blocking skeleton and split (ADR 0042).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job): ``pytest tests_gui/spec/test_gui_blocking.py``
reports 6 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 6 xfailed.

Interface under test
--------------------
- ``Session.skeleton()``: every blocking edge of every component once, as
  ``(SkeletonKey(component, block, edge), nodes (n, 3))``, its *n* nodes
  equally spaced on the line between its skeleton vertices, in the edge's
  direction in ``block``.
- The viewport draws them over the blocks (``viewport.skeleton_keys``, the
  same order; ``uploaded_positions("skeleton")``: the nodes, edges separated
  by NaN rows) and picks them: ``SkeletonPick(key, node, nodes, xyz)``, the
  node nearest the pick. The Display dock's ``"skeleton"`` visibility hides
  them (``set_visible("skeleton", False)``).
- A pick selects the edge (``selected_skeleton``); ``split_selection()``
  runs ``Project.split_blocking(block, axis, node / (nodes - 1))`` for it,
  making its component current first if needed; errors go to the log and
  nothing changes.
- Forms: ``create_blocking``, ``split_blocking``, ``delete_blocks``,
  ``associate_blocking``, ``premesh``.

Cases
-----
A blocking box (2 x 1 x 1, 9 x 5 x 5 nodes) made from its form.
- G1 the skeleton drawn: 12 edges, their nodes uploaded.
- G2 a picked node and the split through it: journal line, 20 edges.
- G3 (x 2) refusals: nothing picked; an end node picked.
- G4 premesh: the blocks drawn, the skeleton picked over them, hidden.
- G5 the delete and associate forms.

Out of scope of this spec: mouse events, the blocking itself
(``tests/spec/test_blocking.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

BOX = {"ni": "9", "nj": "5", "nk": "5", "vmin": "(0, 0, 0)", "vmax": "(2, 1, 1)"}


def window(qapp):
    from shore.gui.window import MainWindow

    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    w.run_form("create_blocking", BOX)
    w.wait_idle()
    w.viewport.show_all()
    w.viewport.draw_now()
    return w


def pick_node(w, block, edge, node):
    nodes = dict(w.session.skeleton())
    key = next(k for k in nodes if (k.block, k.edge) == (block, edge))
    w.viewport.draw_now()
    x, y = w.viewport.project(nodes[key][node])
    return key, w.pick_at(x, y)


class TestBlocking:
    @PENDING
    def test_g1_skeleton(self, qapp):
        w = window(qapp)
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == (
            "p.create_blocking(9, 5, 5, vmin=(0.0, 0.0, 0.0), vmax=(2.0, 1.0, 1.0))"
        )
        skeleton = w.session.skeleton()
        assert len(skeleton) == 12 and w.viewport.skeleton_keys == [k for k, _ in skeleton]
        lengths = sorted(len(pts) for _, pts in skeleton)
        assert lengths == [5] * 8 + [9] * 4
        uploaded = w.viewport.uploaded_positions("skeleton")
        expected = np.concatenate([pts for _, pts in skeleton])
        np.testing.assert_allclose(uploaded[np.isfinite(uploaded).all(axis=1)], expected,
                                   atol=1e-6)  # fmt: skip

    @PENDING
    def test_g2_split_at_pick(self, qapp):
        from shore.gui.viewport import SkeletonPick

        w = window(qapp)
        key, pick = pick_node(w, "b0", "j_lo_k_hi", 2)
        assert isinstance(pick, SkeletonPick) and (pick.key, pick.node, pick.nodes) == (key, 2, 9)
        np.testing.assert_allclose(pick.xyz, [0.5, 0.0, 1.0], atol=1e-12)
        w.split_selection()
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == "p.split_blocking('b0', 'i', 0.25)"
        assert len(w.viewport.skeleton_keys) == 20

    @PENDING
    def test_g3_nothing_picked(self, qapp):
        w = window(qapp)
        lines = w.session.journal_lines()
        w.split_selection()
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert w.session.journal_lines() == lines

    @PENDING
    def test_g3_end_node(self, qapp):
        w = window(qapp)
        lines = w.session.journal_lines()
        _, pick = pick_node(w, "b0", "j_lo_k_hi", 0)  # a corner: any of its three edges
        assert pick.node in (0, pick.nodes - 1)
        w.split_selection()
        assert w.last_error is not None and w.session.journal_lines() == lines

    @PENDING
    def test_g4_premesh(self, qapp):
        from shore.gui.viewport import PickResult, SkeletonPick

        w = window(qapp)
        pick_node(w, "b0", "j_lo_k_hi", 4)
        w.split_selection()
        w.run_form("premesh", {})
        w.wait_idle()
        assert w.last_error is None and w.session.journal_lines()[-1] == "p.premesh()"
        assert [k.label for k in w.viewport.keys] == ["b0", "b1"]
        _, pick = pick_node(w, "b1", "j_lo_k_hi", 2)
        assert isinstance(pick, SkeletonPick)
        w.set_visible("skeleton", False)
        _, pick = pick_node(w, "b1", "j_lo_k_hi", 2)
        assert isinstance(pick, PickResult) and pick.key.label == "b1"

    @PENDING
    def test_g5_forms(self, qapp, sphere_stl):
        w = window(qapp)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.run_form("split_blocking", {"block": "b0", "axis": "i", "fraction": "0.5"})
        w.run_form("associate_blocking", {"block": "b1", "entity": "i_hi", "geometry": "body",
                                          "kind": "body", "indices": "()"})  # fmt: skip
        w.run_form("delete_blocks", {"blocks": "['b0']"})
        assert w.last_error is None
        assert w.session.journal_lines()[-3:] == [
            "p.split_blocking('b0', 'i', 0.5)",
            "p.associate_blocking('b1', 'i_hi', 'body', 'body')",
            "p.delete_blocks(('b0',))",
        ]
        assert w.session.project.blocking.labels() == ["b1"]
