"""FROZEN acceptance spec: GUI spacing editor (ADR 0041).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job): ``pytest tests_gui/spec/test_gui_spacing.py``
reports 5 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 5 xfailed.

Interface under test
--------------------
- ``MainWindow.edge_info()``: for an edge selected with the selection mode
  ``"edge"`` (ADR 0039), ``{"key", "edge", "axis", "nodes", "length",
  "first", "last"}`` (the edge's arc length and its first and last cell),
  else ``None``.
- ``MainWindow.redistribute_selection(law, **params)`` runs
  ``Project.redistribute`` on it, making its component current first if
  needed; errors go to the log and nothing changes.
- The Spacing dock's widgets (``spacing_law``, ``spacing_beta``,
  ``spacing_cluster``, ``spacing_first``, ``spacing_last``,
  ``spacing_interior``) and its Redistribute button
  (``_redistribute_from_panel``) build the same command.

Cases
-----
A box primitive (2 x 1 x 1, 9 x 5 x 5 nodes).
- P1 the selected edge's info.
- P2 redistributing it: the journal line, the nodes, the updated info.
- P3 (x 2) refusals: no edge selected (a face), a bad number in the panel.
- P4 the panel's widgets.

Out of scope of this spec: mouse events, the redistribution itself
(``tests/spec/test_redistribute.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

BOX = {"kind": "box",
       "params": "{'vmin': (0, 0, 0), 'vmax': (2, 1, 1), 'ni': 9, 'nj': 5, 'nk': 5}"}  # fmt: skip


def window(qapp, mode="edge"):
    from shore.gui.viewport import quad_centre
    from shore.gui.window import MainWindow

    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    w.run_form("build_primitive", BOX)
    w.wait_idle()
    key, block = w.session.blocks()[0]
    w.viewport.show_all()
    w.viewport.draw_now()
    # a k_hi quad next to the j_lo edge (along i), at its middle
    x, y = w.viewport.project(quad_centre(block.nodes, "k_hi", (4, 0, 3)))
    w.selection_mode.setCurrentText(mode)
    w.pick_at(x, y)
    return w, key, block


class TestSpacing:
    @PENDING
    def test_p1_edge_info(self, qapp):
        w, key, _ = window(qapp)
        info = w.edge_info()
        assert info is not None and (info["key"], info["edge"], info["axis"]) == (
            key, "j_lo_k_hi", "i",
        )  # fmt: skip
        assert info["nodes"] == 9
        assert info["length"] == pytest.approx(2.0, abs=1e-12)
        assert info["first"] == pytest.approx(0.25, abs=1e-12)
        assert info["last"] == pytest.approx(0.25, abs=1e-12)

    @PENDING
    def test_p2_redistribute(self, qapp):
        from shore.mesh.spacing import Spacing1D

        w, _, _ = window(qapp)
        w.redistribute_selection("tanh", beta=3.0)
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == "p.redistribute('box', 'j_lo_k_hi', 'tanh')"
        x = w.session.blocks()[0][1].nodes[:, 0, -1, 0]
        np.testing.assert_allclose(x, Spacing1D("tanh", beta=3.0).build(9, 0.0, 2.0), atol=1e-12)
        assert w.edge_info()["first"] == pytest.approx(x[1] - x[0], abs=1e-12)

    @PENDING
    def test_p3_no_edge(self, qapp):
        w, _, _ = window(qapp, mode="face")
        lines = w.session.journal_lines()
        assert w.edge_info() is None
        w.redistribute_selection("tanh")
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert w.session.journal_lines() == lines

    @PENDING
    def test_p3_bad_number(self, qapp):
        w, _, _ = window(qapp)
        lines = w.session.journal_lines()
        w.spacing_law.setCurrentText("tanh")
        w.spacing_beta.setText("three")
        w._redistribute_from_panel()
        assert w.last_error is not None and w.session.journal_lines() == lines

    @PENDING
    def test_p4_panel(self, qapp):
        w, _, _ = window(qapp)
        w.spacing_law.setCurrentText("tanh2")
        w.spacing_beta.setText("2")
        w.spacing_cluster.setText("0.3")
        w._redistribute_from_panel()
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == (
            "p.redistribute('box', 'j_lo_k_hi', 'tanh2', beta=2.0, cluster_at=0.3)"
        )
