"""FROZEN acceptance spec: the GUI's O-grid form (ADR 0043).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job): ``pytest tests_gui/spec/test_gui_ogrid.py``
reports 3 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 3 xfailed.

Interface under test
--------------------
The ``ogrid`` form: ``blocks`` (a literal list of labels, positional),
``faces`` and ``exclude`` (literal lists of ``(block, face)`` pairs),
``factor``, ``nodes``. It runs ``Project.ogrid``; the skeleton (ADR 0042) is
redrawn. Errors go to the log and nothing changes.

Cases
-----
A blocking box (2 x 2 x 2, 9 nodes per edge) made from its form.
- W1 an O-grid of the block: the journal line, 32 skeleton edges.
- W2 a half-O (an open face): the journal line, 32 skeleton edges (every
  corner lies on a layered side face).
- W3 a bad factor: refused.

Out of scope of this spec: picking blocks and faces (a later slice), the
O-grid itself (``tests/spec/test_ogrid.py``).
"""

from __future__ import annotations

import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

BOX = {"ni": "9", "nj": "9", "nk": "9", "vmin": "(-1, -1, -1)", "vmax": "(1, 1, 1)"}


def window(qapp):
    from shore.gui.window import MainWindow

    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    w.run_form("create_blocking", BOX)
    w.wait_idle()
    return w


class TestOgridForm:
    @PENDING
    def test_w1_ogrid(self, qapp):
        w = window(qapp)
        w.run_form("ogrid", {"blocks": "['b0']", "factor": "0.6", "nodes": "4"})
        w.wait_idle()
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == "p.ogrid(('b0',), factor=0.6, nodes=4)"
        assert len(w.viewport.skeleton_keys) == 32

    @PENDING
    def test_w2_half(self, qapp):
        w = window(qapp)
        w.run_form("ogrid", {"blocks": "['b0']", "exclude": "[('b0', 'k_lo')]"})
        w.wait_idle()
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == "p.ogrid(('b0',), exclude=(('b0', 'k_lo'),))"
        assert len(w.viewport.skeleton_keys) == 32

    @PENDING
    def test_w3_refused(self, qapp):
        w = window(qapp)
        lines = w.session.journal_lines()
        w.run_form("ogrid", {"blocks": "['b0']", "factor": "1.5"})
        w.wait_idle()
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert w.session.journal_lines() == lines and len(w.viewport.skeleton_keys) == 12
