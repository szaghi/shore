"""FROZEN acceptance spec: the GUI's orthogonal-wall-lines form (ADR 0048).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job):
``pytest tests_gui/spec/test_gui_orthogonalize.py`` reports 2 passed, 0
skipped, 0 xfailed. At freeze: 0 passed, 2 xfailed.

Interface under test
--------------------
The ``orthogonalize`` form (``strength``) runs ``Project.orthogonalize`` on a
worker thread; the blocks are redrawn; errors go to the log and nothing
changes.

Cases
-----
A box blocking around the sphere cut in three along i, the middle block's
i faces associated with the sphere and the block deleted, premeshed.
- G1 the form: the journal line, the blocks' nodes changed.
- G2 a strength of 0: refused.

Out of scope of this spec: the reshaping itself
(``tests/spec/test_orthogonalize.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")


def window(qapp, sphere_stl):
    from shore.gui.window import MainWindow

    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    run = w.session.run
    run("load_geometry", str(sphere_stl))
    run("create_blocking", 7, 7, 7, vmin=(-3, -3, -3), vmax=(3, 3, 3))
    run("split_blocking", "b0", "i", 1 / 3)
    run("split_blocking", "b1", "i", 0.5)
    blocking = w.session.project.blocking
    centres = {lab: blocking.points[list(blocking.block(lab))].mean(axis=0)
               for lab in blocking.labels()}  # fmt: skip
    middle = next(lab for lab, c in centres.items() if np.isclose(c[0], 0.0))
    run("associate_blocking", middle, "i_lo")
    run("associate_blocking", middle, "i_hi")
    run("delete_blocks", [middle])
    run("premesh")
    return w


class TestOrthogonalizeForm:
    @PENDING
    def test_g1_form(self, qapp, sphere_stl):
        w = window(qapp, sphere_stl)
        before = [b.nodes.copy() for _, b in w.session.blocks()]
        w.run_form("orthogonalize", {"strength": "0.8"})
        w.wait_idle()
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == "p.orthogonalize(strength=0.8)"
        after = [b.nodes for _, b in w.session.blocks()]
        assert any(not np.array_equal(a, b) for a, b in zip(after, before, strict=True))

    @PENDING
    def test_g2_refused(self, qapp, sphere_stl):
        w = window(qapp, sphere_stl)
        lines = w.session.journal_lines()
        w.run_form("orthogonalize", {"strength": "0"})
        w.wait_idle()
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert w.session.journal_lines() == lines
