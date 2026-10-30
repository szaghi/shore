"""FROZEN acceptance spec: the GUI's wall-layer form (ADR 0050).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job):
``pytest tests_gui/spec/test_gui_wall_layer.py`` reports 2 passed, 0
skipped, 0 xfailed. At freeze: 0 passed, 2 xfailed.

Interface under test
--------------------
The ``wall_layer`` form (``first_step``, positional) runs
``Project.wall_layer`` on a worker thread; errors go to the log and nothing
changes.

Cases
-----
A 3 x 3 x 3 lattice around the sphere, its centre deleted with its faces
associated, the hole layered (ADR 0043), premeshed.
- W1 the form: the journal line, the first cells.
- W2 a first step too large for the layer: refused.

Out of scope of this spec: the layer itself (``tests/spec/test_wall_layer.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")


def window(qapp, sphere_stl):
    from shore.gui.window import MainWindow

    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    run = w.session.run
    run("load_geometry", str(sphere_stl))
    run("create_blocking", 13, 13, 13, vmin=(-3, -3, -3), vmax=(3, 3, 3))
    splits = (("b0", "i", 1 / 3), ("b1", "i", 0.5), ("b0", "j", 1 / 3), ("b3", "j", 0.5),
              ("b0", "k", 1 / 3), ("b9", "k", 0.5))  # fmt: skip
    for label, axis, fraction in splits:
        run("split_blocking", label, axis, fraction)
    blocking = w.session.project.blocking
    centres = {lab: blocking.points[list(blocking.block(lab))].mean(axis=0)
               for lab in blocking.labels()}  # fmt: skip
    centre = next(lab for lab, c in centres.items() if np.allclose(c, 0.0))
    for f in FACES:
        run("associate_blocking", centre, f)
    hole = {blocking.entity(centre, f)[0] for f in FACES}
    run("delete_blocks", [centre])
    blocking = w.session.project.blocking
    faces = [(lab, f) for lab in blocking.labels() for f in FACES
             if blocking.entity(lab, f)[0] in hole]  # fmt: skip
    run("ogrid", faces=faces, factor=0.5, nodes=9)
    run("premesh")
    return w


class TestWallLayerForm:
    def test_w1_form(self, qapp, sphere_stl):
        w = window(qapp, sphere_stl)
        w.run_form("wall_layer", {"first_step": "0.01"})
        w.wait_idle()
        assert w.last_error is None
        assert w.session.journal_lines()[-1] == "p.wall_layer(0.01)"
        for _, b in w.session.blocks():
            walls = [f.name for f in b.faces() if f.bc.name == "WALL"]
            if len(walls) == 1:
                axis = "ijk".index(walls[0][0])
                lo = walls[0].endswith("lo")
                first = np.take(b.nodes, 1 if lo else -2, axis=axis) - np.take(
                    b.nodes, 0 if lo else -1, axis=axis
                )
                np.testing.assert_allclose(np.linalg.norm(first, axis=-1), 0.01, rtol=1e-9)

    def test_w2_refused(self, qapp, sphere_stl):
        w = window(qapp, sphere_stl)
        lines = w.session.journal_lines()
        w.run_form("wall_layer", {"first_step": "1.0"})
        w.wait_idle()
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert w.session.journal_lines() == lines
