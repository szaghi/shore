"""FROZEN acceptance spec: picking STL feature curves in the GUI (ADR 0049).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job):
``pytest tests_gui/spec/test_gui_stl_curves.py`` reports 2 passed, 0
skipped, 0 xfailed. At freeze: 0 passed, 2 xfailed.

Interface under test
--------------------
An STL's feature curves and points (``StlGeometry.features``, ADR 0049) are
drawn and picked like a CAD model's (ADR 0038, 0039): a pick on a curve is a
``GeometryPick`` of ``EntityRef(CURVE, k)``, and the Association dock links
it to a block edge.

Cases
-----
A blunt-TE wing STL (tip cap, no root cap).
- U1 the curves drawn; a pick on the tip outline picks its curve.
- U2 that curve associated with a box edge from the Association dock: the
  journal line.

Out of scope of this spec: the features themselves
(``tests/spec/test_stl_curves.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")


def wing(tmp_path):
    from shore.profiles import naca
    from shore.shapes import extrude_section

    w = extrude_section(naca("0012", 60, closed_te=False), span=2.0, n_stations=21,
                        span_axis="z", end_caps=True)  # fmt: skip
    root = np.all(np.abs(w.vertices[w.faces][:, :, 2]) < 1e-12, axis=1)
    w.update_faces(~root)
    w.remove_unreferenced_vertices()
    w.export(tmp_path / "wing.stl")
    return tmp_path / "wing.stl"


def tip_outline(curves):
    """The tip cap's airfoil outline (z = 2, longer than a TE base segment)."""
    return next(k for k, c in enumerate(curves) if len(c) > 50 and np.allclose(c[:, 2], 2.0))


def pick_tip(qapp, tmp_path):
    from shore.gui.window import MainWindow

    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    w.run_form("load_geometry", {"path": str(wing(tmp_path)), "skip_hygiene": True})
    w.wait_idle()
    curves = w.session.project.geometries["body"].geometry.features()[0]
    k = tip_outline(curves)
    w.viewport.look_at_point(np.array([0.5, 0.0, 2.0]), 2.0, view_dir=np.array([0.0, 0.0, -1.0]))
    w.viewport.draw_now()
    target = curves[k][np.argmin(np.abs(curves[k][:, 0] - 0.5) + 10 * (curves[k][:, 1] < 0))]
    x, y = w.viewport.project(target)
    return w, k, w.pick_at(x, y)


class TestStlCurves:
    def test_u1_pick_curve(self, qapp, tmp_path):
        from shore.geometry.base import EntityKind, EntityRef
        from shore.gui.viewport import GeometryPick

        w, k, pick = pick_tip(qapp, tmp_path)
        d = w.session.project.geometries["body"].geometry.display()
        assert len(d.curves) == 6 and len(d.vertices) == 4
        assert isinstance(pick, GeometryPick)
        assert pick.entity == EntityRef(EntityKind.CURVE, k)

    def test_u2_associate(self, qapp, tmp_path):
        w, k, _ = pick_tip(qapp, tmp_path)
        box = "{'vmin': (1.0, 0.0, 0.1), 'vmax': (1.5, 0.5, 1.9), 'ni': 3, 'nj': 3, 'nk': 5}"
        w.run_form("build_primitive", {"kind": "box", "params": box})
        w.wait_idle()
        key = w.session.blocks()[0][0]
        w.selected_block = (key, "i_lo_j_lo")
        w.associate_selection()
        assert w.last_error is None
        assert (
            w.session.journal_lines()[-1]
            == f"p.associate('box', 'i_lo_j_lo', 'body', 'curve', ({k},))"
        )
