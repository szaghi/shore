"""FROZEN acceptance spec: GUI slice 1, viewer and command shell (ADR 0037).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The
ONE permitted edit is deleting an ``@PENDING`` line once that case passes
(strict xfail turns a passing pending case into a failure).
``tests/spec/SHA256SUMS`` holds the hashes, taken with ``@PENDING`` lines
removed; ``scripts/check_spec.py`` checks them.

Verdict for acceptance (GUI CI job, ``QT_QPA_PLATFORM=offscreen``, software
Vulkan adapter): ``pytest tests_gui/spec/test_gui_shell.py`` reports 16
passed, 0 skipped, 0 xfailed. At freeze: 1 passed (G11b, the core imports
without the GUI stack: a guard), 15 xfailed.

Interface under test
--------------------
- ``shore.gui.session.Session``: one ``Project``. ``run(method, *args,
  **kwargs)`` calls ``Project.<method>`` (one command, one journal line);
  ``undo``, ``redo``, ``open(dir, replay=)``, ``save(dir)``, ``new``;
  ``journal_lines()``; ``blocks()`` gives ``(BlockKey(component, label),
  block)`` for every component with a mesh; ``tree()``; ``subscribe``.
- ``shore.gui.forms``: ``FORMS``, one ``Form`` per command, whose
  ``arguments(values)`` builds ``Project.<method>``'s arguments.
- ``shore.gui.window.MainWindow(session=None, *, ask_path=..., offscreen_
  viewport=False)``:
  - ``run_form(method, values)``: what the Run button does;
  - long forms run on a worker thread; ``wait_idle()``; ``cancel()``;
  - ``undo()``, ``redo()``, ``open_project(dir, replay=)``,
    ``save_project(dir)``;
  - ``last_error``; ``journal`` (console widget), ``log``, ``tree``;
  - ``viewport``; ``quality_action``, with ``refresh_quality()``.
- ``shore.gui.viewport``: ``FACES``, ``PickResult(key, face, cell, xyz)``,
  ``quad_centre``, ``boundary_values``, ``colormap``.
  ``pygfx_backend.PygfxViewport``: ``set_blocks``, ``set_face_scalars``,
  ``show_all``, ``look_at_point``, ``project``, ``pick``, ``keys``,
  ``uploaded_vertices``, ``face_colours``, ``draw_now``; and
  ``create_canvas(offscreen=True)``.
- ``shore gui [DIR] [--no-replay] [--smoke]``: ``--smoke`` builds the
  window, prints ``shore gui smoke: N blocks, journal M lines`` and exits 0.

Cases
-----
- G1 journal fidelity: a session driven through the forms has one journal
  line per command, the console shows exactly the journal, and the project
  saved from the GUI replays headless to identical grids.
- G2 undo and redo through the window.
- G3 the model tree is the project's: components (current marked), blocks,
  shapes, face BCs.
- G4 the viewport holds every block of every component, its boundary nodes
  uploaded.
- G5 picks on camera-facing boundary cells return the block, face and cell,
  and the float64 centre.
- G6 the float32 trap: wall cells of 1e-6 at an offset of 1e3, picked after
  ``look_at_point``, are exact.
- G7 quality colouring: per boundary quad, the colour map of the scaled
  Jacobian of its cell.
- G8 (x 2) a march on the worker: UI lateness <= 100 ms; cancelling leaves
  mesh and journal unchanged.
- G9 (x 2) errors: a bad form value and a failing command are reported in
  the log; project and journal unchanged; the window still works.
- G10 (x 2) opening a saved two-component project, with and without
  replay (its journal names the saved input copies); ``shore gui DIR
  --smoke``.
- G11 (x 2) isolation: only ``shore/gui/viewport/pygfx_backend.py`` imports
  pygfx, wgpu or rendercanvas, and only ``shore/gui`` imports PySide6; the
  core (``shore.cli``, ``shore.gui.session``, ``shore.gui.forms``) imports
  without the GUI stack.
- G12 every form calls an existing Project method with arguments it accepts.

Out of scope of this spec: performance budgets (measured on native
Windows, ADR 0037), the look of the window, mouse interaction (orbit, double
click), macOS and native Linux, the forms' layout.
"""

from __future__ import annotations

import ast
import inspect
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

SRC = Path(__file__).resolve().parents[2] / "src" / "shore"
CUBED = {"topology": "cubed_sphere", "nk": "5", "ds": "0.05", "growth": "1.1",
         "params": "{'ni': 9, 'nj': 16}"}  # fmt: skip
FLAT_CAPS = {"kind": "flat_caps",
             "params": "{'r_out': 3.0, 'z0': -1.0, 'z1': 4.0, 'ni': 5, 'nc': 6, 'nk': 5, "
                       "'inner_frac': 0.4}"}  # fmt: skip


def window(qapp, **kw):
    from shore.gui.window import MainWindow

    return MainWindow(offscreen_viewport=True, ask_path=lambda *a: None, **kw)


def two_components(w, sphere_stl) -> list[str]:
    """Drive the forms through a two-component session; the journal after each step."""
    steps = []
    w.run_form("load_geometry", {"path": str(sphere_stl)})
    steps.append(list(w.session.journal_lines()))
    w.run_form("build_topology", CUBED)
    w.wait_idle()
    steps.append(list(w.session.journal_lines()))
    w.run_form("march", {})
    w.wait_idle()
    steps.append(list(w.session.journal_lines()))
    w.run_form("use_component", {"name": "background"})
    steps.append(list(w.session.journal_lines()))
    w.run_form("build_primitive", FLAT_CAPS)
    w.wait_idle()
    steps.append(list(w.session.journal_lines()))
    assert w.last_error is None, w.last_error
    return steps


def grids(project) -> dict[tuple[str, str], np.ndarray]:
    return {
        (name, b.label): b.nodes
        for name, mesh in project.components.items()
        if mesh is not None
        for b in mesh.blocks
    }


def same_grids(a, b) -> bool:
    ga, gb = grids(a), grids(b)
    return ga.keys() == gb.keys() and all(np.array_equal(ga[k], gb[k]) for k in ga)


# ---------------------------------------------------------------------------
# G1-G4: the session through the window
# ---------------------------------------------------------------------------


class TestShell:
    @PENDING
    def test_g1_journal_fidelity(self, qapp, sphere_stl, tmp_path):
        from shore.project import Project, run_journal

        w = window(qapp)
        steps = two_components(w, sphere_stl)
        for n, lines in enumerate(steps, start=1):
            assert len(lines) == n
            assert lines[: n - 1] == steps[n - 2] if n > 1 else True
        assert w.journal.toPlainText().splitlines() == w.session.journal_lines()
        w.save_project(str(tmp_path / "saved"))
        assert w.last_error is None
        replayed = Project.load(tmp_path / "saved")
        assert same_grids(w.session.project, replayed)
        journal = run_journal(tmp_path / "saved" / "journal.py", input_dir=tmp_path / "saved")
        assert same_grids(w.session.project, journal)

    @PENDING
    def test_g2_undo_redo(self, qapp, sphere_stl):
        w = window(qapp)
        two_components(w, sphere_stl)
        lines = w.session.journal_lines()
        keys = list(w.viewport.keys)
        w.undo()
        assert w.session.journal_lines() == lines[:-1]
        assert w.journal.toPlainText().splitlines() == lines[:-1]
        assert all(k.component == "main" for k in w.viewport.keys)
        w.redo()
        assert w.session.journal_lines() == lines
        assert w.viewport.keys == keys

    @PENDING
    def test_g3_model_tree(self, qapp, sphere_stl):
        w = window(qapp)
        two_components(w, sphere_stl)
        project = w.session.project
        top = [w.tree.topLevelItem(i) for i in range(w.tree.topLevelItemCount())]
        assert [t.text(0).split(" ")[0] for t in top] == list(project.components)
        assert [("(current)" in t.text(0)) for t in top] == [
            name == project.current_component for name in project.components
        ]
        for item, (_, mesh) in zip(top, project.components.items(), strict=True):
            blocks = [item.child(i) for i in range(item.childCount())]
            assert [b.text(0) for b in blocks] == [b.label for b in mesh.blocks]
            for child, block in zip(blocks, mesh.blocks, strict=True):
                assert child.text(1) == "x".join(map(str, block.nodes.shape[:3]))
                faces = {child.child(i).text(0): child.child(i).text(1) for i in range(6)}
                assert faces == {f.name: f.bc.name for f in block.faces()}

    @PENDING
    def test_g4_viewport_content(self, qapp, sphere_stl):
        from shore.gui.session import BlockKey

        w = window(qapp)
        two_components(w, sphere_stl)
        expected = [
            BlockKey(name, b.label)
            for name, mesh in w.session.project.components.items()
            if mesh is not None
            for b in mesh.blocks
        ]
        assert w.viewport.keys == expected
        for key, block in w.session.blocks():
            ni, nj, nk = block.nodes.shape[:3]
            assert w.viewport.uploaded_vertices(key) == 2 * (nj * nk + ni * nk + ni * nj)


# ---------------------------------------------------------------------------
# G5-G7: the viewport
# ---------------------------------------------------------------------------


def facing_cells(viewport, key, nodes, face, size):
    """Cells of *face* facing the camera squarely, inside the image."""
    from shore.gui.viewport import FACES, face_grids, quad_centre

    grid = face_grids(nodes)[FACES.index(face)]
    eye = np.asarray(viewport.camera.local.position, np.float64) + viewport._origin
    out = []
    for a in range(grid.shape[0] - 1):
        for b in range(grid.shape[1] - 1):
            quad = grid[a : a + 2, b : b + 2].reshape(-1, 3)
            n = np.cross(quad[3] - quad[0], quad[1] - quad[0])
            centre = quad.mean(axis=0)
            view = eye - centre
            cos = abs(n @ view) / (np.linalg.norm(n) * np.linalg.norm(view))
            cell = {"k_hi": (a, b, nodes.shape[2] - 2), "i_lo": (0, a, b)}[face]
            x, y = viewport.project(quad_centre(nodes, face, cell))
            if cos > 0.7 and 10 <= x < size[0] - 10 and 10 <= y < size[1] - 10:
                out.append((cell, x, y, centre))
    return out


class TestViewport:
    @PENDING
    def test_g5_picks(self, qapp, sphere_stl):
        w = window(qapp)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.run_form("build_topology", CUBED)
        w.wait_idle()
        w.run_form("march", {})
        w.wait_idle()
        vp = w.viewport
        vp.draw_now()
        size = tuple(int(v) for v in vp.renderer.logical_size)
        checked = 0
        eye = np.asarray(vp.camera.local.position, np.float64) + vp._origin
        for key, block in w.session.blocks():
            centre = block.nodes.reshape(-1, 3).mean(axis=0)
            if (eye - centre) @ centre <= 0:  # blocks on the far side are hidden
                continue
            for cell, x, y, xyz in facing_cells(vp, key, block.nodes, "k_hi", size):
                r = vp.pick(x, y)
                assert r is not None and (r.key, r.face, r.cell) == (key, "k_hi", cell)
                assert np.array_equal(r.xyz, xyz)
                checked += 1
        assert checked >= 6

    @PENDING
    def test_g6_float32_trap(self, qapp):
        from shore.gui.viewport.pygfx_backend import PygfxViewport, create_canvas

        # a slab at 1e3 with wall cells of 1e-6 along i (geometric growth 1.5)
        x = 1e3 + np.concatenate([[0.0], np.cumsum(1e-6 * 1.5 ** np.arange(30))])
        y = 1e3 + np.linspace(0.0, 1.0, 11)
        z = 1e3 + np.linspace(0.0, 1.0, 11)
        nodes = np.stack(np.meshgrid(x, y, z, indexing="ij"), axis=-1)
        vp = PygfxViewport(create_canvas(offscreen=True, size=(640, 480)))
        vp.set_blocks([("slab", nodes)])
        vp.show_all()
        vp.draw_now()
        for k in (4, 2):
            for i in (0, 1, 2):
                # j_lo face (y = 1e3): cells (i, 0, k), 1e-6 x 0.1 wide
                cell = (i, 0, k)
                quad = nodes[i : i + 2, 0, k : k + 2].reshape(-1, 3)
                centre = quad.mean(axis=0)
                vp.look_at_point(centre, distance=2e-5, view_dir=np.array([0.0, 1.0, 0.0]))
                vp.draw_now()
                px, py = vp.project(centre)
                r = vp.pick(px, py)
                assert r is not None and (r.face, r.cell) == ("j_lo", cell), (i, k, r)
                assert np.array_equal(r.xyz, centre)
        assert vp.rebase_count > 0

    @PENDING
    def test_g7_quality_colours(self, qapp, sphere_stl):
        from shore.gui.viewport import boundary_values, colormap
        from shore.quality import scaled_jacobian

        w = window(qapp)
        two_components(w, sphere_stl)
        w.quality_action.setChecked(True)
        w.refresh_quality()
        for key, block in w.session.blocks():
            expected = colormap(boundary_values(scaled_jacobian(block.nodes)), (-1.0, 1.0))
            assert np.array_equal(w.viewport.face_colours(key), expected)
        w.quality_action.setChecked(False)
        w.refresh_quality()
        key, _ = w.session.blocks()[0]
        assert len(np.unique(w.viewport.face_colours(key), axis=0)) == 1


# ---------------------------------------------------------------------------
# G8-G9: the worker and errors
# ---------------------------------------------------------------------------

LONG = {"topology": "cubed_sphere", "nk": "20", "ds": "0.01", "growth": "1.02",
        "params": "{'ni': 65, 'nj': 128}"}  # fmt: skip


class TestWorker:
    @PENDING
    def test_g8_responsive_march(self, qapp, sphere_stl):
        from PySide6.QtCore import QTimer

        w = window(qapp)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.run_form("build_topology", LONG)
        w.wait_idle()
        late: list[float] = []
        last = [time.perf_counter()]

        def beat():
            now = time.perf_counter()
            late.append(now - last[0] - 0.010)
            last[0] = now

        timer = QTimer()
        timer.setInterval(10)
        timer.timeout.connect(beat)
        timer.start()
        t0 = time.perf_counter()
        w.run_form("march", {"march_scheme": "hyperbolic"})
        w.wait_idle()
        elapsed = time.perf_counter() - t0
        timer.stop()
        assert w.last_error is None
        assert elapsed > 0.5  # a real load on the worker
        assert max(late) <= 0.100, max(late)
        assert w.session.journal_lines()[-1].startswith("p.march(")

    @PENDING
    def test_g8_cancel(self, qapp, sphere_stl):
        from shore.project.commands import Cancelled

        w = window(qapp)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.run_form("build_topology", LONG)
        w.wait_idle()
        before = {k: b.nodes.copy() for k, b in w.session.blocks()}
        lines = w.session.journal_lines()
        w.run_form("march", {"march_scheme": "hyperbolic"})
        w.cancel()
        w.wait_idle()
        assert isinstance(w.last_error, Cancelled)
        assert w.session.journal_lines() == lines
        assert all(np.array_equal(before[k], b.nodes) for k, b in w.session.blocks())


class TestErrors:
    @PENDING
    def test_g9_bad_value(self, qapp, sphere_stl):
        w = window(qapp)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.run_form("build_topology", {**CUBED, "nk": "five"})
        w.wait_idle()
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert len(w.session.journal_lines()) == 1
        w.run_form("build_topology", CUBED)
        w.wait_idle()
        assert w.last_error is None and len(w.session.journal_lines()) == 2

    @PENDING
    def test_g9_failing_command(self, qapp, sphere_stl):
        w = window(qapp)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.run_form("build_topology", {**CUBED, "geometry": "nope"})
        w.wait_idle()
        assert w.last_error is not None and "Error" in w.log.toPlainText()
        assert len(w.session.journal_lines()) == 1 and w.session.project.mesh is None
        w.run_form("build_topology", CUBED)
        w.wait_idle()
        assert w.last_error is None and len(w.viewport.keys) == 6


# ---------------------------------------------------------------------------
# G10-G12: opening, isolation, forms
# ---------------------------------------------------------------------------


class TestOpen:
    @PENDING
    @pytest.mark.parametrize("replay", [True, False])
    def test_g10_open_saved(self, qapp, sphere_stl, tmp_path, replay):
        first = window(qapp)
        two_components(first, sphere_stl)
        first.save_project(str(tmp_path / "saved"))
        w = window(qapp)
        w.open_project(str(tmp_path / "saved"), replay=replay)
        assert w.last_error is None
        lines, saved = w.session.journal_lines(), first.session.journal_lines()
        assert lines[0].startswith("p.load_geometry('inputs/")  # the saved copy
        assert lines[1:] == saved[1:]
        assert w.viewport.keys == first.viewport.keys
        assert same_grids(w.session.project, first.session.project)
        r = subprocess.run(
            [sys.executable, "-c", "from shore.cli import app; app()", "gui",
             str(tmp_path / "saved"), "--smoke", *([] if replay else ["--no-replay"])],
            capture_output=True, text=True, timeout=300,
            env={**__import__("os").environ, "QT_QPA_PLATFORM": "offscreen"},
        )  # fmt: skip
        assert r.returncode == 0, r.stdout + r.stderr
        assert f"shore gui smoke: {len(first.viewport.keys)} blocks, journal 5 lines" in r.stdout


GUI_STACK = ("pygfx", "wgpu", "rendercanvas")


def imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text())
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            out.add(node.module.split(".")[0])
    return out


class TestIsolation:
    @PENDING
    def test_g11_imports(self):
        for path in SRC.rglob("*.py"):
            rel = path.relative_to(SRC).as_posix()
            found = imports(path)
            if rel != "gui/viewport/pygfx_backend.py":
                assert not found & set(GUI_STACK), rel
            if not rel.startswith("gui/"):
                assert "PySide6" not in found, rel
        assert imports(SRC / "gui" / "viewport" / "pygfx_backend.py") & {"pygfx"}

    def test_g11b_core_without_gui_stack(self):  # already holds: a guard
        code = (
            "import sys\n"
            "for m in ('PySide6', 'pygfx', 'wgpu', 'rendercanvas'): sys.modules[m] = None\n"
            "import shore.cli, shore.project\n"
            "try:\n"
            "    import shore.gui.session, shore.gui.forms\n"
            "except ModuleNotFoundError as e:\n"
            "    if not e.name.startswith('shore.gui'): raise\n"
            "print('ok')\n"
        )
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        assert r.returncode == 0 and "ok" in r.stdout, r.stderr


class TestForms:
    @PENDING
    def test_g12_forms_call_project(self):
        from shore.gui.forms import FORMS
        from shore.project import Project

        for f in FORMS:
            method = getattr(Project, f.method)
            args, kwargs = f.arguments({})
            params = inspect.signature(method).parameters
            names = {n for n, p in params.items() if p.kind is not p.VAR_KEYWORD}
            assert set(kwargs) <= names or any(p.kind is p.VAR_KEYWORD for p in params.values()), (
                f.method
            )
            inspect.signature(method).bind(None, *args, **kwargs)
