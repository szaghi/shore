"""FROZEN acceptance spec: every CLI capability in the GUI (ADR 0062).

Frozen 2026-11-01, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job):
``pytest tests_gui/spec/test_gui_parity.py`` reports 16 passed, 0 skipped,
0 xfailed. At freeze: 0 passed, 16 xfailed.

Interface under test
--------------------
``shore.gui.forms.ALL_FORMS`` holds a form for every leaf command of the CLI,
keyed ``"cli:<command path>"`` (``"cli:mesh"``, ``"cli:primitive box"``):
its fields are the command's options and arguments (named as the CLI's
Python parameters, with the CLI's defaults), less the ones the GUI owns
(``project``, ``component``, ``no_replay``, ``journal``). Running one
(``MainWindow.run_form(key, values)``) runs the CLI command in-process on
the session's project: the same defaults, the same files, the same journal
lines as the command line. A failed command sets ``last_error`` and leaves
the project as it was. ``MAPPED`` maps the CLI commands that open a project
directory or a viewer to the GUI's own forms (``write_xall``,
``write_cgns``, ``write_plot3d``, ``run_journal``, ``view_file``);
``EXCLUDED`` holds ``gui``. Every Project command has a form too; ``FORMS``
keeps only the Project-method forms (ADR 0037's G12), ``ALL_FORMS`` adds the
session's and the CLI's. ``form(key)`` finds any of them. Every form is
reachable from ``MainWindow.commands_menu``.

Cases
-----
- P1 every CLI leaf command has a form, or is mapped or excluded.
- P2 every Project command (``shore.project.commands``) has a form.
- P3 the h-grid-tube form's fields carry the CLI's defaults, without the
  project options.
- P4 (x 5) each primitive from its form with the CLI's defaults: the journal
  line of ``shore primitive <kind> --journal``.
- P5 ``mesh`` from its form: the CLI's journal lines and grids.
- P6 ``body sphere`` and ``profile naca`` write the CLI's bytes.
- P7 ``grd``, ``proc-input`` and ``cc-par`` write the CLI's bytes.
- P8 a failing form: ``last_error`` set, the project unchanged.
- P9 ``write_xall`` with the options TOML writes what ``shore xall
  --options`` writes.
- P10 ``run_journal`` replays a journal into the session.
- P11 ``view_file`` shows a ``.geo`` file's block (one journal line,
  ``load_blocks``) and a ``.npz`` surface (display only, no journal line).
- P12 every form is an action of the Commands menu.

Out of scope of this spec: ``shore gui`` itself; the look of the forms.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-11-01; implementation pending")

PRIMITIVES = {
    "box": ["--ni", "5", "--nj", "5", "--nk", "5"],
    "annulus": [],
    "flat-caps": [],
    "h-grid-rect": [],
    "h-grid-tube": [],
}  # fmt: skip
MESH = ["--ni", "9", "--nj", "16", "--nk", "6", "--ds", "0.05"]


def window():
    from shore.gui.window import MainWindow

    return MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)


def cli(*args):
    from shore.cli import app

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = CliRunner().invoke(app, [str(a) for a in args])
    assert r.exit_code == 0, r.output
    return r


def as_values(argv):
    """``--name value`` pairs as form values (the CLI parameter names)."""
    values, it = {}, iter(argv)
    for flag in it:
        values[flag.lstrip("-").replace("-", "_")] = next(it)
    return values


def run(w, key, values):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w.run_form(key, values)
    w._pool.waitForDone()
    from PySide6.QtWidgets import QApplication

    QApplication.processEvents()
    assert w.last_error is None, w.last_error


def leaves():
    from typer.main import get_command

    from shore.cli import app

    out = []

    def walk(c, prefix):
        for name, sub in c.commands.items():
            if hasattr(sub, "commands"):
                walk(sub, f"{prefix}{name} ")
            else:
                out.append(f"{prefix}{name}")

    walk(get_command(app), "")
    return out


def journal_line(path, method):
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.startswith(f"p.{method}(")]
    assert lines, method
    return lines[-1]


class TestCoverage:
    @PENDING
    def test_p1_every_cli_command(self, qapp):
        from shore.gui.forms import EXCLUDED, MAPPED, form

        missing = []
        for path in leaves():
            if path in EXCLUDED:
                continue
            key = MAPPED.get(path, f"cli:{path}")
            try:
                form(key)
            except KeyError:
                missing.append(path)
        assert missing == []
        assert set(MAPPED) | set(EXCLUDED) <= set(leaves())
        assert "gui" in EXCLUDED

    @PENDING
    def test_p2_every_project_command(self, qapp):
        import dataclasses

        from shore.gui.forms import form
        from shore.project import commands

        methods = {
            c.method
            for c in vars(commands).values()
            if isinstance(c, type)
            and dataclasses.is_dataclass(c)
            and isinstance(getattr(c, "method", None), str)
        }
        missing = []
        for m in sorted(methods):
            try:
                form(m)
            except KeyError:
                missing.append(m)
        assert missing == []

    @PENDING
    def test_p3_fields_carry_cli_defaults(self, qapp):
        from shore.gui.forms import form

        f = form("cli:primitive h-grid-tube")
        fields = {x.name: x for x in f.fields}
        assert fields["n_circ"].default == 32
        assert fields["core_frac"].default == 0.5
        assert fields["streamwise_axis"].default == "x"
        assert not {"project", "component", "no_replay", "journal"} & set(fields)


class TestProjectCommands:
    @pytest.mark.parametrize("kind", list(PRIMITIVES))
    @PENDING
    def test_p4_primitives(self, qapp, tmp_path, kind):
        argv = PRIMITIVES[kind]
        cli("primitive", kind, *argv, "-o", tmp_path / "cli", "--journal", tmp_path / "j.py")
        w = window()
        run(w, f"cli:primitive {kind}", {**as_values(argv), "output": str(tmp_path / "gui")})
        assert w.session.project.mesh is not None
        line = journal_line(tmp_path / "j.py", "build_primitive")
        assert line in w.session.journal_lines()

    @PENDING
    def test_p5_mesh(self, qapp, tmp_path, sphere_stl):
        from shore.io.geo import read_geo

        cli("mesh", sphere_stl, *MESH, "-o", tmp_path / "cli" / "s", "--journal", tmp_path / "j.py")
        w = window()
        values = {
            **as_values(MESH),
            "input": str(sphere_stl),
            "output": str(tmp_path / "gui" / "s"),
        }
        run(w, "cli:mesh", values)
        lines = w.session.journal_lines()
        for method in ("build_topology", "march"):
            assert journal_line(tmp_path / "j.py", method) in lines
        for f in sorted((tmp_path / "cli").glob("*.geo")):
            np.testing.assert_array_equal(read_geo(f), read_geo(tmp_path / "gui" / f.name))

    @PENDING
    def test_p6_generators(self, qapp, tmp_path):
        (tmp_path / "gui").mkdir()
        cli("body", "sphere", "-o", tmp_path / "s.stl")
        cli("profile", "naca", "2412", "-o", tmp_path / "n.dat")
        w = window()
        run(w, "cli:body sphere", {"output": str(tmp_path / "gui" / "s.stl")})
        run(
            w,
            "cli:profile naca",
            {"designation": "2412", "output": str(tmp_path / "gui" / "n.dat")},
        )
        for name in ("s.stl", "n.dat"):
            assert (tmp_path / name).read_bytes() == (tmp_path / "gui" / name).read_bytes()

    @PENDING
    def test_p7_xall_file_tools(self, qapp, tmp_path):
        box = PRIMITIVES["box"]
        cli("primitive", "box", *box, "-o", tmp_path / "b", "--adjacency-json", tmp_path / "b.json")
        geo = tmp_path / "b.geo"
        cli("grd", geo, "-o", tmp_path / "cli.grd")
        cli("proc-input", tmp_path / "cli.grd", tmp_path / "cli.proc", "--np", "2")
        cli("cc-par", tmp_path / "b.json", tmp_path / "cli.par", "--grd", "b.grd")
        w = window()
        run(w, "cli:grd", {"inputs": [str(geo)], "output": str(tmp_path / "gui.grd")})
        run(
            w,
            "cli:proc-input",
            {"grd": str(tmp_path / "gui.grd"), "np": "2", "output": str(tmp_path / "gui.proc")},
        )
        run(
            w,
            "cli:cc-par",
            {
                "adjacency": str(tmp_path / "b.json"),
                "output": str(tmp_path / "gui.par"),
                "grd_filename": "b.grd",
            },
        )
        for ext in ("grd", "proc", "par"):
            assert (tmp_path / f"cli.{ext}").read_bytes() == (
                tmp_path / f"gui.{ext}"
            ).read_bytes(), ext

    @PENDING
    def test_p8_failure_changes_nothing(self, qapp, tmp_path):
        w = window()
        run(w, "cli:primitive box", {**as_values(PRIMITIVES["box"]), "output": str(tmp_path / "b")})
        before = list(w.session.journal_lines())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w.run_form("cli:primitive h-grid-tube", {"n_circ": "10", "output": str(tmp_path / "t")})
        w._pool.waitForDone()
        from PySide6.QtWidgets import QApplication

        QApplication.processEvents()
        assert w.last_error is not None
        assert w.session.journal_lines() == before
        assert len(w.session.project.mesh.blocks) == 1


class TestMapped:
    @PENDING
    def test_p9_write_xall_with_options(self, qapp, tmp_path):
        toml = tmp_path / "x.toml"
        toml.write_text("[run]\nnp = 2\n\n[components.main]\nfree_face_bc = 9\n")
        w = window()
        run(w, "cli:primitive box", {**as_values(PRIMITIVES["box"]), "output": str(tmp_path / "b")})
        w.session.project.save(tmp_path / "proj")
        cli("xall", tmp_path / "proj", "--out", tmp_path / "cli", "--options", toml)
        run(w, "write_xall", {"directory": str(tmp_path / "gui"), "options": str(toml)})

        def files(d):
            return sorted(p.relative_to(d) for p in d.rglob("*") if p.is_file())

        names = files(tmp_path / "cli")
        assert names
        assert names == files(tmp_path / "gui")
        for name in names:
            assert (tmp_path / "cli" / name).read_bytes() == (
                tmp_path / "gui" / name
            ).read_bytes(), name

    @PENDING
    def test_p10_run_journal(self, qapp, tmp_path):
        cli(
            "primitive",
            "box",
            *PRIMITIVES["box"],
            "-o",
            tmp_path / "b",
            "--journal",
            tmp_path / "j.py",
        )
        w = window()
        run(w, "run_journal", {"path": str(tmp_path / "j.py")})
        assert len(w.session.project.mesh.blocks) == 1
        assert journal_line(tmp_path / "j.py", "build_primitive") in w.session.journal_lines()

    @PENDING
    def test_p11_view_file(self, qapp, tmp_path, sphere_stl):
        cli("primitive", "box", *PRIMITIVES["box"], "-o", tmp_path / "b")
        cli("project", sphere_stl, "--ni", "9", "--nj", "16", "-o", tmp_path / "surf.npz")
        w = window()
        run(w, "view_file", {"path": str(tmp_path / "b.geo")})
        assert len(w.session.project.mesh.blocks) == 1
        assert any(ln.startswith("p.load_blocks(") for ln in w.session.journal_lines())
        before = list(w.session.journal_lines())
        run(w, "view_file", {"path": str(tmp_path / "surf.npz")})
        assert w.session.journal_lines() == before
        assert "surf" in w.viewport.geometry_names


class TestMenu:
    @PENDING
    def test_p12_menu(self, qapp):
        from shore.gui.forms import ALL_FORMS

        w = window()
        found = []

        def walk(menu):
            for a in menu.actions():
                if a.menu() is not None:
                    walk(a.menu())
                elif a.data():
                    found.append(a.data())

        walk(w.commands_menu)
        assert sorted(found) == sorted(f.method for f in ALL_FORMS)
