"""Tests for the command layer: Project, commands, journal (ADR 0004, phase 3)."""

from __future__ import annotations

import hashlib
import warnings
from pathlib import Path

import numpy as np
import pytest
import trimesh
from typer.testing import CliRunner

from shore.cli import app
from shore.exceptions import MarchingError, MeshInputError, ParameterError
from shore.project import (
    BuildTopology,
    Cancelled,
    CancelToken,
    LoadGeometry,
    March,
    Project,
    run_journal,
)

runner = CliRunner()
SPHERE = {"nk": 6, "ds": 0.05, "growth": 1.1, "ni": 10, "nj": 12, "theta_cap": 0.05 * np.pi}


def _sphere_project(sphere_stl: Path) -> Project:
    p = Project()
    p.load_geometry(sphere_stl)
    p.build_topology("cubed_sphere", **SPHERE)
    return p


def _nodes(p: Project) -> list[np.ndarray]:
    assert p.mesh is not None
    return [b.nodes.copy() for b in p.mesh.blocks]


def _digest(paths: list[Path]) -> dict[str, str]:
    return {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(paths)}


# ---------------------------------------------------------------------------
# Journal rendering
# ---------------------------------------------------------------------------


class TestToSource:
    def test_defaults_omitted_positional_first(self):
        assert LoadGeometry("s.stl").to_source() == "p.load_geometry('s.stl')"
        assert March().to_source() == "p.march()"
        assert March(smooth_iters=0).to_source() == "p.march(smooth_iters=0)"

    def test_build_params_as_keywords(self):
        cmd = BuildTopology("oh", 5, options={"nk_channel": 4}, params={"ni": 16, "nj": 16})
        assert cmd.to_source() == (
            "p.build_topology('oh', nk=5, options={'nk_channel': 4}, ni=16, nj=16)"
        )

    def test_floats_round_trip(self, sphere_stl):
        p = _sphere_project(sphere_stl)
        line = p.history[-1].to_source()
        assert repr(0.05 * np.pi) in line  # repr is exact, so replay is bit-identical


# ---------------------------------------------------------------------------
# Undo / redo
# ---------------------------------------------------------------------------


class TestUndoRedo:
    def test_undo_march_restores_bits(self, sphere_stl):
        p = _sphere_project(sphere_stl)
        before = _nodes(p)
        p.march()
        after = _nodes(p)
        assert p.undo().method == "march"
        for a, b in zip(_nodes(p), before, strict=True):
            np.testing.assert_array_equal(a, b)
        assert p.redo().method == "march"
        for a, b in zip(_nodes(p), after, strict=True):
            np.testing.assert_array_equal(a, b)

    def test_undo_skips_outputs(self, sphere_stl, tmp_path):
        p = _sphere_project(sphere_stl)
        p.march()
        p.write_geo(tmp_path / "s")
        assert p.undo().method == "march"  # the trailing output goes with it
        assert [c.method for c in p.history] == ["load_geometry", "build_topology"]

    def test_undo_build_and_load(self, sphere_stl):
        p = _sphere_project(sphere_stl)
        p.undo()
        assert p.mesh is None and p.build is None and p.body is not None
        p.undo()
        assert p.body is None
        with pytest.raises(ParameterError, match="Nothing to undo"):
            p.undo()

    def test_new_command_clears_redo(self, sphere_stl):
        p = _sphere_project(sphere_stl)
        p.march()
        p.undo()
        p.snap()
        with pytest.raises(ParameterError, match="Nothing to redo"):
            p.redo()

    def test_undo_cache_removed_on_close(self, sphere_stl):
        p = _sphere_project(sphere_stl)
        p.march()
        cache = p.undo_cache._dir
        assert any(cache.iterdir())
        p.close()
        assert not cache.exists()


# ---------------------------------------------------------------------------
# Long commands: progress, cancellation, failure
# ---------------------------------------------------------------------------


class TestLongCommands:
    def test_progress_reported_per_layer(self, sphere_stl):
        p = _sphere_project(sphere_stl)
        seen: list[float] = []
        p.march(progress=lambda f, _m: seen.append(f))
        assert len(seen) == SPHERE["nk"] - 1
        assert seen[-1] == pytest.approx(1.0)

    def test_cancel_restores_and_records_nothing(self, sphere_stl):
        p = _sphere_project(sphere_stl)
        before = _nodes(p)
        token = CancelToken()

        def stop_after_two(fraction: float, _m: str) -> None:
            if fraction >= 2 / (SPHERE["nk"] - 1):
                token.cancel()

        with pytest.raises(Cancelled):
            p.march(progress=stop_after_two, cancel=token)
        for a, b in zip(_nodes(p), before, strict=True):
            np.testing.assert_array_equal(a, b)
        assert [c.method for c in p.history] == ["load_geometry", "build_topology"]

    def test_failed_march_restores(self, sphere_stl):
        p = Project()
        p.load_geometry(sphere_stl)
        p.build_topology("cubed_sphere", nk=6, ni=10, nj=12)
        before = _nodes(p)
        with pytest.raises(MarchingError):  # unlimited smoothing folds ds=1e-3 here
            p.march(smooth_mode="legacy")
        for a, b in zip(_nodes(p), before, strict=True):
            np.testing.assert_array_equal(a, b)


class TestValidation:
    def test_build_needs_body(self):
        with pytest.raises(ParameterError, match="load_geometry"):
            Project().build_topology("ogrid", nk=4, ni=8, nj=8)

    def test_march_needs_mesh(self, sphere_stl):
        p = Project()
        p.load_geometry(sphere_stl)
        with pytest.raises(ParameterError, match="build_topology"):
            p.march()

    def test_unknown_topology(self, sphere_stl):
        p = Project()
        p.load_geometry(sphere_stl)
        with pytest.raises(ParameterError, match="Unknown topology"):
            p.build_topology("tetra", nk=4)

    def test_misspelt_parameter_rejected(self, sphere_stl):
        p = Project()
        p.load_geometry(sphere_stl)
        with pytest.raises(ParameterError, match="nii"):
            p.build_topology("cubed_sphere", nk=4, nii=10, nj=12)


# ---------------------------------------------------------------------------
# Journal replay and persistence
# ---------------------------------------------------------------------------


class TestJournal:
    def test_replay_is_bit_identical(self, sphere_stl, tmp_path):
        p = _sphere_project(sphere_stl)
        p.march()
        journal = p.save_journal(tmp_path / "j.py")
        q = run_journal(journal)
        for a, b in zip(_nodes(q), _nodes(p), strict=True):
            np.testing.assert_array_equal(a, b)

    def test_journal_without_project(self, tmp_path):
        bad = tmp_path / "bad.py"
        bad.write_text("x = 1\n")
        with pytest.raises(MeshInputError, match="Project named 'p'"):
            run_journal(bad)


class TestPersistence:
    def test_save_load_round_trip(self, sphere_stl, tmp_path):
        p = _sphere_project(sphere_stl)
        p.march()
        p.save(tmp_path / "proj")
        assert (tmp_path / "proj" / "adjacency.json").exists()
        q = Project.load(tmp_path / "proj")
        # the saved journal reads the project's own copy of the STL (ADR 0016)
        saved = (tmp_path / "proj" / "journal.py").read_text()
        assert "p.load_geometry('inputs/00_sphere.stl')" in saved
        assert q.journal() == saved

    def test_load_detects_tampered_grid(self, sphere_stl, tmp_path):
        p = _sphere_project(sphere_stl)
        p.march()
        d = p.save(tmp_path / "proj")
        grid = d / "grids" / "main" / "0000.npz"
        with np.load(grid) as g:
            nodes = g["nodes"].copy()
        nodes[0, 0, 0, 0] += 1e-12
        np.savez(grid, nodes=nodes)
        with pytest.raises(MeshInputError, match="differs from the saved grid"):
            Project.load(d)

    def test_not_a_project(self, tmp_path):
        with pytest.raises(MeshInputError, match="not a SHORE project"):
            Project.load(tmp_path)


# ---------------------------------------------------------------------------
# CLI: shore mesh --journal, shore run
# ---------------------------------------------------------------------------

_CGRID = [
    "--topology", "cgrid", "--ni-body", "60", "--ni-wake", "15", "--n-stations", "4",
    "--wake-length", "10.0", "--nk", "8", "--ds", "1e-3", "--smooth", "0.0",
    "--smooth-iters", "0",
]  # fmt: skip
_CH = [
    *(_CGRID[:1]), "ch", *_CGRID[2:],
    "--ch-chord-min", "-2.0", "--ch-chord-max", "15.0",
    "--ch-normal-min", "-2.5", "--ch-normal-max", "2.5", "--nk-channel", "4",
]  # fmt: skip
_SMALL = ["--ni", "10", "--nj", "12", "--nk", "6", "--ds", "0.05", "--growth", "1.1"]


class TestCliJournal:
    @pytest.mark.parametrize(
        ("fixture", "args"),
        [("sphere_stl", _SMALL), ("naca0012_stl", _CGRID), ("naca0012_stl", _CH)],
        ids=["cubed_sphere", "cgrid", "ch"],
    )
    def test_mesh_journal_replays_byte_identical(self, request, fixture, args, tmp_path):
        src = request.getfixturevalue(fixture)
        out = tmp_path / "out"
        out.mkdir()
        journal = tmp_path / "run.py"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = runner.invoke(
                app,
                [
                    "mesh", str(src), *args, "-o", str(out / "m"),
                    "--adjacency-json", str(out / "m.adjacency.json"),
                    "--journal", str(journal),
                ],
            )  # fmt: skip
        assert r.exit_code == 0, r.output
        first = _digest(list(out.iterdir()))
        for f in out.iterdir():
            f.unlink()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = runner.invoke(app, ["run", str(journal), "--save-project", str(tmp_path / "p")])
        assert r.exit_code == 0, r.output
        assert _digest(list(out.iterdir())) == first
        # the saved project keeps the state commands, reading its own input copies
        saved = (tmp_path / "p" / "journal.py").read_text()
        assert "inputs/00_" in saved and "p.write_" not in saved
        state = [ln for ln in journal.read_text().splitlines() if not ln.startswith("p.write_")]
        assert len(saved.splitlines()) == len(state)
        Project.load(tmp_path / "p")  # replays bit-identically from the copies

    def test_run_missing_journal(self, tmp_path):
        r = runner.invoke(app, ["run", str(tmp_path / "nope.py")])
        assert r.exit_code == 1


class TestComponentsAndInputs:
    """Named geometries, components and portable inputs (ADR 0016) beyond the spec."""

    def test_journal_lines(self):
        from shore.project import UseComponent, WriteXall

        assert LoadGeometry("s.stl").to_source() == "p.load_geometry('s.stl')"
        assert LoadGeometry("s.stl", name="b").to_source() == "p.load_geometry('s.stl', name='b')"
        assert UseComponent("collar").to_source() == "p.use_component('collar')"
        assert WriteXall("run", np=4).to_source() == "p.write_xall('run', np=4)"
        line = BuildTopology("ogrid", 5, params={"ni": 9}, geometry="b").to_source()
        assert line == "p.build_topology('ogrid', nk=5, geometry='b', ni=9)"

    def test_save_copies_a_collar_file_once_and_drops_outputs(self, tmp_path):
        box = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
        box.apply_translation((0.0, 0.0, -0.5))
        box.export(str(tmp_path / "box.stl"))
        post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=48)
        post.apply_translation((0.0, 0.0, 0.2))
        post.export(str(tmp_path / "post.stl"))
        p = Project()
        p.load_geometry(tmp_path / "box.stl")
        collar = {
            "seed": (0.55, 0.0, 0.02), "n_loop": 16, "ds_surface": 0.02, "growth_surface": 1.15,
            "m_a": 5, "m_b": 5, "ds_wall": 1e-3, "growth_wall": 1.25, "nk": 8,
        }  # fmt: skip
        p.build_collar(str(tmp_path / "post.stl"), **collar)
        p.write_geo(tmp_path / "out" / "c")
        p.use_component("again")
        p.build_collar(str(tmp_path / "post.stl"), **collar)
        d = p.save(tmp_path / "proj")
        assert sorted(f.name for f in (d / "inputs").iterdir()) == ["00_box.stl", "01_post.stl"]
        saved = (d / "journal.py").read_text()
        assert saved.count("p.build_collar('inputs/01_post.stl'") == 2
        assert "p.write_geo" not in saved
        q = Project.load(d)
        assert list(q.components) == ["main", "again"]

    def test_resolve_input(self, tmp_path):
        from shore.project.commands import input_root, resolve_input

        assert resolve_input("a.stl") == "a.stl"
        token = input_root.set(tmp_path)
        try:
            assert resolve_input("inputs/a.stl") == str(tmp_path / "inputs" / "a.stl")
            assert resolve_input("/abs/a.stl") == "/abs/a.stl"
        finally:
            input_root.reset(token)
