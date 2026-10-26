"""Command layer beyond `shore mesh`: primitives, split, tip caps, Xall writers.

ADR 0004, phase 3 second slice. The CLI commands run over a Project; one
journal can now record the whole Xall pipeline and replay it byte for byte.
"""

from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.cli import app
from shore.exceptions import ParameterError
from shore.project import Project, run_journal

runner = CliRunner()
SPHERE = {"nk": 6, "ds": 0.05, "growth": 1.1, "ni": 10, "nj": 12, "theta_cap": 0.05 * np.pi}
LABELS = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]


def _digest(directory: Path) -> dict[str, str]:
    return {
        f.name: hashlib.sha256(f.read_bytes()).hexdigest()
        for f in sorted(directory.iterdir())
        if f.is_file() and f.suffix != ".py"
    }


def _invoke(*args: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return runner.invoke(app, list(args))


# ---------------------------------------------------------------------------
# One journal for the whole Xall pipeline
# ---------------------------------------------------------------------------


class TestXallPipelineJournal:
    def test_full_pipeline_replays_byte_identical(self, sphere_stl, tmp_path):
        out = tmp_path / "out"
        p = Project()
        p.load_geometry(sphere_stl)
        p.build_topology("cubed_sphere", **SPHERE)
        p.march()
        p.write_geo(out / "wall")
        p.write_adjacency(out / "wall.adjacency.json")
        p.split([(lbl, "k", (3,)) for lbl in LABELS])
        p.write_geo(out / "split", per_block=True)
        p.write_adjacency(out / "split.adjacency.json")
        geos = sorted(out.glob("split_*.geo"))
        p.write_grd(out / "split.grd", geos)
        p.write_cc_par(
            out / "split.adjacency.json",
            out / "split.cc.par",
            grd_filename="split.grd",
            meta={"sub0__k0": {"level": 2}},
        )
        p.write_proc_input(out / "split.grd", out / "split.proc.input", np=3)
        journal = p.save_journal(tmp_path / "pipeline.py")
        first = _digest(out)
        assert {"split.grd", "split.cc.par", "split.proc.input"} <= set(first)

        for f in out.iterdir():
            f.unlink()
        run_journal(journal)
        assert _digest(out) == first

    def test_xall_writers_match_the_cli(self, sphere_stl, tmp_path):
        """The Project writers and the stand-alone CLI commands agree byte for byte."""
        out = tmp_path / "out"
        p = Project()
        p.load_geometry(sphere_stl)
        p.build_topology("cubed_sphere", **SPHERE)
        p.march()
        p.write_geo(out / "w")
        p.write_adjacency(out / "w.adjacency.json")
        geos = [out / f"w_{lbl}.geo" for lbl in LABELS]
        p.write_grd(out / "api.grd", geos)
        p.write_cc_par(out / "w.adjacency.json", out / "api.cc.par", grd_filename="x.grd")
        p.write_proc_input(out / "api.grd", out / "api.proc.input", np=2)

        assert _invoke("grd", *map(str, geos), "-o", str(out / "cli.grd")).exit_code == 0
        assert (
            _invoke(
                "cc-par", str(out / "w.adjacency.json"), str(out / "cli.cc.par"), "--grd", "x.grd"
            ).exit_code
            == 0
        )
        r = _invoke("proc-input", str(out / "cli.grd"), str(out / "cli.proc.input"), "--np", "2")
        assert r.exit_code == 0
        for ext in ("grd", "cc.par", "proc.input"):
            assert (out / f"api.{ext}").read_bytes() == (out / f"cli.{ext}").read_bytes(), ext


# ---------------------------------------------------------------------------
# CLI --journal on the re-expressed commands
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def wall(tmp_path_factory, sphere_stl) -> Path:
    """A cubed-sphere wall component written by `shore mesh` (geo + sidecar)."""
    d = tmp_path_factory.mktemp("wall")
    r = _invoke(
        "mesh", str(sphere_stl), "--ni", "10", "--nj", "12", "--nk", "6",
        "--ds", "0.05", "--growth", "1.1",
        "-o", str(d / "wall"), "--adjacency-json", str(d / "wall.adjacency.json"),
    )  # fmt: skip
    assert r.exit_code == 0, r.output
    return d


def _replay_matches(tmp_path: Path, args: list[str], out: Path) -> None:
    journal = tmp_path / "j.py"
    r = _invoke(*args, "--journal", str(journal))
    assert r.exit_code == 0, r.output
    first = _digest(out)
    for f in out.iterdir():
        if f.is_file():
            f.unlink()
    r = _invoke("run", str(journal))
    assert r.exit_code == 0, r.output
    assert _digest(out) == first


class TestCliJournals:
    @pytest.mark.parametrize(
        "args",
        [
            ["primitive", "box", "--min", "0", "0", "0", "--max", "1", "2", "3",
             "--ni", "4", "--nj", "5", "--nk", "6"],
            ["primitive", "annulus", "--r-in", "0.5", "--r-out", "2", "--z0", "0", "--z1", "1",
             "--ni", "5", "--nj", "8", "--nk", "4"],
            ["primitive", "flat-caps", "--r-out", "3", "--z0", "-1", "--z1", "4",
             "--ni", "5", "--nc", "6", "--nk", "5"],
            ["primitive", "h-grid-rect", "--length", "4", "--width", "1", "--height", "0.5",
             "--ni", "6", "--nj", "4", "--nk", "3"],
            ["primitive", "h-grid-tube", "--length", "3", "--radius", "0.5",
             "--n-axial", "6", "--n-circ", "8", "--n-radial", "4"],
        ],
        ids=["box", "annulus", "flat-caps", "h-grid-rect", "h-grid-tube"],
    )  # fmt: skip
    def test_primitive(self, tmp_path, args):
        out = tmp_path / "out"
        out.mkdir()
        _replay_matches(
            tmp_path,
            [*args, "-o", str(out / "p.geo"), "--adjacency-json", str(out / "p.adjacency.json")],
            out,
        )

    def test_split(self, tmp_path, wall):
        cfg = tmp_path / "split.toml"
        cfg.write_text(
            'version = 1\nkind = "split"\n'
            + "".join(f'\n[[splits]]\nlabel = "{lb}"\naxis = "k"\nat = [3]\n' for lb in LABELS)
        )
        out = tmp_path / "out"
        out.mkdir()
        _replay_matches(
            tmp_path,
            ["split", *map(str, sorted(wall.glob("wall_*.geo"))), "-c", str(cfg),
             "--adjacency", str(wall / "wall.adjacency.json"), "-o", str(out / "s")],
            out,
        )  # fmt: skip

    def test_xall_writers(self, tmp_path, wall):
        geos = [str(wall / f"wall_{lbl}.geo") for lbl in LABELS]
        steps = {
            "grd": ["grd", *geos, "-o", "{out}/w.grd"],
            "cc-par": [
                "cc-par", str(wall / "wall.adjacency.json"), "{out}/w.cc.par", "--grd", "w.grd",
                "--meta", json.dumps({"sub0": {"level": 2}}), "--free-bc", "-9",
            ],
            "proc-input": [
                "proc-input", "{grd}", "{out}/w.proc.input", "--np", "3",
                "--meta", json.dumps({"1": {"body": 1}}),
            ],
        }  # fmt: skip
        grd = tmp_path / "w.grd"
        assert _invoke("grd", *geos, "-o", str(grd)).exit_code == 0
        for name, args in steps.items():
            work = tmp_path / name
            out = work / "out"
            out.mkdir(parents=True)
            args = [a.replace("{out}", str(out)).replace("{grd}", str(grd)) for a in args]
            _replay_matches(work, args, out)

    def test_cap_tip(self, tmp_path, naca0012_stl):
        body = tmp_path / "body"
        body.mkdir()
        r = _invoke(
            "mesh", str(naca0012_stl), "--topology", "cgrid", "--ni-body", "60",
            "--ni-wake", "15", "--n-stations", "4", "--wake-length", "10.0", "--nk", "8",
            "--ds", "1e-3", "--smooth", "0", "--smooth-iters", "0", "-o", str(body / "w"),
        )  # fmt: skip
        assert r.exit_code == 0, r.output
        out = tmp_path / "out"
        out.mkdir()
        _replay_matches(
            tmp_path,
            ["cap-tip", str(body / "w_body.geo"), "--tip", "jlo", "-o", str(out / "c")],
            out,
        )


# ---------------------------------------------------------------------------
# Command semantics
# ---------------------------------------------------------------------------


class TestCommands:
    def test_split_undo_and_no_march_after(self, sphere_stl):
        p = Project()
        p.load_geometry(sphere_stl)
        p.build_topology("cubed_sphere", **SPHERE)
        p.split([("sub0", "k", (3,))])
        assert p.mesh is not None and "sub0__k0" in [b.label for b in p.mesh.blocks]
        with pytest.raises(ParameterError, match="cannot be marched"):
            p.march()
        assert p.undo().method == "split"
        assert [b.label for b in p.mesh.blocks] == LABELS

    def test_split_unknown_label(self, sphere_stl):
        p = Project()
        p.load_geometry(sphere_stl)
        p.build_topology("cubed_sphere", **SPHERE)
        with pytest.raises(ParameterError, match="no block named"):
            p.split([("nope", "k", (3,))])

    def test_load_blocks_labels(self, wall):
        p = Project()
        p.load_blocks(sorted(wall.glob("wall_*.geo")), [wall / "wall.adjacency.json"])
        assert p.mesh is not None
        assert sorted(b.label for b in p.mesh.blocks) == sorted(LABELS)
        assert p.mesh.blocks[0].face("j_lo").partner is not None  # seams wired

    def test_primitive_validation(self):
        p = Project()
        with pytest.raises(ParameterError, match="Unknown primitive"):
            p.build_primitive("sphere")
        with pytest.raises(ParameterError, match="nii"):
            p.build_primitive("box", vmin=(0, 0, 0), vmax=(1, 1, 1), nii=3, nj=3, nk=3)

    def test_primitive_naming_and_per_block(self, tmp_path):
        p = Project()
        p.build_primitive("flat_caps", r_out=3.0, z0=-1.0, z1=4.0, ni=5, nc=6, nk=5, inner_frac=0.4)
        p.write_geo(tmp_path / "fc.geo")  # primitive rule: strips .geo, per-block names
        assert (tmp_path / "fc_center.geo").exists()
        p.write_geo(tmp_path / "x", per_block=True)
        assert (tmp_path / "x_center.geo").exists()

    def test_proc_input_bad_meta(self, wall, tmp_path):
        p = Project()
        geos = [wall / f"wall_{lbl}.geo" for lbl in LABELS]
        p.write_grd(tmp_path / "w.grd", geos)
        with pytest.raises(ParameterError, match="block index"):
            p.write_proc_input(tmp_path / "w.grd", tmp_path / "p.in", meta={"x": {"body": 1}})
