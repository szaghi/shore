"""End-to-end tests for ``shore primitive {box,annulus,flat-caps}``."""

from __future__ import annotations

import textwrap
from pathlib import Path

from typer.testing import CliRunner

from shore.cli import app
from shore.io.geo import read_geo

runner = CliRunner()


# ---------------------------------------------------------------------------
# Box
# ---------------------------------------------------------------------------


class TestPrimitiveBoxCLI:
    def test_basic_invocation(self, tmp_path):
        out = tmp_path / "box.geo"
        r = runner.invoke(
            app,
            [
                "primitive",
                "box",
                "--min",
                "0",
                "0",
                "0",
                "--max",
                "10",
                "5",
                "3",
                "--ni",
                "12",
                "--nj",
                "8",
                "--nk",
                "6",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        nodes = read_geo(out)
        assert nodes.shape == (12, 8, 6, 3)

    def test_invalid_min_max_errors(self, tmp_path):
        r = runner.invoke(
            app,
            [
                "primitive",
                "box",
                "--min",
                "5",
                "0",
                "0",
                "--max",
                "5",
                "5",
                "3",
                "--ni",
                "4",
                "--nj",
                "4",
                "--nk",
                "4",
                "-o",
                str(tmp_path / "bad.geo"),
            ],
        )
        assert r.exit_code != 0
        assert "vmax must exceed vmin" in r.output


# ---------------------------------------------------------------------------
# Annulus
# ---------------------------------------------------------------------------


class TestPrimitiveAnnulusCLI:
    def test_full_annulus(self, tmp_path):
        out = tmp_path / "ann.geo"
        r = runner.invoke(
            app,
            [
                "primitive",
                "annulus",
                "--r-in",
                "0.5",
                "--r-out",
                "2.0",
                "--z0",
                "0",
                "--z1",
                "4",
                "--ni",
                "8",
                "--nj",
                "16",
                "--nk",
                "5",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        assert read_geo(out).shape == (8, 16, 5, 3)

    def test_partial_sector(self, tmp_path):
        out = tmp_path / "ann_sec.geo"
        r = runner.invoke(
            app,
            [
                "primitive",
                "annulus",
                "--r-in",
                "0.5",
                "--r-out",
                "2.0",
                "--z0",
                "0",
                "--z1",
                "4",
                "--theta-start",
                "0",
                "--theta-end",
                "3.14159",
                "--ni",
                "6",
                "--nj",
                "12",
                "--nk",
                "4",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        assert read_geo(out).shape == (6, 12, 4, 3)

    def test_axis_x(self, tmp_path):
        import numpy as np

        out = tmp_path / "ann_x.geo"
        r = runner.invoke(
            app,
            [
                "primitive",
                "annulus",
                "--r-in",
                "0.5",
                "--r-out",
                "2.0",
                "--z0",
                "0",
                "--z1",
                "4",
                "--axis",
                "x",
                "--ni",
                "4",
                "--nj",
                "12",
                "--nk",
                "5",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        nodes = read_geo(out)
        # Axial direction should now run along x.
        axial = nodes[0, 0, :, :]
        assert np.isclose(axial[0, 0], 0.0)
        assert np.isclose(axial[-1, 0], 4.0)


# ---------------------------------------------------------------------------
# Flat-caps
# ---------------------------------------------------------------------------


class TestPrimitiveFlatCapsCLI:
    def test_writes_five_blocks(self, tmp_path):
        stem = tmp_path / "fc"
        r = runner.invoke(
            app,
            [
                "primitive",
                "flat-caps",
                "--r-out",
                "2.0",
                "--z0",
                "0",
                "--z1",
                "4",
                "--ni",
                "6",
                "--nc",
                "8",
                "--nk",
                "4",
                "-o",
                str(stem),
            ],
        )
        assert r.exit_code == 0, r.output
        for label in ("center", "sub_e", "sub_n", "sub_w", "sub_s"):
            assert (tmp_path / f"fc_{label}.geo").exists()


# ---------------------------------------------------------------------------
# Config-file driven invocation
# ---------------------------------------------------------------------------


def _write_box_cfg(tmp_path: Path) -> Path:
    body = textwrap.dedent(
        """
        version = 1
        kind    = "primitive"

        [geometry]
        type = "box"
        vmin = [0.0, 0.0, 0.0]
        vmax = [10.0, 5.0, 3.0]

        [grid]
        ni = 14
        nj = 10
        nk = 6
        """
    )
    p = tmp_path / "cfg.toml"
    p.write_text(body)
    return p


def _write_annulus_cfg(tmp_path: Path) -> Path:
    body = textwrap.dedent(
        """
        version = 1
        kind    = "primitive"

        [geometry]
        type   = "annulus"
        r_in   = 0.5
        r_out  = 2.0
        z0     = 0.0
        z1     = 4.0
        axis   = "z"

        [grid]
        ni = 9
        nj = 16
        nk = 5
        """
    )
    p = tmp_path / "ann.toml"
    p.write_text(body)
    return p


class TestPrimitiveConfigFlag:
    def test_box_from_config(self, tmp_path):
        cfg = _write_box_cfg(tmp_path)
        out = tmp_path / "from_cfg.geo"
        r = runner.invoke(app, ["primitive", "box", "-c", str(cfg), "-o", str(out)])
        assert r.exit_code == 0, r.output
        assert read_geo(out).shape == (14, 10, 6, 3)

    def test_box_cli_overrides_config(self, tmp_path):
        cfg = _write_box_cfg(tmp_path)
        out = tmp_path / "ovr.geo"
        r = runner.invoke(
            app,
            ["primitive", "box", "-c", str(cfg), "--ni", "5", "-o", str(out)],
        )
        assert r.exit_code == 0, r.output
        assert read_geo(out).shape == (5, 10, 6, 3)

    def test_annulus_from_config(self, tmp_path):
        cfg = _write_annulus_cfg(tmp_path)
        out = tmp_path / "from_cfg.geo"
        r = runner.invoke(app, ["primitive", "annulus", "-c", str(cfg), "-o", str(out)])
        assert r.exit_code == 0, r.output
        assert read_geo(out).shape == (9, 16, 5, 3)

    def test_type_mismatch_rejected(self, tmp_path):
        # Box config fed to annulus subcommand → error.
        cfg = _write_box_cfg(tmp_path)
        r = runner.invoke(
            app, ["primitive", "annulus", "-c", str(cfg), "-o", str(tmp_path / "x.geo")]
        )
        assert r.exit_code != 0
        assert "geometry.type='box'" in r.output

    def test_kind_mismatch_rejected(self, tmp_path):
        # A mesh-style config (no kind) given to a primitive command → error.
        bad = tmp_path / "mesh.toml"
        bad.write_text('version = 1\ninput = "x.stl"\n')
        r = runner.invoke(app, ["primitive", "box", "-c", str(bad), "-o", str(tmp_path / "x.geo")])
        assert r.exit_code != 0
        assert "kind" in r.output
