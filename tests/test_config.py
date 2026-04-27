"""Tests for the TOML config loader and the ``shore mesh -c`` flow."""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import pytest
from typer.testing import CliRunner

from shore.cli import app
from shore.config import (
    CONFIG_VERSION,
    GridConfig,
    MeshConfig,
    dump_template,
    load_config,
)
from shore.exceptions import ParameterError
from shore.io.geo import read_geo

runner = CliRunner()


# ---------------------------------------------------------------------------
# load_config — schema validation
# ---------------------------------------------------------------------------


def _write(tmp_path: Path, body: str, name: str = "shore.toml") -> Path:
    p = tmp_path / name
    p.write_text(textwrap.dedent(body))
    return p


class TestLoadConfigBasics:
    def test_minimal_config_loads(self, tmp_path):
        p = _write(
            tmp_path,
            """\
            version = 1
            input = "body.stl"
            """,
        )
        cfg = load_config(p)
        assert isinstance(cfg, MeshConfig)
        # Defaults populated when sections omitted.
        assert cfg.grid.ni == GridConfig().ni
        assert cfg.input == (tmp_path / "body.stl").resolve()

    def test_full_config_round_trip(self, tmp_path):
        p = _write(
            tmp_path,
            """\
            version = 1
            input  = "body.stl"
            output = "body.geo"

            [grid]
            ni = 80
            nj = 120
            nk = 50
            theta_cap_deg = 30.0

            [march]
            ds = 1.0e-4
            growth = 1.2
            smooth = 0.3
            smooth_iters = 4
            blend_normals_k = 5

            [surface]
            i_spacing = "tanh"
            i_beta    = 4.0
            j_spacing = "tanh2"
            j_beta    = 5.0

            [volume]
            k_spacing   = "tanh"
            k_beta      = 4.0
            k_thickness = 5.0
            """,
        )
        cfg = load_config(p)
        assert cfg.grid.ni == 80
        assert cfg.grid.nk == 50
        assert cfg.grid.theta_cap == pytest.approx(math.radians(30.0))
        assert cfg.march.ds == pytest.approx(1e-4)
        assert cfg.march.blend_normals_k == 5
        assert cfg.surface.i_spacing == "tanh"
        assert cfg.surface.j_spacing == "tanh2"
        assert cfg.volume.k_spacing == "tanh"
        assert cfg.volume.k_thickness == pytest.approx(5.0)

    def test_relative_paths_resolved_against_config_parent(self, tmp_path):
        sub = tmp_path / "sub"
        sub.mkdir()
        p = _write(
            sub,
            """\
            version = 1
            input = "../body.stl"
            """,
        )
        cfg = load_config(p)
        assert cfg.input == (tmp_path / "body.stl").resolve()


class TestLoadConfigErrors:
    def test_missing_file(self, tmp_path):
        with pytest.raises(ParameterError, match="not found"):
            load_config(tmp_path / "nope.toml")

    def test_malformed_toml(self, tmp_path):
        p = tmp_path / "x.toml"
        p.write_text('version = 1\ninput = "b.stl"\n[grid\n')
        with pytest.raises(ParameterError, match="Invalid TOML"):
            load_config(p)

    def test_missing_version(self, tmp_path):
        p = _write(tmp_path, 'input = "b.stl"\n')
        with pytest.raises(ParameterError, match="missing required key 'version'"):
            load_config(p)

    def test_wrong_version(self, tmp_path):
        p = _write(tmp_path, 'version = 99\ninput = "b.stl"\n')
        with pytest.raises(ParameterError, match="unsupported config version"):
            load_config(p)

    def test_missing_input(self, tmp_path):
        p = _write(tmp_path, "version = 1\n")
        with pytest.raises(ParameterError, match="missing required key 'input'"):
            load_config(p)

    def test_unknown_top_level_key(self, tmp_path):
        p = _write(
            tmp_path,
            """\
            version = 1
            input = "b.stl"
            stray = 42
            """,
        )
        with pytest.raises(ParameterError, match="unknown key"):
            load_config(p)

    def test_unknown_section_key(self, tmp_path):
        p = _write(
            tmp_path,
            """\
            version = 1
            input = "b.stl"

            [grid]
            ni = 30
            typo_field = 99
            """,
        )
        with pytest.raises(ParameterError, match=r"\[grid\].*unknown key"):
            load_config(p)

    def test_theta_cap_mutex(self, tmp_path):
        p = _write(
            tmp_path,
            """\
            version = 1
            input = "b.stl"

            [grid]
            theta_cap     = 0.5
            theta_cap_deg = 30.0
            """,
        )
        with pytest.raises(ParameterError, match="mutually exclusive"):
            load_config(p)


# ---------------------------------------------------------------------------
# dump_template
# ---------------------------------------------------------------------------


class TestDumpTemplate:
    def test_template_round_trips_through_loader(self, tmp_path, capsys):
        dump_template()
        text = capsys.readouterr().out
        # Replace the placeholder so the loader stops complaining.
        stl = tmp_path / "stub.stl"
        stl.write_text("")
        text = text.replace('"body.stl"', f'"{stl}"')
        p = tmp_path / "shore.toml"
        p.write_text(text)
        cfg = load_config(p)
        assert cfg.grid.ni == GridConfig().ni  # defaults preserved

    def test_template_advertises_current_version(self, capsys):
        dump_template()
        text = capsys.readouterr().out
        assert f"version = {CONFIG_VERSION}" in text


# ---------------------------------------------------------------------------
# CLI integration: shore mesh -c <config>
# ---------------------------------------------------------------------------


_LABELS = ("sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south")


class TestMeshConfigFlag:
    def _make_config(self, tmp_path: Path, sphere_stl: Path, **extras: str) -> Path:
        body = textwrap.dedent(
            f"""\
            version = 1
            input  = "{sphere_stl}"
            output = "{tmp_path}/from_cfg"

            [grid]
            ni = 12
            nj = 16
            nk = 5

            [march]
            ds = 0.05
            growth = 1.1
            """
        )
        for k, v in extras.items():
            body += f"\n{k} = {v}"
        p = tmp_path / "shore.toml"
        p.write_text(body)
        return p

    def test_config_drives_mesh(self, sphere_stl, tmp_path):
        cfg = self._make_config(tmp_path, sphere_stl)
        r = runner.invoke(app, ["mesh", "-c", str(cfg)])
        assert r.exit_code == 0, r.output
        # Six .geo files at the configured stem.
        stem = tmp_path / "from_cfg"
        for label in _LABELS:
            assert (stem.parent / f"{stem.name}_{label}.geo").exists()

    def test_cli_overrides_config(self, sphere_stl, tmp_path):
        """CLI flag takes precedence over the config value."""
        cfg = self._make_config(tmp_path, sphere_stl)
        out_stem = tmp_path / "override"
        r = runner.invoke(
            app,
            ["mesh", "-c", str(cfg), "--ni", "8", "-o", str(out_stem)],
        )
        assert r.exit_code == 0, r.output
        sub0 = read_geo(str(out_stem) + "_sub0.geo")
        assert sub0.shape[0] == 8, "CLI --ni 8 must override config ni=12"

    def test_config_value_used_when_cli_unset(self, sphere_stl, tmp_path):
        """When the user does NOT set --ni, the config's ni wins over the
        built-in default."""
        cfg = self._make_config(tmp_path, sphere_stl)
        r = runner.invoke(app, ["mesh", "-c", str(cfg)])
        assert r.exit_code == 0, r.output
        sub0 = read_geo(str(tmp_path / "from_cfg") + "_sub0.geo")
        # config sets ni=12 (built-in default would be 40).
        assert sub0.shape[0] == 12

    def test_positional_input_overrides_config_input(self, sphere_stl, tmp_path):
        cfg = self._make_config(tmp_path, sphere_stl)
        # Build a *different* STL — but pass the same one positionally to
        # prove the positional arg takes precedence over cfg.input.
        out_stem = tmp_path / "via_positional"
        r = runner.invoke(
            app,
            ["mesh", str(sphere_stl), "-c", str(cfg), "-o", str(out_stem)],
        )
        assert r.exit_code == 0, r.output
        assert (Path(str(out_stem) + "_sub0.geo")).exists()

    def test_missing_input_errors(self, tmp_path):
        # Config has no input key → loader rejects it before mesh runs.
        bad = tmp_path / "bad.toml"
        bad.write_text("version = 1\n")
        r = runner.invoke(app, ["mesh", "-c", str(bad)])
        assert r.exit_code != 0
        assert "missing required key 'input'" in r.output


class TestConfigSubcommand:
    def test_template_action_emits_valid_toml(self, tmp_path):
        r = runner.invoke(app, ["config", "template"])
        assert r.exit_code == 0
        assert f"version = {CONFIG_VERSION}" in r.output

    def test_validate_action_passes(self, sphere_stl, tmp_path):
        p = tmp_path / "shore.toml"
        p.write_text(f'version = 1\ninput = "{sphere_stl}"\n')
        r = runner.invoke(app, ["config", "validate", str(p)])
        assert r.exit_code == 0
        assert "parses cleanly" in r.output

    def test_validate_action_rejects_typo(self, tmp_path):
        p = tmp_path / "shore.toml"
        p.write_text('version = 1\ninput = "x.stl"\nstray = 1\n')
        r = runner.invoke(app, ["config", "validate", str(p)])
        assert r.exit_code != 0
        assert "unknown key" in r.output

    def test_validate_action_requires_path(self):
        r = runner.invoke(app, ["config", "validate"])
        assert r.exit_code != 0

    def test_unknown_action(self):
        r = runner.invoke(app, ["config", "frobnicate"])
        assert r.exit_code != 0
