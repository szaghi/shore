"""Tests for the single-block O-grid topology.

Covers:
- ``OGridTopology.build`` from a Body.
- ``Mesh.march`` dispatching to the O-grid marcher.
- Geometric sanity (j-periodicity, k=0 on the sphere, outer > inner).
- Config-driven invocation through the CLI.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.cli import app
from shore.io.geo import read_geo
from shore.mesh import Body, Mesh, OGridTopology, Spacing1D
from shore.volume.jacobian import jacobian_field

runner = CliRunner()


# ---------------------------------------------------------------------------
# OGridTopology.build
# ---------------------------------------------------------------------------


class TestOGridBuild:
    def test_returns_single_block(self, sphere_stl):
        body = Body.from_stl(sphere_stl)
        blocks = OGridTopology().build(
            surface=None,
            nk=5,
            ni=10,
            nj=12,
            spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
            projector=body.projector,
        )
        assert len(blocks) == 1
        assert blocks[0].label == "ogrid"

    def test_node_shape(self, sphere_stl):
        body = Body.from_stl(sphere_stl)
        (block,) = OGridTopology().build(
            surface=None,
            nk=5,
            ni=10,
            nj=12,
            spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
            projector=body.projector,
        )
        assert block.nodes.shape == (10, 12, 5, 3)

    def test_periodicity_flags(self, sphere_stl):
        body = Body.from_stl(sphere_stl)
        (block,) = OGridTopology().build(surface=None, nk=4, ni=8, nj=12, projector=body.projector)
        assert block.periodic_i is False
        assert block.periodic_j is True

    def test_k0_on_unit_sphere(self, sphere_stl):
        body = Body.from_stl(sphere_stl)
        (block,) = OGridTopology().build(
            surface=None,
            nk=4,
            ni=10,
            nj=16,
            theta_cap=float(np.pi / 12),
            projector=body.projector,
        )
        radii = np.linalg.norm(block.nodes[:, :, 0, :], axis=-1)
        # Icosphere is not perfectly spherical; 3 % tolerance.
        assert radii.min() > 0.95
        assert radii.max() < 1.05

    def test_accepts_pre_projected_surface(self, sphere_stl):
        body = Body.from_stl(sphere_stl)
        surface = body.project(8, 12, theta_cap=float(np.pi / 12))
        (block,) = OGridTopology().build(
            surface=surface,
            nk=4,
            spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
        )
        assert block.nodes.shape == (8, 12, 4, 3)
        np.testing.assert_array_equal(block.nodes[:, :, 0, :], surface)

    def test_missing_projector_errors(self):
        with pytest.raises(ValueError, match="requires a Projector"):
            OGridTopology().build(surface=None, nk=4, ni=8, nj=12)

    def test_missing_ni_nj_errors(self, sphere_stl):
        body = Body.from_stl(sphere_stl)
        with pytest.raises(ValueError, match="needs ni and nj"):
            OGridTopology().build(surface=None, nk=4, projector=body.projector)


# ---------------------------------------------------------------------------
# Mesh.march dispatch
# ---------------------------------------------------------------------------


class TestOGridMarch:
    def test_march_runs_to_completion(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=OGridTopology(),
            ni=20,
            nj=32,
            nk=5,
            theta_cap=float(np.pi / 6),
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        m.march(ds=0.05, growth=1.1)
        # All k layers populated.
        assert not np.any(np.isnan(m.blocks[0].nodes))

    def test_outer_layer_strictly_outside_inner(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=OGridTopology(),
            ni=20,
            nj=32,
            nk=5,
            theta_cap=float(np.pi / 6),
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        m.march(ds=0.05, growth=1.1)
        block = m.blocks[0]
        r_inner = np.linalg.norm(block.nodes[:, :, 0, :], axis=-1)
        r_outer = np.linalg.norm(block.nodes[:, :, -1, :], axis=-1)
        # Every node on the outer layer must be outside the corresponding
        # k=0 node (radial extrusion grows monotonically).
        assert np.all(r_outer > r_inner)

    def test_jacobians_positive_after_march(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=OGridTopology(),
            ni=20,
            nj=32,
            nk=5,
            theta_cap=float(np.pi / 6),
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        m.march(ds=0.05, growth=1.1)
        # j-periodic Jacobian field — every cell strictly positive (modulo
        # FP noise at the polar singular axis rows).
        jac = jacobian_field(m.blocks[0].nodes)
        assert jac.min() > -1e-6


# ---------------------------------------------------------------------------
# CLI: --topology ogrid
# ---------------------------------------------------------------------------


_OGRID_FLAGS = [
    "--topology",
    "ogrid",
    "--ni",
    "20",
    "--nj",
    "32",
    "--nk",
    "5",
    "--ds",
    "0.05",
    "--growth",
    "1.1",
    "--theta-cap-deg",
    "30",
]


class TestOGridCLI:
    def test_writes_single_geo(self, sphere_stl, tmp_path):
        out = tmp_path / "og"
        r = runner.invoke(app, ["mesh", str(sphere_stl), *_OGRID_FLAGS, "-o", str(out)])
        assert r.exit_code == 0, r.output
        assert (tmp_path / "og.geo").exists()
        # And does NOT write any sub*/cap* files.
        for ghost in ("og_sub0.geo", "og_cap_north.geo", "og_ogrid.geo"):
            assert not (tmp_path / ghost).exists()

    def test_geo_has_expected_shape(self, sphere_stl, tmp_path):
        out = tmp_path / "og"
        runner.invoke(app, ["mesh", str(sphere_stl), *_OGRID_FLAGS, "-o", str(out)])
        nodes = read_geo(tmp_path / "og.geo")
        assert nodes.shape == (20, 32, 5, 3)

    def test_default_topology_is_cubed_sphere(self, sphere_stl, tmp_path):
        # No --topology flag → cubed_sphere; 6 files, no <stem>.geo.
        out = tmp_path / "cs"
        cs_flags = [
            "--ni",
            "20",
            "--nj",
            "32",
            "--nk",
            "5",
            "--ds",
            "0.05",
            "--growth",
            "1.1",
            "--theta-cap-deg",
            "30",
        ]
        r = runner.invoke(app, ["mesh", str(sphere_stl), *cs_flags, "-o", str(out)])
        assert r.exit_code == 0, r.output
        assert (tmp_path / "cs_sub0.geo").exists()
        assert not (tmp_path / "cs.geo").exists()

    def test_unknown_topology_rejected(self, sphere_stl, tmp_path):
        out = tmp_path / "bad"
        r = runner.invoke(
            app,
            [
                "mesh",
                str(sphere_stl),
                "--topology",
                "wibble",
                "--ni",
                "10",
                "--nj",
                "12",
                "--nk",
                "4",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code != 0
        assert "unknown --topology" in r.output


# ---------------------------------------------------------------------------
# Config-driven invocation
# ---------------------------------------------------------------------------


class TestOGridConfig:
    def _write_config(self, tmp_path: Path, sphere_stl: Path) -> Path:
        body = textwrap.dedent(
            f"""
            version = 1
            input = "{sphere_stl}"

            [grid]
            ni = 20
            nj = 32
            nk = 5
            theta_cap_deg = 30.0
            topology = "ogrid"

            [march]
            ds = 0.05
            growth = 1.1
            """
        )
        p = tmp_path / "ogrid.toml"
        p.write_text(body)
        return p

    def test_config_drives_ogrid_mesh(self, sphere_stl, tmp_path):
        cfg = self._write_config(tmp_path, sphere_stl)
        out = tmp_path / "og"
        r = runner.invoke(app, ["mesh", "-c", str(cfg), "-o", str(out)])
        assert r.exit_code == 0, r.output
        assert (tmp_path / "og.geo").exists()
        assert read_geo(tmp_path / "og.geo").shape == (20, 32, 5, 3)

    def test_cli_topology_overrides_config(self, sphere_stl, tmp_path):
        # Config says ogrid; CLI says cubed_sphere — CLI wins.
        cfg = self._write_config(tmp_path, sphere_stl)
        out = tmp_path / "cs"
        r = runner.invoke(
            app,
            ["mesh", "-c", str(cfg), "--topology", "cubed_sphere", "-o", str(out)],
        )
        assert r.exit_code == 0, r.output
        assert (tmp_path / "cs_sub0.geo").exists()
        assert not (tmp_path / "cs.geo").exists()

    def test_unknown_topology_in_config_rejected(self, sphere_stl, tmp_path):
        bad = tmp_path / "bad.toml"
        bad.write_text(
            textwrap.dedent(
                f"""
                version = 1
                input = "{sphere_stl}"
                [grid]
                topology = "wibble"
                """
            )
        )
        r = runner.invoke(app, ["mesh", "-c", str(bad), "-o", str(tmp_path / "x")])
        assert r.exit_code != 0
        assert "unknown topology" in r.output
