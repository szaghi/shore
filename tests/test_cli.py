"""End-to-end CLI tests using a synthetic sphere STL.

All commands are exercised via typer.testing.CliRunner so no subprocess
is spawned and coverage is collected normally.

Grid parameters are kept small (ni=10, nj=12, nk=6) to keep the suite fast.
The mesh command produces 6 .geo files (cubed-sphere topology):
  <stem>_sub0.geo .. <stem>_sub3.geo, <stem>_cap_north.geo, <stem>_cap_south.geo
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.cli import app
from shore.io.geo import read_geo, write_geo
from tests.conftest import sphere_layer

runner = CliRunner()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SMALL = ["--ni", "10", "--nj", "12", "--nk", "6", "--ds", "0.05", "--growth", "1.1"]
_PROJ = ["--ni", "10", "--nj", "12"]

_LABELS = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]


def _invoke(*args: str):
    return runner.invoke(app, list(args))


def _geo_files(stem: Path) -> list[Path]:
    """Return the 6 .geo files written by `shore mesh` for a given stem."""
    return [stem.parent / f"{stem.stem}_{label}.geo" for label in _LABELS]


# ---------------------------------------------------------------------------
# shore mesh
# ---------------------------------------------------------------------------


class TestMesh:
    def test_exit_code_ok(self, sphere_stl, tmp_path):
        out = tmp_path / "sphere"
        r = _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_creates_six_geo_files(self, sphere_stl, tmp_path):
        out = tmp_path / "sphere"
        _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(out))
        for path in _geo_files(out):
            assert path.exists(), f"missing {path.name}"

    def test_sub_block_dimensions(self, sphere_stl, tmp_path):
        out = tmp_path / "sphere"
        _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(out))
        # nj=12 → nc = 12//4+1 = 4; sub-blocks: (10, 4, 6, 3)
        for label in ["sub0", "sub1", "sub2", "sub3"]:
            grid = read_geo(tmp_path / f"sphere_{label}.geo")
            assert grid.shape == (10, 4, 6, 3), f"{label}: {grid.shape}"

    def test_cap_block_dimensions(self, sphere_stl, tmp_path):
        out = tmp_path / "sphere"
        _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(out))
        # caps: (nc, nc, nk, 3) = (4, 4, 6, 3)
        for label in ["cap_north", "cap_south"]:
            grid = read_geo(tmp_path / f"sphere_{label}.geo")
            assert grid.shape == (4, 4, 6, 3), f"{label}: {grid.shape}"

    def test_outer_shell_larger_than_inner(self, sphere_stl, tmp_path):
        out = tmp_path / "sphere"
        _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(out))
        for label in _LABELS:
            grid = read_geo(tmp_path / f"sphere_{label}.geo")
            r_inner = np.linalg.norm(grid[:, :, 0, :], axis=-1)
            r_outer = np.linalg.norm(grid[:, :, -1, :], axis=-1)
            assert np.all(r_outer > r_inner), f"{label}: outer not larger than inner"

    def test_missing_input_exits_1(self, tmp_path):
        r = _invoke("mesh", str(tmp_path / "ghost.stl"), *_SMALL)
        assert r.exit_code == 1

    def test_verbose_flag(self, sphere_stl, tmp_path):
        out = tmp_path / "sphere"
        r = _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(out), "-v")
        assert r.exit_code == 0, r.output


class TestMeshCGrid:
    """`shore mesh --topology cgrid` end-to-end."""

    _CGRID: ClassVar[list[str]] = [
        "--topology",
        "cgrid",
        "--ni-body",
        "60",
        "--ni-wake",
        "15",
        "--n-stations",
        "4",
        "--wake-length",
        "10.0",
        "--nk",
        "8",
        "--ds",
        "1e-3",
        "--smooth",
        "0.0",
        "--smooth-iters",
        "0",
    ]

    def test_exit_code_ok(self, naca0012_stl, tmp_path):
        out = tmp_path / "cgrid"
        r = _invoke("mesh", str(naca0012_stl), *self._CGRID, "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_writes_single_geo(self, naca0012_stl, tmp_path):
        out = tmp_path / "cgrid"
        _invoke("mesh", str(naca0012_stl), *self._CGRID, "-o", str(out))
        path = tmp_path / "cgrid.geo"
        assert path.exists()

    def test_block_shape(self, naca0012_stl, tmp_path):
        out = tmp_path / "cgrid"
        _invoke("mesh", str(naca0012_stl), *self._CGRID, "-o", str(out))
        grid = read_geo(tmp_path / "cgrid.geo")
        # 2 * (ni_wake - 1) + ni_body = 28 + 60 = 88; nj=4; nk=8
        assert grid.shape == (88, 4, 8, 3)

    def test_missing_required_flag_exits_1(self, naca0012_stl, tmp_path):
        # Missing --wake-length.
        r = _invoke(
            "mesh",
            str(naca0012_stl),
            "--topology",
            "cgrid",
            "--ni-body",
            "60",
            "--ni-wake",
            "15",
            "--n-stations",
            "4",
        )
        assert r.exit_code == 1

    def test_unknown_topology_exits_1(self, sphere_stl, tmp_path):
        r = _invoke("mesh", str(sphere_stl), "--topology", "garbage")
        assert r.exit_code == 1


# ---------------------------------------------------------------------------
# shore project
# ---------------------------------------------------------------------------


class TestProject:
    def test_exit_code_ok(self, sphere_stl, tmp_path):
        out = tmp_path / "k0.npz"
        r = _invoke("project", str(sphere_stl), *_PROJ, "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_creates_npz(self, sphere_stl, tmp_path):
        out = tmp_path / "k0.npz"
        _invoke("project", str(sphere_stl), *_PROJ, "-o", str(out))
        assert out.exists()

    def test_npz_shape(self, sphere_stl, tmp_path):
        out = tmp_path / "k0.npz"
        _invoke("project", str(sphere_stl), *_PROJ, "-o", str(out))
        data = np.load(out)
        assert "k0" in data
        assert data["k0"].shape == (10, 12, 3)

    def test_points_on_sphere_surface(self, sphere_stl, tmp_path):
        out = tmp_path / "k0.npz"
        _invoke("project", str(sphere_stl), *_PROJ, "-o", str(out))
        k0 = np.load(out)["k0"]
        radii = np.linalg.norm(k0, axis=-1)
        assert np.allclose(radii, 1.0, atol=0.05)

    def test_missing_input_exits_1(self, tmp_path):
        r = _invoke("project", str(tmp_path / "ghost.stl"), *_PROJ)
        assert r.exit_code == 1


# ---------------------------------------------------------------------------
# shore info
# ---------------------------------------------------------------------------


class TestInfo:
    @pytest.fixture
    def small_geo(self, tmp_path) -> Path:
        ni, nj, nk = 4, 5, 3
        i = np.arange(ni, dtype=float)
        j = np.arange(nj, dtype=float)
        k = np.arange(nk, dtype=float)
        ig, jg, kg = np.meshgrid(i, j, k, indexing="ij")
        grid = np.stack([ig, jg, kg], axis=-1)
        p = tmp_path / "small.geo"
        write_geo(p, grid)
        return p

    @pytest.fixture
    def small_npz(self, tmp_path) -> Path:
        k0 = sphere_layer(6, 8)
        p = tmp_path / "surf.npz"
        np.savez_compressed(p, k0=k0)
        return p

    def test_geo_exit_code_ok(self, small_geo):
        r = _invoke("info", str(small_geo))
        assert r.exit_code == 0, r.output

    def test_geo_shows_dimensions(self, small_geo):
        r = _invoke("info", str(small_geo))
        assert "ni=4" in r.output
        assert "nj=5" in r.output
        assert "nk=3" in r.output

    def test_geo_shows_point_count(self, small_geo):
        r = _invoke("info", str(small_geo))
        assert "60" in r.output  # 4*5*3

    def test_npz_exit_code_ok(self, small_npz):
        r = _invoke("info", str(small_npz))
        assert r.exit_code == 0, r.output

    def test_npz_shows_dimensions(self, small_npz):
        r = _invoke("info", str(small_npz))
        assert "ni=6" in r.output
        assert "nj=8" in r.output

    def test_missing_input_exits_1(self, tmp_path):
        r = _invoke("info", str(tmp_path / "ghost.geo"))
        assert r.exit_code == 1

    def test_bad_extension_exits_1(self, tmp_path):
        p = tmp_path / "file.txt"
        p.write_text("nothing")
        r = _invoke("info", str(p))
        assert r.exit_code == 1


# ---------------------------------------------------------------------------
# shore check
# ---------------------------------------------------------------------------


class TestCheck:
    @pytest.fixture
    def valid_geo(self, sphere_stl, tmp_path) -> Path:
        """One sub-block geo file from a valid mesh."""
        _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(tmp_path / "sphere"))
        return tmp_path / "sphere_sub0.geo"

    @pytest.fixture
    def inverted_geo(self, tmp_path) -> Path:
        ni, nj, nk = 6, 8, 3
        surface = sphere_layer(ni, nj, radius=1.0)
        grid = np.empty((ni, nj, nk, 3))
        grid[:, :, 0, :] = surface
        grid[:, :, 1, :] = surface * 1.1
        grid[:, :, 2, :] = surface * 0.5  # collapsed inward → inversion
        p = tmp_path / "inverted.geo"
        write_geo(p, grid)
        return p

    def test_valid_exits_0(self, valid_geo):
        r = _invoke("check", str(valid_geo))
        assert r.exit_code == 0, r.output

    def test_valid_all_ok_in_output(self, valid_geo):
        r = _invoke("check", str(valid_geo))
        assert "ok" in r.output

    def test_inverted_exits_1(self, inverted_geo):
        r = _invoke("check", str(inverted_geo))
        assert r.exit_code == 1

    def test_inverted_shows_inverted_label(self, inverted_geo):
        r = _invoke("check", str(inverted_geo))
        assert "INVERTED" in r.output

    def test_missing_input_exits_1(self, tmp_path):
        r = _invoke("check", str(tmp_path / "ghost.geo"))
        assert r.exit_code == 1


# ---------------------------------------------------------------------------
# shore export
# ---------------------------------------------------------------------------


class TestExport:
    @pytest.fixture
    def valid_geo(self, sphere_stl, tmp_path) -> Path:
        _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(tmp_path / "sphere"))
        return tmp_path / "sphere_sub0.geo"

    @pytest.fixture
    def valid_npz(self, sphere_stl, tmp_path) -> Path:
        out = tmp_path / "sphere_k0.npz"
        _invoke("project", str(sphere_stl), *_PROJ, "-o", str(out))
        return out

    def test_geo_creates_vts(self, valid_geo, tmp_path):
        out = tmp_path / "sphere.vts"
        r = _invoke("export", str(valid_geo), "-o", str(out))
        assert r.exit_code == 0, r.output
        assert out.exists()
        assert out.stat().st_size > 0

    def test_npz_creates_vts(self, valid_npz, tmp_path):
        out = tmp_path / "sphere_k0.vts"
        r = _invoke("export", str(valid_npz), "-o", str(out))
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_default_output_name_geo(self, valid_geo):
        r = _invoke("export", str(valid_geo))
        assert r.exit_code == 0, r.output
        expected = valid_geo.with_suffix(".vts")
        assert expected.exists()

    def test_missing_input_exits_1(self, tmp_path):
        r = _invoke("export", str(tmp_path / "ghost.geo"))
        assert r.exit_code == 1

    def test_bad_extension_exits_1(self, tmp_path):
        p = tmp_path / "file.txt"
        p.write_text("nothing")
        r = _invoke("export", str(p))
        assert r.exit_code == 1

    @pytest.fixture
    def valid_geo_set(self, sphere_stl, tmp_path) -> list[Path]:
        """Six .geo files from a single ``shore mesh`` run."""
        _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(tmp_path / "sphere"))
        return [tmp_path / f"sphere_{lbl}.geo" for lbl in _LABELS]

    def test_multi_geo_creates_vtm(self, valid_geo_set, tmp_path):
        out = tmp_path / "sphere.vtm"
        r = _invoke("export", *[str(p) for p in valid_geo_set], "-o", str(out))
        assert r.exit_code == 0, r.output
        assert out.exists()
        assert out.stat().st_size > 0
        # The user-visible message should list one block per input.
        for lbl in _LABELS:
            assert f"sphere_{lbl}" in r.output

    def test_multi_geo_default_output_is_first_stem_dot_vtm(self, valid_geo_set, tmp_path):
        r = _invoke("export", *[str(p) for p in valid_geo_set])
        assert r.exit_code == 0, r.output
        # Default .vtm path: first input's stem with .vtm suffix.
        expected = valid_geo_set[0].with_suffix(".vtm")
        assert expected.exists()

    def test_multi_input_rejects_npz(self, valid_geo_set, valid_npz, tmp_path):
        r = _invoke("export", str(valid_geo_set[0]), str(valid_npz), "-o", str(tmp_path / "x.vtm"))
        assert r.exit_code != 0
        assert ".geo only" in r.output

    def test_multi_input_missing_file_exits_1(self, valid_geo_set, tmp_path):
        ghost = tmp_path / "missing.geo"
        r = _invoke("export", str(valid_geo_set[0]), str(ghost), "-o", str(tmp_path / "x.vtm"))
        assert r.exit_code == 1
        assert "file not found" in r.output


# ---------------------------------------------------------------------------
# shore profile
# ---------------------------------------------------------------------------


class TestProfileNaca:
    def test_4digit_writes_dat(self, tmp_path):
        out = tmp_path / "naca0012.dat"
        r = _invoke("profile", "naca", "0012", "-o", str(out))
        assert r.exit_code == 0, r.output
        assert out.exists()
        # Header line + data lines; 159 = 2*80 - 1 chord-point coords.
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 1 + 159  # name + coords

    def test_5digit_writes_dat(self, tmp_path):
        out = tmp_path / "naca23012.dat"
        r = _invoke("profile", "naca", "23012", "--n-pts", "40", "-o", str(out))
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_open_te_flag(self, tmp_path):
        out = tmp_path / "naca0012-open.dat"
        r = _invoke("profile", "naca", "0012", "--open-te", "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_invalid_designation_exits_1(self, tmp_path):
        r = _invoke("profile", "naca", "abcd", "-o", str(tmp_path / "x.dat"))
        assert r.exit_code == 1
        assert "digit string" in r.output


# ---------------------------------------------------------------------------
# shore body
# ---------------------------------------------------------------------------


class TestBody:
    def test_sphere(self, tmp_path):
        out = tmp_path / "sphere.stl"
        r = _invoke("body", "sphere", "--radius", "1.5", "-o", str(out))
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_cylinder(self, tmp_path):
        out = tmp_path / "cyl.stl"
        r = _invoke("body", "cylinder", "--radius", "0.3", "--height", "2.0", "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_torus(self, tmp_path):
        out = tmp_path / "torus.stl"
        r = _invoke("body", "torus", "--r-major", "2", "--r-minor", "0.5", "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_box(self, tmp_path):
        out = tmp_path / "box.stl"
        r = _invoke("body", "box", "--extents", "1", "2", "3", "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_capsule(self, tmp_path):
        out = tmp_path / "cap.stl"
        r = _invoke("body", "capsule", "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_invalid_radius_exits_1(self, tmp_path):
        r = _invoke("body", "sphere", "--radius", "-1", "-o", str(tmp_path / "x.stl"))
        assert r.exit_code == 1

    def test_extrude_naca_pipeline(self, tmp_path):
        # 1) generate Selig .dat via shore profile
        dat = tmp_path / "naca0012.dat"
        r1 = _invoke("profile", "naca", "0012", "-o", str(dat))
        assert r1.exit_code == 0, r1.output
        # 2) extrude into a wing STL via shore body extrude
        stl = tmp_path / "wing.stl"
        r2 = _invoke(
            "body",
            "extrude",
            "--section",
            str(dat),
            "--span",
            "1.0",
            "--n-stations",
            "6",
            "--end-caps",
            "-o",
            str(stl),
        )
        assert r2.exit_code == 0, r2.output
        assert stl.exists()

    def test_extrude_missing_section_exits_1(self, tmp_path):
        ghost = tmp_path / "ghost.dat"
        r = _invoke(
            "body",
            "extrude",
            "--section",
            str(ghost),
            "--span",
            "1.0",
            "--n-stations",
            "5",
            "-o",
            str(tmp_path / "x.stl"),
        )
        assert r.exit_code == 1
        assert "file not found" in r.output
