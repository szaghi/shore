"""End-to-end CLI tests using a synthetic sphere STL.

All commands are exercised via typer.testing.CliRunner so no subprocess
is spawned and coverage is collected normally.

Grid parameters are kept small (ni=10, nj=12, nk=6) to keep the suite fast.
The mesh command produces 6 .geo files (caps topology):
  <stem>_sub0.geo .. <stem>_sub3.geo, <stem>_cap_north.geo, <stem>_cap_south.geo
"""

from __future__ import annotations

from pathlib import Path

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
_PROJ  = ["--ni", "10", "--nj", "12"]

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
            r_inner = np.linalg.norm(grid[:, :, 0,  :], axis=-1)
            r_outer = np.linalg.norm(grid[:, :, -1, :], axis=-1)
            assert np.all(r_outer > r_inner), f"{label}: outer not larger than inner"

    def test_missing_input_exits_1(self, tmp_path):
        r = _invoke("mesh", str(tmp_path / "ghost.stl"), *_SMALL)
        assert r.exit_code == 1

    def test_verbose_flag(self, sphere_stl, tmp_path):
        out = tmp_path / "sphere"
        r = _invoke("mesh", str(sphere_stl), *_SMALL, "-o", str(out), "-v")
        assert r.exit_code == 0, r.output


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
        assert "60" in r.output   # 4*5*3

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
        grid[:, :, 2, :] = surface * 0.5   # collapsed inward → inversion
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
