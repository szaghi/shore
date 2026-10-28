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

    def test_writes_three_geo(self, naca0012_stl, tmp_path):
        # Default --topology cgrid emits 3 blocks (upper / body / lower).
        out = tmp_path / "cgrid"
        _invoke("mesh", str(naca0012_stl), *self._CGRID, "-o", str(out))
        for label in ("upper", "body", "lower"):
            assert (tmp_path / f"cgrid_{label}.geo").exists(), f"missing {label}.geo"

    def test_block_shapes(self, naca0012_stl, tmp_path):
        # ni_body=60, ni_wake=15, n_stations=4, nk=8.
        # upper / lower: (15, 4, 8, 3); body: (60, 4, 8, 3).
        out = tmp_path / "cgrid"
        _invoke("mesh", str(naca0012_stl), *self._CGRID, "-o", str(out))
        upper = read_geo(tmp_path / "cgrid_upper.geo")
        body = read_geo(tmp_path / "cgrid_body.geo")
        lower = read_geo(tmp_path / "cgrid_lower.geo")
        assert upper.shape == (15, 4, 8, 3)
        assert body.shape == (60, 4, 8, 3)
        assert lower.shape == (15, 4, 8, 3)

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

    def test_cgrid_blocks_1_writes_single_geo(self, naca0012_stl, tmp_path):
        # Legacy single-block C-grid: --cgrid-blocks 1 emits one .geo.
        out = tmp_path / "legacy"
        r = _invoke(
            "mesh",
            str(naca0012_stl),
            *self._CGRID,
            "--cgrid-blocks",
            "1",
            "-o",
            str(out),
        )
        assert r.exit_code == 0, r.output
        assert (tmp_path / "legacy.geo").exists()
        for label in ("upper", "body", "lower"):
            assert not (tmp_path / f"legacy_{label}.geo").exists()

    def test_cgrid_blocks_4_rejected_for_sharp_te(self, naca0012_stl, tmp_path):
        # 4-block layout requires a blunt TE; the sharp NACA 0012 fixture
        # must be rejected with a clear error.
        r = _invoke(
            "mesh",
            str(naca0012_stl),
            *self._CGRID,
            "--cgrid-blocks",
            "4",
            "-o",
            str(tmp_path / "x"),
        )
        assert r.exit_code == 1
        assert "blunt" in r.output.lower() or "sharp" in r.output.lower()

    def test_cgrid_blocks_invalid_rejected(self, naca0012_stl, tmp_path):
        r = _invoke(
            "mesh",
            str(naca0012_stl),
            *self._CGRID,
            "--cgrid-blocks",
            "2",
            "-o",
            str(tmp_path / "x"),
        )
        assert r.exit_code == 1


class TestMeshCH:
    """`shore mesh --topology ch` end-to-end."""

    _CH: ClassVar[list[str]] = [
        "--topology",
        "ch",
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
        "--ch-chord-min",
        "-2.0",
        "--ch-chord-max",
        "15.0",
        "--ch-normal-min",
        "-2.5",
        "--ch-normal-max",
        "2.5",
        "--nk-channel",
        "4",
    ]

    def test_exit_code_ok(self, naca0012_stl, tmp_path):
        out = tmp_path / "ch"
        r = _invoke("mesh", str(naca0012_stl), *self._CH, "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_writes_12_geo_for_sharp_te(self, naca0012_stl, tmp_path):
        out = tmp_path / "ch"
        _invoke("mesh", str(naca0012_stl), *self._CH, "-o", str(out))
        for label in (
            "upper",
            "body",
            "lower",
            "chan_upper",
            "chan_body_top",
            "chan_body_front",
            "chan_body_bot",
            "chan_lower",
            "chan_outflow_top_upper",
            "chan_outflow_wake_upper",
            "chan_outflow_wake_lower",
            "chan_outflow_top_lower",
        ):
            assert (tmp_path / f"ch_{label}.geo").exists(), f"missing {label}.geo"

    def test_missing_box_extent_exits_1(self, naca0012_stl, tmp_path):
        # Missing --ch-normal-max.
        r = _invoke(
            "mesh",
            str(naca0012_stl),
            "--topology",
            "ch",
            "--ni-body",
            "60",
            "--ni-wake",
            "15",
            "--n-stations",
            "4",
            "--wake-length",
            "10.0",
            "--ch-chord-min",
            "-2.0",
            "--ch-chord-max",
            "15.0",
            "--ch-normal-min",
            "-2.5",
        )
        assert r.exit_code == 1

    def test_missing_cgrid_flag_exits_1(self, naca0012_stl, tmp_path):
        # Missing --wake-length (a C-grid wrap requirement).
        r = _invoke(
            "mesh",
            str(naca0012_stl),
            "--topology",
            "ch",
            "--ni-body",
            "60",
            "--ni-wake",
            "15",
            "--n-stations",
            "4",
            "--ch-chord-min",
            "-2.0",
            "--ch-chord-max",
            "15.0",
            "--ch-normal-min",
            "-2.5",
            "--ch-normal-max",
            "2.5",
        )
        assert r.exit_code == 1

    def test_ni_outflow_knob_changes_outflow_shape(self, naca0012_stl, tmp_path):
        """--ch-ni-outflow flows through to the outflow block i-dimension."""
        out = tmp_path / "ch"
        r = _invoke("mesh", str(naca0012_stl), *self._CH, "--ch-ni-outflow", "12", "-o", str(out))
        assert r.exit_code == 0, r.output
        arr = read_geo(str(tmp_path / "ch_chan_outflow_top_lower.geo"))
        assert arr.shape[0] == 12


class TestMeshOH:
    """`shore mesh --topology oh` end-to-end (cubed-sphere wrap + H-channel)."""

    _OH: ClassVar[list[str]] = [
        "--topology",
        "oh",
        "--ni",
        "16",
        "--nj",
        "16",
        "--nk",
        "5",
        "--ds",
        "0.05",
        "--growth",
        "1.1",
        "--oh-box-xmin",
        "-3.0",
        "--oh-box-xmax",
        "3.0",
        "--oh-box-ymin",
        "-3.0",
        "--oh-box-ymax",
        "3.0",
        "--oh-box-zmin",
        "-3.0",
        "--oh-box-zmax",
        "3.0",
        "--nk-channel",
        "4",
    ]

    _OH_LABELS: ClassVar[list[str]] = [
        "sub0",
        "sub1",
        "sub2",
        "sub3",
        "cap_north",
        "cap_south",
        "chan_sub0",
        "chan_sub1",
        "chan_sub2",
        "chan_sub3",
        "chan_north",
        "chan_south",
    ]

    def test_exit_code_ok(self, sphere_stl, tmp_path):
        out = tmp_path / "oh"
        r = _invoke("mesh", str(sphere_stl), *self._OH, "-o", str(out))
        assert r.exit_code == 0, r.output

    def test_writes_12_geo(self, sphere_stl, tmp_path):
        out = tmp_path / "oh"
        _invoke("mesh", str(sphere_stl), *self._OH, "-o", str(out))
        for label in self._OH_LABELS:
            assert (tmp_path / f"oh_{label}.geo").exists(), f"missing {label}.geo"

    def test_missing_box_extent_exits_1(self, sphere_stl, tmp_path):
        # Drop --oh-box-zmax.
        args = list(self._OH)
        i = args.index("--oh-box-zmax")
        del args[i : i + 2]
        r = _invoke("mesh", str(sphere_stl), *args, "-o", str(tmp_path / "oh"))
        assert r.exit_code == 1

    def test_inverted_box_exits_1(self, sphere_stl, tmp_path):
        # xmin >= xmax.
        args = list(self._OH)
        args[args.index("--oh-box-xmin") + 1] = "3.0"
        args[args.index("--oh-box-xmax") + 1] = "-3.0"
        r = _invoke("mesh", str(sphere_stl), *args, "-o", str(tmp_path / "oh"))
        assert r.exit_code == 1

    def test_theta_cap_warns_and_is_ignored(self, sphere_stl, tmp_path):
        out = tmp_path / "oh"
        with pytest.warns(UserWarning, match="ignored for --topology oh"):
            r = _invoke("mesh", str(sphere_stl), *self._OH, "--theta-cap", "0.3", "-o", str(out))
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

    def test_single_geo_with_vtm_output_creates_vtm(self, valid_geo, tmp_path):
        # Single-block topologies (O-grid, C-grid) sometimes want a .vtm
        # for consistency with multi-block exports.  --output X.vtm on a
        # single-input call should write a 1-block MultiBlock.
        out = tmp_path / "single.vtm"
        r = _invoke("export", str(valid_geo), "-o", str(out))
        assert r.exit_code == 0, r.output
        assert out.exists()


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


# ---------------------------------------------------------------------------
# --adjacency-json sidecar emission
# ---------------------------------------------------------------------------


class TestAdjacencyJsonSidecar:
    """`shore mesh` and `shore primitive *` emit a .adjacency.json on demand."""

    def test_mesh_cubed_sphere_emits_sidecar(self, sphere_stl, tmp_path):
        out = tmp_path / "wall"
        adj = tmp_path / "wall.adjacency.json"
        r = _invoke(
            "mesh",
            str(sphere_stl),
            *_SMALL,
            "-o",
            str(out),
            "--adjacency-json",
            str(adj),
        )
        assert r.exit_code == 0, r.output
        assert adj.exists()

    def test_mesh_sidecar_has_expected_block_count(self, sphere_stl, tmp_path):
        import json

        out = tmp_path / "wall"
        adj = tmp_path / "wall.adjacency.json"
        _invoke(
            "mesh",
            str(sphere_stl),
            *_SMALL,
            "-o",
            str(out),
            "--adjacency-json",
            str(adj),
        )
        payload = json.loads(adj.read_text())
        # Cubed sphere = 4 sub + 2 caps.
        assert len(payload["blocks"]) == 6
        labels = {b["label"] for b in payload["blocks"]}
        assert labels == {"sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"}

    def test_primitive_flat_caps_emits_sidecar(self, tmp_path):
        out = tmp_path / "bg"
        adj = tmp_path / "bg.adjacency.json"
        r = _invoke(
            "primitive",
            "flat-caps",
            "--r-out",
            "5",
            "--z0",
            "-2",
            "--z1",
            "5",
            "--ni",
            "8",
            "--nc",
            "10",
            "--nk",
            "6",
            "-o",
            str(out),
            "--adjacency-json",
            str(adj),
        )
        assert r.exit_code == 0, r.output
        assert adj.exists()

    def test_primitive_flat_caps_sidecar_wires_seams(self, tmp_path):
        import json

        out = tmp_path / "bg"
        adj = tmp_path / "bg.adjacency.json"
        _invoke(
            "primitive",
            "flat-caps",
            "--r-out",
            "5",
            "--z0",
            "-2",
            "--z1",
            "5",
            "--ni",
            "8",
            "--nc",
            "10",
            "--nk",
            "6",
            "-o",
            str(out),
            "--adjacency-json",
            str(adj),
        )
        payload = json.loads(adj.read_text())
        # 1 centre + 4 sub-blocks; SHARED face references between
        # the centre and each sub-block plus sub-sub corner pairs.
        assert len(payload["blocks"]) == 5
        shared = sum(
            1 for b in payload["blocks"] for face in b["faces"].values() if face["bc"] == "SHARED"
        )
        assert shared >= 8

    def test_primitive_box_emits_single_block_sidecar(self, tmp_path):
        import json

        out = tmp_path / "box"
        adj = tmp_path / "box.adjacency.json"
        r = _invoke(
            "primitive",
            "box",
            "--ni",
            "5",
            "--nj",
            "5",
            "--nk",
            "5",
            "-o",
            str(out),
            "--adjacency-json",
            str(adj),
        )
        assert r.exit_code == 0, r.output
        payload = json.loads(adj.read_text())
        assert len(payload["blocks"]) == 1
        # Box is a single block — no SHARED faces.
        no_shared = all(face["bc"] != "SHARED" for face in payload["blocks"][0]["faces"].values())
        assert no_shared


# ---------------------------------------------------------------------------
# shore primitive h-grid-rect / h-grid-tube
# ---------------------------------------------------------------------------


class TestPrimitiveHGridRect:
    def test_writes_geo(self, tmp_path):
        out = tmp_path / "duct"
        r = _invoke(
            "primitive",
            "h-grid-rect",
            "--length",
            "1.0",
            "--width",
            "0.1",
            "--height",
            "0.2",
            "--ni",
            "20",
            "--nj",
            "8",
            "--nk",
            "6",
            "-o",
            str(out),
        )
        assert r.exit_code == 0, r.output
        assert (tmp_path / "duct.geo").exists()

    def test_streamwise_z_block_dimensions(self, tmp_path):
        out = tmp_path / "duct"
        _invoke(
            "primitive",
            "h-grid-rect",
            "--length",
            "1.0",
            "--width",
            "0.1",
            "--height",
            "0.2",
            "--ni",
            "20",
            "--nj",
            "8",
            "--nk",
            "6",
            "--streamwise-axis",
            "z",
            "-o",
            str(out),
        )
        grid = read_geo(tmp_path / "duct.geo")
        assert grid.shape == (20, 8, 6, 3)
        assert abs(grid[..., 2].max() - 1.0) < 1e-12

    def test_adjacency_json_emitted(self, tmp_path):
        import json

        out = tmp_path / "duct"
        adj = tmp_path / "duct.adjacency.json"
        r = _invoke(
            "primitive",
            "h-grid-rect",
            "--length",
            "1.0",
            "--width",
            "0.1",
            "--height",
            "0.2",
            "--ni",
            "10",
            "--nj",
            "5",
            "--nk",
            "5",
            "-o",
            str(out),
            "--adjacency-json",
            str(adj),
        )
        assert r.exit_code == 0, r.output
        payload = json.loads(adj.read_text())
        assert len(payload["blocks"]) == 1
        # H-grid rect has FREE on i-faces, WALL on the four cross-section faces.
        faces = payload["blocks"][0]["faces"]
        assert faces["i_lo"]["bc"] == "FREE"
        assert faces["i_hi"]["bc"] == "FREE"
        assert faces["j_lo"]["bc"] == "WALL"
        assert faces["j_hi"]["bc"] == "WALL"
        assert faces["k_lo"]["bc"] == "WALL"
        assert faces["k_hi"]["bc"] == "WALL"

    def test_invalid_streamwise_axis_exits_1(self, tmp_path):
        r = _invoke(
            "primitive",
            "h-grid-rect",
            "--length",
            "1.0",
            "--width",
            "0.1",
            "--height",
            "0.2",
            "--ni",
            "10",
            "--nj",
            "5",
            "--nk",
            "5",
            "--streamwise-axis",
            "w",
            "-o",
            str(tmp_path / "x"),
        )
        assert r.exit_code == 1


class TestPrimitiveHGridTube:
    def test_writes_geo(self, tmp_path):
        out = tmp_path / "tube"
        r = _invoke(
            "primitive",
            "h-grid-tube",
            "--length",
            "1.0",
            "--radius",
            "0.05",
            "--n-axial",
            "20",
            "--n-circ",
            "16",
            "--n-radial",
            "6",
            "-o",
            str(out),
        )
        assert r.exit_code == 0, r.output
        assert (tmp_path / "tube.geo").exists()

    def test_block_shape(self, tmp_path):
        out = tmp_path / "tube"
        _invoke(
            "primitive",
            "h-grid-tube",
            "--length",
            "1.0",
            "--radius",
            "0.05",
            "--n-axial",
            "20",
            "--n-circ",
            "16",
            "--n-radial",
            "6",
            "-o",
            str(out),
        )
        grid = read_geo(tmp_path / "tube.geo")
        assert grid.shape == (20, 16, 6, 3)

    def test_adjacency_json_emitted(self, tmp_path):
        import json

        out = tmp_path / "tube"
        adj = tmp_path / "tube.adjacency.json"
        r = _invoke(
            "primitive",
            "h-grid-tube",
            "--length",
            "1.0",
            "--radius",
            "0.05",
            "--n-axial",
            "10",
            "--n-circ",
            "8",
            "--n-radial",
            "4",
            "-o",
            str(out),
            "--adjacency-json",
            str(adj),
        )
        assert r.exit_code == 0, r.output
        payload = json.loads(adj.read_text())
        assert len(payload["blocks"]) == 1
        faces = payload["blocks"][0]["faces"]
        # i_lo / i_hi = FREE (inflow/outflow); k_lo (centerline) and k_hi
        # (outer wall) = WALL.
        assert faces["i_lo"]["bc"] == "FREE"
        assert faces["i_hi"]["bc"] == "FREE"
        assert faces["k_lo"]["bc"] == "WALL"
        assert faces["k_hi"]["bc"] == "WALL"

    def test_invalid_radius_exits_1(self, tmp_path):
        r = _invoke(
            "primitive",
            "h-grid-tube",
            "--length",
            "1.0",
            "--radius",
            "0.0",
            "--n-axial",
            "10",
            "--n-circ",
            "8",
            "--n-radial",
            "4",
            "-o",
            str(tmp_path / "x"),
        )
        assert r.exit_code == 1


class TestCollarProjects:
    """CLI project helpers beyond the frozen spec (ADR 0018)."""

    def test_stem_collision_needs_a_name(self, tmp_path):
        import trimesh

        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()
        trimesh.creation.box(extents=(4.0, 4.0, 1.0)).export(str(tmp_path / "a" / "body.stl"))
        trimesh.creation.icosphere(2).export(str(tmp_path / "b" / "body.stl"))
        args = ["--seed", "0.55", "0", "0.02", "--n-loop", "16", "--ds-surface", "0.02",
                "--growth-surface", "1.15", "--m-a", "5", "--m-b", "5", "--ds-wall", "1e-3",
                "--growth-wall", "1.25", "--nk", "8", "--project", str(tmp_path / "p")]  # fmt: skip
        bodies = [str(tmp_path / "a" / "body.stl"), str(tmp_path / "b" / "body.stl")]
        r = runner.invoke(app, ["collar", *bodies, *args])
        assert r.exit_code == 1
        assert "explicit name" in r.output

    def test_open_body_loads_and_fails_only_when_ray_cast(self, tmp_path):
        import trimesh

        from shore.exceptions import MeshInputError
        from shore.mesh.body import Body

        post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=48)
        post.update_faces(post.triangles_center[:, 1] > 0.0)
        post.remove_unreferenced_vertices()
        post.export(str(tmp_path / "half.stl"))
        body = Body.from_stl(tmp_path / "half.stl")
        assert body.anchor_error is not None
        with pytest.raises(MeshInputError, match="anchor"):
            body.projector  # noqa: B018


class TestSplitCollar:
    """``shore collar --split`` (ADR 0024): the collar cut at a sharp trailing edge."""

    @staticmethod
    def bodies(tmp_path):
        import trimesh

        from shore.profiles import naca
        from shore.shapes import extrude_section

        wing = extrude_section(naca("0012", n_pts=40), span=1.0, n_stations=2, end_caps=True)
        wing.apply_translation((0.0, 0.0, -0.3))
        plate = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
        plate.apply_translation((0.5, 0.0, -0.5))
        wing.export(str(tmp_path / "wing.stl"))
        plate.export(str(tmp_path / "plate.stl"))
        return [str(tmp_path / "plate.stl"), str(tmp_path / "wing.stl")]

    ARGS = ("--seed", "0.5", "0.06", "0", "--n-loop", "64", "--ds-surface", "0.01",
            "--growth-surface", "1.15", "--m-a", "6", "--m-b", "6", "--ds-wall", "1e-4",
            "--growth-wall", "1.3", "--nk", "8")  # fmt: skip

    def test_split_writes_six_blocks_and_a_project(self, tmp_path):
        from shore.io.adjacency import read_adjacency_json
        from shore.project import Project

        out = tmp_path / "out" / "c"
        r = runner.invoke(
            app,
            ["collar", *self.bodies(tmp_path), *self.ARGS, "--split", "-o", str(out),
             "--adjacency-json", str(tmp_path / "c.adjacency.json"),
             "--project", str(tmp_path / "p")],
        )  # fmt: skip
        assert r.exit_code == 0, r.output
        assert "6 blocks" in r.output
        assert len(list((tmp_path / "out").glob("c_collar_0_*.geo"))) == 6
        payload = read_adjacency_json(tmp_path / "c.adjacency.json")
        assert len(payload["blocks"]) == 6
        p = Project.load(tmp_path / "p")
        assert len(p.mesh.blocks) == 6
        assert "split=True" in p.journal()

    def test_collar_volume_option(self, tmp_path):
        r = runner.invoke(
            app,
            ["collar", *self.bodies(tmp_path), *self.ARGS, "--split", "--collar-volume",
             "smoothed", "-o", str(tmp_path / "c")],
        )  # fmt: skip
        assert r.exit_code == 0, r.output
        r = runner.invoke(
            app,
            ["collar", *self.bodies(tmp_path), *self.ARGS, "--collar-volume", "fancy",
             "-o", str(tmp_path / "d")],
        )  # fmt: skip
        assert r.exit_code == 1
        assert "volume" in r.output

    def test_no_replay_adds_a_component(self, tmp_path):
        from shore.project import Project

        bodies = self.bodies(tmp_path)
        base = ["collar", *bodies, *self.ARGS, "--project", str(tmp_path / "p")]
        r = runner.invoke(app, base)
        assert r.exit_code == 0, r.output
        r = runner.invoke(app, [*base, "--split", "--component", "split", "--no-replay"])
        assert r.exit_code == 0, r.output
        p = Project.load(tmp_path / "p")  # replays both: the fast-loaded save is consistent
        assert len(p.components["collar"].blocks) == 1
        assert len(p.components["split"].blocks) == 6

    def test_wake_length_needs_split(self, tmp_path):
        r = runner.invoke(
            app,
            ["collar", *self.bodies(tmp_path), *self.ARGS, "--wake-length", "0.1",
             "-o", str(tmp_path / "c")],
        )  # fmt: skip
        assert r.exit_code == 1
        assert "split" in r.output

    def test_project_api(self, tmp_path):
        from shore.project import Project

        bodies = self.bodies(tmp_path)
        plain, split = Project(), Project()
        kw = {"seed": (0.5, 0.06, 0.0), "n_loop": 64, "ds_surface": 0.01,
              "growth_surface": 1.15, "m_a": 6, "m_b": 6, "ds_wall": 1e-4,
              "growth_wall": 1.3, "nk": 8}  # fmt: skip
        for proj in (plain, split):
            proj.load_geometry(bodies[0], name="plate")
        plain.build_collar(bodies[1], geometry="plate", **kw)
        split.build_collar(bodies[1], geometry="plate", split=True, **kw)
        # the sharp TE is cut: six blocks; the plain collar is one
        assert len(plain.mesh.blocks) == 1
        assert len(split.mesh.blocks) == 6


class TestBlendDeprecation:
    """blend_normals_k on CH (ADR 0019): refused like the C-grid it wraps."""

    def test_ch_blend_refused(self, naca0012_stl, tmp_path):
        args = [
            "mesh", str(naca0012_stl), "--topology", "ch", "--ni-body", "60", "--ni-wake", "15",
            "--n-stations", "4", "--wake-length", "10.0", "--nk", "8", "--ds", "1e-3",
            "--ch-chord-min", "-2.0", "--ch-chord-max", "15.0", "--ch-normal-min", "-2.5",
            "--ch-normal-max", "2.5", "--nk-channel", "4", "-o", str(tmp_path / "w"),
            "--blend-normals-k", "3",
        ]  # fmt: skip
        r = runner.invoke(app, args)
        assert r.exit_code == 1
        assert "C-grid" in r.output
