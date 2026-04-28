"""Tests for the Xall binary .grd writer (shore.io.grd) and shore grd CLI."""

from __future__ import annotations

import struct

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.cli import app
from shore.io.geo import write_geo
from shore.io.grd import write_grd

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_grid(ni: int, nj: int, nk: int) -> np.ndarray:
    """Simple (ni, nj, nk, 3) grid: coordinates equal index values."""
    i = np.arange(ni, dtype=float)
    j = np.arange(nj, dtype=float)
    k = np.arange(nk, dtype=float)
    ig, jg, kg = np.meshgrid(i, j, k, indexing="ij")
    return np.stack([ig, jg, kg], axis=-1)


def _read_fortran_record(f) -> bytes:
    """Read one Fortran unformatted record (4-byte length markers)."""
    (nbytes,) = struct.unpack("=i", f.read(4))
    data = f.read(nbytes)
    (nbytes2,) = struct.unpack("=i", f.read(4))
    assert nbytes == nbytes2, "mismatched Fortran record markers"
    return data


# ---------------------------------------------------------------------------
# Unit tests — write_grd
# ---------------------------------------------------------------------------


class TestWriteGrd:
    def test_single_block_file_structure(self, tmp_path):
        grid = _make_grid(4, 5, 6)
        geo = tmp_path / "blk.geo"
        write_geo(geo, grid)
        out = write_grd(tmp_path / "out.grd", [geo])

        assert out.exists()
        assert out.suffix == ".grd"

        with out.open("rb") as f:
            # Record 1: nb = 1
            rec = _read_fortran_record(f)
            (nb,) = struct.unpack("=i", rec)
            assert nb == 1

            # Record 2: cell counts (node counts - 1)
            rec = _read_fortran_record(f)
            ni_c, nj_c, nk_c = struct.unpack("=3i", rec)
            assert (ni_c, nj_c, nk_c) == (3, 4, 5)  # 4-1, 5-1, 6-1

            # Three coordinate records
            for _axis in range(3):
                rec = _read_fortran_record(f)
                coords = np.frombuffer(rec, dtype=np.float64)
                assert coords.shape == (4 * 5 * 6,)

    def test_coordinates_fortran_order(self, tmp_path):
        grid = _make_grid(3, 4, 2)
        geo = tmp_path / "blk.geo"
        write_geo(geo, grid)
        out = write_grd(tmp_path / "out.grd", [geo])

        with out.open("rb") as f:
            _read_fortran_record(f)  # nb
            _read_fortran_record(f)  # dims
            x_rec = _read_fortran_record(f)

        x_flat = np.frombuffer(x_rec, dtype=np.float64)
        expected = grid[:, :, :, 0].ravel(order="F")
        np.testing.assert_allclose(x_flat, expected)

    def test_multi_block_block_count(self, tmp_path):
        grids = [_make_grid(3, 4, 5), _make_grid(2, 3, 4)]
        geo_paths = []
        for i, g in enumerate(grids):
            p = tmp_path / f"blk{i}.geo"
            write_geo(p, g)
            geo_paths.append(p)

        out = write_grd(tmp_path / "multi.grd", geo_paths)
        with out.open("rb") as f:
            rec = _read_fortran_record(f)
            (nb,) = struct.unpack("=i", rec)
        assert nb == 2

    def test_multi_block_cell_counts(self, tmp_path):
        shapes = [(3, 4, 5), (6, 7, 8)]
        geo_paths = []
        for i, (ni, nj, nk) in enumerate(shapes):
            p = tmp_path / f"blk{i}.geo"
            write_geo(p, _make_grid(ni, nj, nk))
            geo_paths.append(p)

        out = write_grd(tmp_path / "multi.grd", geo_paths)
        with out.open("rb") as f:
            _read_fortran_record(f)  # nb
            for ni, nj, nk in shapes:
                rec = _read_fortran_record(f)
                ci, cj, ck = struct.unpack("=3i", rec)
                assert (ci, cj, ck) == (ni - 1, nj - 1, nk - 1)

    def test_grd_extension_added_if_absent(self, tmp_path):
        geo = tmp_path / "blk.geo"
        write_geo(geo, _make_grid(2, 2, 2))
        out = write_grd(tmp_path / "no_ext", [geo])
        assert out.suffix == ".grd"
        assert out.exists()

    def test_empty_list_raises(self, tmp_path):
        with pytest.raises(ValueError, match="empty"):
            write_grd(tmp_path / "out.grd", [])


# ---------------------------------------------------------------------------
# CLI tests — shore grd
# ---------------------------------------------------------------------------


class TestCLIGrd:
    def test_basic_round_trip(self, tmp_path):
        geo = tmp_path / "blk.geo"
        write_geo(geo, _make_grid(4, 5, 6))
        runner = CliRunner()
        r = runner.invoke(app, ["grd", str(geo), "--output", str(tmp_path / "out.grd")])
        assert r.exit_code == 0, r.output
        assert (tmp_path / "out.grd").exists()

    def test_multi_block(self, tmp_path):
        for i, (ni, nj, nk) in enumerate([(3, 4, 5), (2, 3, 4)]):
            write_geo(tmp_path / f"blk{i}.geo", _make_grid(ni, nj, nk))
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "grd",
                str(tmp_path / "blk0.geo"),
                str(tmp_path / "blk1.geo"),
                "--output",
                str(tmp_path / "multi.grd"),
            ],
        )
        assert r.exit_code == 0, r.output
        out = tmp_path / "multi.grd"
        assert out.exists()
        with out.open("rb") as f:
            rec = _read_fortran_record(f)
            (nb,) = struct.unpack("=i", rec)
        assert nb == 2

    def test_missing_input_exits_1(self, tmp_path):
        runner = CliRunner()
        r = runner.invoke(
            app,
            ["grd", str(tmp_path / "nope.geo"), "--output", str(tmp_path / "out.grd")],
        )
        assert r.exit_code == 1
        assert "not found" in r.output.lower()

    def test_non_geo_input_exits_1(self, tmp_path):
        bad = tmp_path / "file.txt"
        bad.write_text("hello")
        runner = CliRunner()
        r = runner.invoke(app, ["grd", str(bad), "--output", str(tmp_path / "out.grd")])
        assert r.exit_code == 1
        assert ".geo" in r.output
