"""Round-trip tests for the geo format writer/reader."""

from __future__ import annotations

import numpy as np
import pytest

from shore.io.geo import read_geo, write_geo


def _make_grid(ni: int, nj: int, nk: int) -> np.ndarray:
    """Simple structured grid: x=i, y=j, z=k."""
    i = np.arange(ni, dtype=float)
    j = np.arange(nj, dtype=float)
    k = np.arange(nk, dtype=float)
    ig, jg, kg = np.meshgrid(i, j, k, indexing="ij")
    return np.stack([ig, jg, kg], axis=-1)


def test_round_trip(tmp_path):
    grid = _make_grid(4, 5, 6)
    p = tmp_path / "test.geo"
    write_geo(p, grid)
    loaded = read_geo(p)
    np.testing.assert_allclose(loaded, grid)


def test_header(tmp_path):
    grid = _make_grid(3, 7, 2)
    p = tmp_path / "test.geo"
    write_geo(p, grid)
    first_line = p.read_text().splitlines()[0]
    assert first_line.strip() == "3 7 2"


def test_wrong_size_raises(tmp_path):
    p = tmp_path / "bad.geo"
    p.write_text("3 3 3\n" + " ".join(["1.0"] * 10) + "\n")
    with pytest.raises(ValueError, match="Expected"):
        read_geo(p)
