"""Tests for VTK export (shore.io.vtk)."""

from __future__ import annotations

import numpy as np
import pytest

pv = pytest.importorskip("pyvista", reason="pyvista not installed")

from shore.io.vtk import geo_to_pyvista, surface_to_pyvista, write_surface_vts, write_vts
from tests.conftest import sphere_layer


def _make_volume(ni: int = 6, nj: int = 8, nk: int = 4) -> np.ndarray:
    from shore.volume.hyperbolic import march
    surface = sphere_layer(ni, nj, radius=1.0)
    return march(surface, nk=nk, ds=0.05, growth=1.1, smooth_strength=0.0)


class TestGeotoPyvista:
    def test_dimensions(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid)
        ni, nj, nk, _ = grid.shape
        assert sg.dimensions == (ni, nj + 1, nk)   # j is closed

    def test_point_count(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid)
        ni, nj, nk, _ = grid.shape
        assert sg.n_points == ni * (nj + 1) * nk

    def test_jacobian_cell_data_present(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid, with_jacobian=True)
        assert "Jacobian" in sg.cell_data

    def test_jacobian_cell_count(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid, with_jacobian=True)
        ni, nj, nk, _ = grid.shape
        expected = (ni - 1) * nj * (nk - 1)
        assert sg.cell_data["Jacobian"].shape[0] == expected

    def test_jacobian_all_positive(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid, with_jacobian=True)
        assert sg.cell_data["Jacobian"].min() > 0.0

    def test_no_jacobian(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid, with_jacobian=False)
        assert "Jacobian" not in sg.cell_data


class TestSurfaceToPyvista:
    def test_dimensions(self):
        k0 = sphere_layer(8, 10, radius=1.0)
        sg = surface_to_pyvista(k0)
        ni, nj, _ = k0.shape
        assert sg.dimensions == (ni, nj + 1, 1)

    def test_point_count(self):
        k0 = sphere_layer(8, 10, radius=1.0)
        sg = surface_to_pyvista(k0)
        ni, nj, _ = k0.shape
        assert sg.n_points == ni * (nj + 1) * 1


class TestWriteVts:
    def test_creates_file(self, tmp_path):
        grid = _make_volume(6, 8, 4)
        out = tmp_path / "test.vts"
        write_vts(out, grid)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_appends_extension(self, tmp_path):
        grid = _make_volume(6, 8, 4)
        out = tmp_path / "test"
        write_vts(out, grid)
        assert (tmp_path / "test.vts").exists()

    def test_readable_by_pyvista(self, tmp_path):
        grid = _make_volume(6, 8, 4)
        out = tmp_path / "test.vts"
        write_vts(out, grid)
        loaded = pv.read(str(out))
        assert loaded.n_points > 0
        assert "Jacobian" in loaded.cell_data


class TestWriteSurfaceVts:
    def test_creates_file(self, tmp_path):
        k0 = sphere_layer(8, 10, radius=1.0)
        out = tmp_path / "surf.vts"
        write_surface_vts(out, k0)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_readable_by_pyvista(self, tmp_path):
        k0 = sphere_layer(8, 10, radius=1.0)
        out = tmp_path / "surf.vts"
        write_surface_vts(out, k0)
        loaded = pv.read(str(out))
        assert loaded.n_points > 0
