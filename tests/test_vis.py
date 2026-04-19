"""Tests for shore.vis — off-screen rendering only (no display required)."""

from __future__ import annotations

import numpy as np
import pytest

pv = pytest.importorskip("pyvista", reason="pyvista not installed")

from shore.vis import save_screenshot, view_geo, view_surface  # noqa: E402
from tests.conftest import sphere_layer  # noqa: E402


def _make_volume(ni: int = 6, nj: int = 8, nk: int = 4) -> np.ndarray:
    """Synthetic O-grid volume: concentric sphere shells at geometric radii."""
    radii = 1.0 * 1.1 ** np.arange(nk)
    layers = [sphere_layer(ni, nj, radius=float(r)) for r in radii]
    return np.stack(layers, axis=2)


class TestViewGeo:
    def test_returns_plotter(self):
        grid = _make_volume()
        pl = view_geo(grid, off_screen=True)
        assert isinstance(pl, pv.Plotter)
        pl.close()

    def test_colorby_none(self):
        grid = _make_volume()
        pl = view_geo(grid, off_screen=True, colorby="none")
        pl.close()

    def test_no_edges(self):
        grid = _make_volume()
        pl = view_geo(grid, off_screen=True, show_edges=False)
        pl.close()


class TestViewSurface:
    def test_returns_plotter(self):
        k0 = sphere_layer(8, 10, radius=1.0)
        pl = view_surface(k0, off_screen=True)
        assert isinstance(pl, pv.Plotter)
        pl.close()

    def test_colorby_none(self):
        k0 = sphere_layer(8, 10, radius=1.0)
        pl = view_surface(k0, off_screen=True, colorby="none")
        pl.close()


class TestSaveScreenshot:
    def test_creates_png(self, tmp_path):
        grid = _make_volume()
        out = tmp_path / "shot.png"
        save_screenshot(grid, str(out))
        assert out.exists()
        # Check PNG magic bytes
        with out.open("rb") as f:
            assert f.read(8) == b"\x89PNG\r\n\x1a\n"
