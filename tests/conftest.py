"""Shared fixtures and helpers for the SHORE test suite."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import trimesh


def sphere_layer(ni: int, nj: int, radius: float = 1.0) -> np.ndarray:
    """Structured lat-lon sphere layer, shape (ni, nj, 3).

    Shared by test_hyperbolic.py and test_cli.py.
    """
    lats = np.linspace(0.05 * np.pi, 0.95 * np.pi, ni)
    lons = np.linspace(0.0, 2.0 * np.pi, nj, endpoint=False)
    la, lo = np.meshgrid(lats, lons, indexing="ij")
    x = radius * np.sin(la) * np.cos(lo)
    y = radius * np.sin(la) * np.sin(lo)
    z = radius * np.cos(la)
    return np.stack([x, y, z], axis=-1)


@pytest.fixture(scope="session")
def sphere_stl(tmp_path_factory) -> Path:
    """Export a unit icosphere as an STL file — created once per session."""
    tmp = tmp_path_factory.mktemp("stl")
    sphere = trimesh.creation.icosphere(subdivisions=3, radius=1.0)
    path = tmp / "sphere.stl"
    sphere.export(str(path))
    return path
