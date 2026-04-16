"""Shared fixtures and helpers for the SHORE test suite."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import trimesh


def sphere_layer(
    ni: int, nj: int, radius: float = 1.0, theta_margin: float | None = None
) -> np.ndarray:
    """Structured lat-lon sphere layer, shape (ni, nj, 3).

    Shared by test_hyperbolic.py and test_cli.py.

    Parameters
    ----------
    ni :
        Number of latitude points.
    nj :
        Number of longitude points.
    radius :
        Sphere radius.
    theta_margin :
        Polar-cap exclusion angle. If None, uses legacy formula.
    """
    if theta_margin is None:
        lats = np.linspace(0.05 * np.pi, 0.95 * np.pi, ni)
    else:
        lats = np.linspace(theta_margin, np.pi - theta_margin, ni)
    lons = np.linspace(0.0, 2.0 * np.pi, nj, endpoint=False)
    la, lo = np.meshgrid(lats, lons, indexing="ij")
    x = radius * np.sin(la) * np.cos(lo)
    y = radius * np.sin(la) * np.sin(lo)
    z = radius * np.cos(la)
    return np.stack([x, y, z], axis=-1)


def yang_layer(
    ni: int, nj: int, radius: float = 1.0, theta_margin: float = np.pi / 6
) -> np.ndarray:
    """Structured Yang O-grid layer, shape (ni, nj, 3).

    Rotated version of the sphere layer by 90 degrees around the y-axis.

    Parameters
    ----------
    ni :
        Number of latitude points.
    nj :
        Number of longitude points.
    radius :
        Sphere radius.
    theta_margin :
        Polar-cap exclusion angle in radians (default pi/6).
    """
    from shore.surface.yinyang import YANG_ROTATION

    if theta_margin is None:
        lats = np.linspace(0.05 * np.pi, 0.95 * np.pi, ni)
    else:
        lats = np.linspace(theta_margin, np.pi - theta_margin, ni)
    lons = np.linspace(0.0, 2.0 * np.pi, nj, endpoint=False)
    la, lo = np.meshgrid(lats, lons, indexing="ij")
    d = np.stack(
        [
            np.sin(la) * np.cos(lo),
            np.sin(la) * np.sin(lo),
            np.cos(la),
        ],
        axis=-1,
    )  # (ni, nj, 3)
    flat = d.reshape(-1, 3)
    rotated = (YANG_ROTATION @ flat.T).T
    return radius * rotated.reshape(ni, nj, 3)


@pytest.fixture(scope="session")
def sphere_stl(tmp_path_factory) -> Path:
    """Export a unit icosphere as an STL file — created once per session."""
    tmp = tmp_path_factory.mktemp("stl")
    sphere = trimesh.creation.icosphere(subdivisions=3, radius=1.0)
    path = tmp / "sphere.stl"
    sphere.export(str(path))
    return path
