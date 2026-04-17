"""Shared fixtures and helpers for the SHORE test suite."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import trimesh


def sphere_layer(
    ni: int,
    nj: int,
    radius: float = 1.0,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    theta_margin: float | None = None,
) -> np.ndarray:
    """Structured lat-lon sphere layer, shape (ni, nj, 3).

    Shared by test_hyperbolic.py, test_sphere_cap.py, and test_cli.py.

    Parameters
    ----------
    ni :
        Number of latitude points.
    nj :
        Number of longitude points.
    radius :
        Sphere radius.
    theta_cap :
        Polar exclusion angle in radians: the band covers
        ``[theta_cap, pi - theta_cap]``.  Default ``0.05*pi`` (~9°).
    theta_cap_deg :
        Convenience alias for ``theta_cap`` in degrees.
    theta_margin :
        Deprecated alias for ``theta_cap`` (radians).
    """
    n_supplied = sum(x is not None for x in (theta_cap, theta_cap_deg, theta_margin))
    if n_supplied > 1:
        raise ValueError("Supply at most one of theta_cap, theta_cap_deg, theta_margin.")
    if theta_cap_deg is not None:
        theta_cap = np.deg2rad(theta_cap_deg)
    elif theta_margin is not None:
        theta_cap = theta_margin
    if theta_cap is None:
        theta_cap = 0.05 * np.pi

    if not (0.0 < theta_cap < np.pi / 2):
        raise ValueError(
            f"theta_cap={theta_cap:.6f} must be in (0, pi/2). "
            "0 would collapse the cap to a point; pi/2 would put the "
            "cap edge at the equator."
        )

    lats = np.linspace(theta_cap, np.pi - theta_cap, ni)
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
