"""Shared fixtures and helpers for the SHORE test suite."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest
import trimesh

# Headless CI: force PyVista to render off-screen.  Without this, any test
# that calls ``Plotter.render()`` segfaults on GitHub Actions runners (no
# X server, no software OpenGL fallback).  The ``pyvista`` import is in a
# try/except because the dependency is optional (``shore[vis]``).
if os.environ.get("CI") == "true":
    try:
        import pyvista

        pyvista.OFF_SCREEN = True
    except ImportError:
        pass

_PUBLIC = Path(__file__).parent.parent / "docs" / "public"


def pytest_sessionfinish(session, exitstatus):
    """After the test run: write coverage.json and coverage-badge.json to docs/public/.

    The VitePress site picks up ``coverage-badge.json`` for a shields.io
    endpoint badge.  The CI ``docs.yml`` workflow uploads ``coverage*.json``
    as an artifact and re-injects it into ``docs/public/`` before building
    the site, so the badge is always in sync with the suite that just
    passed.
    """
    _PUBLIC.mkdir(parents=True, exist_ok=True)
    try:
        import coverage as coverage_lib

        cov = coverage_lib.Coverage()
        cov.load()
        cov.json_report(outfile=str(_PUBLIC / "coverage.json"), pretty_print=True)
        data = json.loads((_PUBLIC / "coverage.json").read_text())
        pct = float(data["totals"]["percent_covered_display"])
    except Exception:
        return

    if pct >= 90:
        color = "brightgreen"
    elif pct >= 75:
        color = "green"
    elif pct >= 60:
        color = "yellow"
    elif pct >= 40:
        color = "orange"
    else:
        color = "red"

    badge = {
        "schemaVersion": 1,
        "label": "coverage",
        "message": f"{pct:.0f}%",
        "color": color,
    }
    (_PUBLIC / "coverage-badge.json").write_text(json.dumps(badge, indent=2))


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


@pytest.fixture(scope="session")
def sphere_stl(tmp_path_factory) -> Path:
    """Export a unit icosphere as an STL file — created once per session."""
    tmp = tmp_path_factory.mktemp("stl")
    sphere = trimesh.creation.icosphere(subdivisions=3, radius=1.0)
    path = tmp / "sphere.stl"
    sphere.export(str(path))
    return path
