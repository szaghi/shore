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


# ---------------------------------------------------------------------------
# NACA 0012 STL fixture (for the C-grid topology tests)
# ---------------------------------------------------------------------------


def _naca0012_section(n_pts: int, blunt_te_thickness: float = 0.0) -> np.ndarray:
    """NACA 0012 cross-section as a ``(n_pts - 1, 2)`` Selig-ordered polyline.

    Thin wrapper over :func:`shore.profiles.naca` for the test suite.
    The legacy fixture used ``n_pts`` to mean "node count along the
    closed curve", with one shared LE node — i.e. ``n_pts // 2`` chord
    samples per side and ``n_pts - 1`` total points after the LE
    merge.  ``shore.profiles.naca`` takes ``n_pts`` as chord samples
    per side instead, so we pass ``n_pts // 2`` to keep the legacy
    point-count contract.

    For blunt-TE sections (``blunt_te_thickness > 0``) the upper- and
    lower-TE vertices are shifted by ``blunt_te_thickness / 2`` in
    opposite y directions.  ``shore.profiles`` doesn't expose a
    blunt knob; this small post-process keeps the test invariant
    that ``te_thickness == blunt_te_thickness``.
    """
    if n_pts < 8 or n_pts % 2 != 0:
        raise ValueError(f"_naca0012_section: n_pts must be even >= 8, got {n_pts}")
    from shore.profiles import naca

    polyline = naca("0012", n_pts=n_pts // 2, closed_te=True)
    if blunt_te_thickness > 0.0:
        polyline = polyline.copy()
        polyline[0, 1] += 0.5 * blunt_te_thickness  # upper TE
        polyline[-1, 1] -= 0.5 * blunt_te_thickness  # lower TE
    return polyline


def _make_naca0012_stl(
    path: Path,
    *,
    n_section_pts: int = 80,
    n_stations: int = 8,
    chord: float = 1.0,
    span: float = 0.5,
    span_axis: str = "z",
    blunt_te_thickness: float = 0.0,
) -> Path:
    """Build a NACA 0012 STL via :func:`shore.shapes.extrude_section`.

    Thin wrapper for the test suite.  Returns a watertight wing
    with centroid-fan end caps; for blunt TE the upper / lower TE
    vertices are pre-shifted in 2D so that ``extrude_section``'s
    side-surface stripping naturally bridges the TE gap.
    """
    from shore.shapes import extrude_section

    section = _naca0012_section(n_section_pts, blunt_te_thickness=blunt_te_thickness) * chord
    mesh = extrude_section(
        section,
        span=span,
        n_stations=n_stations,
        span_axis=span_axis,
        end_caps=True,
    )
    mesh.export(str(path))
    return path


@pytest.fixture(scope="session")
def naca0012_stl(tmp_path_factory) -> Path:
    """NACA 0012 sharp-TE STL — created once per session."""
    tmp = tmp_path_factory.mktemp("naca")
    return _make_naca0012_stl(tmp / "naca0012.stl")


@pytest.fixture(scope="session")
def naca0012_blunt_stl(tmp_path_factory) -> Path:
    """NACA 0012 blunt-TE (1% chord thickness) STL — created once per session."""
    tmp = tmp_path_factory.mktemp("naca_blunt")
    return _make_naca0012_stl(tmp / "naca0012_blunt.stl", blunt_te_thickness=0.01)
