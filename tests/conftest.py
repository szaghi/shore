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
    """NACA 0012 cross-section as a closed (n_pts, 2) polyline.

    Parameters
    ----------
    n_pts :
        Total node count along the closed curve.  Must be even; half on
        suction side, half on pressure side.  The two TE endpoints are
        coincident for a sharp TE (``blunt_te_thickness == 0``) or
        offset by ``blunt_te_thickness * chord`` perpendicular to the
        chord for a blunt TE.
    blunt_te_thickness :
        TE gap as a fraction of chord (= 1.0).  Default ``0.0`` (sharp).

    Returns
    -------
    numpy.ndarray, shape ``(n_pts, 2)``
        Closed polyline starting at the upper TE, going forward over
        the suction surface to the LE, back over the pressure surface
        to the lower TE.  ``polyline[0]`` is the upper TE point;
        ``polyline[-1]`` is the lower TE point.
    """
    if n_pts < 8 or n_pts % 2 != 0:
        raise ValueError(f"_naca0012_section: n_pts must be even >= 8, got {n_pts}")
    n_half = n_pts // 2
    # Cosine-spaced x along chord — denser at LE and TE, the natural NACA
    # parameterisation.
    beta = np.linspace(0.0, np.pi, n_half)
    x = 0.5 * (1.0 - np.cos(beta))  # 0 → 1
    # NACA 4-digit thickness equation (12% chord), with the TE-closing
    # ``-0.1015 * x^4`` term (sharp TE convention).
    t = 0.12
    y_t = (
        5.0 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    )
    # Force the analytic TE ordinate to exactly 0 (the polynomial
    # leaves a ~0.1% residual; we want a clean closed curve).
    y_t[-1] = 0.5 * blunt_te_thickness
    # Suction (upper) surface: TE → LE
    suction = np.column_stack([x[::-1], y_t[::-1]])
    # Pressure (lower) surface: LE → TE
    pressure = np.column_stack([x[1:], -y_t[1:]])
    polyline = np.vstack([suction, pressure])
    assert polyline.shape == (n_pts - 1, 2), polyline.shape
    # We've produced n_pts - 1 points (LE shared between branches).
    # For sharp TE the polyline is open; we don't close it because the
    # consumer code will treat polyline[0] / polyline[-1] as the upper /
    # lower TE points and decide sharp vs blunt.
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
    """Build a NACA 0012 STL by extruding the analytic section along *span_axis*.

    The result is a watertight triangulation: two surface meshes
    (upper + lower) joined along the LE and TE feature lines, capped
    at both span ends.  Suitable input for both the C-grid TE
    detector (sharp dihedral at the TE) and the cubed-sphere O-grid
    path (though it would fail star-shape there — that's fine, this
    fixture is for C-grid tests).
    """
    section = _naca0012_section(n_section_pts, blunt_te_thickness=blunt_te_thickness) * chord
    # Extrude along span_axis at n_stations stations.
    ax = {"x": 0, "y": 1, "z": 2}[span_axis]
    other_axes = [a for a in (0, 1, 2) if a != ax]
    span_coords = np.linspace(0.0, span, n_stations)

    # vertices[s, i] = (3,) world point at station s, section node i.
    n_sect = section.shape[0]
    vertices = np.zeros((n_stations, n_sect, 3))
    for s, sc in enumerate(span_coords):
        vertices[s, :, other_axes[0]] = section[:, 0]
        vertices[s, :, other_axes[1]] = section[:, 1]
        vertices[s, :, ax] = sc

    # Triangulate the side surface (between consecutive stations + section nodes).
    # The section is open (not closed), so we connect it as a strip; the
    # last → first edge is the TE, which we leave as a sharp feature.
    tris: list[tuple[int, int, int]] = []

    def vid(s: int, i: int) -> int:
        return s * n_sect + i

    for s in range(n_stations - 1):
        for i in range(n_sect - 1):
            v00 = vid(s, i)
            v10 = vid(s + 1, i)
            v01 = vid(s, i + 1)
            v11 = vid(s + 1, i + 1)
            tris.append((v00, v10, v11))
            tris.append((v00, v11, v01))

    # End-cap triangulation: fan from the section centroid at each span end.
    # NACA 0012 is convex (good enough — only the wake-side of the chord
    # is missing for a sharp section; a blunt section is a closed polygon).
    # For sharp TE, the section's last and first nodes coincide (TE),
    # so we fan from the chord-midpoint centroid covering all interior
    # edges.
    section_2d = section  # (n_sect, 2)
    centroid_2d = section_2d.mean(axis=0)
    # Add the two cap-centroid vertices at the start- and end-span planes.
    n_total_pts = n_stations * n_sect
    cap_lo_3d = np.zeros(3)
    cap_lo_3d[other_axes[0]] = centroid_2d[0]
    cap_lo_3d[other_axes[1]] = centroid_2d[1]
    cap_lo_3d[ax] = span_coords[0]
    cap_hi_3d = cap_lo_3d.copy()
    cap_hi_3d[ax] = span_coords[-1]
    centroid_lo_vid = n_total_pts
    centroid_hi_vid = n_total_pts + 1
    all_vertices = np.vstack([vertices.reshape(-1, 3), cap_lo_3d, cap_hi_3d])

    # Cap at low-span end: fan from centroid around the section.
    for i in range(n_sect - 1):
        tris.append((centroid_lo_vid, vid(0, i + 1), vid(0, i)))
    # Closing edge of the polyline (lower-TE → upper-TE) only exists
    # for blunt sections; for sharp TE the two endpoints coincide
    # after merge_vertices(), so the fan triangle is degenerate
    # and gets dropped harmlessly.
    if blunt_te_thickness > 0.0:
        tris.append((centroid_lo_vid, vid(0, 0), vid(0, n_sect - 1)))

    # Cap at high-span end: same, opposite winding.
    for i in range(n_sect - 1):
        tris.append((centroid_hi_vid, vid(n_stations - 1, i), vid(n_stations - 1, i + 1)))
    if blunt_te_thickness > 0.0:
        tris.append((centroid_hi_vid, vid(n_stations - 1, n_sect - 1), vid(n_stations - 1, 0)))

    # Side-surface bridging triangles for the blunt-TE wake gap (one
    # quad per span step, between (s, 0) — upper TE — and (s, n_sect - 1)
    # — lower TE).  Sharp TE: those two vertices coincide after
    # merge_vertices(), so the bridging triangles collapse to degenerates.
    if blunt_te_thickness > 0.0:
        for s in range(n_stations - 1):
            v00 = vid(s, 0)
            v01 = vid(s, n_sect - 1)
            v10 = vid(s + 1, 0)
            v11 = vid(s + 1, n_sect - 1)
            tris.append((v00, v01, v11))
            tris.append((v00, v11, v10))

    mesh = trimesh.Trimesh(vertices=all_vertices, faces=np.asarray(tris, dtype=np.int64))
    # Merge coincident vertices so the TE feature line and the LE
    # cap-fan rays share their endpoints with the side surface.
    # ``process()`` without ``validate`` skips winding-fix (which needs
    # networkx); the analytic construction above already produces a
    # consistent winding.
    mesh.merge_vertices()
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
