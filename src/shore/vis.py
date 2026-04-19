"""Visualisation helpers for SHORE grids.

Provides two backends:

- **pyvista** (default): opens an interactive 3D window.  Requires a display.
  Works on Linux/X11, Wayland (via XWayland), macOS, and Windows.

- **notebook** (``jupyter=True``): uses PyVista's HTML backend, which renders
  inline in a Jupyter notebook via a WebGL canvas — no display required.

Usage
-----
Interactive window::

    from shore.vis import view_geo
    from shore.io.geo import read_geo

    grid = read_geo("sphere.geo")
    view_geo(grid)

Jupyter notebook::

    from shore.vis import view_geo
    from shore.io.geo import read_geo

    grid = read_geo("sphere.geo")
    view_geo(grid, jupyter=True)

"""

from __future__ import annotations

import numpy as np

try:
    import pyvista as pv
except ImportError as e:
    raise ImportError(
        "pyvista is required for visualisation.  "
        "Install it with: pip install shore[vis]  or  pip install pyvista"
    ) from e

from shore.io.vtk import geo_to_pyvista, surface_to_pyvista
from shore.volume.jacobian import jacobian_field

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _make_plotter(jupyter: bool, off_screen: bool = False) -> pv.Plotter:
    """Create a plotter configured for the chosen backend."""
    if jupyter:
        pv.set_jupyter_backend("html")
    return pv.Plotter(off_screen=off_screen)


def _jacobian_clim(grid: np.ndarray) -> tuple[float, float]:
    """Colour limits for the Jacobian: 5th-95th percentile, floored at 0."""
    jac = jacobian_field(grid)
    lo = max(0.0, float(np.percentile(jac, 5)))
    hi = float(np.percentile(jac, 95))
    return lo, hi


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def view_geo(
    grid: np.ndarray,
    *,
    jupyter: bool = False,
    off_screen: bool = False,
    show_edges: bool = True,
    opacity: float = 0.6,
    colorby: str = "Jacobian",
    cmap: str = "RdYlGn",
    wall_color: str = "steelblue",
) -> pv.Plotter:
    """Visualise a structured volume grid.

    Renders the full O-grid block coloured by the Jacobian field.  The k=0
    wall layer is highlighted as a solid surface.

    Parameters
    ----------
    grid :
        Volume grid, shape (ni, nj, nk, 3).
    jupyter :
        If True, use the HTML backend for inline Jupyter display.
    off_screen :
        If True, render off-screen (for testing or image export).
    show_edges :
        Draw cell edges (default True).
    opacity :
        Grid block opacity [0-1] (default 0.6).
    colorby :
        Cell data array to colour by.  ``"Jacobian"`` (default) or ``"none"``.
    cmap :
        Matplotlib colourmap name (default ``"RdYlGn"`` — red=bad, green=good).
    wall_color :
        Colour for the k=0 wall surface (default ``"steelblue"``).

    Returns
    -------
    pyvista.Plotter
        Call ``.show()`` to open the window (already called when
        ``jupyter=False`` and ``off_screen=False``).
    """
    pl = _make_plotter(jupyter, off_screen)

    sg = geo_to_pyvista(grid, with_jacobian=True)

    clim = _jacobian_clim(grid) if colorby == "Jacobian" else None
    pl.add_mesh(
        sg,
        scalars=colorby if colorby != "none" else None,
        cmap=cmap,
        clim=clim,
        show_edges=show_edges,
        opacity=opacity,
        scalar_bar_args={"title": colorby, "vertical": True},
    )

    # Highlight wall surface (k=0 layer)
    wall = surface_to_pyvista(grid[:, :, 0, :])
    pl.add_mesh(wall, color=wall_color, opacity=1.0, show_edges=False)

    pl.add_axes()
    pl.add_bounding_box()

    if not jupyter and not off_screen:
        pl.show()

    return pl


def view_surface(
    k0: np.ndarray,
    *,
    jupyter: bool = False,
    off_screen: bool = False,
    show_edges: bool = True,
    cmap: str = "viridis",
    colorby: str = "radius",
) -> pv.Plotter:
    """Visualise a k=0 surface layer.

    Colours nodes by distance from the origin (``"radius"``) or by a flat
    colour when ``colorby="none"``.

    Parameters
    ----------
    k0 :
        Surface layer, shape (ni, nj, 3).
    jupyter :
        If True, use the HTML backend for inline Jupyter display.
    off_screen :
        Render off-screen (for testing).
    show_edges :
        Draw cell edges.
    cmap :
        Colourmap name.
    colorby :
        ``"radius"`` (default) or ``"none"``.

    Returns
    -------
    pyvista.Plotter
    """
    pl = _make_plotter(jupyter, off_screen)

    sg = surface_to_pyvista(k0)

    if colorby == "radius":
        pts = sg.points
        radii = np.linalg.norm(pts, axis=1)
        sg.point_data["radius"] = radii
        scalars = "radius"
    else:
        scalars = None

    pl.add_mesh(
        sg,
        scalars=scalars,
        cmap=cmap,
        show_edges=show_edges,
        scalar_bar_args={"title": colorby, "vertical": True} if scalars else {},
    )

    pl.add_axes()

    if not jupyter and not off_screen:
        pl.show()

    return pl


def save_screenshot(
    grid: np.ndarray,
    path: str,
    *,
    window_size: tuple[int, int] = (1200, 900),
) -> None:
    """Render a volume grid off-screen and save a PNG screenshot.

    Useful in CI or scripted workflows where no display is available.

    Parameters
    ----------
    grid :
        Volume grid, shape (ni, nj, nk, 3).
    path :
        Output PNG path.
    window_size :
        Render resolution in pixels.
    """
    pl = view_geo(grid, off_screen=True)
    pl.window_size = list(window_size)
    pl.screenshot(path)
    pl.close()
