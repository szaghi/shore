---
title: shore.vis
outline: [2, 3]
---

# `shore.vis`

Visualisation helpers for SHORE grids.

Provides two backends:

- **pyvista** (default): opens an interactive 3D window.  Requires a display.
  Works on Linux/X11, Wayland (via XWayland), macOS, and Windows.

- **notebook** (`jupyter=True`): uses PyVista's HTML backend, which renders
  inline in a Jupyter notebook via a WebGL canvas — no display required.

Usage
-----
Interactive window:

```
from shore.vis import view_geo
from shore.io.geo import read_geo

grid = read_geo("sphere.geo")
view_geo(grid)

```
Jupyter notebook:

```
from shore.vis import view_geo
from shore.io.geo import read_geo

grid = read_geo("sphere.geo")
view_geo(grid, jupyter=True)
```

## `view_geo` <Badge type="info" text="function" /> {#view_geo}

```python
view_geo(
    grid: np.ndarray,
    *,
    jupyter: bool = False,
    off_screen: bool = False,
    show_edges: bool = True,
    opacity: float = 0.6,
    colorby: str = 'Jacobian',
    cmap: str = 'RdYlGn',
    wall_color: str = 'steelblue',
) -> pv.Plotter
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/vis.py#L71">source</a>

Visualise a structured volume grid.

Renders the full O-grid block coloured by the Jacobian field.  The k=0
wall layer is highlighted as a solid surface.

**Parameters**

- **`grid`** — Volume grid, shape (ni, nj, nk, 3).
- **`jupyter`** — If True, use the HTML backend for inline Jupyter display.
- **`off_screen`** — If True, render off-screen (for testing or image export).
- **`show_edges`** — Draw cell edges (default True).
- **`opacity`** — Grid block opacity [0-1] (default 0.6).
- **`colorby`** — Cell data array to colour by.  `"Jacobian"` (default) or `"none"`.
- **`cmap`** — Matplotlib colourmap name (default `"RdYlGn"` — red=bad, green=good).
- **`wall_color`** — Colour for the k=0 wall surface (default `"steelblue"`).

**Returns**

- `pyvista.Plotter` — Call `.show()` to open the window (already called when `jupyter=False` and `off_screen=False`).

## `view_surface` <Badge type="info" text="function" /> {#view_surface}

```python
view_surface(
    k0: np.ndarray,
    *,
    jupyter: bool = False,
    off_screen: bool = False,
    show_edges: bool = True,
    cmap: str = 'viridis',
    colorby: str = 'radius',
) -> pv.Plotter
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/vis.py#L140">source</a>

Visualise a k=0 surface layer.

Colours nodes by distance from the origin (`"radius"`) or by a flat
colour when `colorby="none"`.

**Parameters**

- **`k0`** — Surface layer, shape (ni, nj, 3).
- **`jupyter`** — If True, use the HTML backend for inline Jupyter display.
- **`off_screen`** — Render off-screen (for testing).
- **`show_edges`** — Draw cell edges.
- **`cmap`** — Colourmap name.
- **`colorby`** — `"radius"` (default) or `"none"`.

**Returns**

- `pyvista.Plotter`

## `save_screenshot` <Badge type="info" text="function" /> {#save_screenshot}

```python
save_screenshot(
    grid: np.ndarray,
    path: str,
    *,
    window_size: tuple[int, int] = (1200, 900),
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/vis.py#L201">source</a>

Render a volume grid off-screen and save a PNG screenshot.

Useful in CI or scripted workflows where no display is available.

**Parameters**

- **`grid`** — Volume grid, shape (ni, nj, nk, 3).
- **`path`** — Output PNG path.
- **`window_size`** — Render resolution in pixels.
