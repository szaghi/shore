# shore.vis

Interactive and off-screen visualization using PyVista.

::: warning Optional dependency
This module requires **pyvista** (`pip install shore[vis]`). Import will fail with `ImportError` if pyvista is not installed.
:::

## Functions

### `view_geo`

```python
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
) -> pv.Plotter
```

Render a SHORE volume grid in an interactive 3-D viewer.

The k=0 wall layer is highlighted in `wall_color`. The outer layers are coloured by the chosen scalar field with the colour scale clamped to the 5th–95th percentile to avoid outliers washing out the range.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `grid` | `np.ndarray (ni, nj, nk, 3)` | — | Full O-grid volume |
| `jupyter` | `bool` | `False` | Set the HTML backend and return without calling `pl.show()` |
| `off_screen` | `bool` | `False` | Render off-screen (for testing and screenshots) |
| `show_edges` | `bool` | `True` | Draw cell edges |
| `opacity` | `float` | 0.6 | Mesh opacity (0 = transparent, 1 = opaque) |
| `colorby` | `str` | `"Jacobian"` | Scalar field to colour by; `"none"` for a solid colour |
| `cmap` | `str` | `"RdYlGn"` | Matplotlib colormap name |
| `wall_color` | `str` | `"steelblue"` | Colour for the k=0 wall layer |

**Returns**

`pyvista.Plotter` — the plotter instance. Call `pl.show()` if needed; the function calls it automatically unless `jupyter=True` or `off_screen=True`.

**Example — desktop window**

```python
from shore.io.geo import read_geo
from shore.vis import view_geo

grid = read_geo("sphere.geo")
view_geo(grid)   # opens interactive window
```

**Example — Jupyter notebook**

```python
import pyvista as pv
pv.set_jupyter_backend("html")

from shore.vis import view_geo
from shore.io.geo import read_geo

grid = read_geo("sphere.geo")
pl = view_geo(grid, jupyter=True)
pl.show()
```

---

### `view_surface`

```python
def view_surface(
    k0: np.ndarray,
    *,
    jupyter: bool = False,
    off_screen: bool = False,
    show_edges: bool = True,
    cmap: str = "viridis",
    colorby: str = "radius",
) -> pv.Plotter
```

Render a structured surface layer (k=0) in a 3-D viewer, coloured by the distance from the origin.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `k0` | `np.ndarray (ni, nj, 3)` | — | Structured surface layer |
| `jupyter` | `bool` | `False` | HTML backend; return without `pl.show()` |
| `off_screen` | `bool` | `False` | Off-screen rendering |
| `show_edges` | `bool` | `True` | Draw cell edges |
| `cmap` | `str` | `"viridis"` | Matplotlib colormap name |
| `colorby` | `str` | `"radius"` | Scalar field; `"none"` for solid colour |

**Returns**

`pyvista.Plotter`

**Example**

```python
import numpy as np
from shore.vis import view_surface

k0 = np.load("sphere_k0.npz")["k0"]
view_surface(k0)
```

---

### `save_screenshot`

```python
def save_screenshot(
    grid: np.ndarray,
    path: str,
    *,
    window_size: tuple[int, int] = (1200, 900),
) -> None
```

Render a volume grid off-screen and save a PNG screenshot.

Does not open a window. Suitable for CI pipelines and automated report generation.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `grid` | `np.ndarray (ni, nj, nk, 3)` | — | Full O-grid volume |
| `path` | `str` | — | Output PNG path |
| `window_size` | `tuple[int, int]` | `(1200, 900)` | Pixel dimensions of the rendered image |

**Example**

```python
from shore.io.geo import read_geo
from shore.vis import save_screenshot

grid = read_geo("sphere.geo")
save_screenshot(grid, "sphere_preview.png")
```

## Jupyter workflow

The recommended Jupyter workflow uses the `html` backend, which renders inline in the notebook without requiring a display:

```python
import pyvista as pv
pv.set_jupyter_backend("html")

from shore.vis import view_geo, view_surface
from shore.io.geo import read_geo
import numpy as np

# Volume
grid = read_geo("sphere.geo")
pl = view_geo(grid, jupyter=True)
pl.show()

# Surface layer
k0 = np.load("sphere_k0.npz")["k0"]
pl = view_surface(k0, jupyter=True)
pl.show()
```

See `examples/visualize.ipynb` for a full notebook with Jacobian analysis plots.

## See also

- [`shore view`](view) — CLI wrapper for `view_geo` / `view_surface`
- [`shore.io.vtk`](api-vtk) — write VTK files for ParaView
- [`shore.volume.jacobian`](api-jacobian) — access the raw Jacobian array
