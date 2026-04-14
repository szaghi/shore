# shore.io.vtk

VTK StructuredGrid export for ParaView, VisIt, and PyVista.

::: warning Optional dependency
This module requires **pyvista** (`pip install shore[vis]`). Import will fail with `ImportError` if pyvista is not installed.
:::

## O-grid periodicity and VTK

VTK StructuredGrid has no concept of periodic connectivity. To represent the closed O-grid in VTK, the j-direction is **closed**: the first column (`j=0`) is appended as an extra `j=nj` column, giving dimensions `(ni, nj+1, nk)`. All exported grids follow this convention.

## Functions

### `geo_to_pyvista`

```python
def geo_to_pyvista(
    grid: np.ndarray,
    with_jacobian: bool = True,
) -> pv.StructuredGrid
```

Convert a SHORE volume grid to a PyVista `StructuredGrid`.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `grid` | `np.ndarray (ni, nj, nk, 3)` | — | Full O-grid volume |
| `with_jacobian` | `bool` | `True` | Attach per-cell Jacobian as `cell_data["Jacobian"]` |

**Returns**

`pyvista.StructuredGrid` with dimensions `(ni, nj+1, nk)`.

- `n_points = ni * (nj+1) * nk`
- `cell_data["Jacobian"]` shape: `(ni-1) * nj * (nk-1)` (if `with_jacobian=True`)

**Example**

```python
from shore.io.geo import read_geo
from shore.io.vtk import geo_to_pyvista

grid = read_geo("sphere.geo")
sg = geo_to_pyvista(grid)
print(sg.dimensions)               # (ni, nj+1, nk)
print(sg.cell_data["Jacobian"].min())  # > 0 for a valid grid
```

---

### `surface_to_pyvista`

```python
def surface_to_pyvista(k0: np.ndarray) -> pv.StructuredGrid
```

Convert a SHORE surface layer (k=0) to a PyVista `StructuredGrid`.

**Parameters**

| Name | Type | Description |
|---|---|---|
| `k0` | `np.ndarray (ni, nj, 3)` | Structured surface layer |

**Returns**

`pyvista.StructuredGrid` with dimensions `(ni, nj+1, 1)`.

**Example**

```python
import numpy as np
from shore.io.vtk import surface_to_pyvista

k0 = np.load("sphere_k0.npz")["k0"]
sg = surface_to_pyvista(k0)
sg.plot(show_edges=True)
```

---

### `write_vts`

```python
def write_vts(
    path: str | Path,
    grid: np.ndarray,
    with_jacobian: bool = True,
) -> None
```

Write a volume grid to a VTK StructuredGrid file (`.vts`).

Appends the `.vts` extension if the path has no suffix.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `path` | `str \| Path` | — | Output path (`.vts` appended if needed) |
| `grid` | `np.ndarray (ni, nj, nk, 3)` | — | Full O-grid volume |
| `with_jacobian` | `bool` | `True` | Include Jacobian cell data |

**Example**

```python
from shore.io.geo import read_geo
from shore.io.vtk import write_vts

grid = read_geo("sphere.geo")
write_vts("sphere.vts", grid)
# → sphere.vts ready for ParaView
```

---

### `write_surface_vts`

```python
def write_surface_vts(path: str | Path, k0: np.ndarray) -> None
```

Write a surface layer to a VTK StructuredGrid file (`.vts`).

**Parameters**

| Name | Type | Description |
|---|---|---|
| `path` | `str \| Path` | Output path |
| `k0` | `np.ndarray (ni, nj, 3)` | Structured surface layer |

**Example**

```python
import numpy as np
from shore.io.vtk import write_surface_vts

k0 = np.load("sphere_k0.npz")["k0"]
write_surface_vts("sphere_surface.vts", k0)
```

## Point ordering

PyVista StructuredGrid expects points in Fortran (column-major, i-fastest) order. Internally, `geo_to_pyvista` and `surface_to_pyvista` use `ravel(order='F')` per coordinate after closing the j-direction, matching this requirement exactly.

## See also

- [`shore export`](export) — CLI wrapper for `write_vts` / `write_surface_vts`
- [`shore.vis`](api-vis) — interactive PyVista viewer built on top of this module
- [`shore.io.geo`](api-geo) — read `.geo` files before exporting
