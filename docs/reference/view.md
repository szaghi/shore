# shore view

Open an interactive PyVista viewer for a `.geo` volume or `.npz` surface layer.

::: warning Optional dependency
`shore view` requires **pyvista** (`pip install shore-mesh[vis]`). The command exits with code 1 and a clear error message if pyvista is not installed.
:::

## Synopsis

```
shore view [OPTIONS] INPUT
```

## Arguments

| Argument | Description |
|---|---|
| `INPUT` | Path to a `.geo` volume file or a `.npz` surface layer file |

## Options

| Option | Default | Description |
|---|---|---|
| `--colorby TEXT` | `Jacobian` / `radius` | Scalar field used for colouring (`none` to disable) |
| `--no-edges` | off | Hide mesh edges |
| `--opacity FLOAT` | 0.6 | Volume opacity (0 = transparent, 1 = opaque) |
| `--help` | — | Show this help and exit |

## Output

An interactive 3-D PyVista window. For `.geo` input, the full O-grid volume is rendered with the wall layer (k=0) highlighted in steelblue. The outer layers are coloured by the per-cell Jacobian using a `RdYlGn` colormap (red = low quality, green = high quality) with the colour scale clamped to the 5th–95th percentile.

For `.npz` input, the projected surface layer is rendered and coloured by the distance from the origin (`radius`).

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Window closed normally |
| 1 | Input not found, unrecognised extension, or pyvista not installed |

## Examples

View a volume grid:

```bash
shore view sphere.geo
```

View only the surface layer:

```bash
shore view sphere_k0.npz
```

Colour by nothing (solid grey), hide edges:

```bash
shore view sphere.geo --colorby none --no-edges
```

## Jupyter usage

`shore view` opens a native window and is not suited for Jupyter notebooks. Use the Python API instead:

```python
import pyvista as pv
pv.set_jupyter_backend("html")

from shore.vis import view_geo
import numpy as np
from shore.io.geo import read_geo

grid = read_geo("sphere.geo")
pl = view_geo(grid, jupyter=True)
pl.show()
```

See [Quick Start — Visualization](/guide/quickstart#visualization) and [`shore.vis`](api-vis) for the full Jupyter workflow.

## See also

- [`shore export`](export) — write a `.vts` file for ParaView
- [`shore.vis`](api-vis) — Python API for interactive and off-screen rendering
- [`shore.io.vtk`](api-vtk) — Python API for VTK file export
