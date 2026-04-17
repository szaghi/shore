# Quick Start

## Your first O-grid

Given a body surface `sphere.stl`, generate a structured O-grid with 40 latitude points, 60 longitude points, and 30 wall-normal layers:

```bash
shore mesh sphere.stl --ni 40 --nj 60 --nk 30 --ds 1e-3 --growth 1.15 -o sphere.geo -v
```

SHORE will:

1. Load the STL and count triangles
2. Project a structured `(40×60)` sphere template onto the body surface
3. March 29 layers outward from the surface, growing each step by 15 %
4. Write `sphere.geo` in Xall_GPU format

The `-v` flag shows the live progress bar:

```
⠸ Loading STL     ████████████████████████  1/1   320 triangles         0:00:00
⠸ Projecting      ████████████████████████  1/1   (40×60) surface layer  0:00:00
⠸ Marching        ████████████████████████ 29/29  jmin=8.23e-03  Δs=4.28e-02  0:00:00
✓ Wrote sphere.geo  (40×60×30)
```

## Inspecting the result

```bash
shore info sphere.geo
```

```
┌──────────────────────────────────────────────┐
│ sphere.geo                                   │
├─────────────┬───────────────────────────────┤
│ Dimensions  │ ni=40  nj=60  nk=30           │
│ Points      │ 72,000                         │
│ Cells       │ 47,151                         │
│ x bounds    │ [-1.042,  1.042]               │
│ y bounds    │ [-1.042,  1.042]               │
│ z bounds    │ [-1.042,  1.042]               │
└─────────────┴───────────────────────────────┘
```

```bash
shore check sphere.geo
```

```
✓ All Jacobians positive.
```

## Running stages individually

For debugging, you can run each stage separately:

```bash
# Stage 1: project the surface only
shore project sphere.stl --ni 40 --nj 60 -o sphere_k0.npz -v

# Inspect the projected layer in Python
python -c "
import numpy as np
k0 = np.load('sphere_k0.npz')['k0']
print(k0.shape)        # (40, 60, 3)
import numpy.linalg as la
print(la.norm(k0, axis=-1).mean())  # ≈ 1.0
"

# Stage 2: extrude the saved layer
shore march sphere_k0.npz --nk 30 --ds 1e-3 --growth 1.15 -o sphere.geo -v
```

## Choosing marching parameters

| Parameter | Effect | Typical range |
|---|---|---|
| `--ds` | First-layer thickness (geometric spacing) | 1e-5 – 1e-2 (body units) |
| `--growth` | Layer-to-layer thickness ratio | 1.05 – 1.25 |
| `--smooth` | Laplacian blend per layer | 0.1 – 0.5 |
| `--smooth-iters` | Laplacian sweeps per layer | 1 – 5 |
| `--nk` | Total layers | 20 – 60 |

**Grid quality vs. speed trade-offs:**

- Reducing `--growth` gives a smoother grid but requires more layers for the same far-field extent.
- Increasing `--smooth` improves orthogonality at the cost of some cell-size uniformity.
- If `shore check` reports a near-inversion (jmin < 1e-6) on early layers, the body likely has concavities relative to its centroid — increase `--smooth` or reduce `--ds`.

## Controlling cell distribution with spacing laws

The default geometric spacing (`ds * growth^k`) works well when you only care about the first-layer thickness. When you need more control — a fixed far-field boundary, clustering around a geometric feature, or finer cells at the cap boundaries — switch to a tanh-based law.

### Wall-normal (k-direction)

```bash
# Geometric: fixed first-layer, unbounded outer boundary
shore mesh wing.stl --nk 40 --ds 1e-3 --growth 1.10 -o wing.geo

# Tanh: fixed total thickness, cluster near wall
shore mesh wing.stl --nk 40 --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 5.0 -o wing.geo

# Tanh2: fixed total thickness, cluster at wall AND far-field
shore mesh wing.stl --nk 40 --volume-k-spacing tanh2 --volume-k-beta 3.0 --volume-k-thickness 5.0 -o wing.geo
```

Key difference: with `geometric`, `--ds` controls the first layer and the outer boundary floats. With `tanh`/`tanh2`, `--volume-k-thickness` fixes the outer boundary and the first-layer size is determined by `nk` and `--volume-k-beta`.

### Surface (i- and j-directions)

Surface spacing is set at projection time and cannot be changed without re-projecting. Use `--surface-i-spacing` for latitude (pole-to-pole) and `--surface-j-spacing` for longitude (around the body).

```bash
# Longitude clustering at the leading edge (phi = 0)
shore mesh airfoil.stl \
  --surface-j-spacing tanh2 --surface-j-beta 5.0 --j-cluster-lon 0.0 \
  --ni 40 --nj 80 --nk 30

# Latitude clustering at cap boundaries (for caps topology)
shore mesh body.stl --topology caps \
  --surface-i-spacing tanh --surface-i-beta 4.0 \
  --cap-pull 0.4
```

The `--j-cluster-lon` value is a longitude in radians: `0.0` for the positive x-axis, `pi/2` for the y-axis, `pi` for the negative x-axis. Set it to the azimuthal angle of the geometric feature you want to refine.

See [Spacing laws](/concepts/spacing) for the full explanation of each law, what physical region it clusters toward, and how to choose $\beta$.

::: warning Star-shaped constraint
The surface projection assumes the body is **star-shaped** with respect to its centroid: every ray from the centroid must hit the surface exactly once. This is satisfied by spheres, ellipsoids, fuselages, and most convex CFD bodies. Concave bodies (e.g. re-entrant cavities, forked geometries) will produce a `ValueError` during projection — see [Limitations](/concepts/pipeline#limitations).
:::

## Using SHORE from Python

All CLI commands are thin wrappers over a pure Python API:

```python
from shore.surface.io       import load_stl
from shore.surface.remesh   import project_sphere_to_surface
from shore.volume.hyperbolic import march
from shore.volume.spacing   import build_k_steps
from shore.io.geo            import write_geo

mesh = load_stl("sphere.stl")
k0   = project_sphere_to_surface(mesh, ni=40, nj=60)

# Geometric spacing (default — same as ds=1e-3, growth=1.15)
grid = march(k0, nk=30, ds=1e-3, growth=1.15)

# Tanh spacing — fixed outer boundary, cluster at wall
steps = build_k_steps(30, ds=1e-3, spacing="tanh", beta=4.0, total_thickness=2.0)
grid = march(k0, nk=30, steps=steps)

write_geo("sphere.geo", grid)
```

The `march` function also accepts a `layer_callback` for custom progress reporting:

```python
def my_callback(k: int, jmin: float, step: float) -> None:
    print(f"layer {k:3d}  jmin={jmin:.3e}  step={step:.3e}")

grid = march(k0, nk=30, ds=1e-3, growth=1.15, layer_callback=my_callback)
```

## Visualization

SHORE ships two complementary visualization paths: **ParaView** via VTK export and **PyVista** for interactive Python views.

::: tip Install the vis extras
```bash
pip install shore[vis]   # adds pyvista
```
:::

### ParaView

Export any `.geo` or `.npz` file to VTK StructuredGrid format and open it in ParaView:

```bash
shore export sphere.geo           # → sphere.vts
shore export sphere_k0.npz        # → sphere_k0.vts
paraview sphere.vts
```

The `.vts` file contains the grid coordinates and, for volume grids, a `Jacobian` cell-data array. In ParaView:

1. Set representation to **Surface with Edges** to inspect cell structure.
2. Colour by **Jacobian** to spot low-quality cells.
3. Use **Clip** or **Slice** filters to examine interior layers.

### Interactive desktop viewer

```bash
shore view sphere.geo         # opens a PyVista window
shore view sphere_k0.npz      # view the projected surface
```

The wall layer (k=0) is highlighted in steelblue; outer layers are coloured by the Jacobian on a `RdYlGn` scale. Use `--colorby none` for a plain view.

### Jupyter notebook

```python
import pyvista as pv
pv.set_jupyter_backend("html")   # inline rendering, no display required

from shore.vis import view_geo, view_surface
from shore.io.geo import read_geo
import numpy as np

# Volume grid
grid = read_geo("sphere.geo")
pl = view_geo(grid, jupyter=True)
pl.show()

# Surface layer
k0 = np.load("sphere_k0.npz")["k0"]
pl = view_surface(k0, jupyter=True)
pl.show()
```

See `examples/visualize.ipynb` for a complete notebook that covers projection, marching, Jacobian analysis, and spacing-law comparison with plots.

### Automated screenshots

```python
from shore.vis import save_screenshot
from shore.io.geo import read_geo

grid = read_geo("sphere.geo")
save_screenshot(grid, "sphere_preview.png", window_size=(1600, 1200))
```

### Accessing Jacobian data directly

```python
from shore.volume.jacobian import jacobian_field, jacobian_per_layer
from shore.io.geo import read_geo
import numpy as np

grid = read_geo("sphere.geo")

J      = jacobian_field(grid)      # (ni-1, nj, nk-1) — per-cell
jmins  = jacobian_per_layer(grid)  # (nk-1,)  — per-layer minimum

import matplotlib.pyplot as plt
plt.plot(np.arange(1, len(jmins)+1), jmins, marker="o")
plt.xlabel("Layer")
plt.ylabel("Min Jacobian")
plt.yscale("log")
plt.title("Jacobian profile")
plt.show()
```
