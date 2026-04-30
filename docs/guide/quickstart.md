# Quick Start

::: tip Not installed yet?
The fastest install on Linux / macOS / WSL:

```bash
curl -fsSL https://raw.githubusercontent.com/szaghi/shore/main/install.sh | bash
```

See the full [Installation guide](./installation.md) for `pipx` and
`pip` alternatives.
:::

## Your first O-grid

Given a body surface `sphere.stl`, generate a 6-block cubed-sphere
volume grid with 40 equator latitude **nodes**, 60 longitude **nodes**,
and 30 wall-normal **nodes** (29 cells extruded from the wall):

```bash
shore mesh sphere.stl --ni 40 --nj 60 --nk 30 --ds 1e-3 --growth 1.15 -o sphere -v
```

::: tip Nodes vs. cells
`--ni`, `--nj`, `--nk` are **node counts** — the array shape of the
`(ni, nj, nk, 3)` block SHORE produces.  ``--nk 30`` means 30
wall-normal node rows including the wall row at k=0, which is 29
extruded cells.  The only file SHORE writes that uses cell counts
is the binary ``.grd`` (geogrd Fortran convention); the conversion
is automatic.
:::

SHORE will:

1. Load the STL and count triangles
2. Build the south and north cap k=0 seed faces by gnomonic flat-square
   projection (4 corners on the colatitude-`theta_cap` parallel)
3. Build the equator k=0 from the cap seam rings via great-circle
   meridians (with C1 spacing match at every seam vertex)
4. March 29 cells outward with seam-aware normals, Laplacian smoothing,
   and per-layer Jacobian guards
5. Write 6 `.geo` files: `sphere_sub0.geo` ... `sphere_cap_south.geo`

The `-v` flag shows the live progress bar:

```
⠸ Loading STL     ████████████████  1/1   5,120 triangles  0:00:00
⠸ Marching        ████████████████ 29/29  jmin=2.30e-05    0:00:00
✓ Wrote 6 .geo files  (sub: 40x16x30, caps: 16x16x30)
```

## Inspecting the result

```bash
shore info sphere_sub0.geo
```

```
┌─────────────────────────────────────────────────┐
│ sphere_sub0.geo                                 │
├─────────────┬───────────────────────────────────┤
│ Dimensions  │ ni=40  nj=16  nk=30               │
│ Points      │ 19,200                            │
│ Cells       │ 16,965                            │
│ x bounds    │ [-1.042,  1.042]                  │
│ y bounds    │ [-1.042,  1.042]                  │
│ z bounds    │ [-1.042,  1.042]                  │
└─────────────┴───────────────────────────────────┘
```

```bash
shore check sphere_sub0.geo
```

```
✓ All Jacobians positive.
```

## Choosing marching parameters

| Parameter | Effect | Typical range |
|---|---|---|
| `--ds` | First-layer thickness (geometric spacing) | 1e-5 – 1e-2 (body units) |
| `--growth` | Layer-to-layer thickness ratio | 1.05 – 1.25 |
| `--smooth` | Laplacian blend per layer | 0.1 – 0.5 |
| `--smooth-iters` | Laplacian sweeps per layer | 1 – 5 |
| `--nk` | Wall-normal node count (includes the wall row at k=0); the marcher extrudes `nk - 1` cells | 20 – 60 |

**Grid quality vs. speed trade-offs:**

- Reducing `--growth` gives a smoother grid but requires more layers for
  the same far-field extent.
- Increasing `--smooth` improves orthogonality at the cost of some
  cell-size uniformity.
- If `shore check` reports a near-inversion (jmin < 1e-6) on early
  layers, the body likely has concavities relative to the centroid —
  increase `--smooth` or reduce `--ds`.

## Cap / equator compatibility

The cubed-sphere topology pins the equator's first/last meridional step to the
cap-edge first cell at every seam vertex (C1 spacing match). For the
pinning to be active, the cap radial step and the equator's natural
radial step should be within ~2× of each other:

$$
\frac{\theta_{\text{cap}}}{n_c - 1} \approx \frac{\pi - 2\theta_{\text{cap}}}{n_i - 1},
\quad n_c = \tfrac{n_j}{4} + 1
$$

The defaults `--ni 40 --nj 60 --theta-cap` (= $\pi/6$) land in the
compatibility window. If you push `theta_cap` very small with a coarse
`ni`, SHORE falls back to unpinned spacing and emits one
`UserWarning` per build with guidance.

## Controlling cell distribution with spacing laws

The default geometric spacing (`ds * growth^k`) works well when you only
care about the first-layer thickness. When you need more control —
fixed far-field, clustering at cap seams or cap corners — switch to a
tanh-based law.

### Wall-normal (k-direction)

```bash
# Geometric: fixed first-layer, unbounded outer boundary
shore mesh wing.stl --nk 40 --ds 1e-3 --growth 1.10 -o wing

# Tanh: fixed total thickness, cluster near wall
shore mesh wing.stl --nk 40 \
    --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 5.0 \
    -o wing

# Tanh2: fixed total thickness, cluster at wall AND far-field
shore mesh wing.stl --nk 40 \
    --volume-k-spacing tanh2 --volume-k-beta 3.0 --volume-k-thickness 5.0 \
    -o wing
```

Key difference: with `geometric`, `--ds` controls the first layer and
the outer boundary floats. With `tanh`/`tanh2`, `--volume-k-thickness`
fixes the outer boundary and the first-layer size is determined by
`nk` and `--volume-k-beta`.

### Surface (i- and j-directions)

```bash
# Equator meridional clustering at the cap seams (tanh i-spacing)
shore mesh body.stl --ni 60 --nj 80 \
    --surface-i-spacing tanh --surface-i-beta 4.0 \
    --theta-cap-deg 30 \
    -o body

# Cap-corner clustering (tanh2 j-spacing on the cap edges,
# inherited by the equator seam)
shore mesh body.stl --ni 40 --nj 60 \
    --surface-j-spacing tanh2 --surface-j-beta 4.0 \
    --theta-cap-deg 30 \
    -o body
```

See [Spacing laws](/algorithm/spacing-laws) for the full explanation of each
law, what physical region it clusters toward, and how to choose $\beta$.

::: warning Star-shaped constraint
The surface projection assumes the body is **star-shaped** with respect
to its anchor (centroid by default). Every ray from the anchor must hit
the surface exactly once. This is satisfied by spheres, ellipsoids,
fuselages, and most convex CFD bodies. Concave bodies (re-entrant
cavities, forked geometries) will fail during projection — see
[Limitations](/guide/pipeline#limitations).
:::

## Using SHORE from Python

The CLI is a thin wrapper over the OOP Python API:

```python
import numpy as np
from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

m = Mesh.from_stl(
    "sphere.stl",
    topology=CubedSphereTopology(),
    ni=40, nj=60, nk=30,
    theta_cap=float(np.pi / 6),
    spacing_k=Spacing1D(law="geometric", ds=1e-3, growth=1.15),
)

m.march(ds=1e-3, growth=1.15, smooth_strength=0.2, smooth_iters=2)
m.write_geo("sphere")  # writes 6 sphere_<label>.geo files
```

The `march` method also accepts a `layer_callback` for custom progress
reporting:

```python
def my_callback(k: int, jmin: float) -> None:
    print(f"layer {k:3d}  jmin={jmin:.3e}")

m.march(ds=1e-3, growth=1.15, layer_callback=my_callback)
```

For wall-normal spacing other than geometric, build a `Spacing1D`:

```python
from shore.mesh import Spacing1D

# Tanh spacing — cluster at wall, fixed total thickness
sp_k = Spacing1D(law="tanh", beta=4.0, total_thickness=2.0)
steps = sp_k.build_steps(40)

m = Mesh.from_stl(
    "wing.stl", topology=CubedSphereTopology(),
    ni=60, nj=80, nk=40,
    spacing_k=sp_k,
)
m.march(steps=steps)
```

## Visualization

SHORE ships two complementary visualization paths: **ParaView** via VTK
export and **PyVista** for interactive Python views.

::: tip Install the vis extras
```bash
pip install shore-mesh[vis]   # adds pyvista
```
:::

### ParaView

Export the 6 `.geo` files to a single VTK MultiBlock `.vtm` and open
it in ParaView:

```bash
shore export sphere_sub0.geo  # → sphere_sub0.vts
# ... or via the Python API for MultiBlock:
python -c "
from shore.io.geo import read_geo
from shore.io.vtk import write_caps_vtm
labels = ['sub0', 'sub1', 'sub2', 'sub3', 'cap_north', 'cap_south']
blocks = [read_geo(f'sphere_{lbl}.geo') for lbl in labels]
write_caps_vtm('sphere.vtm', blocks, with_jacobian=True)
"
paraview sphere.vtm
```

### Interactive desktop viewer

```bash
shore view sphere_sub0.geo   # opens a PyVista window
```

The wall layer (k=0) is highlighted; outer layers are coloured by the
Jacobian on a `RdYlGn` scale. Use `--colorby none` for a plain view.

### Jupyter notebook

```python
import pyvista as pv
pv.set_jupyter_backend("html")   # inline rendering, no display required

from shore.io.geo import read_geo
from shore.io.vtk import open_to_pyvista

grid = read_geo("sphere_sub0.geo")
sg = open_to_pyvista(grid, with_jacobian=True)
pl = pv.Plotter()
pl.add_mesh(sg, scalars="Jacobian", cmap="RdYlGn", show_edges=True, opacity=0.7)
pl.show()
```

See `examples/01-body-fitted/tour.ipynb` for a complete notebook that covers
mesh generation, Jacobian analysis, spacing-law comparison, and the C1
seam-matching demo.

### Accessing Jacobian data directly

```python
from shore.volume.jacobian import jacobian_field, jacobian_per_layer
from shore.io.geo import read_geo
import numpy as np
import matplotlib.pyplot as plt

grid = read_geo("sphere_sub0.geo")

J = jacobian_field(grid)             # (ni-1, nj, nk-1) — per-cell
jmins = jacobian_per_layer(grid)     # (nk-1,)  — per-layer minimum

plt.plot(np.arange(1, len(jmins)+1), jmins, marker="o")
plt.xlabel("Layer")
plt.ylabel("Min Jacobian")
plt.yscale("log")
plt.title("Jacobian profile")
plt.show()
```
