# shore.topology.sphere_cap

Polar cap block generation for sphere O-grids. Closes the two polar holes of a
lat-lon O-grid by splitting the equatorial block into four 90° sub-blocks and
extruding two hexahedral cap blocks, one per pole.

See [Sphere polar cap topology](/concepts/sphere-cap) for a conceptual
overview.

## Public API

### `split_ogrid_with_polar_caps`

```python
def split_ogrid_with_polar_caps(
    grid: np.ndarray,
    nk: int,
    ds: float,
    growth: float = 1.15,
    smooth_strength: float = 0.2,
    smooth_iters: int = 2,
    cap_pull: float = 0.0,
    cap_squareness: float = 0.0,
    layer_callback: Callable[[int, float], None] | None = None,
) -> list[np.ndarray]
```

Split a lat-lon O-grid surface into 4 equatorial sub-blocks, build polar cap
seed faces, and march all 6 blocks simultaneously by hyperbolic extrusion.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `grid` | `np.ndarray (ni, nj, 3)` | — | k=0 lat-lon O-grid surface layer |
| `nk` | `int` | — | Total wall-normal layers including k=0 |
| `ds` | `float` | — | First-layer thickness |
| `growth` | `float` | 1.15 | Geometric growth ratio between layers |
| `smooth_strength` | `float` | 0.2 | Laplacian blend factor [0–1] |
| `smooth_iters` | `int` | 2 | Laplacian sweeps per layer |
| `cap_pull` | `float` | 0.0 | Dome blend factor in [0, 1). Pulls the cap interior toward the pole, improving interior cell quality. 0 = flat Laplacian cap (default). Values 0.3–0.5 give best results; above ~0.6 may cause Jacobian failures for aggressive ds/growth. See [cap_pull](/concepts/sphere-cap#cap_pull). |
| `cap_squareness` | `float` | 0.0 | Blend factor in [0, 1] that straightens each cap boundary edge from a circular arc toward the geodesic chord connecting its corner nodes. 0 = parallel-circle arc (~157° corner angle, default). 1 = inscribed-square chord (~90° corner angle). Full range is numerically stable. Applied at every layer. See [cap_squareness](/concepts/sphere-cap#cap_squareness). |
| `layer_callback` | `Callable[[int, float], None] \| None` | None | Called after each accepted layer with `(k, jmin)` where `k` is 1-based and `jmin` is the global minimum Jacobian across all 6 blocks |

**Returns**

`list[np.ndarray]` of length 6, in order:

| Index | Name | Shape |
|---|---|---|
| 0 | `sub0` | `(ni, nc, nk, 3)` |
| 1 | `sub1` | `(ni, nc, nk, 3)` |
| 2 | `sub2` | `(ni, nc, nk, 3)` |
| 3 | `sub3` | `(ni, nc, nk, 3)` |
| 4 | `cap_north` | `(nc, nc, nk, 3)` |
| 5 | `cap_south` | `(nc, nc, nk, 3)` |

where `nc = nj // 4 + 1`.

**Raises**

| Exception | Condition |
|---|---|
| `RuntimeError` | Negative Jacobian detected in any sub-block or cap interior cell. Message includes the block label and layer index. |

::: warning nk limit for unit-sphere grids
With `ds = 0.02` and `growth = 1.10` on a unit sphere, keep `nk ≤ 10`
(total thickness ~0.32). Pushing beyond ~0.6× the sphere radius inverts
equatorial cells near the poles.
:::

**Example**

```python
import numpy as np
from tests.conftest import sphere_layer
from shore.topology.sphere_cap import split_ogrid_with_polar_caps
from shore.io.vtk import write_caps_vtm

surface = sphere_layer(40, 60)

def on_layer(k: int, jmin: float) -> None:
    print(f"k={k:3d}  jmin={jmin:.3e}")

blocks = split_ogrid_with_polar_caps(
    surface, nk=10, ds=0.02, growth=1.10,
    cap_squareness=0.5,    # straighten boundary edges toward inscribed square
    cap_pull=0.4,          # dome cap interior toward pole
    layer_callback=on_layer,
)

# blocks[0..3] — equatorial sub-blocks (ni × nc × nk × 3)
# blocks[4]   — cap_north             (nc × nc × nk × 3)
# blocks[5]   — cap_south             (nc × nc × nk × 3)
write_caps_vtm("sphere_caps.vtm", blocks, with_jacobian=True)
```

---

### `_adjust_nj` *(internal)*

```python
def _adjust_nj(nj: int) -> int
```

Round `nj` down to the nearest multiple of 4. Emits `warnings.warn` if
adjusted; returns `nj` unchanged if already divisible by 4.

---

### `_gnomonic_seed_face` *(internal)*

```python
def _gnomonic_seed_face(corners: np.ndarray, nc: int) -> np.ndarray
```

Build a gnomonic seed face from four corner points on a sphere using bilinear
interpolation followed by radial projection. Returns shape `(nc, nc, 3)`.
Used as a utility; not called in the main marching loop (the seed face is
built by `_build_cap_k0` instead).

---

### `_build_cap_k0` *(internal)*

```python
def _build_cap_k0(
    eq_rows: list[np.ndarray],
    nc: int,
    pole: str = "north",
) -> np.ndarray
```

Build the k=0 cap seed face by solving the discrete Laplacian with Dirichlet
boundary conditions. Boundary edges are set from the four equatorial pole rows;
interior nodes converge from zero to the minimum-energy surface.

---

### `_overwrite_cap_boundaries` *(internal)*

```python
def _overwrite_cap_boundaries(
    cap_layer: np.ndarray,
    eq_layers: list[np.ndarray],
    pole: str = "north",
) -> None
```

Overwrite the four boundary edges of a cap layer with the equatorial
sub-block pole rows. The mapping differs between the north and south caps due
to opposite face windings:

**North cap:**

| Cap edge | Equatorial source |
|---|---|
| `cap[0, :]` | `sub0[ni-1, :]` |
| `cap[:, -1]` | `sub1[ni-1, :]` |
| `cap[-1, :]` | `sub2[ni-1, ::-1]` |
| `cap[:, 0]` | `sub3[ni-1, ::-1]` |

**South cap:**

| Cap edge | Equatorial source |
|---|---|
| `cap[:, 0]` | `sub0[0, :]` |
| `cap[-1, :]` | `sub1[0, :]` |
| `cap[:, -1]` | `sub2[0, ::-1]` |
| `cap[0, :]` | `sub3[0, ::-1]` |

---

### `_squareness_rows` / `_squareness_row` *(internal)*

```python
def _squareness_rows(
    eq_rows: list[np.ndarray],
    squareness: float,
) -> list[np.ndarray]

def _squareness_row(
    row: np.ndarray,
    squareness: float,
    sphere_r: float,
) -> np.ndarray
```

Straighten each cap boundary edge from a circular arc toward a geodesic chord.
`_squareness_row` operates on one edge: for each intermediate node at
parametric position `t`, it slerp-interpolates the unit vector between the two
corner directions, blends with the original arc direction, and re-projects onto
the sphere.  `_squareness_rows` applies `_squareness_row` to all 4 edges,
estimating the sphere radius from the corner node radii.

---

### `_pull_cap_interior` *(internal)*

```python
def _pull_cap_interior(
    seed: np.ndarray,
    eq_rows: list[np.ndarray],
    cap_pull: float,
    pole: str = "north",
) -> np.ndarray
```

Pull the interior nodes of a cap k=0 seed face toward the geographic pole
and re-project onto the sphere. The 4 boundary edges are kept exactly
unchanged (boundary nodes are overwritten after the pull). The pull weight is
a 2-D Hann product that peaks at the face centre and falls to zero at all 4
boundary edges. Returns a new array; the input is not modified.

---

### `_march_cap_layer` *(internal)*

```python
def _march_cap_layer(
    layer: np.ndarray,
    step: float,
    smooth_strength: float,
    smooth_iters: int,
    eq_rows: list[np.ndarray],
    pole: str = "north",
) -> np.ndarray
```

Advance one cap layer by:
1. Normal-based displacement (same as equatorial blocks).
2. Laplacian smoothing (same strength and iteration count).
3. Dirichlet boundary pinning — overwrite the four boundary edges with the
   equatorial pole rows for this layer.

::: tip Internal functions
Functions prefixed with `_` are internal implementation details. Use
`split_ogrid_with_polar_caps` for all application code.
:::

## VTK export

The six blocks can be exported to a ParaView MultiBlock file via
`shore.io.vtk.write_caps_vtm`:

```python
from shore.io.vtk import write_caps_vtm

write_caps_vtm("sphere_caps.vtm", blocks, with_jacobian=True)
```

This writes:
- Six `.vts` structured grid files — one per block.
- One `.vtm` MultiBlock file that groups them under the names
  `sub0`, `sub1`, `sub2`, `sub3`, `cap_north`, `cap_south`.

Each block uses the open (non-periodic) Jacobian `_jacobian_field_open` for
cell data. The four corner cells of each cap are topologically degenerate and
are excluded from the Jacobian check; they may show near-zero values in
ParaView.

## See also

- [Sphere polar cap topology](/concepts/sphere-cap) — conceptual overview
- [`shore.volume.hyperbolic`](api-hyperbolic) — underlying marching kernel
- [`shore.volume.jacobian`](api-jacobian) — Jacobian functions used internally
- [`shore.io.vtk`](api-vtk) — VTK export including `write_caps_vtm`
