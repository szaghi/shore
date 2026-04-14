# shore.volume.jacobian

Per-cell Jacobian computation for structured O-grids.

::: tip O-grid only
All functions in this module assume an **O-grid topology**: the j-index is periodic (`j=nj-1` and `j=0` are adjacent). Passing a plain Cartesian box grid will produce incorrect results because the periodic j-wrap creates a non-physical backward cell across the box boundary.
:::

## Functions

### `jacobian_field`

```python
def jacobian_field(grid: np.ndarray) -> np.ndarray
```

Compute the scalar-triple-product Jacobian for every cell in the volume.

$$J_{i,j,k} = (\mathbf{e}_2 \times \mathbf{e}_1) \cdot \mathbf{e}_3$$

where $\mathbf{e}_1 = P_{i,j+1,k} - P_{i,j,k}$ (j-tangent, periodic), $\mathbf{e}_2 = P_{i+1,j,k} - P_{i,j,k}$ (i-tangent), and $\mathbf{e}_3 = P_{i,j,k+1} - P_{i,j,k}$ (wall-normal).

**Parameters**

| Name | Type | Description |
|---|---|---|
| `grid` | `np.ndarray (ni, nj, nk, 3)` | Full structured O-grid volume |

**Returns**

`numpy.ndarray` of shape `(ni-1, nj, nk-1)` — one Jacobian value per cell. The j-dimension equals `nj` (not `nj-1`) because the periodic wrap closes the last cell.

**Example**

```python
from shore.volume.jacobian import jacobian_field
from shore.io.geo import read_geo

grid = read_geo("sphere.geo")
J = jacobian_field(grid)
print(J.shape)        # (ni-1, nj, nk-1)
print(J.min())        # > 0 for a valid O-grid
```

---

### `jacobian_per_layer`

```python
def jacobian_per_layer(grid: np.ndarray) -> np.ndarray
```

Minimum Jacobian over all `(i, j)` cells for each wall-normal layer interval.

**Parameters**

| Name | Type | Description |
|---|---|---|
| `grid` | `np.ndarray (ni, nj, nk, 3)` | Full structured O-grid volume |

**Returns**

`numpy.ndarray` of shape `(nk-1,)` — one value per layer interval `[k, k+1]`.

**Example**

```python
from shore.volume.jacobian import jacobian_per_layer
import numpy as np

jmins = jacobian_per_layer(grid)
for k, jmin in enumerate(jmins):
    print(f"layer {k+1:3d}  jmin={jmin:.3e}")
```

This is the same quantity that `shore check` displays in its table and that the hyperbolic marching loop monitors at every step.

---

### `min_jacobian_two_layers`

```python
def min_jacobian_two_layers(layer_k: np.ndarray, layer_k1: np.ndarray) -> float
```

Minimum Jacobian over all cells between two consecutive layers.

Equivalent to `jacobian_field(np.stack([layer_k, layer_k1], axis=2)).min()`, but without allocating the full 4-D array when only the scalar result is needed.

**Parameters**

| Name | Type | Description |
|---|---|---|
| `layer_k` | `np.ndarray (ni, nj, 3)` | Inner layer (k) |
| `layer_k1` | `np.ndarray (ni, nj, 3)` | Outer layer (k+1) |

**Returns**

`float` — minimum Jacobian. Negative if any cell is inverted.

**Example**

```python
from shore.volume.jacobian import min_jacobian_two_layers
import numpy as np
from tests.conftest import sphere_layer

inner = sphere_layer(20, 24, radius=1.0)
outer = sphere_layer(20, 24, radius=1.5)

jmin = min_jacobian_two_layers(inner, outer)
assert jmin > 0.0   # valid outward extrusion
```

Used internally by `march()` to check cell quality before accepting each new layer.

## Notes

### Sign convention

Positive Jacobian means the cell is correctly oriented: the outward normal of the cross product $\mathbf{e}_2 \times \mathbf{e}_1$ aligns with the wall-normal direction $\mathbf{e}_3$.

For a sphere O-grid with rays cast outward from the centroid, the lat-lon directions give:

- $\mathbf{e}_1$ (longitude, j) — eastward
- $\mathbf{e}_2$ (latitude, i) — southward
- $\mathbf{e}_2 \times \mathbf{e}_1$ — outward radial (positive Jacobian when $\mathbf{e}_3$ also points outward)

### Relationship to `shore check`

`shore check` calls `jacobian_per_layer` and reports the result per layer. The same function is used in the `march` loop via `min_jacobian_two_layers`.

## See also

- [`shore.volume.hyperbolic`](api-hyperbolic) — marching kernel that uses these functions
- [`shore check`](check) — CLI command built on `jacobian_per_layer`
- [Hyperbolic marching — Jacobian check](/concepts/hyperbolic-marching#jacobian-check) — conceptual explanation
