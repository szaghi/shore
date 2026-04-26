# shore.volume.jacobian

Per-cell Jacobian computation for structured O-grids.

::: tip O-grid vs open topology
`jacobian_field` and `jacobian_per_layer` assume a **periodic** j-direction
(O-grid). For open (non-periodic) blocks — such as the equatorial sub-blocks
and cap blocks of the [sphere-cap topology](/topologies/caps) — use
`_jacobian_field_open` instead. `min_jacobian_two_layers` accepts a
`periodic_j` flag to select the correct behaviour automatically.
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
def min_jacobian_two_layers(
    layer_k: np.ndarray,
    layer_k1: np.ndarray,
    periodic_j: bool = True,
) -> float
```

Minimum Jacobian over all cells between two consecutive layers.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `layer_k` | `np.ndarray (ni, nj, 3)` | — | Inner layer (k) |
| `layer_k1` | `np.ndarray (ni, nj, 3)` | — | Outer layer (k+1) |
| `periodic_j` | `bool` | `True` | When `True` the j-direction wraps (O-grid, `nj` cells). When `False` it is open (`nj-1` cells). |

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

---

### `_jacobian_field_open` *(internal)*

```python
def _jacobian_field_open(grid: np.ndarray) -> np.ndarray
```

Jacobian field for a non-periodic (open) structured volume. Identical to
`jacobian_field` but the j-direction is not wrapped: there are `nj-1` cells
in j rather than `nj`.

**Parameters**

| Name | Type | Description |
|---|---|---|
| `grid` | `np.ndarray (ni, nj, nk, 3)` | Open structured volume |

**Returns**

`numpy.ndarray` of shape `(ni-1, nj-1, nk-1)`.

Used by the [sphere-cap topology](/topologies/caps) for equatorial
sub-blocks and cap blocks, which have open boundaries in both i and j.

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
- [Hyperbolic marching — Jacobian check](/algorithm/hyperbolic-marching#jacobian-check) — conceptual explanation
