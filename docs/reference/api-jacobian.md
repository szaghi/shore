# shore.volume.jacobian

Per-cell Jacobian computation for structured O-grids.

::: tip O-grid vs open topology
`jacobian_field` assumes a **periodic** j-direction (O-grid). For open
(non-periodic) blocks — such as the equatorial sub-blocks and cap blocks of the
[cubed-sphere topology](/topologies/cubed-sphere) — use `_jacobian_field_open`
instead. The per-layer minimum-Jacobian guard used by the marcher and by
`Mesh.check` is built on these two field functions (see
[Relationship to `shore check`](#relationship-to-shore-check)).
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

::: tip Per-layer minimum Jacobian
There is no standalone `jacobian_per_layer` helper. The per-layer minimum that
`shore check` tabulates is computed inside [`Mesh.check`](./api-mesh#meshcheck)
and `shore.mesh.inspect`: for each layer interval `[k, k+1]` they slice the two
layers, stack them into a 2-layer volume, call `jacobian_field` (periodic-j) or
`_jacobian_field_open` (open block), and take `jac.min()`. The marcher applies
the same per-layer guard through `shore.topology.sphere_cap._check_jacobian`
(see [`api-sphere-cap`](./api-sphere-cap)).
:::

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

Used by the [cubed-sphere topology](/topologies/cubed-sphere) for equatorial
sub-blocks and cap blocks, which have open boundaries in both i and j.

## Notes

### Sign convention

Positive Jacobian means the cell is correctly oriented: the outward normal of the cross product $\mathbf{e}_2 \times \mathbf{e}_1$ aligns with the wall-normal direction $\mathbf{e}_3$.

For a sphere O-grid with rays cast outward from the centroid, the lat-lon directions give:

- $\mathbf{e}_1$ (longitude, j) — eastward
- $\mathbf{e}_2$ (latitude, i) — southward
- $\mathbf{e}_2 \times \mathbf{e}_1$ — outward radial (positive Jacobian when $\mathbf{e}_3$ also points outward)

### Relationship to `shore check`

`shore check` calls [`Mesh.check`](./api-mesh#meshcheck), which computes the
per-layer minimum Jacobian per block by slicing each consecutive layer pair into
a 2-layer volume and reducing `jacobian_field` / `_jacobian_field_open` with
`.min()`. The marcher applies the same per-layer guard at every accepted layer
via `shore.topology.sphere_cap._check_jacobian`.

## See also

- [`shore.volume.hyperbolic`](api-hyperbolic) — marching kernel that uses these functions
- [`shore.topology.sphere_cap`](./api-sphere-cap) — `_check_jacobian`, the per-layer guard
- [`shore check`](check) — CLI command built on `Mesh.check`
- [Hyperbolic marching — Jacobian check](/algorithm/hyperbolic-marching#jacobian-check) — conceptual explanation
