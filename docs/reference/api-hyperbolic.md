# shore.volume.hyperbolic

Hyperbolic-marching primitives — per-layer normal estimation,
Laplacian smoothing, and STL-normal query — used by the OOP mesh
marcher. **The user-facing entry point is
[`Mesh.march()`](./api-mesh)**, not the functions on this page.

::: tip Internal helpers
Every function in `shore.volume.hyperbolic` is underscore-prefixed
and intended for internal use by `shore.mesh.Mesh._march_caps` /
`Mesh._march_ogrid` and the cubed-sphere marching primitives in
[`shore.topology.sphere_cap`](./api-sphere-cap). They are documented
here for users who want to read the kernel or build a new
topology.
:::

## `_compute_normals`

```python
def _compute_normals(
    layer: np.ndarray,
    periodic: tuple[bool, bool] = (False, True),
    ghost_rows: dict[str, np.ndarray] | None = None,
) -> np.ndarray
```

Estimate outward unit normals for every node of a structured
`(ni, nj, 3)` layer via the cross product of the in-plane tangents.
Central differences in both directions; one-sided fallback at
non-periodic boundaries unless a ghost row is supplied for that side.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `layer` | `np.ndarray (ni, nj, 3)` | — | Node coordinates of one wall-normal layer. |
| `periodic` | `tuple[bool, bool]` | `(False, True)` | `(periodic_i, periodic_j)`. When `True` the corresponding direction wraps around for end nodes; when `False` the boundary uses one-sided differences (or a ghost row, if supplied). Default matches the legacy single-block O-grid (i open, j periodic). Cubed-sphere blocks use `(False, False)`. |
| `ghost_rows` | `dict[str, np.ndarray] \| None` | `None` | Virtual neighbour rows for any of the 4 faces. Keys: `'i_lo'` / `'i_hi'` give shape `(nj, 3)`; `'j_lo'` / `'j_hi'` give shape `(ni, 3)`. When a key is present the corresponding boundary's tangent becomes the centred difference between the next interior row and the ghost row — the seam-aware kernel that keeps multi-block spacing continuous during marching. |

**Returns**

`numpy.ndarray (ni, nj, 3)` — unit outward normals.

**Example — single-block O-grid normals**

```python
import numpy as np
from shore.volume.hyperbolic import _compute_normals

layer = np.load("sphere_k0.npz")["k0"]   # (ni, nj, 3)
n = _compute_normals(layer, periodic=(False, True))
print(np.linalg.norm(n, axis=-1).max())  # ≈ 1.0
```

**Example — cubed-sphere block with seam-aware ghost rows**

```python
ghost = {"i_hi": cap_north_first_row, "j_lo": sub3_pole_row}
n = _compute_normals(sub0_layer, periodic=(False, False), ghost_rows=ghost)
```

See [Seam-aware normals](/algorithm/seam-aware-normals) for the
algorithm and [Cubed-sphere topology](/topologies/cubed-sphere) for
the wiring.

---

## `_laplacian_smooth`

```python
def _laplacian_smooth(
    layer: np.ndarray,
    strength: float,
    periodic: tuple[bool, bool] = (False, True),
) -> np.ndarray
```

One Jacobi sweep of in-plane 4-neighbour averaging:

```
P_smoothed = P + strength * (P_laplacian - P)
```

with `P_laplacian` the average of the 4 in-plane neighbours.
Non-periodic boundary rows / columns are clamped (never moved); the
k=0 wall row is always clamped regardless of periodicity.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `layer` | `np.ndarray (ni, nj, 3)` | — | Node coordinates of one layer. |
| `strength` | `float` | — | Blend factor in `[0, 1]`. 0 = no change; 1 = full average. Typical: 0.2. |
| `periodic` | `tuple[bool, bool]` | `(False, True)` | Same semantics as `_compute_normals`. Non-periodic boundaries are clamped. |

**Returns**

`numpy.ndarray (ni, nj, 3)` — the smoothed layer.

`Mesh.march()` calls this `smooth_iters` times per accepted layer
(default 2 sweeps).

---

## `_stl_normals`

```python
def _stl_normals(
    trimesh_obj: object,
    points: np.ndarray,
    interpolate: bool = True,
) -> np.ndarray
```

Query surface normals at arbitrary world-space points by
closest-point projection onto a `trimesh.Trimesh`. Used by the
optional `--blend-normals-k` feature of `shore mesh` to mix
grid-computed normals with true STL surface normals near the wall.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `trimesh_obj` | `trimesh.Trimesh` | — | The body STL. |
| `points` | `np.ndarray (N, 3)` | — | World-space query points. |
| `interpolate` | `bool` | `True` | When `True`, returns the barycentric blend of the closest triangle's three vertex normals at the closest point — C0-continuous across triangle edges, suppresses facet-aligned banding on coarse STLs. When `False`, returns the flat face normal of the closest triangle (legacy behaviour, kept for regression debugging). |

**Returns**

`numpy.ndarray (N, 3)` — unit surface normals.

::: warning Mesh hygiene
Vertex normals are area-weighted by trimesh and require a clean
mesh (consistent face orientation, no degenerate faces, merged
vertices). [`shore.surface.io.load_stl`](./api-surface-io) enforces
all three. Bypassing it on an unclean mesh produces wrong vertex
normals and a banded blended march.
:::

## See also

- [`shore.mesh.Mesh`](./api-mesh) — the public OOP entrypoint that
  drives the marcher.
- [`shore.topology.sphere_cap`](./api-sphere-cap) — single-layer
  cubed-sphere marching primitives that orchestrate
  `_compute_normals` + `_laplacian_smooth` per block.
- [`shore.volume.jacobian`](./api-jacobian) — per-cell Jacobian
  field; the marcher's per-layer guard against inversion.
- [Hyperbolic marching](/algorithm/hyperbolic-marching) — the
  algorithm in detail.
- [Seam-aware normals](/algorithm/seam-aware-normals) — what
  `ghost_rows` enables.
