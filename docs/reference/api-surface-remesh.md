# shore.surface.remesh

Unstructured STL surface → structured `(ni, nj)` quad grid via lat-lon
sphere projection (single-block path).

::: tip Cubed-sphere topology uses different modules
The 6-block cubed-sphere topology built by `shore mesh` does **not** use the
functions on this page. It uses `shore.surface.cap.build_cap_k0`
(gnomonic flat-square cap seed) and `shore.surface.equator.build_equator`
(great-circle meridian equator with C1 seam pinning). The functions
below build a single-block lat-lon grid and are used by `shore project`,
the `Body.project()` method, and downstream tools that consume a
single-block surface array.
:::

## Functions

### `sphere_latlon_directions`

```python
def sphere_latlon_directions(
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    theta_margin: float | None = None,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
) -> np.ndarray
```

Build `(ni, nj)` unit direction vectors on a lat-lon sphere.

Poles are excluded (no singular points). The j-direction is periodic.
The equatorial band spans `[theta_cap, pi - theta_cap]` in colatitude.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `ni` | `int` | — | Latitude points (exclusive of poles) |
| `nj` | `int` | — | Longitude points (full circle, periodic) |
| `theta_cap` | `float \| None` | `0.05*pi` | Polar exclusion angle in radians. Must be in `(0, pi/2)`. The cap boundary parallel has radius `sin(theta_cap)` on the unit sphere. |
| `theta_cap_deg` | `float \| None` | None | Convenience alias for `theta_cap` in degrees. Mutually exclusive with `theta_cap`. |
| `theta_margin` | `float \| None` | None | Deprecated alias for `theta_cap` (radians). |
| `i_spacing` | `str` | `"uniform"` | Latitude spacing law: `"uniform"` / `"tanh"` / `"tanh2"`. |
| `i_beta` | `float` | `3.0` | Clustering strength for tanh / tanh2 on i. |
| `j_spacing` | `str` | `"uniform"` | Longitude spacing law: `"uniform"` or `"tanh2"`. |
| `j_beta` | `float` | `3.0` | Clustering strength for tanh2 on j. |
| `j_cluster_lon` | `float \| None` | None | Longitude (radians) the tanh2 j-distribution clusters toward. |

**Returns**

`numpy.ndarray` of shape `(ni, nj, 3)` — unit direction vectors.

**Raises**

| Exception | Condition |
|---|---|
| `ValueError` | More than one of `theta_cap`, `theta_cap_deg`, `theta_margin` supplied, or resolved angle outside `(0, pi/2)`. |

---

### `project_sphere_to_surface`

```python
def project_sphere_to_surface(
    mesh: trimesh.Trimesh,
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    theta_margin: float | None = None,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
    projector: object | None = None,
) -> np.ndarray
```

Project a structured sphere grid onto a triangulated surface.

**Algorithm:**
1. Build `(ni, nj)` lat-lon direction vectors centred at `mesh.centroid`
2. Cast `ni × nj` rays from the centroid outward
3. Find the first intersection of each ray with the STL triangulation (BVH, `rtree` backend)
4. Return the intersection points as the k=0 structured surface layer

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `mesh` | `trimesh.Trimesh` | — | The body surface (raw STL) |
| `ni` | `int` | — | Latitude resolution of the output structured grid |
| `nj` | `int` | — | Longitude resolution of the output structured grid |
| `theta_cap` | `float \| None` | `0.05*pi` | Polar exclusion angle in radians. Controls how far the equatorial band extends from the poles and therefore the physical size of the polar cap faces. Must be in `(0, pi/2)`. |
| `theta_cap_deg` | `float \| None` | None | Convenience alias for `theta_cap` in degrees. Mutually exclusive with `theta_cap`. |
| `theta_margin` | `float \| None` | None | Deprecated alias for `theta_cap` (radians). |
| `i_spacing` / `i_beta` | `str` / `float` | `"uniform"` / `3.0` | Latitude spacing law and clustering strength. |
| `j_spacing` / `j_beta` | `str` / `float` | `"uniform"` / `3.0` | Longitude spacing law (`"uniform"` / `"tanh2"`) and clustering strength. |
| `j_cluster_lon` | `float \| None` | None | Longitude (radians) the tanh2 j-distribution clusters toward. |
| `projector` | `object \| None` | None | A pre-built [`Projector`](./api-projector) with a cached BVH to reuse; when `None`, a fresh one is built from `mesh`. |

**Returns**

`numpy.ndarray` of shape `(ni, nj, 3)` — coordinates of the k=0 structured surface layer.

**Raises**

| Exception | Condition |
|---|---|
| `ValueError` | More than one angle parameter supplied, angle outside `(0, pi/2)`, or any ray misses the body. |

**Example — default cap size (~9°):**

```python
from shore.surface.io     import load_stl
from shore.surface.remesh import project_sphere_to_surface

mesh = load_stl("sphere.stl")
k0   = project_sphere_to_surface(mesh, ni=40, nj=60)

print(k0.shape)                            # (40, 60, 3)
print(np.linalg.norm(k0, axis=-1).mean()) # ≈ 1.0 for unit sphere
```

**Example — larger cap (20°), degrees alias:**

```python
k0 = project_sphere_to_surface(mesh, ni=40, nj=60, theta_cap_deg=20.0)
# Cap boundary parallel radius ≈ sin(20°) ≈ 0.342 × sphere radius
```

## theta_cap reference table

| `theta_cap` (°) | `theta_cap` (rad) | Cap boundary radius / sphere radius |
|---|---|---|
| 5° | 0.087 | 0.087 |
| 9° (default) | 0.157 | 0.156 |
| 15° | 0.262 | 0.259 |
| 20° | 0.349 | 0.342 |
| 30° | 0.524 | 0.500 |

The cap boundary parallel radius equals `sin(theta_cap)`. Larger `theta_cap`
gives a physically larger cap face and a shorter equatorial band.

## Notes

- **Star-shaped constraint**: the algorithm requires every ray from the centroid to hit the body exactly once. Bodies with concavities relative to their centroid will trigger the `ValueError`. See [Pipeline limitations](/guide/pipeline#limitations).
- **Multiple hits**: `trimesh` is called with `multiple_hits=False` — only the nearest intersection per ray is used.
- **Performance**: for `ni=60, nj=80` (4,800 rays) on a mesh with ~10,000 triangles, projection takes ~50 ms on a modern CPU.

## See also

- [Cubed-sphere topology](/topologies/cubed-sphere) — the production 6-block path
  (uses `shore.surface.cap.build_cap_k0` and
  `shore.surface.equator.build_equator`, not the functions on this page)
- [`shore.topology.sphere_cap`](/reference/api-sphere-cap) — marching
  primitives shared by the cubed-sphere topology
- `shore.surface.cap.build_cap_k0` — gnomonic flat-square cap k=0 seed
- `shore.surface.equator.build_equator` — equator k=0 from cap seam rings
