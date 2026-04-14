# shore.surface.remesh

Unstructured STL surface → structured `(ni, nj)` quad grid via sphere projection.

## Functions

### `sphere_latlon_directions`

```python
def sphere_latlon_directions(ni: int, nj: int) -> np.ndarray
```

Build `(ni, nj)` unit direction vectors on a lat-lon sphere.

Poles are excluded (no singular points). The j-direction is periodic.

**Parameters**

| Name | Type | Description |
|---|---|---|
| `ni` | `int` | Latitude points (exclusive of poles) |
| `nj` | `int` | Longitude points (full circle, periodic) |

**Returns**

`numpy.ndarray` of shape `(ni, nj, 3)` — unit direction vectors.

---

### `project_sphere_to_surface`

```python
def project_sphere_to_surface(
    mesh: trimesh.Trimesh,
    ni: int,
    nj: int,
) -> np.ndarray
```

Project a structured sphere grid onto a triangulated surface.

**Algorithm:**
1. Build `(ni, nj)` lat-lon direction vectors centred at `mesh.centroid`
2. Cast `ni × nj` rays from the centroid outward
3. Find the first intersection of each ray with the STL triangulation (BVH, `rtree` backend)
4. Return the intersection points as the k=0 structured surface layer

**Parameters**

| Name | Type | Description |
|---|---|---|
| `mesh` | `trimesh.Trimesh` | The body surface (raw STL) |
| `ni` | `int` | Latitude resolution of the output structured grid |
| `nj` | `int` | Longitude resolution of the output structured grid |

**Returns**

`numpy.ndarray` of shape `(ni, nj, 3)` — coordinates of the k=0 structured surface layer.

**Raises**

| Exception | Condition |
|---|---|
| `ValueError` | Any ray misses the body. Reports the number of missed rays and suggests checking star-shaped property or repairing the STL. |

**Example**

```python
from shore.surface.io     import load_stl
from shore.surface.remesh import project_sphere_to_surface

mesh = load_stl("sphere.stl")
k0   = project_sphere_to_surface(mesh, ni=40, nj=60)

print(k0.shape)                                      # (40, 60, 3)
print(np.linalg.norm(k0, axis=-1).mean())            # ≈ 1.0 for unit sphere
```

## Notes

- **Star-shaped constraint**: the algorithm requires every ray from the centroid to hit the body exactly once. Bodies with concavities relative to their centroid will trigger the `ValueError`. See [Pipeline limitations](/concepts/pipeline#limitations).
- **Multiple hits**: `trimesh` is called with `multiple_hits=False` — only the nearest intersection per ray is used.
- **Performance**: for `ni=60, nj=80` (4,800 rays) on a mesh with ~10,000 triangles, projection takes ~50 ms on a modern CPU.
