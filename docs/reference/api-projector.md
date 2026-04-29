# shore.surface.projector

Ray-mesh intersector with a single cached BVH, used by every SHORE
surface-projection step (cubed-sphere caps, equator, single-block
O-grid). Wraps `trimesh.ray.ray_triangle.RayMeshIntersector` with a
batch entry point that returns rich per-ray status, plus a typed
result dataclass.

For the algorithm, see [Surface projection](/algorithm/projection).

## `Projector`

```python
class Projector:
    mesh:    trimesh.Trimesh
    anchor:  np.ndarray   # shape (3,)

    def __init__(self, mesh: trimesh.Trimesh, anchor: np.ndarray) -> None: ...
    def project(
        self,
        dirs: np.ndarray,
        on_failure: str = "raise",
    ) -> ProjectionResult: ...
```

The constructor builds the BVH (an `rtree`-backed
`RayMeshIntersector`) once and caches it; subsequent
`Projector.project(dirs)` calls reuse the same acceleration
structure, so per-ray cost is `O(log T)` in the triangle count.

### `Projector.project`

```python
def project(self, dirs: np.ndarray, on_failure: str = "raise") -> ProjectionResult
```

Cast rays from `self.anchor` in directions `dirs` against the
cached body, vectorised across all rays in one call.

| Parameter | Description |
|---|---|
| `dirs` | Unit direction vectors, shape `(..., 3)`. Flattened internally. |
| `on_failure` | `"raise"` (default) — raise `MeshInputError` if any ray is not `HIT`. `"nearest"` — fall back to `trimesh.proximity.closest_point` for non-`HIT` rays and flag them in the status array. `"skip"` — leave non-`HIT` locations as `NaN`. |

**Returns** a [`ProjectionResult`](#projectionresult).

**Per-ray validation.** Internally each ray's hit is classified by
priority `DEGENERATE > BACKFACE > TOO_NEAR > TOO_FAR > HIT`:

- `BACKFACE` — `dot(face_normal, ray_direction) <= 0`. Casting from
  the body's anchor outward, a valid hit must show positive dot
  product with the outward face normal; a negative dot means the
  ray hit a triangle from the wrong side (typically a non-watertight
  mesh).
- `TOO_NEAR` — hit within `0.05 * bbox_diagonal` of the anchor;
  usually means the anchor grazes a face.
- `TOO_FAR` — hit further than `5 * bbox_diagonal`; usually a numeric
  edge case.
- `DEGENERATE` — `NaN` in the hit coordinates.

The `on_failure="raise"` default surfaces these conditions with a
detailed `MeshInputError` (counts of each non-`HIT` status) so that
the caller can take the appropriate corrective action (repair the
STL, change the anchor strategy).

## `ProjectionResult`

```python
@dataclass
class ProjectionResult:
    locations:  np.ndarray   # (n, 3) — hit positions; NaN for non-HIT rays
    status:     np.ndarray   # (n,)   — ProjectionStatus per ray
    distances:  np.ndarray   # (n,)   — anchor-to-hit distance; NaN for non-HIT
```

| Property | Returns |
|---|---|
| `result.hit_mask` | Boolean array, `True` where `status == HIT`. |
| `result.all_hit` | `True` iff every ray has `status == HIT`. |

## `ProjectionStatus`

```python
class ProjectionStatus(Enum):
    HIT         = "hit"
    MISSED      = "missed"
    BACKFACE    = "backface"
    DEGENERATE  = "degenerate"
    TOO_NEAR    = "too_near"
    TOO_FAR     = "too_far"
```

## Example

```python
import numpy as np
from shore.surface.io import load_stl
from shore.surface.projector import Projector
from shore.surface.anchor import resolve_anchor

mesh = load_stl("body.stl")
anchor = resolve_anchor(mesh)

proj = Projector(mesh, anchor)               # builds BVH once
dirs = np.random.randn(1000, 3)              # 1000 random directions
dirs /= np.linalg.norm(dirs, axis=-1, keepdims=True)

result = proj.project(dirs, on_failure="skip")
print(result.all_hit)                        # likely False for random dirs
print(result.hit_mask.sum(), "out of", len(result.status))
```

## See also

- [`shore.surface.io`](./api-surface-io) — `load_stl` for the
  watertight-checked input.
- [`shore.surface.anchor`](#) — `resolve_anchor` (centroid-based
  anchor selection inside the body).
- [`shore.surface.cap`, `shore.surface.equator`](./api-cap-equator) —
  the cubed-sphere k=0 builders that consume `Projector.project`.
- [`shore.surface.remesh`](./api-surface-remesh) — the legacy
  single-block lat-lon path (also drives `Projector.project`).
- [Surface projection](/algorithm/projection) — the algorithm.
