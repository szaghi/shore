# Surface projection (k=0 ray-casting)

The first algorithmic step turns an unstructured triangulated STL into
a **structured array of surface nodes** by ray-casting from a single
anchor point through a set of pre-chosen direction vectors onto the
body. The output is the k=0 layer that all subsequent stages build on.

## Inputs and outputs

| | |
|---|---|
| **Input** | A watertight `trimesh.Trimesh` and a set of unit direction vectors `dirs` of shape `(N, 3)` from a single anchor inside the body. |
| **Output** | A `(N, 3)` array of intersection points on the body surface, one per ray. |
| **Module** | `shore.surface.projector.Projector` |

## Anchor selection

The **anchor** is the common origin of all rays. SHORE places it at the
body centroid by default (`Body.from_stl`), with a `mesh.contains`
check that the chosen point sits inside the closed surface. The body
must be **star-shaped from the anchor**: every ray from the anchor must
hit the body exactly once. This rules out concave geometries with
re-entrant features visible from the centroid; for those, you'd need
to pre-decompose the body into star-shaped pieces.

![Anchor and ray-cast](/figures/projection-raycast.png)

The diagram shows a 2D cross-section: the anchor (red) is inside a
slightly distorted blob; rays fan outward; each ray's first
intersection with the body is a k=0 node (blue).

## Ray-cast machinery

`Projector.__init__` builds a `trimesh.RayMeshIntersector` once and
caches it. Subsequent `Projector.project(dirs)` calls reuse the same
acceleration structure (BVH), so the per-ray cost is `O(log T)` in the
triangle count.

```python
from shore.surface.projector import Projector

proj = Projector(mesh, anchor)        # builds BVH once
result = proj.project(dirs)           # ray-cast, status array per ray
hits   = result.locations             # (N, 3) hit points
status = result.status                # 'hit' / 'miss' / 'tooNear' / ...
```

The `ProjectionResult.status` field flags rays that miss the body, hit
too close to the anchor (likely the anchor itself if it grazes a face),
or hit too far. By default `Projector.project(dirs, on_failure="raise")`
raises `MeshInputError` if any ray is non-`HIT`; pass `"warn"` or
`"none"` to relax.

## Direction vectors

The projector itself is topology-agnostic — it casts whatever
directions you give it. The directions for the **caps topology** come
from two places:

1. The cap k=0 builder (`shore.surface.cap.build_cap_k0`) computes
   directions for a flat tangent-plane diamond at each pole and feeds
   them through the projector. See
   [k=0 mesh construction](k0-construction).
2. The equator builder (`shore.surface.equator.build_equator`) computes
   per-meridian slerp-interpolated directions from south seam to north
   seam and feeds the interior nodes through the projector. The seam
   nodes themselves are bit-exact copies of the cap edges, so they
   skip the projector round-trip.

For the legacy single-block path, `shore.surface.remesh.sphere_latlon_directions`
generates a uniform lat-lon grid of directions that
`Body.project()` ray-casts as a single batch.

## Performance

A single `Projector.project(dirs)` call vectorises over all rays. For a
20k-triangle STL with a 60×60 cap (3600 rays) the projection takes
single-digit milliseconds on a modern CPU. The BVH build is the
dominant cost on first use, ~30 ms for the same triangle count.

## Failure modes

- **Anchor outside the body**: `mesh.contains([anchor])` returns False;
  raises `MeshInputError` at `Body.from_stl`.
- **Body not star-shaped from anchor**: some rays miss, some hit twice
  (only the first hit is recorded). `Projector.project` raises with the
  `MISS` status by default.
- **Mesh not watertight**: rays may pass through gaps in the surface;
  raises with the `TOO_FAR` status.

## See also

- [k=0 mesh construction](k0-construction) — how the projected nodes
  are arranged into the 6-block caps topology
- [`shore.surface.projector`](/reference/api-surface-remesh) (Python API)
- [Caps topology](/topologies/caps) — full pipeline view
