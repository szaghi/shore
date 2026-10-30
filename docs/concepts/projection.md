# Surface projection

Surface projection turns an unstructured body surface into a **structured
array of surface nodes**. Rays leave one anchor point inside the body
along a structured set of directions; each ray's hit on the body is a k=0
node. The projector knows nothing about blocks: the
[k=0 construction](/concepts/k0-construction) decides which directions to
cast.

<Figure name="projection-raycast" ext="svg" caption="Rays from the anchor inside a star-shaped section; each first hit on the body becomes a k=0 node." />

## Inputs and outputs

| | |
|---|---|
| Input | the body, an anchor point inside it, and unit directions of shape `(N, 3)` |
| Output | `ProjectionResult`: `locations` `(N, 3)` (NaN where a ray failed), `status` per ray, `distances` from the anchor |
| Code | [`shore.surface.projector.Projector`](/reference/api/shore.surface.projector#Projector) (STL); `AnchoredCaster` for CAD bodies |

## Anchor

All rays share one origin, so the body must be **star-shaped from the
anchor**: every ray must cross the surface exactly once.
[`resolve_anchor`](/reference/api/shore.surface.anchor#resolve_anchor)
picks the anchor; `Body.from_stl` / `Body.from_file` use the `auto`
strategy:

| Strategy | Anchor |
|---|---|
| `auto` (default) | `incenter`, else `bbox`, each validated; raises if both fail |
| `incenter` | centre of the minimum enclosing sphere of the vertices |
| `bbox` | bounding-box centre |
| `centroid` | volumetric centroid (can lie outside a non-convex body) |
| `user` | a given point, validated |

Validation is a `mesh.contains` test. A body whose anchor cannot be
resolved still loads: the error is kept and raised when the projector is
first used, so open half-model bodies can still serve collar grids
([ADR 0018](/adr/0018-cli-collars-and-projects)).

## Ray casting

For an STL the `Projector` builds a trimesh BVH (`RayMeshIntersector`) once
per body and reuses it for every call, so a batch costs one vectorised
query. Each hit is classified, in this priority:

| Status | Condition |
|---|---|
| `DEGENERATE` | NaN in the hit |
| `BACKFACE` | the hit triangle's normal does not face along the ray (inverted or wrong-side facet) |
| `TOO_NEAR` | hit closer than 5 % of the body's bounding-box diagonal |
| `TOO_FAR` | hit farther than 5 bounding-box diagonals |
| `MISSED` | no hit |
| `HIT` | everything else |

Only the first hit along each ray is kept. `project(dirs,
on_failure=...)` decides what a non-`HIT` ray does: `"raise"` (default,
`MeshInputError` with the counts per status), `"nearest"` (fall back to the
closest surface point and keep the status flag) or `"skip"` (leave NaN).

For a CAD body the caster ray-casts the merged tessellation and then snaps
every hit to the exact surface, so the k=0 nodes lie on the CAD faces, not
on the facets ([ADR 0002](/adr/0002-geometry-backend)).

## Who supplies the directions

- **Cubed sphere**: `build_cap_k0` casts a gnomonic diamond at each pole,
  then `build_equator` casts slerped great-circle meridians between the two
  cap seam rings; the seam rows are copied from the caps, not re-cast
  ([k=0 construction](/concepts/k0-construction)).
- **O-grid**: `sphere_latlon_directions` builds a lat-lon grid of
  directions over colatitudes $[\theta_{cap},\, \pi - \theta_{cap}]$, cast
  in one batch.
- **OH hybrid**: the cubed-sphere directions with `theta_cap` and the
  azimuth offset both pinned to $\pi/4$.

## Failure modes

| Symptom | Cause |
|---|---|
| `MeshInputError` naming the anchor | anchor outside the body: the body is not closed, or not star-shaped from any of the tried points |
| `MeshInputError` with `missed` / `backface` counts | some rays leave through a gap, or hit an inverted facet; repair the STL |
| several hits along a ray, only the first kept | the body is not star-shaped from the anchor: the grid skips the hidden surface; nothing warns, so inspect the k=0 layer |

For bodies that cannot be star-shaped, see
[Pipeline limits](/concepts/pipeline-limits#star-shaped-bodies).

## See also

- [k=0 construction](/concepts/k0-construction): how the projected nodes become the 6 cubed-sphere blocks
- [Geometry](/concepts/geometry): the STL and CAD backends behind the caster
- [`shore.surface.projector`](/reference/api/shore.surface.projector), [`shore.surface.anchor`](/reference/api/shore.surface.anchor), [`shore.surface.remesh`](/reference/api/shore.surface.remesh)
