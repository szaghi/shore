---
title: Geometry
---

# Geometry

The mesh core never talks to trimesh or OpenCASCADE directly. It talks to a
`Geometry`: a small protocol of batched numpy queries, implemented by two
backends, one for triangulated STL and one for exact CAD (STEP, IGES,
BREP) ([ADR 0002](/adr/0002-geometry-backend)). Everything that touches the
body (k=0 ray casting, association snaps, surface marching, intersection
tracing, the GUI's drawing and picking) goes through it.

## The protocol

`shore.geometry.base.Geometry` asks five things of a backend:

| Member | Answer |
|---|---|
| `bounds` | axis-aligned box, shape `(2, 3)` |
| `entities(kind)` | every entity of a kind, as `EntityRef(kind, index)` |
| `ray_cast(anchor, dirs)` | hits of rays from one interior anchor (the star-shaped templates; see [Surface projection](/concepts/projection)) |
| `project(points, entity=None, *, exact=False)` | closest points, distances, outward normals and, where the backend has them, parameters |
| `evaluate(entity, params)` | positions at `(u, v)` on a surface or `t` on a curve |
| `fingerprint(entity)` | a persistent `EntityKey` (below) |

All arrays are batched and keep their leading dimensions: projecting a
`(ni, nj, 3)` face returns `(ni, nj, 3)` points. `entity=None` means the
whole body. A `Projection` holds `points`, `distances`, `normals` (NaN on a
curve or point), `params` (`(..., 2)` on a surface, `(...)` on a curve,
`None` on STL and for the body) and, for a body projection on a CAD model,
`entity_index`: the face each point landed on.

`load_geometry(path)` picks the backend from the suffix: `.stl` gives
`StlGeometry`; `.step`/`.stp`/`.iges`/`.igs`/`.brep`/`.brp` give
`OcctGeometry`, which needs the `cad` extra (`pip install 'shore-mesh[cad]'`;
without it the load raises `MeshInputError` naming the extra). In a project,
`p.load_geometry(path, name=...)` does the same under a name.

## Entity kinds

| Kind | STL | CAD |
|---|---|---|
| `BODY` | the triangulation | the whole shape |
| `SURFACE` | none | each B-rep face |
| `CURVE` | each feature curve (polyline) | each B-rep edge |
| `POINT` | each feature point | each B-rep vertex |

An `EntityRef` index is valid only for the geometry instance that issued it.
Anything that must survive a reload stores a fingerprint instead.

## STL geometry

`StlGeometry(mesh, crease_deg=30.0)` wraps a `trimesh.Trimesh`. The facets
*are* the geometry, so `project` is exact for either value of `exact`, and
there are no parameters (`evaluate` is refused).

**Crease-aware normals** ([ADR 0014](/adr/0014-crease-aware-stl-normals)).
Each triangle corner gets the corner-angle-weighted mean of the normals of
the faces around its vertex that lie within `crease_deg` of the triangle's
own normal (Thürrner–Wüthrich weighting, restricted to one side of every
sharper edge). A projected point's normal is the barycentric blend of its
triangle's corner normals. Flat faces bounded by sharp edges get their exact
normal; curved regions stay smooth across facets. Corners on an open boundary
that lies in a plane (the cut of a half model) have their normal projected
into that plane, so a half model's normals mirror the whole model's. On a
box the old vertex-normal blend was up to 54.5° off; the crease-aware normal
is exact.

These normals answer "what is the surface normal here", which is what
surface marching, intersection tracing and collars need. `Body.normals`
deliberately keeps trimesh's plain vertex-normal blend for STL bodies
([ADR 0015](/adr/0015-stl-normals-for-marching-blend)); since the removal of
normal blending in the marcher ([ADR 0034](/adr/0034-remove-normal-blending))
it is used only as the reference of the wall-angle metric.

**Feature curves and points** ([ADR 0049](/adr/0049-stl-feature-curves)).
A mesh edge is a feature edge if its two facets' normals differ by more than
`crease_deg`, or if it is an open-boundary edge. Feature edges are chained
into polylines (the `CURVE` entities), split at every vertex where other
than two of them meet (junctions and open ends: the `POINT` entities).
Curves are ordered by their lowest vertex index, so the order is
deterministic. A smooth sphere has no curves; a box has 12 curves and 8
points; a blunt-trailing-edge wing has its TE corner lines, root and tip
outlines and the TE corners. Projection onto a curve is the closest point of
the polyline; onto a point, the point.

<Figure name="stl-features" caption="Feature curves of an STL, coloured per curve, and its feature points." />

```python
import numpy as np

from shore.geometry import EntityKind, StlGeometry
from shore.geometry.keys import match_keys
from shore.shapes import box

g = StlGeometry(box((2.0, 1.0, 0.5)))  # crease_deg=30 by default
curves = g.entities(EntityKind.CURVE)
points = g.entities(EntityKind.POINT)
print(len(curves), len(points))  # 12 8: the box edges and corners

proj = g.project(np.array([[0.3, 0.9, 0.9]]), curves[0])
print(proj.normals)  # [[nan nan nan]]: curves have no normal

key = g.fingerprint(curves[0])  # kind, "polyline", centroid, length, extent
stretched = StlGeometry(box((2.1, 1.0, 0.5)))  # an edited model
print(match_keys([key], stretched).summary())  # 1 keys: 0 same, 1 changed, 0 lost
```

Each query point's closest point is computed once and remembered by its
exact bytes (up to 200 000 points per geometry), which returns bit-identical
answers and roughly halves the time of iterative tracers
([ADR 0032](/adr/0032-stl-projection-cache)).

## CAD geometry

`OcctGeometry` reads STEP and BREP through build123d and IGES through OCP's
`IGESControl_Reader`; IGES faces are then sewn (tolerance $10^{-6}$), because
IGES stores faces as independent trimmed surfaces and association needs
shared edges. The shape is healed with `ShapeFix_Shape` on load. Faces,
edges and vertices are indexed in `TopExp.MapShapes` order.

At load every face is tessellated once (`BRepMesh_IncrementalMesh`, linear
deflection $10^{-3}$ × the bounding-box diagonal by default, angular
deflection 0.2 rad), keeping `(u, v)` and the oriented normal at every
tessellation node; every edge is sampled to a polyline. Projection then has
two tiers:

| | `exact=False` (tier 1) | `exact=True` (tier 2) |
|---|---|---|
| Method | closest point on the tessellation (trimesh, R-trees built at load); `(u, v)` and normals interpolated barycentrically | per-point OCCT Newton projection seeded from tier 1, with a trimming check; a foot outside the trimmed face falls back to an exact shape-distance query and lands on the face boundary; curves by clamped Newton |
| Accuracy | near the surface ≤ about 2 × deflection; grows with distance (0.045 at distance 1 on the unit sphere) | exact: sphere radius error $2 \cdot 10^{-16}$ |
| Throughput (measured) | $2.5 \cdot 10^4$–$1.3 \cdot 10^5$ points/s | $1.4 \cdot 10^3$–$4 \cdot 10^4$ points/s |
| Used for | smoothing loops, interactive picks, intermediate passes | the final snap of each surface algorithm; distant queries |

At singular parameters (a sphere pole, a cone apex) the arbitrary parameter
value at the singular node is left out of the interpolation. `ray_cast`
always answers on the tessellation; for a CAD body the k=0 builders snap the
hits to the exact surface afterwards (`AnchoredCaster(snap=True)`).

Algorithms that need one closed triangulation (anchor resolution, plane
sections and trailing-edge detection of the C-grid) read `merged_proxy`:
the tessellation with coincident vertices merged and degenerate triangles
dropped, watertight for a closed solid. They run at tessellation accuracy.

<Figure name="cad-entities" caption="A CAD wing–body: faces coloured per face, edges drawn as dark curves." />

## Fingerprints and matching

OpenCASCADE gives shapes no identity that survives re-reading an edited file
(the "topological naming problem"). SHORE therefore names entities by what
they look like ([ADR 0003](/adr/0003-association)). An `EntityKey`
(`shore.geometry.keys`) holds:

- `kind` and `geom_type`: the OCCT surface or curve type (`Plane`,
  `Cylinder`, `BSplineSurface`, `Line`, …), `"polyline"` for an STL curve,
  `"Point"`, `"mesh"` / `"brep"` for a body; these must match exactly;
- `centroid`, `measure` (area, length, 0 for a point) and `extent`
  (bounding-box size), normalised by the model's bounding box: coordinates
  as $(x - x_{\min})/d$, areas over $d^2$, lengths over $d$, with $d$ the
  bounding-box diagonal. A uniform rescale leaves keys unchanged.

`match_keys(keys, geometry)` re-finds stored keys in a geometry. Within each
`(kind, geom_type)` group it solves a one-to-one assignment
(`scipy.optimize.linear_sum_assignment`) minimising the summed distance
$|\Delta c| + |\Delta e| + |\Delta m| / \max(m)$, so two similar faces cannot
both claim one entity. Each key is reported:

| Status | Cost |
|---|---|
| `same` | ≤ `tol_same` = $10^{-9}$ |
| `changed` | between the two: re-found after an edit |
| `lost` | > `tol_match` = 0.25, or no candidate of its type |

Because keys are normalised by the model box, an edit that changes the box
(a span stretch) marks every entity `changed`, including ones that did not
move: `changed` means "re-found", not "edited". Large topological edits (a
face split in CAD) lose keys; the association that holds them then fails to
resolve and asks to be re-associated. Measured on a polygonal NACA 0012 wing
(40 faces, 114 edges): every association resolves after a +5 % span edit;
some are lost after +30 %. The names slot exists but STEP labels are not read
yet.

This is what lets a journal survive a CAD edit: replay the same journal on
the edited file and every association re-finds its entity
([Projects and journals](/concepts/projects)).

## Display

`display()` returns how a geometry is drawn and picked in the GUI
([ADR 0038](/adr/0038-gui-geometry-and-slices)): float64 triangles with the
entity owning each (the body of an STL, a face of a CAD model), the curve
polylines and the point positions. For an STL the curves and points are its
feature curves and points; for CAD, the tier-1 tessellation, edge polylines
and vertices.

## Limits

- `ray_cast` templates (cubed sphere, O-grid, OH) need a star-shaped body
  around one anchor; association-driven projection has no such restriction.
- OCCT calls hold the GIL and run one point at a time; tier 2 is for final
  passes, not inner loops.
- Tier 1 is inaccurate far from the surface: use `exact=True` there.
- STL curves come from one crease threshold (the same as the normals'); curves
  from cut planes are not provided.

## See also

- [ADR 0002](/adr/0002-geometry-backend), [ADR 0003](/adr/0003-association),
  [ADR 0014](/adr/0014-crease-aware-stl-normals),
  [ADR 0015](/adr/0015-stl-normals-for-marching-blend),
  [ADR 0032](/adr/0032-stl-projection-cache),
  [ADR 0038](/adr/0038-gui-geometry-and-slices), [ADR 0049](/adr/0049-stl-feature-curves)
- API: [`shore.geometry`](/reference/api/shore.geometry),
  [`shore.geometry.base`](/reference/api/shore.geometry.base#Geometry),
  [`shore.geometry.stl`](/reference/api/shore.geometry.stl#StlGeometry),
  [`shore.geometry.occt`](/reference/api/shore.geometry.occt#OcctGeometry),
  [`shore.geometry.keys`](/reference/api/shore.geometry.keys#match_keys)
- [Load STL and CAD](/cookbook/load-geometry), [Feature curves of an STL](/cookbook/stl-curves)
- [Association](/concepts/association)
