---
title: shore.geometry.occt
outline: [2, 3]
---

# `shore.geometry.occt`

OCCT backend of the `Geometry` protocol (ADR 0002).

Exact B-rep geometry from STEP, BREP or IGES via OpenCASCADE (OCP, installed
with build123d: `pip install 'shore-mesh[cad]'`). Entities are the shape's
faces (`SURFACE`), edges (`CURVE`) and vertices (`POINT`), indexed in
`TopExp.MapShapes` order; indices are stable for one load only (persistent
identity is ADR 0003).

Projection has two tiers:

- `exact=False` (tier 1, vectorised): closest point on a tessellation that
  is built once at load and keeps (u, v) and the oriented surface normal at
  every node; parameters and normals are interpolated barycentrically.
  Closest points come from trimesh on the tessellation (R-trees built once at
  load, ~2.5e4-1.3e5 points/s). BRepMesh enforces its deflection on edges;
  triangle interiors sag up to ~2x further. For a query at distance `h`
  from the surface the error is about ``2 \* deflection + h \* (facet normal
  error)``: a couple of deflections for points near the surface (smoothing,
  dragging), more for distant points, which should use `exact=True`.
  Where a parameter is singular (a sphere pole, a cone apex: the surface
  derivative along it vanishes) the arbitrary value at the singular node is
  left out of that parameter's interpolation.
- `exact=True` (tier 2, per point): Newton projection onto the exact face
  surface, seeded from tier 1, with a trimming check; points whose foot lies
  outside the trimmed face fall back to an exact shape-distance query, so the
  answer stays on the face (or its boundary). Curves use a clamped Newton
  iteration on the exact edge curve.

`ray_cast` answers on the tessellation (tier 1), through the same
`Projector` the STL backend uses.

These calls run in the caller's thread and hold the GIL; a GUI must run long
loads in a subprocess (ADR 0006).

## `OcctGeometry` <Badge type="tip" text="class" /> {#OcctGeometry}

```python
OcctGeometry(
    shape: TopoDS_Shape,
    *,
    deflection: float | None = None,
    rel_deflection: float = 0.001,
    angular_deflection: float = 0.2,
    heal: bool = True,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L119">source</a>

B-rep body as a `Geometry`.

**Parameters**

- **`shape`** — OCCT shape (`TopoDS_Shape`), e.g. `build123d_object.wrapped`.
- **`deflection`** — Absolute linear deflection of the tier-1 tessellation. Default: `rel_deflection` times the bounding-box diagonal.
- **`rel_deflection`** — Linear deflection relative to the bounding-box diagonal.
- **`angular_deflection`** — Angular deflection of the tessellation, radians.
- **`heal`** — Run `ShapeFix_Shape` on the input first.

### `OcctGeometry.from_file` <Badge type="info" text="classmethod" /> {#OcctGeometry-from_file}

```python
from_file(path: str | Path, **kwargs: object) -> OcctGeometry
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L178">source</a>

Load `.step`/`.stp` (build123d), `.brep`/`.brp` (build123d)
or `.iges`/`.igs` (OCP `IGESControl_Reader`).

### `OcctGeometry.bounds` <Badge type="info" text="property" /> {#OcctGeometry-bounds}

```python
bounds: np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L283">source</a>

### `OcctGeometry.display` <Badge type="info" text="method" /> {#OcctGeometry-display}

```python
display() -> Display
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L286">source</a>

The model as drawn (ADR 0038): the tier-1 tessellation, each
triangle owned by its face; edges as polylines; vertices.

### `OcctGeometry.entities` <Badge type="info" text="method" /> {#OcctGeometry-entities}

```python
entities(kind: EntityKind) -> tuple[EntityRef, ...]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L299">source</a>

### `OcctGeometry.ray_cast` <Badge type="info" text="method" /> {#OcctGeometry-ray_cast}

```python
ray_cast(
    anchor: np.ndarray,
    dirs: np.ndarray,
    on_failure: str = 'raise',
) -> ProjectionResult
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L308">source</a>

### `OcctGeometry.project` <Badge type="info" text="method" /> {#OcctGeometry-project}

```python
project(
    points: np.ndarray,
    entity: EntityRef | None = None,
    *,
    exact: bool = False,
) -> Projection
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L313">source</a>

### `OcctGeometry.evaluate` <Badge type="info" text="method" /> {#OcctGeometry-evaluate}

```python
evaluate(entity: EntityRef, params: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L337">source</a>

### `OcctGeometry.fingerprint` <Badge type="info" text="method" /> {#OcctGeometry-fingerprint}

```python
fingerprint(entity: EntityRef) -> EntityKey
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L355">source</a>

Fingerprint: OCCT geometric type, centroid, area/length, bbox (cached).

### `OcctGeometry.proxy` <Badge type="info" text="property" /> {#OcctGeometry-proxy}

```python
proxy: trimesh.Trimesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L396">source</a>

The tier-1 tessellation (faces not merged across edges).

### `OcctGeometry.merged_proxy` <Badge type="info" text="method" /> {#OcctGeometry-merged_proxy}

```python
merged_proxy() -> trimesh.Trimesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/occt.py#L401">source</a>

The tessellation with coincident vertices merged and degenerate
triangles dropped: watertight for a closed solid.

For mesh-based consumers that need a closed triangulation
(inside/outside tests, plane sections, dihedral-angle TE detection).
Its vertex normals are smoothed across face boundaries, so do not
read normals from it: use `project`.
