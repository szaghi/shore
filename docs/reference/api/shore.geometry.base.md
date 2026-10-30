---
title: shore.geometry.base
outline: [2, 3]
---

# `shore.geometry.base`

Geometry protocol: the only geometry API the mesh core sees (ADR 0002).

A `Geometry` answers batched numpy queries about a body: anchored ray casts
(the existing templates), closest-point projection onto the body or onto one
of its entities, and parametric evaluation where the backend has parameters.
Backends: `StlGeometry` (triangulated STL, via trimesh) and, later,
`OcctGeometry` (STEP/IGES/BREP, via build123d/OCP).

## `EntityKind` <Badge type="tip" text="class" /> {#EntityKind}

```python
class EntityKind(Enum)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L24">source</a>

Topological kind of a geometry entity.

## `EntityRef` <Badge type="tip" text="dataclass" /> {#EntityRef}

```python
EntityRef(kind: EntityKind, index: int = 0)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L34">source</a>

Backend-local handle to one geometry entity.

`index` is only meaningful within the backend instance that issued the
reference; persistent identity across re-imports is the job of the
association layer (ADR 0003), not of this handle.

## `Display` <Badge type="tip" text="dataclass" /> {#Display}

```python
Display(
    points: np.ndarray,
    triangles: np.ndarray,
    triangle_entity: np.ndarray,
    surface_kind: EntityKind,
    curves: tuple[np.ndarray, ...] = (),
    vertices: np.ndarray = field(default_factory=lambda: np.empty((0, 3))),
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L51">source</a>

A geometry as it is drawn (ADR 0038): triangles, curves and points.

Each triangle belongs to entity ``EntityRef(surface_kind,
triangle_entity[t])``: the body of an STL, a face of a CAD model. Curve
`c` is the polyline of `EntityRef(CURVE, c)`, point `p` the
position of `EntityRef(POINT, p)`. Coordinates are float64.

## `Projection` <Badge type="tip" text="dataclass" /> {#Projection}

```python
Projection(
    points: np.ndarray,
    distances: np.ndarray,
    normals: np.ndarray,
    params: np.ndarray | None = None,
    entity_index: np.ndarray | None = None,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L69">source</a>

Closest-point projection of a batch of points.

**Parameters**

- **`points`** (np.ndarray, shape (..., 3)) — Projected positions.
- **`distances`** (np.ndarray, shape (...)) — Distance from each query point to its projection.
- **`normals`** (np.ndarray, shape (..., 3)) — Unit outward surface normals at the projected points (for a curve entity, NaN).
- **`params`** (np.ndarray or None) — Parameters of the projected points on the target entity, shape (..., 2) for a surface and (...) for a curve; *None* when the backend has no parametrisation (STL), and for body projections (parameters only make sense per entity).
- **`entity_index`** (np.ndarray or None) — For body projections on a multi-entity backend: index of the SURFACE entity each point landed on, shape (...); *None* otherwise.

## `Geometry` <Badge type="tip" text="class" /> {#Geometry}

```python
class Geometry(Protocol)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L99">source</a>

Batched geometry queries used by the mesh core.

All array arguments and results are numpy arrays; leading dimensions of
query arrays are preserved in the results.

### `Geometry.bounds` <Badge type="info" text="property" /> {#Geometry-bounds}

```python
bounds: np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L107">source</a>

Axis-aligned bounding box, shape (2, 3): `[min, max]`.

### `Geometry.entities` <Badge type="info" text="method" /> {#Geometry-entities}

```python
entities(kind: EntityKind) -> tuple[EntityRef, ...]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L111">source</a>

All entities of *kind* (empty when the backend has none).

### `Geometry.ray_cast` <Badge type="info" text="method" /> {#Geometry-ray_cast}

```python
ray_cast(
    anchor: np.ndarray,
    dirs: np.ndarray,
    on_failure: str = 'raise',
) -> ProjectionResult
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L115">source</a>

Cast rays from one interior *anchor* along unit *dirs* (..., 3).

Semantics and per-ray status follow `Projector.project`; this is the
projection the star-shaped templates (cubed-sphere, O-grid) use.

### `Geometry.project` <Badge type="info" text="method" /> {#Geometry-project}

```python
project(
    points: np.ndarray,
    entity: EntityRef | None = None,
    *,
    exact: bool = False,
) -> Projection
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L125">source</a>

Closest-point projection of *points* (..., 3) onto *entity*.

`entity=None` projects onto the whole body. `exact=False` allows a
backend to answer from a tessellated proxy (ADR 0002 tier 1);
`exact=True` requires the exact geometry (tier 2).

### `Geometry.evaluate` <Badge type="info" text="method" /> {#Geometry-evaluate}

```python
evaluate(entity: EntityRef, params: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L140">source</a>

Positions of *entity* at *params*: (..., 2) on a surface, (...) on a curve.

### `Geometry.fingerprint` <Badge type="info" text="method" /> {#Geometry-fingerprint}

```python
fingerprint(entity: EntityRef) -> EntityKey
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L144">source</a>

Persistent geometric fingerprint of *entity* (ADR 0003).

## `AnchoredCaster` <Badge type="tip" text="class" /> {#AnchoredCaster}

```python
AnchoredCaster(geometry: Geometry, anchor: np.ndarray, snap: bool = False)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L149">source</a>

`Projector`-compatible anchored ray caster over any `Geometry`.

Exposes `anchor` and `project(dirs, on_failure)` (the interface the
k=0 builders call: `surface.cap`, `surface.equator`,
`surface.remesh`), backed by `geometry.ray_cast`. With `snap=True`
every HIT location is then moved to the exact closest point on the body
(`geometry.project(..., exact=True)`): for a CAD body the ray cast runs
on the tessellation, and snapping puts k=0 on the exact surface. For an
STL body snapping would be a no-op, so it is left off.

**Parameters**

- **`geometry`** — Body geometry.
- **`anchor`** — Ray origin inside the body, shape (3,).
- **`snap`** — Snap hits to the exact surface.

### `AnchoredCaster.project` <Badge type="info" text="method" /> {#AnchoredCaster-project}

```python
project(dirs: np.ndarray, on_failure: str = 'raise') -> ProjectionResult
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/base.py#L175">source</a>

Cast from the anchor along *dirs* (..., 3); see `Projector.project`.
