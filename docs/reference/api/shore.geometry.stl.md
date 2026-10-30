---
title: shore.geometry.stl
outline: [2, 3]
---

# `shore.geometry.stl`

STL backend of the `Geometry` protocol (ADR 0002).

Wraps a `trimesh.Trimesh`. The body is the only entity; there is no
parametrisation. Ray casts go through the existing `Projector` (one cached
BVH per anchor), so points, distances and ray casts are bit-identical to the
pre-protocol code.

Normals are crease-aware (ADR 0014): each triangle corner gets the
corner-angle-weighted mean of the normals of the faces around its vertex
that lie within `crease_deg` of the triangle's own, and the corners of an
open boundary that lies in a plane (a half model's cut) get their normal
projected into that plane. They differ from the legacy
`volume.hyperbolic._stl_normals` wherever faces meet at a sharp edge.

## `StlGeometry` <Badge type="tip" text="class" /> {#StlGeometry}

```python
StlGeometry(mesh: trimesh.Trimesh, crease_deg: float = 30.0)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L36">source</a>

Triangulated body surface as a `Geometry`.

**Parameters**

- **`mesh`** — Triangulated body surface; `load_stl` (used by `from_file`) cleans it so that face orientation and normals are consistent.
- **`crease_deg`** — Faces whose normals differ by more than this angle are not averaged into each other's corner normals (sharp edges). 180 on a closed mesh gives trimesh's vertex normals.

**Raises**

- `ParameterError` — *crease_deg* outside (0, 180].

### `StlGeometry.from_file` <Badge type="info" text="classmethod" /> {#StlGeometry-from_file}

```python
from_file(
    path: str | Path,
    skip_hygiene: bool = False,
    crease_deg: float = 30.0,
) -> StlGeometry
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L69">source</a>

Load an STL file (see `shore.surface.io.load_stl`).

### `StlGeometry.bounds` <Badge type="info" text="property" /> {#StlGeometry-bounds}

```python
bounds: np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L78">source</a>

### `StlGeometry.display` <Badge type="info" text="method" /> {#StlGeometry-display}

```python
display() -> Display
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L81">source</a>

The mesh as drawn (ADR 0038): every triangle belongs to the body;
the feature curves and their junctions (ADR 0049) drawn over it.

### `StlGeometry.entities` <Badge type="info" text="method" /> {#StlGeometry-entities}

```python
entities(kind: EntityKind) -> tuple[EntityRef, ...]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L95">source</a>

### `StlGeometry.features` <Badge type="info" text="method" /> {#StlGeometry-features}

```python
features() -> tuple[list[np.ndarray], np.ndarray]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L105">source</a>

Feature curves and points (ADR 0049), computed once.

Feature edges are the mesh edges whose two facets' normals differ by
more than `crease_deg`, and the open-boundary edges. They are
chained into polylines (the CURVE entities) split at every vertex
where other than two feature edges meet; those vertices are the POINT
entities. Chains, and the points, are ordered by their lowest vertex
index, closed loops starting there.

### `StlGeometry.ray_cast` <Badge type="info" text="method" /> {#StlGeometry-ray_cast}

```python
ray_cast(
    anchor: np.ndarray,
    dirs: np.ndarray,
    on_failure: str = 'raise',
) -> ProjectionResult
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L119">source</a>

### `StlGeometry.project` <Badge type="info" text="method" /> {#StlGeometry-project}

```python
project(
    points: np.ndarray,
    entity: EntityRef | None = None,
    *,
    exact: bool = False,
) -> Projection
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L124">source</a>

Closest point on the triangulation.

The facets *are* the geometry of an STL, so the answer is exact for
either value of *exact*. Normals are barycentric blends of the
crease-aware corner normals of the closest triangle (ADR 0014); on a
crease edge the closest triangle, and so the normal, is either
side's.

### `StlGeometry.evaluate` <Badge type="info" text="method" /> {#StlGeometry-evaluate}

```python
evaluate(entity: EntityRef, params: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L165">source</a>

### `StlGeometry.fingerprint` <Badge type="info" text="method" /> {#StlGeometry-fingerprint}

```python
fingerprint(entity: EntityRef) -> EntityKey
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L168">source</a>

Fingerprint of the body (area-weighted facet centroid, area, bbox),
a feature curve (length-weighted centroid, length, bbox) or a point.

### `StlGeometry.projector` <Badge type="info" text="method" /> {#StlGeometry-projector}

```python
projector(anchor: np.ndarray) -> Projector
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/stl.py#L214">source</a>

The cached `Projector` for *anchor* (one BVH build per anchor).
