---
title: shore.mesh.association
outline: [2, 3]
---

# `shore.mesh.association`

Block-to-geometry association (ADR 0003).

Every vertex, edge and face of a block may be *associated* with geometry: the
whole body, a set of CAD surfaces, a chain of CAD curves, or a CAD point. The
nodes of an associated entity are projected onto its geometry after any
operation that moves them, with priority vertex -&gt; edge -&gt; face (ICEM's
projection order): face nodes first, then edge nodes, then the corners, so
the lower-dimensional constraint wins where they meet.

Associations name geometry persistently: `geometry` is a key into the
mesh's geometries (`"body"` for the meshed body) and `keys` are
`EntityKey` fingerprints, re-found in the current geometry by
`match_keys`. An association of kind `BODY` has no keys.

## `Association` <Badge type="tip" text="dataclass" /> {#Association}

```python
Association(geometry: str, kind: EntityKind, keys: tuple[EntityKey, ...] = ())
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L35">source</a>

Persistent link from a block entity to geometry.

**Parameters**

- **`geometry`** — Key of the geometry in the mesh's geometry map (`"body"` for the meshed body).
- **`kind`** — `BODY` (whole body), `SURFACE` (surface set), `CURVE` (curve chain) or `POINT`.
- **`keys`** — Fingerprints of the associated entities; empty for `BODY`.

### `Association.body` <Badge type="info" text="classmethod" /> {#Association-body}

```python
body(geometry: str = BODY_GEOMETRY) -> Association
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L63">source</a>

Association with the whole body.

### `Association.of` <Badge type="info" text="classmethod" /> {#Association-of}

```python
of(geometry_name: str, geometry: Geometry, refs: list[EntityRef]) -> Association
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L68">source</a>

Association with entities *refs* of *geometry* (all of one kind).

### `Association.to_json` <Badge type="info" text="method" /> {#Association-to_json}

```python
to_json() -> dict[str, Any]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L77">source</a>

### `Association.from_json` <Badge type="info" text="classmethod" /> {#Association-from_json}

```python
from_json(data: dict[str, Any], *, source: str = '<association>') -> Association
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L85">source</a>

## `AssociationResolver` <Badge type="tip" text="class" /> {#AssociationResolver}

```python
AssociationResolver(geometries: dict[str, Geometry], *, exact: bool = False)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L93">source</a>

Resolves associations to current `EntityRef`s and projects nodes.

**Parameters**

- **`geometries`** — Geometry map, keyed like `Association.geometry`.
- **`exact`** — Project with `exact=True` (tier 2) instead of the tessellation.

### `AssociationResolver.refs` <Badge type="info" text="method" /> {#AssociationResolver-refs}

```python
refs(assoc: Association) -> list[EntityRef]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L110">source</a>

Current entities of *assoc*; raises if any key is lost.

### `AssociationResolver.project` <Badge type="info" text="method" /> {#AssociationResolver-project}

```python
project(assoc: Association, points: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L127">source</a>

Project *points* (..., 3) onto *assoc*; nearest entity of a set wins.

## `snap_block` <Badge type="info" text="function" /> {#snap_block}

```python
snap_block(block: HexBlock, resolver: AssociationResolver) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L149">source</a>

Project every associated node of *block* onto its geometry, in place.

Faces first, then edges, then vertices, so that the lower-dimensional
association wins at shared nodes (vertex -&gt; edge -&gt; face priority).

## `associate_wall` <Badge type="info" text="function" /> {#associate_wall}

```python
associate_wall(block: HexBlock, assoc: Association) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/association.py#L176">source</a>

Associate *block*'s k_lo face, its 4 edges and 4 corners with *assoc*.
