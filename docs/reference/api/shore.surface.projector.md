---
title: shore.surface.projector
outline: [2, 3]
---

# `shore.surface.projector`

Ray-mesh projector with cached BVH and per-ray status tracking.

Phase 0: Projector class + ProjectionResult decoupling.
Phase 2: Full status array (BACKFACE, TOO_NEAR, TOO_FAR, DEGENERATE) + on_failure.

## `ProjectionStatus` <Badge type="tip" text="class" /> {#ProjectionStatus}

```python
class ProjectionStatus(Enum)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/projector.py#L22">source</a>

Per-ray outcome of a Projector.project call.

## `ProjectionResult` <Badge type="tip" text="dataclass" /> {#ProjectionResult}

```python
ProjectionResult(locations: np.ndarray, status: np.ndarray, distances: np.ndarray)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/projector.py#L34">source</a>

Outcome of a batch ray-cast from a single anchor point.

**Parameters**

- **`locations`** (np.ndarray, shape (n, 3)) — Hit positions for HIT rays; NaN for non-HIT rays.
- **`status`** (np.ndarray, shape (n,), dtype=object) — Per-ray ProjectionStatus value.
- **`distances`** (np.ndarray, shape (n,)) — Anchor-to-hit distance; NaN for non-HIT rays.

### `ProjectionResult.hit_mask` <Badge type="info" text="property" /> {#ProjectionResult-hit_mask}

```python
hit_mask: np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/projector.py#L52">source</a>

Boolean mask where True means HIT.

### `ProjectionResult.all_hit` <Badge type="info" text="property" /> {#ProjectionResult-all_hit}

```python
all_hit: bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/projector.py#L57">source</a>

True iff every ray has status HIT.

## `Projector` <Badge type="tip" text="class" /> {#Projector}

```python
Projector(mesh: trimesh.Trimesh, anchor: np.ndarray)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/projector.py#L62">source</a>

Ray-mesh intersector with a single cached BVH.

Holds a `RayMeshIntersector` built once at construction time so that
multiple `project` calls share the same acceleration structure.

**Parameters**

- **`mesh`** — Triangulated body surface.
- **`anchor`** — Ray origin, shape (3,).

### `Projector.project` <Badge type="info" text="method" /> {#Projector-project}

```python
project(dirs: np.ndarray, on_failure: str = 'raise') -> ProjectionResult
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/projector.py#L83">source</a>

Cast rays from anchor in directions *dirs*.

**Parameters**

- **`dirs`** — Unit direction vectors, shape (..., 3).  Flattened internally.
- **`on_failure`** — `'raise'` — raise MeshInputError if any ray is not HIT. `'nearest'` — fall back to trimesh closest-point for non-HIT rays and flag them in the status array. `'skip'` — leave non-HIT locations as NaN.

**Returns**

- `ProjectionResult`
