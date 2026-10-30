---
title: shore.mesh.vertex
outline: [2, 3]
---

# `shore.mesh.vertex`

Vertex: corner descriptor for a structured block.

## `Vertex` <Badge type="tip" text="dataclass" /> {#Vertex}

```python
Vertex(index: int, ijk: tuple[int, int, int], association: Association | None = None)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/vertex.py#L16">source</a>

One corner of a structured block.

**Parameters**

- **`index`** — 0-based index within the block's vertex list.
- **`ijk`** — `(i, j, k)` corner indices into `block.nodes`.  Use `-1` for the last index along an axis (standard Python convention).
- **`association`** — Geometry this corner is projected onto (ADR 0003), or `None`.

### `Vertex.body_projected` <Badge type="info" text="property" /> {#Vertex-body_projected}

```python
body_projected: bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/vertex.py#L35">source</a>

True when associated with a whole body (legacy flag, read-only).

### `Vertex.coords` <Badge type="info" text="method" /> {#Vertex-coords}

```python
coords(block: Block) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/vertex.py#L39">source</a>

Return the (3,) coordinate vector for this vertex.
