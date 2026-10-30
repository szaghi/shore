---
title: shore.mesh.edge
outline: [2, 3]
---

# `shore.mesh.edge`

Edge: arc descriptor for a structured block edge.

## `Edge` <Badge type="tip" text="dataclass" /> {#Edge}

```python
Edge(
    index: int,
    name: EdgeName,
    axis: Literal['i', 'j', 'k'] | None,
    vertex_indices: tuple[int, int],
    n_nodes: int,
    spacing: Spacing1D,
    association: Association | None = None,
    first_step: float | None = None,
    last_step: float | None = None,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/edge.py#L22">source</a>

One edge of a structured block.

**Parameters**

- **`index`** — 0-based index within the block's edge list.
- **`name`** — Canonical name encoding the two constant axes, e.g. `"j_lo_k_lo"`.
- **`axis`** — The free axis along which this edge runs (`'i'`, `'j'`, or `'k'`), or `None` for non-axis-aligned edges (future use).
- **`vertex_indices`** — Indices of the two endpoint vertices in the block's vertex list.
- **`n_nodes`** — Number of nodes along this edge.
- **`spacing`** — 1D spacing law for node distribution along this edge.
- **`association`** — Geometry this edge's nodes are projected onto (ADR 0003), or `None`.
- **`first_step, last_step`** — Optional per-vertex first / last interval lengths.  When set, the topology / primitive builder treats them as a C1 matching request at that endpoint of the edge: the first (or last) cell along the edge will have this physical length, with the interior shaped by `spacing`.  When unset (`None`), the edge is built purely from `spacing` with no endpoint pinning. Setting these does *not* require `spacing.law == 'pinned'`: the consumer is free to read the steps and switch laws as needed.  At corners where multiple incident edges request different first/last steps, the builder applies a precedence rule (cap / body geometry wins) and emits a `UserWarning` naming the overridden edges.

### `Edge.body_projected` <Badge type="info" text="property" /> {#Edge-body_projected}

```python
body_projected: bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/edge.py#L68">source</a>

True when associated with a whole body (legacy flag, read-only).

### `Edge.nodes` <Badge type="info" text="method" /> {#Edge-nodes}

```python
nodes(block: Block) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/edge.py#L72">source</a>

Return a (n_nodes, 3) view of this edge's node coordinates.
