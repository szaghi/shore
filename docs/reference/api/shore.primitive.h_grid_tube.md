---
title: shore.primitive.h_grid_tube
outline: [2, 3]
---

# `shore.primitive.h_grid_tube`

Circular-bore H-grid primitive (ADR 0060).

A straight circular tube filled by five hexahedral blocks in an O-H
("butterfly") cross-section:

- `core`: a square in the middle, its corners at `core_frac` of the
  radius, a Cartesian grid;
- `ring_0` .. `ring_3`: four blocks joining the core's sides to the
  outer wall, each a quarter of the circle. Their lines run straight from
  the core's side to the wall, along the 45 deg diagonals at the corners.

No singular axis: every cell is a proper hexahedron, and every block is
right-handed whatever the streamwise axis. Faces:

- `i_lo` / `i_hi` (streamwise inflow / outflow) of every block: `FREE`;
- the rings' `k_hi` (the outer wall): `WALL`;
- every face between two blocks: `SHARED`, with its axis map.

Before ADR 0060 the tube was one block on a polar grid: its `k = 0` row
collapsed onto the centreline, which was marked `WALL`; its first ring
of cells was degenerate and the block was left-handed for the x and z
axes.

## `HGridTubeMesh` <Badge type="tip" text="class" /> {#HGridTubeMesh}

```python
HGridTubeMesh(blocks: list[HexBlock])
```
Bases: `PrimitiveMesh`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/h_grid_tube.py#L50">source</a>

Circular-bore H-grid: a straight tube, O-H cross-section, five blocks.

Node arrays: `ring_q` is `(n_axial, n_circ // 4 + 1, n_radial, 3)`,
i streamwise, j along its quarter of the wall, k from the core's side
(k = 0) to the wall (k = n_radial - 1); `core` is ``(n_axial,
n_circ // 4 + 1, n_circ // 4 + 1, 3)`. `ring_q`` covers the angles
-45 + 90 q to 45 + 90 q degrees in the cross-section's plane.

**Examples**

```python
>>> # 1m-long tube of radius 0.05m with 60 cells round the wall.
>>> ht = HGridTubeMesh.build(
...     length=1.0, radius=0.05,
...     n_axial=80, n_circ=60, n_radial=12,
...     streamwise_axis="z",
...     k_law="tanh", k_beta=3.0,  # cluster the ring's nodes at the wall
... )
>>> [b.label for b in ht.blocks]
# ['ring_0', 'ring_1', 'ring_2', 'ring_3', 'core']
>>> ht.blocks[0].nodes.shape
# (80, 16, 12, 3)
```

### `HGridTubeMesh.build` <Badge type="info" text="classmethod" /> {#HGridTubeMesh-build}

```python
build(
    length: float,
    radius: float,
    n_axial: int,
    n_circ: int,
    n_radial: int,
    streamwise_axis: str = 'x',
    i_law: str = 'uniform',
    k_law: str = 'uniform',
    i_beta: float = 3.0,
    k_beta: float = 3.0,
    label: str | None = None,
    core_frac: float = 0.5,
) -> HGridTubeMesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/h_grid_tube.py#L81">source</a>

Build a circular-bore H-grid.

**Parameters**

- **`length`** — Tube length along the streamwise axis.  Must be positive.
- **`radius`** — Tube outer radius.  Must be positive.
- **`n_axial`** — Streamwise node count, `>= 2`.
- **`n_circ`** — Cells round the wall, a multiple of 4 and `>= 4`: each ring block takes a quarter, and the core has `n_circ // 4` cells a side.
- **`n_radial`** — The ring blocks' node count from the core to the wall, `>= 2`.
- **`streamwise_axis`** — World axis the tube is aligned with: `"x"` (default), `"y"`, or `"z"`.  The cross-section spans the other two.
- **`i_law, k_law`** — Spacing laws (`"uniform"` / `"tanh"` / `"tanh2"`) along the stream (`i`) and across the rings, from the core to the wall (`k`). `"tanh"` clusters at the start of the stream, and at the wall across the rings.
- **`i_beta, k_beta`** — Clustering strengths for tanh / tanh2.
- **`label`** — Ignored since ADR 0060 (a `DeprecationWarning` when given): the blocks are `TUBE_BLOCKS`.
- **`core_frac`** — Radius of the core's corners as a fraction of the tube's, in (0, 1). Default `0.5`.

**Returns**

- `HGridTubeMesh` — Five right-handed blocks (`TUBE_BLOCKS`).

**Raises**

- `ParameterError` — On invalid dimensions or node counts.
