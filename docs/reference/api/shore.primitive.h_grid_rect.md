---
title: shore.primitive.h_grid_rect
outline: [2, 3]
---

# `shore.primitive.h_grid_rect`

Rectangular-duct H-grid primitive.

A single hex block filling a Cartesian box, with default BCs that
match the typical internal-flow / duct configuration:

- `i_lo` / `i_hi` (streamwise inflow / outflow) → `FREE`.
- `j_lo` / `j_hi` / `k_lo` / `k_hi` (the four duct walls) → `WALL`.

Geometrically the mesh is bit-identical to :class:`shore.primitive.BoxMesh`;
the only differences are the BC defaults and the user-facing
`streamwise_axis` parameter.  Internally :func:`HGridRectMesh.build`
delegates to `BoxMesh.build` for the node array, then permutes the
axes so `i` becomes the streamwise direction and overrides the
default face BCs.

For more general duct geometries (variable cross-section, non-trivial
inner shapes), see the planned STL-driven H-grid topology
(`shore mesh --topology hgrid`).

## `HGridRectMesh` <Badge type="tip" text="class" /> {#HGridRectMesh}

```python
HGridRectMesh(blocks: list[HexBlock])
```
Bases: `PrimitiveMesh`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/h_grid_rect.py#L42">source</a>

Single-block rectangular H-grid (duct without an internal body).

Node array shape: `(ni, nj, nk, 3)` where `i` is the streamwise
direction.  All three axes can carry independent spacing laws
(uniform / tanh / tanh2); tanh / tanh2 on the `j` and `k` axes
typically clusters near the duct walls.

**Examples**

```python
>>> # 1m-long square duct, 0.1m x 0.1m cross-section, streamwise = z.
>>> hg = HGridRectMesh.build(
...     length=1.0, width=0.1, height=0.1,
...     ni=80, nj=20, nk=20,
...     streamwise_axis="z",
...     j_law="tanh2", j_beta=4.0,
...     k_law="tanh2", k_beta=4.0,
... )
>>> hg.blocks[0].nodes.shape
# (80, 20, 20, 3)
>>> hg.blocks[0].face("k_lo").bc.value
# 'wall'
```

### `HGridRectMesh.build` <Badge type="info" text="classmethod" /> {#HGridRectMesh-build}

```python
build(
    length: float,
    width: float,
    height: float,
    ni: int,
    nj: int,
    nk: int,
    streamwise_axis: str = 'x',
    i_law: str = 'uniform',
    j_law: str = 'uniform',
    k_law: str = 'uniform',
    i_beta: float = 3.0,
    j_beta: float = 3.0,
    k_beta: float = 3.0,
    edge_pins: dict[str, Spacing1D] | None = None,
    label: str = 'h_grid_rect',
) -> HGridRectMesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/h_grid_rect.py#L73">source</a>

Build a rectangular-duct H-grid.

**Parameters**

- **`length`** — Duct length along the streamwise axis.  Must be positive.
- **`width, height`** — Cross-section dimensions, perpendicular to the streamwise axis.  Both must be positive.  `width` is along the block's `j` axis; `height` is along `k`.
- **`ni, nj, nk`** — Node counts along the (streamwise, j, k) directions of the block.  All must be `>= 2`.
- **`streamwise_axis`** — World axis the duct is aligned with: `"x"`, `"y"`, or `"z"`.  The cross-section spans the other two axes. Default `"x"`.
- **`i_law, j_law, k_law`** — Spacing laws for the three block axes (`"uniform"` / `"tanh"` / `"tanh2"`).  Defaults uniform.  Tanh / tanh2 on `j` and `k` cluster near the duct walls.
- **`i_beta, j_beta, k_beta`** — Clustering strengths for tanh / tanh2.
- **`edge_pins`** — Optional sparse per-edge spacing override; see :func:`shore.primitive.box.BoxMesh.build`.
- **`label`** — Block label written into the `.geo` filename when used in multi-block contexts.

**Returns**

- `HGridRectMesh` — Single-block primitive with default face BCs: `i_lo` / `i_hi` = `FREE` (inflow / outflow), `j_lo` / `j_hi` / `k_lo` / `k_hi` = `WALL` (duct walls).

**Raises**

- `ParameterError` — On invalid dimensions, node counts, or streamwise axis.
