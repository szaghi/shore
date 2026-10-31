---
title: shore.primitive.box
outline: [2, 3]
---

# `shore.primitive.box`

Cartesian box primitive — a single axis-aligned hexahedral block.

## `BoxExtent` <Badge type="tip" text="dataclass" /> {#BoxExtent}

```python
BoxExtent(vmin: tuple[float, float, float], vmax: tuple[float, float, float])
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/box.py#L22">source</a>

Min/max corners of an axis-aligned box.

The min corner is the (i=0, j=0, k=0) node; the max corner is the
(ni-1, nj-1, nk-1) node.  All three components of *vmax* must
strictly exceed *vmin*.

## `BoxMesh` <Badge type="tip" text="class" /> {#BoxMesh}

```python
BoxMesh(blocks: list[HexBlock])
```
Bases: `PrimitiveMesh`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/box.py#L56">source</a>

Single-block axis-aligned Cartesian box.

Node array shape: `(ni, nj, nk, 3)`.  Each axis has its own
spacing law (uniform / tanh / tanh2) by default.  Per-edge
overrides via `edge_pins` switch the construction to transfinite
interpolation between the 4 parallel edges of the affected axis.
No projection, no marching.

**Examples**

Tensor-product (legacy) construction with global per-axis laws:

```python
>>> box = BoxMesh.build(
...     vmin=(0.0, 0.0, 0.0),
...     vmax=(10.0, 5.0, 3.0),
...     ni=40, nj=30, nk=20,
... )
>>> box.blocks[0].nodes.shape
# (40, 30, 20, 3)
```

Per-edge override (one i-edge gets tanh clustering, others stay
uniform — interior nodes blend the four edge distributions):

```python
>>> from shore.mesh.spacing import Spacing1D
>>> box = BoxMesh.build(
...     vmin=(0.0, 0.0, 0.0), vmax=(10.0, 5.0, 3.0),
...     ni=40, nj=30, nk=20,
...     edge_pins={"j_lo_k_lo": Spacing1D(law="tanh", beta=4.0)},
... )
```

### `BoxMesh.build` <Badge type="info" text="classmethod" /> {#BoxMesh-build}

```python
build(
    vmin: tuple[float, float, float],
    vmax: tuple[float, float, float],
    ni: int,
    nj: int,
    nk: int,
    i_law: str = 'uniform',
    j_law: str = 'uniform',
    k_law: str = 'uniform',
    i_beta: float = 3.0,
    j_beta: float = 3.0,
    k_beta: float = 3.0,
    edge_pins: dict[str, Spacing1D] | None = None,
    label: str = 'box',
) -> BoxMesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/box.py#L94">source</a>

Build a box primitive.

**Parameters**

- **`vmin, vmax`** — 3-tuples giving the min and max corners (world coordinates). Each component of *vmax* must strictly exceed the corresponding component of *vmin*.
- **`ni, nj, nk`** — Node counts along the three axes.  All must be &gt;= 2.
- **`i_law, j_law, k_law`** — Default spacing laws for each axis (`"uniform"`, `"tanh"`, `"tanh2"`).  Used for any of the 4 parallel edges of an axis that aren't explicitly overridden in *edge_pins*.
- **`i_beta, j_beta, k_beta`** — Clustering strengths for tanh / tanh2.
- **`edge_pins`** — Optional sparse override mapping edge name → `Spacing1D`. Edge names are the canonical 12 in :mod:`shore.mesh.edge` (e.g. `"j_lo_k_lo"` for the i-axis edge at j=0, k=0). When any edge is pinned, construction switches to transfinite interpolation between the 4 parallel edges of that edge's axis.  Edges from other axes that aren't mentioned use the global `*_law` default.
- **`label`** — Block label written into the `.geo` filename when used in multi-block contexts.

**Returns**

- `BoxMesh`
