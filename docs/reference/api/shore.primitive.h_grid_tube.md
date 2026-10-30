---
title: shore.primitive.h_grid_tube
outline: [2, 3]
---

# `shore.primitive.h_grid_tube`

Circular-bore H-grid primitive.

A single hex block filling a straight circular tube, with default BCs
that match the typical internal-flow / duct configuration:

- `i_lo` / `i_hi` (streamwise inflow / outflow) → `FREE`.
- `j_lo` / `j_hi` (circumferential) → SHARED via `periodic_j`.
- `k_lo` (centerline singular axis) → `WALL` (degenerate).
- `k_hi` (outer wall) → `WALL`.

Cross-section is a 2D polar grid (lat-lon style): `j` is the
circumferential angle and `k` is the radial coordinate.  The
`k = 0` row is the **polar singular axis** at the duct centerline —
all `n_circ` angular nodes collapse to a single point.  This is the
analogue of the polar axis row on
:class:`shore.mesh.OGridTopology` and is acceptable because Xall's
sphere-grid code handles polar singularities cleanly.

For non-trivial duct geometries (variable cross-section, non-circular
bores), see the planned STL-driven H-grid topology
(`shore mesh --topology hgrid`).

## `HGridTubeMesh` <Badge type="tip" text="class" /> {#HGridTubeMesh}

```python
HGridTubeMesh(blocks: list[HexBlock])
```
Bases: `PrimitiveMesh`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/h_grid_tube.py#L43">source</a>

Single-block circular-bore H-grid (straight tube, no internal body).

Node array shape: `(n_axial, n_circ, n_radial, 3)`.  The block's
`i` axis is streamwise, `j` is the circumferential angle
(periodic by construction), and `k` is the radial coordinate
from the centerline (`k = 0`) to the outer wall (`k = n_radial - 1`).

**Examples**

```python
>>> # 1m-long tube of radius 0.05m with 60 circumferential cells.
>>> ht = HGridTubeMesh.build(
...     length=1.0, radius=0.05,
...     n_axial=80, n_circ=60, n_radial=12,
...     streamwise_axis="z",
...     k_law="tanh", k_beta=3.0,  # cluster radial nodes near the wall
... )
>>> ht.blocks[0].nodes.shape
# (80, 60, 12, 3)
>>> ht.blocks[0].periodic_j
# True
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
    label: str = 'h_grid_tube',
) -> HGridTubeMesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/h_grid_tube.py#L73">source</a>

Build a circular-bore H-grid.

**Parameters**

- **`length`** — Tube length along the streamwise axis.  Must be positive.
- **`radius`** — Tube outer radius.  Must be positive.
- **`n_axial, n_circ, n_radial`** — Node counts along the streamwise, circumferential, and radial directions.  `n_axial >= 2`, `n_circ >= 4`, `n_radial >= 2`.
- **`streamwise_axis`** — World axis the tube is aligned with: `"x"` (default), `"y"`, or `"z"`.  The cross-section spans the other two axes.
- **`i_law, k_law`** — Spacing laws for the streamwise (`i`) and radial (`k`) axes (`"uniform"` / `"tanh"` / `"tanh2"`). Circumferential (`j`) is always uniform — clustering on a periodic axis is unsupported.
- **`i_beta, k_beta`** — Clustering strengths for tanh / tanh2.
- **`label`** — Block label.

**Returns**

- `HGridTubeMesh` — Single-block primitive with default face BCs: `i_lo` / `i_hi` = `FREE` (inflow / outflow), `j` direction = SHARED (periodic), `k_lo` = `WALL` (centerline singular axis), `k_hi` = `WALL` (outer wall).

**Raises**

- `ParameterError` — On invalid dimensions or node counts.
