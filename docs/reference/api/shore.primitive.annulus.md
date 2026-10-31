---
title: shore.primitive.annulus
outline: [2, 3]
---

# `shore.primitive.annulus`

Annular cylinder primitive — a single block with optional theta-periodicity.

Geometry: an annular cylinder bounded by an inner radius, an outer
radius, and two axial planes.  The structured grid is parametrised by
`(r, theta, z)`; the *theta* direction can be either fully periodic
(`2*pi`) or a partial sector `[theta_start, theta_end]`.

Index convention (i,j,k) → (radial, azimuthal, axial).  This matches
the legacy single-block O-grid: the j direction is the periodic /
sector direction.

## `AnnulusExtent` <Badge type="tip" text="dataclass" /> {#AnnulusExtent}

```python
AnnulusExtent(r_in: float, r_out: float, z0: float, z1: float)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/annulus.py#L34">source</a>

Inner radius, outer radius, axial extent of an annular cylinder.

## `AnnularCylinderMesh` <Badge type="tip" text="class" /> {#AnnularCylinderMesh}

```python
AnnularCylinderMesh(blocks: list[HexBlock])
```
Bases: `PrimitiveMesh`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/annulus.py#L65">source</a>

Single-block annular cylinder.

Index meaning
-------------
- `i` (axis 0): radial direction, `r_in` → `r_out`.
- `j` (axis 1): azimuthal direction, `theta_start` → `theta_end`.
- `k` (axis 2): axial direction along *axis* unit vector.

For a full annulus (`theta_end - theta_start == 2*pi`) the j
direction is periodic and the last azimuthal node coincides with
the first.  For a partial sector j is non-periodic.

The cylinder axis defaults to z.  Pass `axis="x"`, `axis="y"`,
or a 3-vector to rotate the geometry into world frame.

Per-edge spacing
----------------
The `edge_pins` argument lets individual edges carry a different
`Spacing1D` than the global `i_law` / `j_law` / `k_law`.
When any edge of an axis is pinned, that axis's parametric
coordinate becomes a transfinite-interpolation blend of the 4
parallel edges' distributions.

Periodic-j constraint: for a full annulus, the j-direction is
azimuthally periodic and pinning any j-edge (`i_lo_k_lo`,
`i_hi_k_lo`, `i_lo_k_hi`, `i_hi_k_hi`) is rejected — the
"first" and "last" cells of a periodic edge are the same cell, so
per-edge endpoint pinning is meaningless there.  Use a partial
sector if you need azimuthal clustering control.

**Examples**

```python
>>> ann = AnnularCylinderMesh.build(
...     r_in=0.5, r_out=2.0, z0=0.0, z1=4.0,
...     ni=20, nj=64, nk=30,
... )
>>> ann.blocks[0].periodic_j
# True
```

### `AnnularCylinderMesh.build` <Badge type="info" text="classmethod" /> {#AnnularCylinderMesh-build}

```python
build(
    r_in: float,
    r_out: float,
    z0: float,
    z1: float,
    ni: int,
    nj: int,
    nk: int,
    theta_start: float = 0.0,
    theta_end: float | None = None,
    axis: str | tuple[float, float, float] = 'z',
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    i_law: str = 'uniform',
    j_law: str = 'uniform',
    k_law: str = 'uniform',
    i_beta: float = 3.0,
    j_beta: float = 3.0,
    k_beta: float = 3.0,
    edge_pins: dict[str, Spacing1D] | None = None,
    label: str = 'annulus',
) -> AnnularCylinderMesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/annulus.py#L112">source</a>

Build an annular cylinder primitive.

**Parameters**

- **`r_in, r_out`** — Inner and outer radii (`r_out > r_in >= 0`).
- **`z0, z1`** — Axial extent along *axis* (`z1 > z0`).
- **`ni, nj, nk`** — Node counts (radial, azimuthal, axial).  All &gt;= 2.  These are the array shape of the produced `(ni, nj, nk, 3)` block, not cell counts.  `ni = 20` means 20 radial nodes and 19 cells.
- **`theta_start, theta_end`** — Azimuthal sector in radians.  When *theta_end* is None, the cylinder is a full annulus (`theta_start` to `theta_start + 2*pi`) and the j direction is periodic. When given explicitly, the sector is partial and j is non-periodic; `theta_end > theta_start` is required.
- **`axis`** — `"x"`, `"y"`, `"z"` (case-insensitive), or a 3-component vector giving an arbitrary axis direction.
- **`center`** — World-space anchor for the cylinder axis.  The axis line passes through *center* with direction *axis*; `z=0` of the local frame corresponds to *center*.
- **`i_law, j_law, k_law, *beta`** — Per-axis default spacing laws and clustering strengths.
- **`edge_pins`** — Sparse per-edge overrides; see class docstring.

**Returns**

- `AnnularCylinderMesh`
