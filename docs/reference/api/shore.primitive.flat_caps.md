---
title: shore.primitive.flat_caps
outline: [2, 3]
---

# `shore.primitive.flat_caps`

Flat-caps cylinder primitive — 5-block conforming alternative to a
hollow Chimera background.

Geometry
--------

At every axial slice the disk of radius `r_out` is decomposed into:

- 1 *central square* block, axis-aligned in the (e1, e2) plane and
  centred on the axis.  Its 4 corners sit at radius ``r_inner = inner_frac
  \* r_out` (default `inner_frac = 0.5``) at azimuths
  `pi/4 + k*pi/2` (k = 0..3) — i.e. the "diamond" orientation.
  Side length `s = sqrt(2) * r_inner`.
- 4 *trapezoidal sub-blocks* around the central square.  Each sub-block
  has:
    - inner edge: one side of the central square (a straight line of
      length `s`, both endpoints at radius `r_inner`);
    - outer edge: the matching 90° arc on the outer circle at
      radius `r_out`;
    - two side edges: straight radial segments connecting a corner of
      the central square at radius `r_inner` to the same azimuth on
      the outer circle at radius `r_out`.

The decomposition is then extruded along *axis* with `nk` axial
layers.  This produces 5 blocks total — no hole at the centre, no
overset fringes between blocks (all seams are SHARED structured
boundaries with bit-exact node coincidence).

C1 first-step continuity
------------------------

The radial direction inside each sub-block has its first cell
**pinned** to the central square's first cell (in the corresponding
direction), so that the radial cell-size is C1-continuous across the
central / sub-block seam:

    central_first_step = sqrt(2) \* r_inner / (nc - 1)   (for c_law='uniform')
    sub_first_radial   = central_first_step             (always pinned)

The remaining `ni - 2` cells of the sub-block radial distribution
absorb the difference, distributed by the chosen `i_law`.  When the
pin is incompatible with `i_law`/`ni` (would require a non-monotone
or trivially-empty interior), SHORE falls back to the unpinned
distribution and emits a `UserWarning`.

Use case
--------

When you need a cylindrical background mesh but want to avoid Chimera
fringes (high-order accuracy, conservation, GPU memory access), this
primitive replaces the
`annular cylinder + hole` overset combination with a single
conforming structured assembly.

Block labels
------------

- `center`           — central square block.
- `sub_e`, `sub_n`, `sub_w`, `sub_s` — surrounding
  sub-blocks (east / north / west / south in the (e1, e2) plane).

Each sub-block uses index `i` for the radial direction (square edge
at `i = 0`, circle arc at `i = ni - 1`) and `j` for the
along-edge / along-arc direction.  `k` is the axial direction in
every block.

Index counts
------------

- Center block:    `(nc, nc, nk, 3)`
- Each sub-block:  `(ni, nc, nk, 3)` — `nc` is shared with the
  central block to guarantee SHARED conformity at every seam.

## `FlatCapsExtent` <Badge type="tip" text="dataclass" /> {#FlatCapsExtent}

```python
FlatCapsExtent(r_out: float, z0: float, z1: float)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/flat_caps.py#L97">source</a>

Outer radius and axial extent of a flat-caps cylinder.

## `FlatCapsCylinder` <Badge type="tip" text="class" /> {#FlatCapsCylinder}

```python
FlatCapsCylinder(blocks: list[HexBlock])
```
Bases: `PrimitiveMesh`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/flat_caps.py#L171">source</a>

5-block flat-caps cylinder primitive.

**Examples**

```python
>>> fc = FlatCapsCylinder.build(r_out=2.0, z0=0.0, z1=4.0,
...                             ni=10, nc=12, nk=20)
>>> [b.label for b in fc.blocks]
# ['center', 'sub_e', 'sub_n', 'sub_w', 'sub_s']
```

### `FlatCapsCylinder.build` <Badge type="info" text="classmethod" /> {#FlatCapsCylinder-build}

```python
build(
    r_out: float,
    z0: float,
    z1: float,
    ni: int,
    nc: int,
    nk: int,
    inner_frac: float = 0.5,
    axis: str | tuple[float, float, float] = 'z',
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    i_law: str = 'uniform',
    c_law: str = 'uniform',
    k_law: str = 'uniform',
    i_beta: float = 3.0,
    c_beta: float = 3.0,
    k_beta: float = 3.0,
    edge_pins: dict[str, Spacing1D] | None = None,
) -> FlatCapsCylinder
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/flat_caps.py#L188">source</a>

Build a flat-caps cylinder primitive.

**Parameters**

- **`r_out`** — Outer radius of the disk.
- **`z0, z1`** — Axial extent along *axis* (`z1 > z0`).
- **`ni`** — Radial node count for each sub-block (square edge → outer circle).  Must be &gt;= 2.
- **`nc`** — Shared node count along each side of the central square. Each sub-block uses *nc* nodes along the arc; the central block is `(nc, nc)` in (e1, e2).  Must be &gt;= 2.
- **`nk`** — Axial node count.  Must be &gt;= 2.
- **`inner_frac`** — Ratio `r_inner / r_out` of the central square's corner-radius to the outer disk radius.  Must lie in `(0, 1)`.  Default 0.5.
- **`axis`** — Cylinder axis (`"x"`, `"y"`, `"z"`, or a 3-vector).
- **`center`** — World-space anchor for the cylinder axis.
- **`i_law, c_law, k_law, *beta`** — Spacing laws.  `i` controls the radial sub-block distribution (with first cell pinned to the central first cell); `c` controls the central-square / sub-block-arc distribution (shared between center and sub-blocks); `k` controls the axial distribution.
- **`edge_pins`** — Optional sparse override mapping k-axis edge name → `Spacing1D`.  When any k-edge is pinned the axial coordinate is built by transfinite interpolation between the 4 parallel k-edges of every block, allowing different axial spacing on (e.g.) the inner radius vs the outer circumference at the same k-layer. Pins on i- or c-axis edges are **rejected**: the radial (i) and tangential (c) distributions on this primitive are shaped by the geometric construction (per-j first-step matching at the central-square seam) and overriding them would require rewriting that geometry.

**Returns**

- `FlatCapsCylinder`
