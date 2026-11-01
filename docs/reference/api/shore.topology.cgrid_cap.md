---
title: shore.topology.cgrid_cap
outline: [2, 3]
---

# `shore.topology.cgrid_cap`

Wing-tip caps of the C-grid (ADR 0059).

The C-grid wraps a wing's sides and leaves its span ends open. A cap closes
one tip as an overset (Chimera) component: a body-fitted grid on the flat
tip face, extruded outwards.

- **Footprint, in the tip plane:** an O-H layout of five blocks. A
  `core` sits inside the section, inset by `core_frac` of the
  half-thickness from both sides. A ring of four blocks joins it to the
  section's wall: `upper` (the C-grid's i = 0 side), `le` round the
  leading edge between the two LE cuts, `lower`, and `te` across the
  trailing-edge base (a blunt or rounded TE). The ring's lines run straight
  from the core to the wall nodes, which are the C-grid's own tip-station
  nodes. The core's corners sit on 45 deg diagonals from the LE cuts and TE
  corners, so the three blocks at each corner share it at about 90 / 135 /
  135 deg.
- **Volume:** k leaves the tip plane along its outward normal, at the
  distances `k_distances` (the wrap's wall-normal schedule). Ring lines
  lean outwards by up to `flare_deg`, growing linearly from the core to
  the section's edge, so the cap overhangs the wing's side.
- **Faces:** every block is right-handed. `k_lo` is the wall (the tip
  face). The faces between cap blocks are SHARED, with axis maps found from
  their corner nodes. Every other face is FREE: an overset fringe.

The caps are built, not marched; `check_cap_blocks` refuses a cap with a
folded cell or a block without cells (ADR 0058).

## `check_cap_blocks` <Badge type="info" text="function" /> {#check_cap_blocks}

```python
check_cap_blocks(blocks: list) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L42">source</a>

Refuse a wing-tip cap with a folded cell or a block without cells (ADR 0058).

The caps are built, not marched, so no marcher guard sees them. Every
block must have at least one cell along each axis and every cell a
minimum corner scaled Jacobian above 0 (`shore.quality`).

**Raises**

- `ParameterError` — Naming the first invalid block and its worst cell.

## `section_normal` <Badge type="info" text="function" /> {#section_normal}

```python
section_normal(cross_section: np.ndarray, outward: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L84">source</a>

Unit normal of the section's plane, on the side of *outward*.

## `build_cap_blocks` <Badge type="info" text="function" /> {#build_cap_blocks}

```python
build_cap_blocks(
    cross_section: np.ndarray,
    *,
    outward: np.ndarray,
    le_cut_upper_idx: int,
    le_cut_lower_idx: int,
    k_distances: np.ndarray,
    n_cap_j: int = 12,
    j_spacing: str = 'tanh',
    j_beta: float = 4.0,
    core_frac: float = 0.5,
    flare_deg: float = 30.0,
    tip_inset: float = 0.0,
) -> dict[str, np.ndarray]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L92">source</a>

The O-H cap of one wing tip (ADR 0059): a core and four ring blocks.

**Parameters**

- **`cross_section`** — `(ni, 3)` wall nodes of the tip station, from the upper TE (i = 0) round the LE to the lower TE, in one plane. The TE must be blunt or rounded (`auto_round_sharp_te`): the `te` block spans its base. These nodes become the cap's edge, bit for bit when `tip_inset` is 0.
- **`outward`** — Any vector pointing away from the wing: it picks the side of the section's plane the cap grows on (`section_normal`).
- **`le_cut_upper_idx, le_cut_lower_idx`** — i-indices of the LE cuts, symmetric about the middle node (`find_symmetric_le_cuts`); they bound the `le` block, and the `te` block gets the same node count across the TE base.
- **`k_distances`** — `(nk,)` distances of the layers from the wall, from 0 and increasing.
- **`n_cap_j`** — Ring blocks' node count from the core to the wall.
- **`j_spacing, j_beta`** — Spacing law of that direction (`"uniform"` / `"tanh"` / `"tanh2"`; `"tanh"` clusters at the wall) and its strength.
- **`core_frac`** — The core's half-thickness as a fraction of the section's, in (0, 1).
- **`flare_deg`** — Outward lean of the ring's k-lines at the section's edge, in [0, 90) degrees, growing linearly from 0 at the core.
- **`tip_inset`** — Distance of the wall from the section, along the outward normal.

**Returns**

- `dict[str, ndarray]` — `CAP_BLOCKS` -&gt; `(ni, nj, nk, 3)` right-handed node arrays, the wall at `k = 0`. In the ring blocks j runs from the core to the section's wall.

**Raises**

- `ParameterError` — On a sharp TE, cuts that are not symmetric or do not straddle the LE, a section that is not single-valued along its chord, or parameters out of range.

## `wire_cap_blocks` <Badge type="info" text="function" /> {#wire_cap_blocks}

```python
wire_cap_blocks(blocks: list[HexBlock]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L295">source</a>

Make the faces two cap blocks share SHARED, with their axis maps.

Two faces are shared when their four corner nodes are the same points
(bit for bit). The axis map follows from where each corner lands on the
partner face, as `premesh` does for a blocking.

## `find_symmetric_le_cuts` <Badge type="info" text="function" /> {#find_symmetric_le_cuts}

```python
find_symmetric_le_cuts(
    cross_section: np.ndarray,
    le_frac: float,
    span_axis_world: int,
) -> tuple[int, int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L334">source</a>

Find the upper / lower LE-cut indices in a body cross-section.

The cross-section sweeps upper-TE-cut → suction → LE → pressure →
lower-TE-cut.  This locates the i-indices on each side of the LE
point that correspond to chord position `le_frac`, snapped to be
**symmetric** around the LE point in i-station count (so
`upper_idx + lower_idx == 2 * le_idx`).  Symmetry is required by
the cap-fill builder for the chord-cut SHARED face.

**Parameters**

- **`cross_section`** — `(ni, 3)` polyline.
- **`le_frac`** — Chord-fraction of the LE-cut, in `(0, 1)`.
- **`span_axis_world`** — Which world axis (0 / 1 / 2) is the span direction; the chord axis is the in-plane axis with the largest extent.

**Returns**

- **`(le_cut_upper_idx, le_cut_lower_idx)`** (tuple of int)

**Raises**

- `ParameterError` — On invalid `le_frac` or when the i-axis is too coarse to resolve the requested cut.

## `detect_te_thickness` <Badge type="info" text="function" /> {#detect_te_thickness}

```python
detect_te_thickness(cross_section: np.ndarray, span_axis_world: int) -> float
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L410">source</a>

Return the TE thickness of a cross-section polyline.

Sharp TE → `0.0` (i=0 and i=ni-1 coincide).  Blunt TE → distance
between i=0 and i=ni-1 in the airfoil's thickness direction.

## `auto_round_sharp_te` <Badge type="info" text="function" /> {#auto_round_sharp_te}

```python
auto_round_sharp_te(
    cross_section: np.ndarray,
    te_thickness: float,
    span_axis_world: int,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L425">source</a>

Synthesise a blunt cross-section by splitting the sharp TE node.

Returns a copy of `cross_section` with the i=0 and i=ni-1 nodes
moved apart by `te_thickness` along the airfoil-thickness
direction (perpendicular to chord, in the airfoil's 2D plane), so
the TE-base segment has finite extent.  Used by `CGridTipCap` to
auto-round sharp-TE inputs *only inside the cap*; the user's
body wrap is left untouched.

The thickness direction is the in-plane axis with the smaller
spread (the same heuristic used by `_detect_thickness_axis`).
The split is anti-symmetric: i=0 moves by +te_thickness/2 and
i=ni-1 by -te_thickness/2 so the TE-base midpoint stays at the
original sharp-TE location.

**Parameters**

- **`cross_section`** — `(ni, 3)` polyline.  Must have `cross_section[0]` and `cross_section[-1]` coincident (sharp TE).
- **`te_thickness`** — Desired TE thickness.  Must be &gt; 0.
- **`span_axis_world`** — World axis (0 / 1 / 2) that is the spanwise direction; used to identify the airfoil-thickness in-plane axis.

**Returns**

- `ndarray, shape ``(ni, 3)``, same dtype as input.`

**Raises**

- `ParameterError` — On invalid inputs or when the cross-section is not actually sharp.
