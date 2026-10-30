---
title: shore.topology.cgrid_cap
outline: [2, 3]
---

# `shore.topology.cgrid_cap`

Cap-fill builder for the C-grid topology's wing-tip closures.

The C-grid topology wraps a wing's outer surface but leaves the spanwise
boundaries (j=0 and j=n_stations-1) FREE — the wing is implicitly assumed
to extend to infinity in span, or to terminate at a Chimera-fringe boundary.
For a finite wing with closed tips, the j_lo / j_hi faces of the body wrap
need a **wall** closure: a thin spanwise slab covering the airfoil interior
cross-section at the tip plane.

This module builds that closure as **3 cap sub-blocks** (sharp TE) or **4
cap sub-blocks** (blunt TE) per tip, in an H-style "T-fill" layout that has
no singular axis (Xall-compatible).

Cap-block axis convention
-------------------------
Each cap sub-block is shaped `(ni_cap, nj_cap, 2, 3)`:

- **i-axis**: along the airfoil profile (chordwise-ish).  For `cap_upper`
  i runs from upper-LE-cut (i=0) to upper-TE-cut (i=ni_cap-1) — *reversed*
  w.r.t. the body's i-axis.  For `cap_lower` i runs from lower-LE-cut to
  lower-TE-cut — forward w.r.t. body i.  For `cap_le` i wraps from
  upper-LE-cut around the LE point to lower-LE-cut.
- **j-axis**: airfoil-thickness direction (perpendicular to chord, in the
  cross-section's 2D plane).  For `cap_upper` j runs from chord-cut at
  y=0 (j=0) to airfoil suction wall (j=nj_cap-1).  `cap_lower` is the
  mirror.  `cap_le` j runs from chord-perp inner line (j=0) to airfoil
  profile wrapping the LE (j=nj_cap-1).
- **k-axis**: spanwise.  **2 nodes, 1 cell** thick.  k=0 is the wing-tip
  wall (offset by `tip_inset` along `span_dir` from the body cross-section),
  k=1 is the body cross-section plane (SHARED with the body wrap's j_lo).

The resulting block is **not marched**.  It is a pure parametric fill,
fully determined by the body's j=0 cross-section polyline plus
`tip_inset` and the cap j-spacing law.

## `build_cap_blocks` <Badge type="info" text="function" /> {#build_cap_blocks}

```python
build_cap_blocks(
    cross_section: np.ndarray,
    le_cut_upper_idx: int,
    le_cut_lower_idx: int,
    span_dir: np.ndarray,
    tip_inset: float,
    n_cap_j: int,
    nk_cap: int = 2,
    j_spacing: str = 'tanh',
    j_beta: float = 4.0,
    is_blunt: bool = False,
    wake_base_cross_section: np.ndarray | None = None,
    te_cut_upper_idx: int | None = None,
    te_cut_lower_idx: int | None = None,
) -> dict[str, np.ndarray]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L45">source</a>

Build the cap sub-blocks at one wing tip of a C-grid.

K-axis convention
-----------------
The cap is a 3D structured slab spanning the wing-tip region.
Its k-axis is **spanwise**:

- **k=0**:        outboard wing-tip wall, offset from
                  `cross_section` by `tip_inset * span_dir`.
                  Wired to `WALL` by the topology layer.
- **k=nk_cap-1**: inboard chimera-overlap plane, coincident with
                  `cross_section`.  Wired to `FREE` (chimera
                  fringe receiving from the C-grid body wrap) by
                  the topology layer.

Interior layers (`1 <= k <= nk_cap - 2`) interpolate linearly
between the two end faces along `span_dir`.

**Parameters**

- **`cross_section`** — `(ni_body_post_ext, 3)` array.  The body's spanwise-extreme cross-section polyline.  i sweeps upper-TE-cut → suction → LE → pressure → lower-TE-cut.
- **`le_cut_upper_idx, le_cut_lower_idx`** — i-indices in `cross_section` at the LE-cut on the suction and pressure sides respectively.  Must satisfy `0 < le_cut_upper_idx < le_cut_lower_idx < ni_body - 1`.
- **`te_cut_upper_idx, te_cut_lower_idx`** — Optional i-indices in `cross_section` at the TE-cut on each side.  When `None` (default), the cap chord-band blocks extend all the way to the body's i=0 / i=ni-1 endpoints — which **collapses to a singular axis at a sharp TE** (where the wall point coincides with the chord cut).  For sharp TE always pass valid TE-cut indices that sit *strictly inboard* of the TE point so the cap_upper.i_hi / cap_lower.i_hi columns are non-degenerate.  The fluid sliver between the TE-cut and the actual TE point is left unmeshed by the cap; chimera-fringe overlap with the C-grid wake blocks closes the gap.  Must satisfy `0 <= te_cut_upper_idx < le_cut_upper_idx < le_cut_lower_idx < te_cut_lower_idx <= ni_body - 1` and the chord-band symmetry constraint `le_cut_upper_idx - te_cut_upper_idx == te_cut_lower_idx - le_cut_lower_idx` (equal i-station counts for cap_upper and cap_lower so the chord-cut SHARED face is 1:1).
- **`span_dir`** — Unit vector pointing **outward** along the span axis at this tip — i.e. away from the body interior.  Used to offset the wing-tip wall plane from the body cross-section by `tip_inset`.
- **`tip_inset`** — Spanwise thickness of the cap slab (k-axis extent in world units). Must be &gt; 0.
- **`n_cap_j`** — Node count along the cap's j-axis (airfoil-thickness direction). Must be &gt;= 2.  All cap sub-blocks share this resolution along their j-axis.
- **`nk_cap`** — Node count along the cap's k-axis (spanwise).  Must be &gt;= 2. Default `2` keeps a single-cell-thick slab; larger values give finer spanwise resolution at the wing tip.
- **`j_spacing`** — Spacing law for the cap j-axis: `'uniform'` / `'tanh'` / `'tanh2'`.  `'tanh'` (default) clusters nodes near the airfoil wall (j_hi side) for boundary-layer resolution.
- **`j_beta`** — Clustering strength for tanh / tanh2.
- **`is_blunt`** — When True, includes a 4th cap sub-block `'wake_base'` for the TE-base region.  `wake_base_cross_section` must be supplied.
- **`wake_base_cross_section`** — Required when `is_blunt=True`.  Shape `(ni_wake, nk_w, 3)`: the wake_base block's j=0 face nodes (the wake_base's spanwise end face, which the cap closes).  `ni_wake` is the wake's node count, `nk_w` is the wake_base's k-axis resolution across the TE thickness.

**Returns**

- `dict[str, np.ndarray]` — Maps cap sub-block label to its node array of shape `(ni_cap, nj_cap, nk_cap, 3)`.  Keys are always `'le'`, `'upper'`, `'lower'`, `'te_upper'`, and `'te_lower'` (the 5-block H-fill that closes the airfoil interior without singular axes — provided the cross-section has a non-zero TE thickness).  When `is_blunt=True` *and* the user's body wrap has a real `wake_base` block, an extra `'wake_base'` cap sub-block is included.

**Raises**

- `ParameterError` — On invalid shapes, indices, or unsupported spacing law.

## `find_symmetric_le_cuts` <Badge type="info" text="function" /> {#find_symmetric_le_cuts}

```python
find_symmetric_le_cuts(
    cross_section: np.ndarray,
    le_frac: float,
    span_axis_world: int,
) -> tuple[int, int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L489">source</a>

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
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L565">source</a>

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
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/cgrid_cap.py#L580">source</a>

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
