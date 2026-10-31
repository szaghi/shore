# CH hybrid topology (C-grid wrap + H-channel)

The CH hybrid extends a [C-grid](c-grid) wrap (3 blocks for a sharp TE,
4 for a blunt one) with channel blocks out to a rectangular frame in the
chord/normal plane, extruded across the span. The result is one grid of
12 (sharp TE) or 14 (blunt TE) blocks from the wing wall to the frame,
with point-matched seams throughout and no chimera fringe between the
near-body grid and the far field. It is the wing counterpart of the
[OH hybrid](oh-hybrid). Built by `shore mesh --topology ch`.

<Figure name="ch-hybrid" caption="CH hybrid in a span section: the C-grid wrap, the five channel blocks to the top, front and bottom of the frame, and the outflow blocks behind the wake." />

## Block decomposition

12 (sharp TE) or 14 (blunt TE) structured blocks in canonical
order:

| Index | Block | Shape | Role |
|---|---|---|---|
| 0 | `upper` | `(ni_wake, n_stations, nk_wrap, 3)` | Wrap — upper wake (i sweeps upper-far-tip → upper-TE) |
| 1 | `body` | `(ni_body, n_stations, nk_wrap, 3)` | Wrap — body (i sweeps upper-TE → LE → lower-TE) |
| 2 | `lower` | `(ni_wake, n_stations, nk_wrap, 3)` | Wrap — lower wake (i sweeps lower-TE → lower-far-tip) |
| 3 | `wake_base` *(blunt only)* | `(ni_wake, n_stations, nk_w, 3)` | Wrap — wake-base block bridging the TE thickness gap; `i_hi` is the bluff-base wall |
| 3/4 | `chan_upper` | `(ni_wake, n_stations, nk_channel, 3)` | k_hi-side channel — bridges upper.k_hi to the top frame edge |
| 4/5 | `chan_body_top` | `(i_split_upper+1, n_stations, nk_channel, 3)` | k_hi-side channel — bridges body's i=0..i_split_upper segment to the top frame edge |
| 5/6 | `chan_body_front` | `(i_split_lower-i_split_upper+1, n_stations, nk_channel, 3)` | k_hi-side channel — bridges body's middle segment around the LE to the front frame edge |
| 6/7 | `chan_body_bot` | `(ni_body-i_split_lower, n_stations, nk_channel, 3)` | k_hi-side channel — bridges body's i_split_lower..ni_body-1 segment to the bottom frame edge |
| 7/8 | `chan_lower` | `(ni_wake, n_stations, nk_channel, 3)` | k_hi-side channel — bridges lower.k_hi to the bottom frame edge |
| 8/9 | `chan_outflow_top_upper` | `(ni_outflow, n_stations, nk_channel, 3)` | Outflow — top half on the upper side (conformal SHARED to `chan_upper.i_lo`) |
| 9/10 | `chan_outflow_wake_upper` | `(ni_outflow, n_stations, nk_wrap, 3)` | Outflow — wake half on the upper side (conformal SHARED to `upper.i_lo`) |
| 10/11 | `chan_outflow_wake_lower` | `(ni_outflow, n_stations, nk_wrap, 3)` | Outflow — wake half on the lower side (conformal SHARED to `lower.i_hi`) |
| 11/12 | `chan_outflow_top_lower` | `(ni_outflow, n_stations, nk_channel, 3)` | Outflow — top half on the lower side (conformal SHARED to `chan_lower.i_hi`) |
| 13 *(blunt only)* | `chan_outflow_wake_base` | `(ni_outflow, n_stations, nk_w, 3)` | Outflow — bridge between wake_lower and wake_upper across the TE-thickness gap (conformal SHARED to `wake_base.i_lo`) |

`i_split_upper` / `i_split_lower` are body-block i-indices at which
`body.k_hi` is partitioned into 3 sub-channels; defaults are
`ni_body // 3` and `2 * ni_body // 3` respectively. `ni_outflow`
sets the chord-direction node count of the four (or five) outflow
sub-blocks; default `8`.

## Far-field frame

CH's frame is **rectangular in the chord/normal plane**, extruded
across span. It is defined by four world-coord extents
(`chord_min`, `chord_max`, `normal_min`, `normal_max`) plus a span
direction inherited from the C-grid wrap (`span_axis`) and a
freestream direction (`freestream_world`) used to pick the chord
axis:

- `chord_axis` = the non-span world axis with the largest absolute
  freestream component (with `freestream_world=None`: the first non-span axis, x unless the span is x, positive direction).
- `normal_axis` = the remaining non-span axis.
- `span_axis_idx` = `"xyz".index(span_axis)`.

Each channel's outer face is a planar quad on the frame. The two body
splits (`i_split_body_upper`, `i_split_body_lower`) are mapped to the
frame's top-front and bottom-front corners, so `chan_body_top`,
`chan_body_front` and `chan_body_bot` end on the top, front and bottom
edges of the frame, and adjacent quads share their corners. The outflow
blocks close the region behind the wrap up to the rear edge
(`chord_max`).

## Memory layout (SHARED-by-view)

The wrap and the **5 k_hi-side channels** share a single contiguous
numpy buffer per wrap segment, of depth `nk_wrap + nk_channel - 1`.
Each k_hi-side channel is a numpy view onto the buffer's deeper
slots:

- `nodes_all[:, :, :nk_wrap, :]` — wrap views (upper / body / lower)
- `nodes_all[:, :, nk_wrap-1:, :]` — corresponding channel views

`wrap.k_hi` (the marched outer envelope) and `channel.k_lo` are
therefore the same buffer slot — bit-identical with no copy. The
**outflow sub-blocks** use a similar two-stack layout: `wake_X` +
`top_X` per side (k=0 of `top_X` aliased to k=last of `wake_X` at
the marched envelope normal). Wake-cut SHARED between the two
sides is enforced by data copy in `fill_channel_blocks`, with the
face-graph correctly stamped.

## Conformal SHARED to the wrap

The layout avoids any non-conformal interface between the
wrap and the outflow channels:

| Outflow sub-block | Conformal SHARED partner | nk constraint |
|---|---|---|
| `chan_outflow_top_upper.i_lo` | `chan_upper.i_lo` | `nk_channel` |
| `chan_outflow_top_lower.i_lo` | `chan_lower.i_hi` | `nk_channel` |
| `chan_outflow_wake_upper.i_lo` | `upper.i_lo` | `nk_wrap` |
| `chan_outflow_wake_lower.i_lo` | `lower.i_hi` | `nk_wrap` |
| `chan_outflow_wake_base.i_lo` *(blunt)* | `wake_base.i_lo` | `nk_w` |

Same-axis i↔i seams. The blocks above the wake cut (`top_upper`,
`wake_upper`) and the blunt `wake_base` would be left-handed with their
channel k pointing up, so their j runs against their partners' and those
seams map j with a flip; the lower-side seams are identity
([ADR 0057](/adr/0057-right-handed-ch-outflow)). Every CH block is
right-handed. The outflow's i_lo face data is copied bit-exact from the
wrap face (j-reversed where flipped) in `fill_channel_blocks`, so the
geometric coincidence at the seam is exact (not a tolerance-based match).

## Channel construction (deferred fill)

Channel block interiors are not built at `topology.build()` time.
The wrap's `k_hi` envelope is unknown until marching completes. The
fill proceeds as follows in `CHTopology.fill_channel_blocks`:

**k_hi-side channels (5)**: same algorithm as
[OH hybrid](oh-hybrid#channel-construction-deferred-fill) — bilinear
sweep of an outer-frame quad, plus optional per-`(i, j)` C1 first-step
pin to the wrap's last marching step.

**Outflow sub-blocks (4 or 5)**:

1. Copy the inner face (channel `i_lo`) bit-exact from the
   conformal-SHARED partner (`chan_X.i_lo` / `wrap.X.i_lo` /
   `wake_base.i_lo`).
2. Construct the outer face (channel `i_hi`) by replacing the chord
   coordinate with `chord_max` (rear frame edge); span and normal
   match the inner face.
3. Sweep the interior `i = 1..ni_outflow-2` columns by linear
   interpolation in chord under `k_law_channel`.

`Mesh.march()` runs this fill automatically when the topology is
`CHTopology`; users who drive marching manually must call
`CHTopology.fill_channel_blocks(blocks)` after marching.

## C1 seam continuity (k_hi-side channels only)

By default (`c1_seam=True`), each k_hi-side channel's first cell at
every `(i, j)` is **pinned to the wrap's last marching step at the
same `(i, j)`** — same algorithm as the OH hybrid. The remaining
`nk_channel - 1` cells grow geometrically out to the frame edge,
filling the gap exactly with smooth growth.

Outflow sub-blocks do not use C1 pinning: their channel-i runs along
chord (not normal), so there is no wrap k-axis to inherit from. The
chord distribution shape is set directly by `k_law_channel` /
`k_beta_channel`.

## Boundary conditions

**Wrap blocks** (same as standalone [C-grid](c-grid)):

| Block | Face | BC | Partner |
|---|---|---|---|
| `body` | `k_lo` | `WALL` | — (body STL) |
| `upper` / `lower` | `k_lo` | `SHARED` | wake-cut partner (or `wake_base` for blunt) |
| `body` | `i_lo` | `SHARED` | `upper.i_hi` (upper-TE column) |
| `body` | `i_hi` | `SHARED` | `lower.i_lo` (lower-TE column) |
| `upper` | `i_lo` | `SHARED` | `chan_outflow_wake_upper.i_lo` |
| `lower` | `i_hi` | `SHARED` | `chan_outflow_wake_lower.i_lo` |
| `wake_base` | `i_lo` *(blunt)* | `SHARED` | `chan_outflow_wake_base.i_lo` |
| `wake_base` | `i_hi` *(blunt)* | `WALL` | — (TE-base bluff wall) |
| Any wrap | `j_lo`, `j_hi` | `FREE` | — (wing tips) |

**k_hi-side channels (5)**:

| Block | Face | BC | Partner |
|---|---|---|---|
| `chan_upper` / `chan_body_*` / `chan_lower` | `k_lo` | `SHARED` | Wrap's `k_hi` (view-aliased) |
| Same | `k_hi` | `FREE` | Top / bottom / front frame edge |
| Same | `i_lo` / `i_hi` | `SHARED` | Adjacent channel (lateral seam) or outflow sub-block |
| Same | `j_lo`, `j_hi` | `FREE` | Wing tips |

**Outflow sub-blocks**:

| Block | Face | BC | Partner |
|---|---|---|---|
| `chan_outflow_top_*` | `k_lo` | `SHARED` | Adjacent `chan_outflow_wake_*` (envelope seam, view-aliased) |
| `chan_outflow_top_upper` | `k_hi` | `FREE` | Top frame edge |
| `chan_outflow_top_lower` | `k_hi` | `FREE` | Bottom frame edge |
| `chan_outflow_wake_*` | `k_lo` | `SHARED` | Wake-cut partner (sharp: across; blunt: via wake_base) |
| `chan_outflow_wake_*` | `k_hi` | `SHARED` | Adjacent `chan_outflow_top_*` |
| Any outflow | `i_lo` | `SHARED` | Conformal partner (see above) |
| Any outflow | `i_hi` | `FREE` | Rear frame edge |
| Any outflow | `j_lo`, `j_hi` | `FREE` | Wing tips |

## Parameters

| `shore mesh` flag | `CHTopology` argument | Default | Meaning |
|---|---|---|---|
| `--ch-chord-min`, `--ch-chord-max`, `--ch-normal-min`, `--ch-normal-max` | `chord_min`, ... | required | frame extents, world coordinates; must contain the marched wrap |
| `--ch-freestream-{x,y,z}` | `freestream_world` | `None` (first non-span axis) | picks the chord axis |
| `--span-axis` | `span_axis` | `z` | span axis |
| `--nk-channel` | `nk_channel` | 5 | channel layers, the shared seam included (≥ 2) |
| `--k-law-channel`, `--k-beta-channel` | `k_law_channel`, `k_beta_channel` | `uniform`, 3.0 | channel law (k_hi-side channels with the C1 pin off, and the outflow chord sweep) |
| `--c1-seam-channel` / `--no-c1-seam-channel` | `c1_seam` | on | pin each channel's first cell to the wrap's last step |
| `--ch-i-split-body-upper`, `--ch-i-split-body-lower` | `i_split_body_upper`, `i_split_body_lower` | `ni_body // 3`, `2 * ni_body // 3` | where `body.k_hi` is split into three channels |
| `--ch-ni-outflow` | `ni_outflow` | 8 | chord-direction nodes of each outflow block |

The wrap takes the C-grid flags (`--ni-body`, `--ni-wake`, `--n-stations`,
`--wake-length`, `--wake-growth`, `--te-dihedral-deg`) and the marching
flags. As for the C-grid, `Mesh.from_stl` does not accept `CHTopology`:
in Python, build the blocks with `CHTopology(...).build(...)` and wrap
them in a `Mesh`. The [wing tutorial](/tutorials/wing) and
[`shore mesh`](/reference/mesh) have runnable commands.

## Limits

- The frame must contain the marched wrap envelope.
- The channel blocks do not march: their quality depends on how far the
  frame is from the envelope relative to `nk_channel`.
- The hyperbolic scheme's acceptance spec does not cover the channel
  blocks or blunt-TE wraps ([ADR 0005](/adr/0005-marching-schemes)).

## See also

- [C-grid topology](/concepts/topologies/c-grid) — the wrap standalone (without
  far-field channel).
- [OH hybrid topology](/concepts/topologies/oh-hybrid) — the cubed-sphere
  analogue (closed envelope wrap + 6-face H-channel).
- [`CHTopology`](/reference/api/shore.mesh.topology#CHTopology) — Python API reference.
- [Marching](/concepts/hyperbolic-marching) — the wrap uses the same
  marcher as the standalone C-grid.
