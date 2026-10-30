# OH hybrid topology (cubed-sphere wrap + H-channel)

The OH hybrid extends a [cubed-sphere](cubed-sphere) wrap, marched from
the body, with six channel blocks out to an axis-aligned box. Wrap and
channels share their seam by memory: each wrap block and its channel are
views of one buffer, so `wrap.k_hi` and `channel.k_lo` are the same
nodes. The result is one 12-block grid from the wall to the box, with no
chimera fringe between a near-body grid and a background. Built by
`shore mesh --topology oh`.

<Figure name="oh-hybrid" caption="OH hybrid, cut away: the marched cubed-sphere wrap and the six channel blocks reaching the faces of the box." />

## Block decomposition

12 structured blocks in canonical order:

| Index | Block | Shape | Role |
|---|---|---|---|
| 0 | `sub0` | `(ni, nc, nk_wrap, 3)` | Wrap equatorial sector facing the **+y** box face |
| 1 | `sub1` | `(ni, nc, nk_wrap, 3)` | Wrap equatorial sector facing the **−x** box face |
| 2 | `sub2` | `(ni, nc, nk_wrap, 3)` | Wrap equatorial sector facing the **−y** box face |
| 3 | `sub3` | `(ni, nc, nk_wrap, 3)` | Wrap equatorial sector facing the **+x** box face |
| 4 | `cap_north` | `(nc, nc, nk_wrap, 3)` | Wrap polar cap facing the **−z** box face |
| 5 | `cap_south` | `(nc, nc, nk_wrap, 3)` | Wrap polar cap facing the **+z** box face |
| 6 | `chan_sub0` | `(ni, nc, nk_channel, 3)` | Channel paired with `sub0`; outer face on **+y** |
| 7 | `chan_sub1` | `(ni, nc, nk_channel, 3)` | Channel paired with `sub1`; outer face on **−x** |
| 8 | `chan_sub2` | `(ni, nc, nk_channel, 3)` | Channel paired with `sub2`; outer face on **−y** |
| 9 | `chan_sub3` | `(ni, nc, nk_channel, 3)` | Channel paired with `sub3`; outer face on **+x** |
| 10 | `chan_north` | `(nc, nc, nk_channel, 3)` | Channel paired with `cap_north`; outer face on **−z** |
| 11 | `chan_south` | `(nc, nc, nk_channel, 3)` | Channel paired with `cap_south`; outer face on **+z** |

with `nc = nj // 4 + 1` (the cubed-sphere convention).

## How the 1:1 wrap-to-channel mapping works

The cubed-sphere wrap is built with **`theta_cap = pi/4`** and
**`phi_offset = pi/4`**, both pinned by `OHTopology` regardless of
what the user passes. With these values:

- The four equatorial sub-block *centres* (rather than corners) face
  the world cardinal axes — sub0→+y, sub1→−x, sub2→−y, sub3→+x.
- The two polar caps face ±z naturally.
- Each of the 6 wrap blocks therefore points at exactly one of the
  bounding box's 6 faces, allowing a 1:1 wrap-to-channel mapping.
- Adjacent wrap envelopes meet on **cube edges** (where two box
  faces meet); their shared lateral seams snap to the cube edge
  bit-exactly under the outer-quad construction.

Any other `theta_cap` would put the cap-equator seam inside a cube
face, forcing the channel outer faces to be L-shaped across cube
edges and giving up structured-block planarity. The pinning is
intrinsic to the OH geometry — not a knob.

## Channel construction (deferred fill)

Channel block interiors are not built at `topology.build()` time.
The wrap's `k_hi` envelope is unknown until marching completes, so
the channel's `k=0` row (which **is** the wrap's `k_hi` by view
aliasing) gets populated after `Mesh.march()` finishes. Fill
proceeds in three steps per channel:

1. **Pick a box face** by the wrap envelope's centroid (full mean
   over all `(i, j)` envelope nodes, not the corners: the corners sit at
   cap-equator junctions where the bulk is not).
2. **Snap the 4 envelope corners to that face's cube vertices**:
   face-normal component → plane value, in-plane components → vmin /
   vmax by sign relative to the box centre. Each corner lands at a
   cube vertex; adjacent channels share their seam corners exactly
   on the common cube edge.
3. **Ruled-surface fill**: bilinear `(a, b)` sweep of the snapped
   outer quad, then per-`k` interpolation between the wrap envelope
   and the outer face.

`Mesh.march()` runs this fill automatically when the topology is
`OHTopology`; manual drivers should call
`OHTopology.fill_channel_blocks(blocks)` after marching.

## C1 seam continuity

By default (`c1_seam=True`), the channel's first cell at every
`(i, j)` is **pinned to the wrap's last marching step at the same
`(i, j)`** — bit-exact spacing match across the seam. The remaining
`nk_channel - 1` cells then grow geometrically with a per-`(i, j)`
ratio solved from
$$
\sum_{k=0}^{n-1} r^k = \frac{\text{gap}_{ij}}{\text{first\_step}_{ij}}
$$
so the channel column fills the gap exactly while preserving smooth
geometric growth from the pinned first step out to the box face.
The channel continues the wrap's growth: no clustering at the seam and
no jump at the second cell.

When `c1_seam=False`, the user-supplied `k_law_channel` /
`k_beta_channel` apply directly with no pinning. Useful only if you
deliberately want a cell-size discontinuity at the wrap-channel
seam (rarely wanted).

## Boundary conditions

| Block | Face | BC | Partner |
|---|---|---|---|
| Any wrap | `k_lo` | `WALL` | — (body STL surface) |
| Any wrap | `k_hi` | `SHARED` | Paired channel's `k_lo` (view-aliased) |
| Any channel | `k_lo` | `SHARED` | Paired wrap's `k_hi` (view-aliased) |
| Any channel | `k_hi` | `FREE` | — (box face, far-field / Chimera fringe) |
| Wrap sub-blocks | `j_lo`, `j_hi` | `SHARED` | Adjacent sub-block's `j_hi`, `j_lo` |
| Channel sub-blocks | `j_lo`, `j_hi` | `SHARED` | Adjacent channel sub-block's `j_hi`, `j_lo` |
| Wrap sub-blocks | `i_lo`, `i_hi` | `DIRICHLET` | Cap blocks (cubed-sphere seam table) |
| Channel sub-blocks | `i_lo`, `i_hi` | `DIRICHLET` | Channel cap blocks (mirrors wrap) |

The lateral seam wiring on the channel side mirrors the wrap exactly
— same partners, same `AxisMap` — so any downstream `cc.par`
generation that handles the cubed-sphere lateral seams handles the
channel lateral seams without changes.

## Parameters

| `shore mesh` flag | `OHTopology` argument | Default | Meaning |
|---|---|---|---|
| `--oh-box-{x,y,z}{min,max}` | `box_vmin`, `box_vmax` | required | the box; must contain the marched wrap |
| `--nk-channel` | `nk_channel` | 5 (CLI) | channel layers, the shared seam included (≥ 2) |
| `--k-law-channel`, `--k-beta-channel` | `k_law_channel`, `k_beta_channel` | `uniform`, 3.0 | channel law, used only with the C1 pin off |
| `--c1-seam-channel` / `--no-c1-seam-channel` | `c1_seam` | on | pin each channel's first cell to the wrap's last step |

The wrap takes the cubed-sphere flags (`--ni`, `--nj`, `--nk`, spacing,
marching). `--theta-cap` / `--theta-cap-deg` are ignored with a
`UserWarning`: the wrap pins $\theta_{cap} = \pi/4$. `Mesh.from_stl`
accepts `OHTopology`, and `Mesh.march` fills the channels after marching
the wrap. Commands: [`shore mesh`](/reference/mesh).

## Limits

- The body must be star-shaped from the anchor, as for the cubed sphere.
- With $\theta_{cap} = \pi/4$ the caps are large; the cap–equator
  spacing pin still needs compatible `(ni, nj)`
  ([k=0 construction](/concepts/k0-construction#compatibility-window-and-fallback)).
- The channel blocks do not march: their quality depends on how far the
  box is from the envelope relative to `nk_channel`. The hyperbolic
  scheme's acceptance does not cover them.

## See also

- [Cubed-sphere topology](/concepts/topologies/cubed-sphere) — the wrap construction
  used by OH (also usable standalone, with `theta_cap` and `phi_offset`
  free instead of pinned to π/4).
- [CH hybrid topology](/concepts/topologies/ch-hybrid) — the airfoil-/hydrofoil
  analogue (C-grid wrap + chord/normal-aligned far-field frame).
- [H-grid topology](/concepts/topologies/h-grid) — body-less channel meshes
  for internal flow (not the same as the OH channel block, which
  bridges a body wrap to a box).
- [`OHTopology`](/reference/api/shore.mesh.topology#OHTopology) — Python API reference.
- [Marching](/concepts/hyperbolic-marching) — the wrap uses the
  cubed-sphere marcher.
