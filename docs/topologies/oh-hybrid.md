# OH hybrid topology (cubed-sphere wrap + H-channel)

The **OH hybrid** is a Style-2 conforming hybrid topology: a 6-block
cubed-sphere wrap (marched outward from the body STL surface) plus a
6-block H-channel (frustum blocks bridging the wrap envelope to a
user-supplied axis-aligned bounding box). Wrap and channel are
**SHARED-by-view** at the seam — they share one contiguous numpy
buffer per pair, so `wrap.k_hi` and `channel.k_lo` are bit-identical
with no copy.

The result is a **12-block, single-component near-body mesh** that
covers everything from the body surface out to a Cartesian box,
without the chimera-fringe handover that would otherwise be needed
between a body-fitted near-body mesh and a separate background. Use
it when you want a contiguous structured mesh from the wall to the
far field in one piece.

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
   over all `(i, j)` envelope nodes — robust because corners sit at
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
This is the natural continuation of the wrap's hyperbolic marching
into the channel — no clustering at the seam, no second-cell jump.

When `c1_seam=False`, the user-supplied `k_law_channel` /
`k_beta_channel` apply directly with no pinning. Useful only if you
deliberately want a cell-size discontinuity at the wrap-channel
seam (rarely the right thing for boundary-layer meshes).

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

## CLI

`shore mesh --topology oh` builds the OH hybrid end-to-end. The
far-field box is given by six scalar extents; channel depth and law
reuse the shared channel flags (also used by `--topology ch`):

```bash
shore mesh sphere.stl -o oh --topology oh \
  --ni 24 --nj 32 --nk 20 --ds 0.02 --growth 1.10 \
  --oh-box-xmin -3 --oh-box-xmax 3 \
  --oh-box-ymin -3 --oh-box-ymax 3 \
  --oh-box-zmin -3 --oh-box-zmax 3 \
  --nk-channel 10
# writes 12 oh_<label>.geo files
```

`--theta-cap` is ignored for `--topology oh` (the wrap pins it to
`pi/4`); passing it emits a `UserWarning`. OH requires a star-shaped
body (centroid ray-cast, like `cubed_sphere`).

## Python API

The Python API exposes the same topology with full control over the
box as 3-vectors:

```python
from shore.mesh import Mesh, OHTopology, Spacing1D

mesh = Mesh.from_stl(
    "sphere.stl",
    topology=OHTopology(
        box_vmin=(-3.0, -3.0, -3.0),
        box_vmax=(+3.0, +3.0, +3.0),
        nk_channel=10,
    ),
    ni=24, nj=32, nk=20,
    spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.10),
)
mesh.march(ds=0.02, growth=1.10)
mesh.write_geo("oh")          # 12 oh_<label>.geo files
mesh.write_vtk("oh")          # oh.vtm MultiBlock for ParaView
```

Constructor parameters:

| Name | Default | Description |
|---|---|---|
| `box_vmin`, `box_vmax` | — | Length-3 tuples; axis-aligned bounding box. Must contain the marched wrap envelope. |
| `nk_channel` | — | Number of channel layers (≥ 2). |
| `k_law_channel` | `"uniform"` | Channel wall-normal law (`"uniform"` / `"tanh"` / `"tanh2"`). Honoured only when `c1_seam=False`. |
| `k_beta_channel` | `3.0` | Clustering strength for tanh / tanh2 (when `c1_seam=False`). |
| `c1_seam` | `True` | When true, pin channel's first cell to the wrap's last marching step per `(i, j)` and grow geometrically out to the box face. |

`build()` accepts the same surface arguments as `CubedSphereTopology`
(`ni`, `nj`, `theta_cap`, `i_spacing` / `j_spacing` etc.). Note that
`theta_cap` is *ignored* — OH pins it to `pi/4`.

## Example

[`examples/08-oh-hybrid/visualize_oh_hybrid.py`](https://github.com/szaghi/shore/tree/main/examples/08-oh-hybrid)
is the canonical demo. It generates a unit-icosphere STL, builds the
12-block OH hybrid, marches the wrap, fills the channels, and writes
`examples_out/oh_<label>.geo` (×12) plus `examples_out/oh.vtm` for
ParaView.

```bash
python examples/08-oh-hybrid/visualize_oh_hybrid.py --generate --vtk --no-view
```

Channel knobs (`--box-half`, `--nk-channel`, `--channel-k-law`,
`--channel-k-beta`) are exposed on the CLI; pass `--help` for the
full list.

## See also

- [Cubed-sphere topology](/topologies/cubed-sphere) — the wrap construction
  used by OH (also usable standalone, with `theta_cap` and `phi_offset`
  free instead of pinned to π/4).
- [CH hybrid topology](/topologies/ch-hybrid) — the airfoil-/hydrofoil
  analogue (C-grid wrap + chord/normal-aligned far-field frame).
- [H-grid topology](/topologies/h-grid) — body-less channel meshes
  for internal flow (not the same as the OH channel block, which
  bridges a body wrap to a box).
- [`shore.mesh.topology.OHTopology`](/reference/api-topology#ohtopology)
  — Python API reference.
- [Hyperbolic marching](/algorithm/hyperbolic-marching) — the wrap
  uses the same cubed-sphere kernel.
