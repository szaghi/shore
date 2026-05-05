# CH hybrid topology (C-grid wrap + H-channel)

The **CH hybrid** is the airfoil/hydrofoil analogue of the
[OH hybrid](oh-hybrid): a C-grid wrap (3-block sharp-TE or 4-block
blunt-TE) plus an H-channel that bridges the wrap to a far-field
rectangular frame in the chord/normal plane (extruded across span).
The result is a **12-block** (sharp TE) or **14-block** (blunt TE),
single-component mesh that covers everything from the body wall to
the far field — body wrap + wake + far-field channel — without a
chimera-fringe handover.

Use it for sharp- or blunt-TE wing / airfoil meshes where you want
a contiguous structured mesh from the wall to the box in one piece.

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
  freestream component (default `+x` when `freestream_world=None`).
- `normal_axis` = the remaining non-span axis.
- `span_axis_idx` = `"xyz".index(span_axis)`.

The frame's six perimeter regions and which channel sub-block sits
on each:

```
                                    normal_max
        ┌───────────────────────────────────────┬──────────────┐
        │       chan_upper                      │              │
        │       (top frame, wake side)          │ chan_outflow │
        │                                       │  _top_upper  │
        ├──────────┬─────────────────┬──────────┼──────────────┤
        │ chan_    │                 │          │ chan_outflow │
        │ body_top │   wrap (3 blk)  │ chan_    │ _wake_upper  │
        │          │                 │ body_    │              │
        │   ┌──────┴──┐         ┌────┤ front    │  wake-cut    │
        │   │ wrap    │  body   │ wrap   ──────┤──────────────┤
        │   │ upper   │         │ lower│       │ chan_outflow │
        │   └──────┬──┘         └────┤  body_  │ _wake_lower  │
        │ chan_    │                 │  bot     │              │
        │ body_bot │                 │          │ chan_outflow │
        ├──────────┴─────────────────┴──────────┼  _top_lower  │
        │       chan_lower                      │              │
        │       (bottom frame, wake side)       │              │
        └───────────────────────────────────────┴──────────────┘
                                                  chord_max
                          normal_min
```

(Schematic — actual block extents depend on body shape and wake
length; only the chord/normal directions are shown, span is the
out-of-page axis.)

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

The Phase-1.2 design avoids any non-conformal interface between the
wrap and the outflow channels:

| Outflow sub-block | Conformal SHARED partner | nk constraint |
|---|---|---|
| `chan_outflow_top_upper.i_lo` | `chan_upper.i_lo` | `nk_channel` |
| `chan_outflow_top_lower.i_lo` | `chan_lower.i_hi` | `nk_channel` |
| `chan_outflow_wake_upper.i_lo` | `upper.i_lo` | `nk_wrap` |
| `chan_outflow_wake_lower.i_lo` | `lower.i_hi` | `nk_wrap` |
| `chan_outflow_wake_base.i_lo` *(blunt)* | `wake_base.i_lo` | `nk_w` |

All identity AxisMap (same-axis i↔i). The outflow's i_lo face data
is copied bit-exact from the wrap face in `fill_channel_blocks`, so
the geometric coincidence at the seam is exact (not a tolerance-based
match).

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
| `upper` | `i_lo` | `SHARED` | `chan_outflow_wake_upper.i_lo` (Phase-1.2 conformal) |
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

## CLI

```bash
shore mesh airfoil.stl --topology ch \
    --ni-body 80 --ni-wake 20 --n-stations 5 --wake-length 15.0 \
    --nk 12 --ds 5e-4 --growth 1.08 \
    --ch-chord-min -2.0 --ch-chord-max 20.0 \
    --ch-normal-min -2.5 --ch-normal-max 2.5 \
    --nk-channel 8 --k-law-channel uniform \
    -o ch
```

Required flags (in addition to the standard `--ni`-style options):

| Flag | Description |
|---|---|
| `--ni-body INTEGER` | C-grid body wrap node count |
| `--ni-wake INTEGER` | Per-branch wake node count |
| `--n-stations INTEGER` | Span-direction station count |
| `--wake-length FLOAT` | Wake length downstream of TE (body units) |
| `--ch-chord-min FLOAT` | Frame minimum chord coordinate |
| `--ch-chord-max FLOAT` | Frame maximum chord coordinate |
| `--ch-normal-min FLOAT` | Frame minimum normal coordinate |
| `--ch-normal-max FLOAT` | Frame maximum normal coordinate |

Channel-specific flags:

| Flag | Default | Description |
|---|---|---|
| `--nk-channel INTEGER` | 5 | Total channel-block layers (incl. shared seam at k=0). |
| `--k-law-channel TEXT` | `uniform` | Channel k-axis spacing law (`uniform` / `tanh` / `tanh2`). |
| `--k-beta-channel FLOAT` | 3.0 | Clustering strength for tanh / tanh2. |
| `--c1-seam-channel/--no-c1-seam-channel` | on | Pin channel first step to wrap last step (k_hi-side channels only). |

The standard wrap flags (`--span-axis`, `--wake-growth`, `--ds`,
`--growth`, `--smooth*`) carry over from `--topology cgrid`.

## Python API

```python
from shore.mesh import Body, CHTopology, Mesh, Spacing1D

body = Body.from_stl("airfoil.stl", skip_hygiene=True)

topo = CHTopology(
    chord_min=-2.0, chord_max=20.0,
    normal_min=-2.5, normal_max=2.5,
    nk_channel=8,
    k_law_channel="uniform",
    c1_seam=True,
)
blocks = topo.build(
    surface=None, nk=12, mesh=body.mesh,
    spacing_k=Spacing1D(law="geometric", ds=5e-4, growth=1.08),
    ni_body=80, ni_wake=20, n_stations=5, wake_length=15.0,
)
mesh = Mesh(body=body, blocks=blocks, topology=topo)
mesh.march(ds=5e-4, growth=1.08)
mesh.write_geo("ch")          # 12 ch_<label>.geo files (sharp TE)
```

Constructor parameters:

| Name | Default | Description |
|---|---|---|
| `chord_min`, `chord_max`, `normal_min`, `normal_max` | — | Frame extents in world coords. Must contain the marched wrap envelope. |
| `nk_channel` | — | k_hi-side channel layers (≥ 2). |
| `freestream_world` | `None` (= +x) | Length-3 vector picking the chord axis. |
| `span_axis` | `"z"` | `"x"` / `"y"` / `"z"` — span world axis. |
| `k_law_channel` | `"uniform"` | k-spacing law for k_hi-side channels (when `c1_seam=False`) and for outflow chord-direction sweeps. |
| `k_beta_channel` | `3.0` | Clustering strength for tanh / tanh2. |
| `c1_seam` | `True` | Pin k_hi-side channel's first step to wrap's last step (per `(i, j)`). |
| `i_split_body_upper`, `i_split_body_lower` | `ni_body // 3`, `2 * ni_body // 3` | Body i-indices for the 3-way `chan_body_*` split. |
| `ni_outflow` | `8` | Chord-direction node count of every outflow sub-block. |

`build()` accepts the same C-grid wrap arguments as
[`CGridTopology.build`](c-grid#python-api): `ni_body`, `ni_wake`,
`n_stations`, `wake_length`, `wake_growth`, `te_dihedral_deg`, plus
`mesh=` (a `trimesh.Trimesh`) when `surface=None`. `Mesh.from_stl`
does **not** accept CHTopology — use `Mesh(body, blocks, topology)`
directly.

## Example

[`examples/09-ch-hybrid/visualize_ch_hybrid.py`](https://github.com/szaghi/shore/tree/main/examples/09-ch-hybrid)
generates a NACA 0012 wing STL, builds the CH hybrid (12 blocks for
sharp TE), marches the wrap, fills the channels, and writes
`examples_out/ch_<label>.geo` (×12) plus `examples_out/ch.vtm` for
ParaView.

```bash
python examples/09-ch-hybrid/visualize_ch_hybrid.py --generate --vtk --no-view
```

A bash sibling
[`examples/09-ch-hybrid/visualize_ch_hybrid.sh`](https://github.com/szaghi/shore/tree/main/examples/09-ch-hybrid)
chains `shore profile naca` → `shore body extrude` → `shore mesh
--topology ch` → `shore export` to produce the same outputs:

```bash
./examples/09-ch-hybrid/visualize_ch_hybrid.sh
./examples/09-ch-hybrid/visualize_ch_hybrid.sh --ni-body 120 --nk 16 --nk-channel 12
```

## See also

- [C-grid topology](/topologies/c-grid) — the wrap standalone (without
  far-field channel).
- [OH hybrid topology](/topologies/oh-hybrid) — the cubed-sphere
  analogue (closed envelope wrap + 6-face H-channel).
- [`shore.mesh.topology.CHTopology`](/reference/api-topology#chtopology)
  — Python API reference.
- [Hyperbolic marching](/algorithm/hyperbolic-marching) — the wrap
  uses the same kernel as standalone C-grid.
