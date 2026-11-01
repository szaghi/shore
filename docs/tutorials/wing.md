---
title: A wing with a C-grid
---

# A wing with a C-grid

You build a NACA 0012 wing from its section, wrap it in a C-grid with wake
branches, and extend the C-grid to a far-field box with the CH hybrid. You
learn how sharp and blunt trailing edges change the block layout, and why
the tips are not capped yet.

::: info You need
- `shore-mesh` (no extras); `shore export` needs the `vis` extra.
- About 10 minutes. Every grid builds in about a second.
- The script: [`examples/tutorials/wing.py`](https://github.com/szaghi/shore/blob/main/examples/tutorials/wing.py).
  Each CLI command below writes the same bytes as the matching region of the
  script.
:::

## 1. Make the section and the wing

`shore profile naca` writes the section as a Selig `.dat` (upper TE → LE →
lower TE, `2 n_pts − 1` points). `shore body extrude` stacks copies of it
along the span; `--end-caps` closes the ends with triangle fans so the STL is
watertight. The STL has more span stations (8) than the grid will (5): the
C-grid slices the body between its span extremes, and the slices must not
land on the end-cap fans.

::: code-group
```bash [CLI]
shore profile naca 0012 --n-pts 40 -o naca0012.dat
shore body extrude --section naca0012.dat --span 0.5 --n-stations 8 \
    --span-axis z --end-caps -o wing.stl
```
<<< @/../examples/tutorials/wing.py#body{python} [Python]
:::

```
✓ Wrote naca0012.dat  (79 points, closed TE)
✓ Wrote wing.stl  (626 vertices, 1248 faces, wing, 8 stations, capped)
```

The chord is 1, the span 0.5 along z (from z = 0 to 0.5). The Python
version reads the section back from the `.dat`, as the CLI does: the file
holds 7 decimals, and the same coordinates give the same STL.

## 2. Wrap it in a C-grid

The C-grid wraps the section from the trailing edge round the leading edge
and back, and continues along two wake branches. With a sharp trailing edge
it is three blocks: `body` (the wrap) and `upper` / `lower` (the wake
branches), joined by SHARED seams at the TE columns and, along the wake cut,
`upper.k_lo` to `lower.k_lo`.

::: code-group
```bash [CLI]
shore mesh wing.stl --topology cgrid --ni-body 120 --ni-wake 30 --n-stations 5 \
    --wake-length 15 --nk 20 --ds 5e-4 --growth 1.15 --smooth 0.1 --smooth-iters 1 \
    -o wing --adjacency-json wing.adjacency.json -v
```
<<< @/../examples/tutorials/wing.py#cgrid{python} [Python]
:::

```
✓ Loaded wing.stl  (1248 triangles)
  Marching       ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 19/19 layer 19/19, min J 9.78e-06 0:00:00
✓ Wrote 3 .geo files (upper, body, lower)
✓ Wrote wing.adjacency.json
```

The script prints the node counts and the range of the scaled Jacobian of
each block:

```
upper                    (30, 5, 20)  sJ [+0.966, +1.000]
body                     (120, 5, 20)  sJ [+0.903, +1.000]
lower                    (30, 5, 20)  sJ [+0.966, +1.000]
```

| Parameter | Meaning |
|---|---|
| `ni_body` / `--ni-body` | nodes round the section, TE → LE → TE |
| `ni_wake` / `--ni-wake` | nodes along each wake branch, the TE node included |
| `n_stations` / `--n-stations` | span stations: the j nodes of every block |
| `wake_length` / `--wake-length` | length of each wake branch from the TE, in body units (here chords) |
| `--wake-growth` | geometric growth of the wake i-spacing (default 1.15) |

The wall-normal thickness here is $5\times10^{-4}\,(1.15^{19}-1)/0.15 \approx 0.044$
chords: a thin near-body grid that a background grid overlaps.

<Figure name="cgrid-wing" caption="C-grid around the NACA 0012 wing: a span slice of the volume and the wall." />

### Sharp and blunt trailing edges

SHORE detects the trailing edge by the dihedral angle of the STL
(`--te-dihedral-deg`, default 30°). A blunt TE needs a fourth block,
`wake_base`, that meshes the gap between the two wake branches. The default
layout (`--cgrid-blocks 3`) refuses a blunt body:

```bash
shore profile naca 0012 --n-pts 40 --open-te -o naca0012_blunt.dat
shore body extrude --section naca0012_blunt.dat --span 0.5 --n-stations 8 \
    --end-caps -o wing_blunt.stl
shore mesh wing_blunt.stl --topology cgrid --ni-body 120 --ni-wake 30 --n-stations 5 \
    --wake-length 15 --nk 20 --ds 5e-4 --growth 1.15 --smooth 0.1 --smooth-iters 1 -o blunt
```

```
Error: CGridTopology.build: 3-block layout requires a sharp TE, but the body's TE was
detected as blunt. Either close the TE (e.g. shore profile naca --closed-te) or pass
blocks=4 to get the wake-base block W that meshes the TE gap.
```

Add `--cgrid-blocks 4` (`blocks=4` in Python) and the run writes
`upper`, `body`, `lower` and `wake_base`; `wake_base` is `30 × 5 × 3`
nodes, its three k nodes spanning the 0.25 % chord TE gap.

## 3. Close the tips

The C-grid is open at the two span ends. SHORE had two ways to close a tip,
`--cap-tips` on `shore mesh` and `shore cap-tip`, but both built folded or
empty blocks (ADR 0057), so both are now refused (ADR 0058) until the cap is
redesigned:

```bash
shore mesh wing.stl --topology cgrid --ni-body 120 --ni-wake 30 --n-stations 5 \
    --wake-length 15 --nk 20 --ds 5e-4 --growth 1.15 --cap-tips both -o capped
```

```
Error: Folded cell (scaled Jacobian -1.000) in tip cap block cap_jlo_le at cell
(0, 0, 0), 1881 of 1881 folded: the cap construction cannot close this section;
no cap was built (ADR 0058).
```

## 4. Extend to a far-field box: the CH hybrid

A C-grid alone ends a few percent of a chord off the wall. The CH hybrid
adds non-marched channel blocks that bridge the wrap's outer layer to an
axis-aligned box in the chord–normal plane, spanning the body's span. The
box is set by `chord_min` / `chord_max` / `normal_min` / `normal_max`; the
chord axis is the largest non-span component of the freestream
(`--ch-freestream-{x,y,z}`, here the default).

::: code-group
```bash [CLI]
shore mesh wing.stl --topology ch --ni-body 80 --ni-wake 20 --n-stations 5 \
    --wake-length 15 --nk 12 --ds 5e-4 --growth 1.08 --smooth 0 --smooth-iters 0 \
    --ch-chord-min -2 --ch-chord-max 20 --ch-normal-min -2.5 --ch-normal-max 2.5 \
    --nk-channel 8 -o ch -v
```
<<< @/../examples/tutorials/wing.py#ch{python} [Python]
:::

```
✓ Wrote 12 .geo files (upper, body, lower, chan_upper, chan_body_top, chan_body_front,
chan_body_bot, chan_lower, chan_outflow_top_upper, chan_outflow_wake_upper,
chan_outflow_wake_lower, chan_outflow_top_lower)
```

```
upper                    (20, 5, 12)  sJ [+1.000, +1.000]
body                     (80, 5, 12)  sJ [+0.878, +1.000]
lower                    (20, 5, 12)  sJ [+1.000, +1.000]
chan_upper               (20, 5, 8)  sJ [+0.479, +0.973]
chan_body_top            (27, 5, 8)  sJ [+0.710, +0.996]
chan_body_front          (28, 5, 8)  sJ [+0.443, +0.999]
chan_body_bot            (27, 5, 8)  sJ [+0.710, +0.996]
chan_lower               (20, 5, 8)  sJ [+0.479, +0.973]
chan_outflow_top_upper   (8, 5, 8)  sJ [+1.000, +1.000]
chan_outflow_wake_upper  (8, 5, 12)  sJ [+1.000, +1.000]
chan_outflow_wake_lower  (8, 5, 12)  sJ [+1.000, +1.000]
chan_outflow_top_lower   (8, 5, 8)  sJ [+1.000, +1.000]
far-field box [-2.   -2.5   0.01] .. [20.    2.5   0.49]
```

Twelve blocks for a sharp TE (14 for a blunt one): the 3-block wrap, five
channels on the wrap's outer face (`chan_upper` .. `chan_lower`, sharing the
wrap's `k_hi` nodes) and four outflow blocks that close the box behind the
wake. `--nk-channel` (8) is the channel depth including the shared layer;
by default each channel's first step equals the wrap's last marching step
(`--c1-seam-channel`). Every block is right-handed: the upper outflow
blocks run j against their partners (ADR 0057). The box spans z = 0.01 to 0.49: the span of the slices,
not of the STL.

<Figure name="ch-hybrid" caption="CH hybrid: the C-grid wrap and the channel blocks to the far-field box, span slice." />

## 5. Look at it

```bash
shore export wing_*.geo -o wing.vtm
shore export ch_*.geo -o ch.vtm
```

```
✓ Wrote wing.vtm (3 blocks: wing_body, wing_lower, wing_upper)
```

## What you learned

- `shore profile naca` + `shore body extrude` build a test wing; the STL
  needs more span stations than the grid.
- `--topology cgrid` builds 3 blocks for a sharp TE and needs
  `--cgrid-blocks 4` for a blunt one.
- The tip caps (`--cap-tips`, `shore cap-tip`) are refused until their
  redesign: their blocks folded (ADR 0058).
- `--topology ch` extends the C-grid to a far-field box with non-marched
  channels.

## Next

- [3. A Chimera assembly](/tutorials/chimera-assembly): a near-body grid
  inside a background grid.
- [C-grid topology](/concepts/topologies/c-grid),
  [CH hybrid](/concepts/topologies/ch-hybrid),
  [C-grid construction](/concepts/c-grid-construction),
  [Airfoil profiles](/concepts/airfoil-profiles).
- [`shore mesh`](/reference/mesh), [`shore cap-tip`](/reference/cap-tip),
  [`shore profile`](/reference/profile), [`shore body`](/reference/body).
