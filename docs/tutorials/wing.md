---
title: A wing with a C-grid
---

# A wing with a C-grid

You build a NACA 0012 wing from its section, wrap it in a C-grid with wake
branches, close its tips with cap blocks, and extend the C-grid to a
far-field box with the CH hybrid. You learn how sharp and blunt trailing
edges change the block layout, the two ways to cap a tip, and how to read
`shore check` on blocks that are mirrored rather than folded.

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

The C-grid is open at the two span ends. SHORE closes a tip in one of two
ways:

| | Conformal caps | Chimera cap |
|---|---|---|
| Command | `shore mesh --topology cgrid --cap-tips jlo\|jhi\|both` | `shore cap-tip <stem>_body.geo` |
| Blocks per tip | 5 cap blocks + `cap_<tip>_wake_base` | 5 (`cap_le`, `cap_upper`, `cap_lower`, `cap_te_upper`, `cap_te_lower`) |
| C-grid | promoted to 4 blocks, TE rounded over the whole span | unchanged |
| Cap–body interface | nodes coincide bit for bit | overset: the cap's `k_hi` is FREE (a fringe) |

### Conformal caps

::: code-group
```bash [CLI]
shore mesh wing.stl --topology cgrid --ni-body 120 --ni-wake 30 --n-stations 5 \
    --wake-length 15 --nk 20 --ds 5e-4 --growth 1.15 --smooth 0.1 --smooth-iters 1 \
    --cap-tips both --cap-tip-inset 0.02 -o capped -v
```
<<< @/../examples/tutorials/wing.py#capped{python} [Python]
:::

```
CGridTopology: cap_tips='both' on a sharp-TE body — auto-promoting blocks=3 to blocks=4
and rounding the body's TE by te_thickness=0.0019952 so the cap and body meshes match
conformally.
✓ Wrote 16 .geo files (upper, body, lower, wake_base, cap_jlo_le, cap_jlo_upper,
cap_jlo_lower, cap_jlo_te_upper, cap_jlo_te_lower, cap_jlo_wake_base, cap_jhi_le, ...)
```

The TE thickness is `--cap-te-thickness-frac` (default 0.002) times the
chord. `--cap-le-frac` and `--cap-te-frac` (0.05 and 0.95) place the cuts
that split the cap into its LE, mid and TE parts; `--cap-n-cap-j` (12) is the
node count across the section thickness. The caps are built, not marched.

### A chimera cap

`shore cap-tip` reads the body block of a C-grid already on disk and builds
an independent cap for one tip. It rounds its own copy of a sharp TE
(`--te-thickness-frac`, default 0.002 chord); the wrap is not touched. A
blunt body block (from `--cgrid-blocks 4`) is capped without rounding.

::: code-group
```bash [CLI]
shore cap-tip wing_body.geo --tip jlo -o cap_jlo --adjacency-json cap_jlo.adjacency.json
```
<<< @/../examples/tutorials/wing.py#tipcap{python} [Python]
:::

```
CGridTipCap: sharp TE detected; auto-rounding cap with te_thickness=0.0019952 (the body wrap is unchanged).
auto tip_inset = 0.000997602
auto nk = 20 (matches body)
wrote cap_jlo_cap_le.geo
wrote cap_jlo_cap_upper.geo
...
```

The cap's `k_lo` face is the tip wall and its `k_hi` face a chimera fringe
that receives from the wrap: `shore cc-par` on its sidecar writes it as an
overset component ([3. A Chimera assembly](/tutorials/chimera-assembly)).

<Figure name="cgrid-tip-cap" caption="Tip cap blocks over the C-grid wall." />

### Reading `shore check` on the caps

`shore check` reports some cap blocks as inverted:

```
shore check cap_jlo_cap_lower.geo
...
✗ Inverted cells detected.
```

The scaled Jacobians the script prints tell two different stories:

```
cap_jlo_le               (10, 12, 20)  sJ [-1.000, -0.000]
cap_jlo_upper            (54, 12, 20)  sJ [+0.973, +1.000]
cap_jlo_lower            (54, 12, 20)  sJ [-1.000, -0.973]
cap_jlo_te_upper         (3, 12, 20)  sJ [+0.990, +1.000]
cap_jlo_te_lower         (3, 12, 20)  sJ [-1.000, -0.990]
cap_jlo_wake_base        (30, 3, 20)  sJ [-1.000, -1.000]
...
cap_jhi_le               (10, 12, 20)  sJ [-1.000, +0.995]
```

- A block whose every cell is negative, such as `cap_jlo_lower`
  (−1.000 to −0.973), is not folded: it is left-handed, the mirror image of
  `cap_jlo_upper` with the same index directions. Its cells are as good as
  its partner's.
- The `le` blocks are different. In `cap_jlo_le` 190 cells of the middle
  i-plane have a zero scaled Jacobian (the block collapses to a line there),
  and `cap_jhi_le` mixes signs: 1672 cells positive, 190 zero, 19 negative.

::: warning
The LE cap blocks contain degenerate cells, and the left-handed blocks are
reported as inverted by `shore check`. Check that your solver accepts both
before you run on a capped wing.
:::

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
chan_outflow_top_upper   (8, 5, 8)  sJ [-1.000, -1.000]
chan_outflow_wake_upper  (8, 5, 12)  sJ [-1.000, -1.000]
chan_outflow_wake_lower  (8, 5, 12)  sJ [+1.000, +1.000]
chan_outflow_top_lower   (8, 5, 8)  sJ [+1.000, +1.000]
far-field box [-2.   -2.5   0.01] .. [20.    2.5   0.49]
```

Twelve blocks for a sharp TE (14 for a blunt one): the 3-block wrap, five
channels on the wrap's outer face (`chan_upper` .. `chan_lower`, sharing the
wrap's `k_hi` nodes) and four outflow blocks that close the box behind the
wake. `--nk-channel` (8) is the channel depth including the shared layer;
by default each channel's first step equals the wrap's last marching step
(`--c1-seam-channel`). The two `*_upper` outflow blocks are left-handed
(scaled Jacobian −1 everywhere) and `shore check` reports them as inverted,
as for the caps. The box spans z = 0.01 to 0.49: the span of the slices,
not of the STL.

<Figure name="ch-hybrid" caption="CH hybrid: the C-grid wrap and the channel blocks to the far-field box, span slice." />

## 5. Look at it

```bash
shore export wing_*.geo -o wing.vtm
shore export capped_*.geo -o capped.vtm
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
- Conformal caps (`--cap-tips`) share nodes with a rounded-TE C-grid;
  `shore cap-tip` builds an overset cap from a body block.
- `--topology ch` extends the C-grid to a far-field box with non-marched
  channels.
- A block with every scaled Jacobian at −1 is mirrored, not folded; the LE
  cap blocks have degenerate cells.

## Next

- [3. A Chimera assembly](/tutorials/chimera-assembly): a near-body grid
  inside a background grid.
- [C-grid topology](/concepts/topologies/c-grid),
  [CH hybrid](/concepts/topologies/ch-hybrid),
  [C-grid construction](/concepts/c-grid-construction),
  [Airfoil profiles](/concepts/airfoil-profiles).
- [`shore mesh`](/reference/mesh), [`shore cap-tip`](/reference/cap-tip),
  [`shore profile`](/reference/profile), [`shore body`](/reference/body).
