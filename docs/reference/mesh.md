# shore mesh

Full pipeline: STL → structured volume → geo files.

Five topologies are supported, selected via `--topology`:

- **[`cubed_sphere`](/topologies/cubed-sphere)** (default): 6-block
  decomposition (4 equatorial sub-blocks + 2 polar caps), gnomonic
  flat-square cap k=0, C1 cap–equator seam pin, seam-aware normals.
  Writes six `<stem>_<label>.geo` files.
- **[`ogrid`](/topologies/o-grid)**: single-block lat-lon grid,
  j-periodic. Writes one `<stem>.geo` file. Polar singular axis at
  `i = 0` and `i = ni - 1`.
- **[`cgrid`](/topologies/c-grid)**: 3-block (sharp TE) or 4-block
  (blunt TE) C-grid for airfoils / hydrofoils with downstream wake
  branches. Optional 1-block legacy layout via `--cgrid-blocks 1`.
  Writes one `<stem>_<label>.geo` per block (or one `<stem>.geo` for
  the legacy 1-block). Requires C-grid wrap flags (see below).
- **[`ch`](/topologies/ch-hybrid)**: 12-block (sharp TE) or 14-block
  (blunt TE) hybrid: C-grid wrap + H-channel out to a chord/normal-
  aligned far-field frame. Requires both C-grid wrap flags and frame
  extents (`--ch-chord-min/max`, `--ch-normal-min/max`).
- **[`oh`](/topologies/oh-hybrid)**: 12-block hybrid: cubed-sphere
  wrap + 6-block H-channel out to an axis-aligned far-field box.
  Requires the box extents (`--oh-box-{x,y,z}{min,max}`) and reuses
  the shared channel flags. Star-shaped bodies only; `--theta-cap` is
  ignored (the wrap pins it to `pi/4`).

## Synopsis

```
shore mesh [OPTIONS] [INPUT]
```

`INPUT` is optional when `--config` supplies it. CLI flags override
config-file values; config-file values override built-in defaults.

## Arguments

| Argument | Description |
|---|---|
| `INPUT` | Path to the input `.stl` file (binary or ASCII). Optional when `--config` is given and the config provides `input`. |

## Options

### Grid resolution

`--ni`, `--nj`, `--nk` are **node counts** — the array shape of the
`(ni, nj, nk, 3)` coordinate block SHORE produces. A block with
`--nk 30` has 30 wall-normal layers including the wall row, i.e. 29
cells extruded from the wall. The only file SHORE writes that uses
cell counts is the binary `.grd` (geogrd convention); the conversion
is automatic.

| Option | Default | Description |
|---|---|---|
| `--ni INTEGER` | 40 | Equator meridional **node count** (south seam → north seam). |
| `--nj INTEGER` | 60 | Longitude **node count**; multiple of 4 (rounded down with a warning if not). |
| `--nk INTEGER` | 30 | Wall-normal **node count**, including the wall row at k=0. With `--nk 30` the marcher extrudes 29 cells from the wall. |

### Cap geometry

| Option | Default | Description |
|---|---|---|
| `--theta-cap FLOAT` | $\pi/6$ | Polar exclusion angle in radians (cap corners sit at this colatitude) |
| `--theta-cap-deg FLOAT` | — | Same as `--theta-cap` but in degrees (mutually exclusive) |

::: tip Cap / equator compatibility
Choose `ni`, `nj`, `theta_cap` so the cap radial step
`theta_cap / (nc − 1)` (with `nc = nj // 4 + 1`) is within ~2× of the
equator's natural radial step `(π − 2*theta_cap) / (ni − 1)`. This keeps
the per-j C1 seam pinning active. When out of range, the equator falls
back to its unpinned `i_spacing` law and emits a single `UserWarning`.
:::

### Surface spacing

| Option | Default | Description |
|---|---|---|
| `--surface-i-spacing TEXT` | `uniform` | Equator meridional spacing: `'uniform'`, `'tanh'`, `'tanh2'` |
| `--surface-i-beta FLOAT` | 3.0 | Clustering strength for tanh/tanh2 i-spacing |
| `--surface-j-spacing TEXT` | `uniform` | Cap-edge tangent-plane spacing: `'uniform'` or `'tanh2'`; inherited by the equator seam |
| `--surface-j-beta FLOAT` | 3.0 | Clustering strength for tanh2 j-spacing |

### Volume spacing

| Option | Default | Description |
|---|---|---|
| `--volume-k-spacing TEXT` | `geometric` | Wall-normal spacing: `'geometric'`, `'tanh'`, `'tanh2'` |
| `--volume-k-beta FLOAT` | 3.0 | Clustering strength for tanh/tanh2 |
| `--volume-k-thickness FLOAT` | — | Total wall-normal thickness (required for tanh/tanh2) |
| `--ds FLOAT` | 0.001 | First-layer thickness (geometric spacing) |
| `--growth FLOAT` | 1.15 | Geometric growth ratio |

### Smoothing & normals

| Option | Default | Description |
|---|---|---|
| `--smooth FLOAT` | 0.2 | Laplacian blend strength per layer [0–1] |
| `--smooth-iters INTEGER` | 2 | Laplacian sweeps applied after each layer |
| `--blend-normals-k INTEGER` | 0 | Number of layers over which to blend STL surface normals into the grid-computed normals (0 = off; useful for sharp-feature bodies) |

### Topology

| Option | Default | Description |
|---|---|---|
| `--topology TEXT` | `cubed_sphere` | One of `cubed_sphere` (6 blocks), `ogrid` (1 block, j-periodic), `cgrid` (3-/4-/1-block C-grid for airfoils — see [C-grid](#c-grid-options-topology-cgrid-or-ch) options), `ch` (12-/14-block CH hybrid — see [CH hybrid](#ch-hybrid-options-topology-ch) options), or `oh` (12-block cubed-sphere + H-channel hybrid — see [OH hybrid](#oh-hybrid-options-topology-oh) options). |

### C-grid options (`--topology cgrid` or `ch`)

These flags are required for both `cgrid` and `ch` topologies. They
configure the C-grid wrap geometry around an airfoil / hydrofoil.

| Option | Default | Description |
|---|---|---|
| `--ni-body INTEGER` | — | Per-station body wrap node count (upper-TE → LE → lower-TE). |
| `--ni-wake INTEGER` | — | Per-branch wake node count (TE node included). |
| `--n-stations INTEGER` | — | Span-direction station count. |
| `--wake-length FLOAT` | — | Length of each wake branch downstream of TE (body units). |
| `--wake-growth FLOAT` | 1.15 | Geometric growth ratio for wake i-spacing. |
| `--span-axis TEXT` | `z` | `x` / `y` / `z` — span world axis. |
| `--te-dihedral-deg FLOAT` | 30.0 | Sharpness threshold for TE detection. |
| `--cgrid-blocks INTEGER` | 3 | `cgrid` only: 3 (Xall-correct sharp-TE), 4 (Xall-correct blunt-TE), or 1 (legacy single-block). Ignored for `--topology ch`. |
| `--cap-tips TEXT` | `none` | `cgrid` only: `none` / `jlo` / `jhi` / `both` — wing-tip cap blocks at the requested spanwise extreme(s). See [`shore cap-tip`](cap-tip) for a standalone CLI. |
| `--cap-te-thickness*`, `--cap-le-frac`, `--cap-te-frac`, `--cap-tip-inset`, `--cap-n-cap-j`, `--cap-j-spacing-cgrid`, `--cap-j-beta-cgrid`, `--cap-nk` | various | `cgrid + caps` only — cap-block geometry. See [C-grid topology](/topologies/c-grid#tip-caps-closed-finite-wing). |

### CH hybrid options (`--topology ch`)

In addition to the C-grid wrap flags above, CH requires the frame
extents in world coordinates and exposes channel knobs.

| Option | Default | Description |
|---|---|---|
| `--ch-chord-min FLOAT` | — | Far-field box minimum chord coordinate. |
| `--ch-chord-max FLOAT` | — | Far-field box maximum chord coordinate. |
| `--ch-normal-min FLOAT` | — | Far-field box minimum normal coordinate. |
| `--ch-normal-max FLOAT` | — | Far-field box maximum normal coordinate. |
| `--nk-channel INTEGER` | 5 | Total channel-block layers (incl. shared seam at k=0). |
| `--k-law-channel TEXT` | `uniform` | Channel k-axis spacing law (`uniform` / `tanh` / `tanh2`). |
| `--k-beta-channel FLOAT` | 3.0 | Clustering strength for tanh / tanh2 channel spacing. |
| `--c1-seam-channel/--no-c1-seam-channel` | on | Pin each k_hi-side channel's first step to the wrap's last marching step (per `(i, j)`). |
| `--ch-freestream-x/y/z FLOAT` | auto | Freestream direction components (world). Selects the chord axis; defaults to the largest non-span component when unset. |
| `--ch-i-split-body-upper/lower INTEGER` | auto | Body-wrap i index where the upper / lower channel segment splits. Auto-placed at frame corners when unset. |
| `--ch-ni-outflow INTEGER` | 8 | Outflow sub-block i resolution. |

### OH hybrid options (`--topology oh`)

OH requires an axis-aligned far-field box (six world-coordinate
extents) and reuses the shared channel knobs (`--nk-channel`,
`--k-law-channel`, `--k-beta-channel`, `--c1-seam-channel`). The wrap
is star-shaped (centroid ray-cast); `--theta-cap` is ignored (the
wrap pins it to `pi/4`, emitting a `UserWarning` if set).

| Option | Default | Description |
|---|---|---|
| `--oh-box-xmin/xmax FLOAT` | — | Far-field box x extents (world). |
| `--oh-box-ymin/ymax FLOAT` | — | Far-field box y extents (world). |
| `--oh-box-zmin/zmax FLOAT` | — | Far-field box z extents (world). Box must contain the marched wrap envelope. |

### I/O and display

| Option | Default | Description |
|---|---|---|
| `-c` / `--config PATH` | — | Read parameter values from a TOML file (see [Config files](config)). Explicit CLI flags still override the file. |
| `-o` / `--output PATH` | `<INPUT stem>` | Output stem. Multi-block topologies (`cubed_sphere`, `cgrid`, `ch`, `oh`) write one `<stem>_<label>.geo` per block; `ogrid` and `cgrid --cgrid-blocks 1` write one `<stem>.geo`. |
| `--adjacency-json PATH` | — | Also write a `.adjacency.json` sidecar describing the SHARED / DIRICHLET face graph (required by [`shore cc-par`](cc-par)). |
| `-v` / `--verbose` | off | Show Rich progress bar with per-layer stats. |
| `--help` | — | Show this help and exit. |

### Deprecated

| Option | Status |
|---|---|
| `--cap-pull FLOAT` | Ignored. Cap interior is now built by gnomonic flat-square projection — square cells by construction. |
| `--cap-squareness FLOAT` | Ignored. Same reason. |
| `--j-cluster-lon FLOAT` | Ignored. The gnomonic cap has 4-fold symmetry; use `--surface-j-spacing tanh2` to cluster nodes near the cap-corner longitudes. |

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Input file not found, projection failed, or negative Jacobian during march |

## Output

The exact file set depends on `--topology`:

**`cubed_sphere`** (default) — 6 `.geo` files:

```
<stem>_sub0.geo        sub-block 0 (phi 0°→90°)         (ni × nc × nk)
<stem>_sub1.geo        sub-block 1 (phi 90°→180°)       (ni × nc × nk)
<stem>_sub2.geo        sub-block 2 (phi 180°→270°)      (ni × nc × nk)
<stem>_sub3.geo        sub-block 3 (phi 270°→360°)      (ni × nc × nk)
<stem>_cap_north.geo   north polar cap                  (nc × nc × nk)
<stem>_cap_south.geo   south polar cap                  (nc × nc × nk)
```

with `nc = nj // 4 + 1`.

**`ogrid`** — 1 `.geo` file: `<stem>.geo` of shape `(ni × nj × nk)`.

**`cgrid`** with `--cgrid-blocks 3` (default sharp-TE) — 3 `.geo`
files: `<stem>_upper.geo`, `<stem>_body.geo`, `<stem>_lower.geo`.
`--cgrid-blocks 4` (blunt-TE) adds `<stem>_wake_base.geo`.
`--cgrid-blocks 1` (legacy) writes a single `<stem>.geo`. The
`--cap-tips` option adds 5–6 cap blocks per requested tip (see
[C-grid](/topologies/c-grid#tip-caps-closed-finite-wing)).

**`ch`** — 12 `.geo` files (sharp TE) or 14 (blunt TE). See
[CH hybrid block decomposition](/topologies/ch-hybrid#block-decomposition)
for the full list.

**`oh`** — 12 `.geo` files: the 6 cubed-sphere wrap blocks
(`<stem>_sub0..sub3`, `<stem>_cap_north`, `<stem>_cap_south`) plus
6 channel blocks (`<stem>_chan_sub0..chan_sub3`,
`<stem>_chan_north`, `<stem>_chan_south`).

## Examples

Minimal invocation — defaults for all grid parameters:

```bash
shore mesh sphere.stl
# → sphere_sub*.geo, sphere_cap_*.geo  (40 x 60 x 30, theta_cap = pi/6)
```

High-resolution boundary layer with verbose output:

```bash
shore mesh fuselage.stl \
  --ni 80 --nj 120 --nk 50 \
  --ds 1e-5 --growth 1.12 --smooth 0.3 \
  -o fuselage_bl -v
```

Tanh-spaced boundary layer with fixed total thickness:

```bash
shore mesh wing.stl \
  --ni 60 --nj 80 --nk 40 \
  --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 2.0 \
  -o wing_tanh
```

Equator meridional clustering at the cap seams:

```bash
shore mesh body.stl \
  --ni 60 --nj 80 --nk 30 \
  --surface-i-spacing tanh --surface-i-beta 4.0 \
  --theta-cap-deg 30 \
  -o body
```

Cap-corner refinement (j-spacing tanh2):

```bash
shore mesh body.stl \
  --ni 40 --nj 60 --nk 30 \
  --surface-j-spacing tanh2 --surface-j-beta 4.0 \
  --theta-cap-deg 30 \
  -o body
```

STL-normal blending for sharp features (first 5 layers blend toward the
true surface normal of the underlying triangulation):

```bash
shore mesh airfoil.stl \
  --ni 40 --nj 60 --nk 30 \
  --blend-normals-k 5 \
  -o airfoil
```

Single-block O-grid (one `.geo` file, j-periodic):

```bash
shore mesh sphere.stl --topology ogrid \
  --ni 40 --nj 60 --nk 30 --ds 0.02 --growth 1.1 \
  -o sphere
# → sphere.geo  (40 x 60 x 30, j-periodic)
```

C-grid for an airfoil (3-block sharp-TE default, 3 `.geo` files):

```bash
shore mesh airfoil.stl --topology cgrid \
  --ni-body 80 --ni-wake 20 --n-stations 5 --wake-length 15.0 \
  --nk 12 --ds 5e-4 --growth 1.08 \
  -o airfoil
# → airfoil_upper.geo, airfoil_body.geo, airfoil_lower.geo
```

CH hybrid (C-grid wrap + far-field H-channel; 12 `.geo` files for sharp TE):

```bash
shore mesh airfoil.stl --topology ch \
  --ni-body 80 --ni-wake 20 --n-stations 5 --wake-length 15.0 \
  --nk 12 --ds 5e-4 --growth 1.08 \
  --ch-chord-min -2.0 --ch-chord-max 20.0 \
  --ch-normal-min -2.5 --ch-normal-max 2.5 \
  --nk-channel 8 \
  -o ch
# → ch_upper.geo, ch_body.geo, ch_lower.geo,
#   ch_chan_upper.geo, ch_chan_body_top.geo, ch_chan_body_front.geo,
#   ch_chan_body_bot.geo, ch_chan_lower.geo,
#   ch_chan_outflow_top_upper.geo, ch_chan_outflow_wake_upper.geo,
#   ch_chan_outflow_wake_lower.geo, ch_chan_outflow_top_lower.geo
```

Drive the full pipeline from a TOML config (see [Config files](config)):

```bash
shore mesh -c fuselage.toml          # all parameters from the file
shore mesh -c fuselage.toml --nk 80  # CLI --nk wins over fuselage.toml
```

## Notes

- The STL must be a single connected, watertight mesh. Multi-body STL files
  (multiple shells) will fail during projection — split them first with a
  CAD tool or trimesh.
- For `cubed_sphere`, `--nj` must be a multiple of 4 (otherwise rounded
  down with a warning, so the equatorial band can be split into 4
  equal sectors).
- For `cgrid` / `ch`, the body STL must be a finite-span wing-like mesh
  with a detectable trailing-edge polyline (sharp or blunt). See
  [C-grid construction](/algorithm/c-grid) for the TE-detection
  algorithm.
- For `ch`, the frame extents (`--ch-chord-min/max` /
  `--ch-normal-min/max`) must contain the marched wrap envelope.
  Validation runs at `fill_channel_blocks` time and raises a
  `ParameterError` on violation.
- Use `shore check` to verify the per-block Jacobians.
- See [Spacing laws](/algorithm/spacing-laws) for a detailed explanation of each
  distribution and when to use it.

## See also

- [Config files](config) — drive `shore mesh` from a TOML file (`-c`)
- [`shore project`](project) — generate the k=0 surface only (single-block,
  legacy `body.project` path; not the cubed-sphere topology)
- [`shore check`](check) — validate Jacobians on existing `.geo` files
- [`shore export`](export) — convert `.geo` to VTK MultiBlock for ParaView
- [Cubed-sphere topology](/topologies/cubed-sphere) — block layout, cap construction,
  seam matching
- [O-grid topology](/topologies/o-grid) — single-block alternative
- [C-grid topology](/topologies/c-grid) — sharp-/blunt-TE wrap for airfoils
- [CH hybrid topology](/topologies/ch-hybrid) — C-grid wrap + far-field channel
- [OH hybrid topology](/topologies/oh-hybrid) — cubed-sphere wrap + 6-face H-channel (`--topology oh`)
- [Spacing laws](/algorithm/spacing-laws) — clustering distributions
- [Hyperbolic marching](/algorithm/hyperbolic-marching) — algorithm details
