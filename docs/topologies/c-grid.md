# C-grid topology

A C-grid is a structured grid topology that **wraps a body with a
sharp or blunt trailing edge** (airfoil, hydrofoil, wing) and
**extends two wake branches downstream**. Looking at a span station,
the grid lines form the letter C: the outer envelope opens
downstream of the TE, leaving a wake-cut where the upper and lower
surfaces meet.

The C-grid is a **first-class topology** in SHORE: build it with
`shore mesh --topology cgrid` (or `topology = "cgrid"` in a TOML
config). It coexists with the [O-grid](o-grid) and
[cubed-sphere](cubed-sphere) topologies as a peer.

## Block layouts

`shore mesh --topology cgrid` supports three layouts via
`--cgrid-blocks`:

| Layout | Blocks | TE | Use case |
|---|---|---|---|
| **3-block** (default) | upper / body / lower | sharp | The Xall-correct sharp-TE C-grid; SHARED seams at the upper-TE column, lower-TE column, and wake cut |
| **4-block** | upper / body / lower / wake_base | blunt | 3-block + a wake-base block W bridging the TE thickness gap; W's `i_hi` face is the bluff-base WALL |
| **1-block** (legacy) | one concatenated block | sharp only | Single hex block with the wake cut encoded as an i-direction overlap; **Xall cannot apply BCs cleanly** on this overlap (the j-counter "flips direction" along the self-overlapping edge).  Kept for opt-in compatibility with downstream tools that don't have the j-flip issue. |

::: warning Migration from v0.1.0
v0.1.0 shipped the 1-block layout as the default.  v0.2.0 makes
**3-block the default**.  If you have a v0.1.0 workflow that relies
on the single-block layout, pass `--cgrid-blocks 1` to opt in.  See
[the v0.2.0 changelog](/guide/changelog) for the migration note.
:::

## 3-block layout (sharp TE, default)

The 3-block layout splits the C-grid into three hex blocks:

```
                      BLOCK upper (upper wake)
upper-wake-far-tip ──────────────────────── upper-TE
                                                │
                                                │ BLOCK body (suction → LE → pressure)
                                                │
                                            lower-TE
lower-wake-far-tip ──────────────────────── lower-TE
                      BLOCK lower (lower wake)
```

Per-block dimensions:

- **`upper`** `(ni_wake, n_stations, nk, 3)` — i sweeps `upper-wake-far-tip → upper-TE`.
- **`body`** `(ni_body, n_stations, nk, 3)` — i sweeps `upper-TE → suction → LE → pressure → lower-TE`.
- **`lower`** `(ni_wake, n_stations, nk, 3)` — i sweeps `lower-TE → lower-wake-far-tip`.

j is the span direction; k is wall-normal extrusion.

### SHARED seams

Three SHARED face pairs:

| Seam | Faces | AxisMap (this side) | Notes |
|---|---|---|---|
| Upper-TE column | `upper.i_hi ↔ body.i_lo` | identity on `(j, k)` | Same-axis i↔i; node coincidence via numpy view aliasing. |
| Lower-TE column | `body.i_hi ↔ lower.i_lo` | identity on `(j, k)` | Same as above. |
| Wake cut       | `upper.k_lo ↔ lower.k_lo` | `((i, True), (j, False))` | Same-axis k↔k; **i flipped** because upper sweeps far-tip→TE while lower sweeps TE→far-tip.  This i-reversal is the manifestation of what was the "j-counter flips direction" pathology in the 1-block layout's self-overlap. |

For sharp TE the wake-cut nodes coincide bit-exactly in space at
k=0 — the SHARED face is a true memory-aliased view.  Xall sees a
clean face-to-face adjacency on all three seams and applies BCs
correctly.

### Face boundary conditions (defaults)

| Block | i_lo | i_hi | j_lo | j_hi | k_lo | k_hi |
|---|---|---|---|---|---|---|
| `upper` | `FREE` (far-tip) | `SHARED → body.i_lo` | `FREE` | `FREE` | `SHARED → lower.k_lo` | `FREE` |
| `body`  | `SHARED → upper.i_hi` | `SHARED → lower.i_lo` | `FREE` | `FREE` | `WALL` (airfoil surface) | `FREE` |
| `lower` | `SHARED → body.i_hi` | `FREE` (far-tip) | `FREE` | `FREE` | `SHARED → upper.k_lo` | `FREE` |

For a chimera assembly, the `FREE` faces become BC code `40`
(chimera fringe) or `-9` (extrapolation) in a downstream `cc.par`.

## 4-block layout (blunt TE)

For a body with a blunt (finite-thickness) TE, the upper and lower
wake k=0 strips are **parallel but separated by the TE thickness**.
The 3-block layout's wake-cut SHARED face would have a non-zero gap
between U.k_lo and L.k_lo — invalid as SHARED — and the TE base
wall of the airfoil would not be meshed.

The 4-block layout adds a **wake-base block W** that fills the gap:

```
y axis up
   ▲    upper-wake-far-tip ────── upper-TE     (Block upper, k_lo at y = +te/2)
   │                                  │
   │                                  │ Block body
   │                                  │
   │                                  │
   │       far-tip W ────── TE base wall      (Block wake_base)
   │                                  │
   │    lower-wake-far-tip ────── lower-TE    (Block lower, k_lo at y = -te/2)
   │
```

- **`wake_base`** `(ni_wake, n_stations, nk_w, 3)` — i along the
  wake (far-tip → TE), j is span, k spans the TE thickness from
  `y = -te/2` to `y = +te/2`.

Two extra SHARED seams (replacing the 3-block wake-cut):

| Seam | Faces | AxisMap (W side) |
|---|---|---|
| W upper face ↔ upper k_lo | `wake_base.k_hi ↔ upper.k_lo` | identity on `(i, j)` (W's i and upper's i match) |
| W lower face ↔ lower k_lo | `wake_base.k_lo ↔ lower.k_lo` | `((i, True), (j, False))` (W's i and lower's i are reversed) |

W's `i_hi` face — the cross-section at the TE — is **`WALL`**: this
is the bluff-base wall of the airfoil that the 3-block layout
cannot represent.  The 1-block legacy layout simply ignores it.

W is not marched; its geometry is fully determined at build time by
linear interpolation in k between the lower- and upper-wake k=0
strips.  The marcher's per-block loop skips it.

## Building from the CLI

```bash
# Sharp-TE 3-block (default):
shore mesh naca0012.stl --topology cgrid \
    --ni-body 120 --ni-wake 30 --n-stations 5 --wake-length 15 \
    --nk 30 --ds 5e-4 --growth 1.15 \
    -o naca0012
# → naca0012_upper.geo + naca0012_body.geo + naca0012_lower.geo

# Blunt-TE 4-block:
shore mesh naca0012_blunt.stl --topology cgrid --cgrid-blocks 4 \
    --ni-body 120 --ni-wake 30 --n-stations 5 --wake-length 15 \
    --nk 30 --ds 5e-4 --growth 1.15 \
    -o naca0012_blunt
# → naca0012_blunt_{upper,body,lower,wake_base}.geo

# Legacy 1-block (sharp only):
shore mesh naca0012.stl --topology cgrid --cgrid-blocks 1 \
    --ni-body 120 --ni-wake 30 --n-stations 5 --wake-length 15 \
    --nk 30 --ds 5e-4 --growth 1.15 \
    -o naca0012
# → naca0012.geo
```

Required flags (any layout):

- `--ni-body` — per-station body wrap node count.
- `--ni-wake` — per-branch wake node count (TE included).
- `--n-stations` — number of span-axis stations (becomes block `nj`).
- `--wake-length` — physical length of each wake branch from the TE,
  in body-coordinate units. A typical value is 15-20 chord lengths.

Optional flags:

- `--cgrid-blocks {1,3,4}` — block layout (default `3`).
- `--span-axis x|y|z` — defaults to `z`. The body's longest non-span
  axis is taken as the chord direction; the freestream defaults to
  `+x_world`-projected-onto-the-station-plane.
- `--wake-growth` — geometric growth ratio for wake i-spacing
  (default `1.15`). Wake cells start at the body's TE-cell length and
  grow toward the far field.
- `--te-dihedral-deg` — sharpness threshold for TE detection (default
  `30°`). Lower for blunter TEs (e.g. `20°`); higher to reject
  spurious tessellation seams.

The output is one `<stem>_<label>.geo` file per block (or one
`<stem>.geo` for `--cgrid-blocks 1`).

## Python API

```python
from shore.mesh import Body, Mesh, CGridTopology, Spacing1D

body = Body.from_stl("naca0012.stl", skip_hygiene=True)
topo = CGridTopology()
blocks = topo.build(
    surface=None,
    nk=30,
    mesh=body.mesh,
    spacing_k=Spacing1D(law="geometric", ds=5e-4, growth=1.15),
    ni_body=120, ni_wake=30, n_stations=5,
    wake_length=15.0,
    span_axis="z",
    blocks=3,    # 3 (default), 4 (blunt-TE), or 1 (legacy)
)
m = Mesh(body=body, blocks=blocks, topology=topo)
m.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)
m.write_geo("naca0012")
# → naca0012_upper.geo, naca0012_body.geo, naca0012_lower.geo
```

`Mesh.from_stl()` does **not** support `CGridTopology` because the
C-grid takes different parameters (`ni_body` / `ni_wake` /
`n_stations` / `wake_length`) that don't fit the legacy
`(ni, nj, nk)` signature. Build the block list directly with
`CGridTopology().build(...)`, then construct the `Mesh` from the
resulting blocks as shown above.

## Pipeline (algorithm summary)

The C-grid k=0 layer is built in four stages — see
[Algorithm — C-grid construction](/algorithm/c-grid) for derivations.

1. **TE detection** (`detect_te_edges`) — find the TE polyline(s)
   via dihedral classification.  Sharp TEs produce one strand;
   blunt TEs produce two (upper + lower).
2. **Span slicing** (`slice_by_span`) — slice the STL with
   `n_stations` planes perpendicular to the span axis.  Each slice
   is reordered to start at the upper-TE and resampled by
   arclength fraction.  For blunt TEs the body wrap terminates at
   the lower-TE node, excluding the TE-base segment.
3. **C-curve construction** (`build_c_curve`) — append upper and
   lower wake branches per station with geometric i-spacing.
4. **Topology assembly** + **hyperbolic march** — split the
   per-station polyline into 1, 3, or 4 blocks (`per_block_layers`
   on `CGridK0Result`); wire SHARED seams; march each block layer
   by layer with seam-aware ghost rows at the upper-TE and
   lower-TE seams.

## Body requirements

- **STL must be watertight enough that
  `trimesh.intersections.mesh_plane` can stitch a closed cross
  section at every span station.** Multi-component slices
  (interior holes, internal cavities) are rejected.
- **TE must be approximately span-aligned**, within ~30° (controls
  the candidate-edge filter inside `detect_te_edges`). Highly swept
  wings need a coarser threshold or a manually specified TE.
- **The body's chord axis must dominate** — the longer of the two
  non-span axes is auto-picked as the chord direction; the
  freestream defaults to chord-axis +1.

Two end-to-end examples run on a NACA 0012 fixture without external
data:

- [`examples/06-c-grid/visualize_naca0012.py`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012.py)
  — Python script with optional PyVista viewer.
- [`examples/06-c-grid/visualize_naca0012.sh`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012.sh)
  — bash equivalent using only `shore` CLI calls (`shore profile naca`
  → `shore body extrude` → `shore mesh --topology cgrid` → `shore export`).

## See also

- [Algorithm — C-grid construction](/algorithm/c-grid) — full derivation
  of TE detection, span slicing, C-curve build, and the SHARED-seam
  wiring for 3-block / 4-block topologies.
- [`shore mesh` reference](/reference/mesh) — full CLI command reference.
- [Pipeline overview](/guide/pipeline) — where the C-grid fits in the
  body-fitted workflow.
