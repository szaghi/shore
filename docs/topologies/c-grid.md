# C-grid topology

A C-grid is a structured grid topology that **wraps a body with a
sharp trailing edge** (airfoil, hydrofoil, wing) and **extends two
wake branches downstream**. Looking at a span station, the grid
lines form the letter C: the outer envelope opens downstream of the
TE, leaving a wake-cut where the upper and lower surfaces meet.

The single-block C-grid is a **first-class topology** in SHORE: build
it with `shore mesh --topology cgrid` (or `topology = "cgrid"` in a
TOML config). It coexists with the [O-grid](o-grid) and
[cubed-sphere](cubed-sphere) topologies as a peer; the three pick
different trade-offs:

| | C-grid | O-grid | Cubed sphere |
|---|---|---|---|
| Blocks | 1 | 1 | 6 |
| Body shape | sharp-TE (airfoil) | star-shaped (sphere-like) | star-shaped |
| Wake handling | wake branches in i-direction | none | none |
| Periodicity | non-periodic in i, j | j-periodic | conforming SHARED seams |
| Best for | airfoils, hydrofoils | Chimera background, smooth bodies | production near-body meshing on smooth bodies |

## i-axis layout

The C-grid uses a **single continuous parametric line** along i,
sweeping from one wake far-field tip, around the body, to the other
wake far-field tip:

```
i = 0                                  i = ni_total - 1
upper-wake-far  →  upper-TE  →  suction  →  LE  →  pressure  →  lower-TE  →  lower-wake-far
                       i = ni_wake - 1                              i = ni_wake + ni_body - 1
```

`ni_total = 2 * (ni_wake - 1) + ni_body`, where `ni_body` is the
per-station body wrap node count (upper-TE → LE → lower-TE) and
`ni_wake` is the per-branch wake node count (the TE node is shared
with the body, so it is counted once but appears twice in the layout
above for clarity).

j is the **span direction** — one node per span station. k is
wall-normal as usual.

The TE sits at two i-indices: `i = ni_wake - 1` (upper-TE) and
`i = ni_wake + ni_body - 1` (lower-TE). For a **sharp** TE the two
TE nodes coincide in space; for a **blunt** TE they are separated by
the TE thickness.

## Building from the CLI

```bash
shore mesh naca0012.stl --topology cgrid \
    --ni-body 120 --ni-wake 30 --n-stations 5 --wake-length 15 \
    --nk 30 --ds 5e-4 --growth 1.15 \
    -o naca0012
```

Required flags (non-default for `cgrid` only):

- `--ni-body` — per-station body wrap node count.
- `--ni-wake` — per-branch wake node count (TE included).
- `--n-stations` — number of span-axis stations (becomes block `nj`).
- `--wake-length` — physical length of each wake branch from the TE,
  in body-coordinate units. A typical value is 15-20 chord lengths.

Optional flags:

- `--span-axis x|y|z` — defaults to `z`. The body's longest non-span
  axis is taken as the chord direction; the freestream defaults to
  `+x_world`-projected-onto-the-station-plane.
- `--wake-growth` — geometric growth ratio for wake i-spacing
  (default `1.15`). Wake cells start at the body's TE-cell length and
  grow toward the far field.
- `--te-dihedral-deg` — sharpness threshold for TE detection (default
  `30°`). Lower for blunter TEs (e.g. `20°`); higher to reject
  spurious tessellation seams.

The output is one `<stem>.geo` file with shape `(ni_total, n_stations, nk, 3)`.

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
)
m = Mesh(body=body, blocks=blocks, topology=topo)
m.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)
m.write_geo("naca0012")  # → naca0012_cgrid.geo
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

1. **TE detection** ([`detect_te_edges`](/algorithm/c-grid#trailing-edge-detection))
   — find the TE polyline on the STL via dihedral classification,
   span-axis alignment filter, and chord-extreme filter. Sharp TEs
   produce one strand; blunt TEs produce two (upper + lower).

2. **Span slicing** ([`slice_by_span`](/algorithm/c-grid#span-slicing))
   — slice the STL with `n_stations` planes perpendicular to the
   span axis. Each slice is stitched into a closed polyline, rotated
   so it starts at the upper-TE and walks the suction surface first,
   and **resampled by arclength fraction** to a common `ni_body`
   node count. This common parametrisation is what makes
   tapered / swept wings line up at the same parametric i-position
   on every station.

3. **C-curve construction** ([`build_c_curve`](/algorithm/c-grid#c-curve))
   — append two wake branches per station, each extending in the
   freestream direction with geometric i-spacing rooted at the
   body's TE-cell length.

4. **Hyperbolic march** — standard SHORE hyperbolic extrusion with
   both axes non-periodic. The wake far-field tips and wing tips
   are clamped during Laplacian smoothing.

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

The example at
[`examples/06-c-grid/visualize_naca0012.py`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012.py)
runs end-to-end on a NACA 0012 fixture without external data.

## Face boundary conditions

| Face | Default BC | Meaning |
|---|---|---|
| `i_lo` | `FREE` | upper-wake far-field tip (chimera fringe / extrapolation) |
| `i_hi` | `FREE` | lower-wake far-field tip |
| `j_lo` | `FREE` | inboard wing tip |
| `j_hi` | `FREE` | outboard wing tip |
| `k_lo` | `WALL` | body wrap (over the body i-range; degenerate over the wake i-range — the wake-cut sits at k=0) |
| `k_hi` | `FREE` | outer / chimera fringe |

For a chimera assembly, replace `i_lo` / `i_hi` / `j_lo` / `j_hi` /
`k_hi` with a downstream `cc.par` `free_face_bc=40` (chimera fringe)
or `-9` (extrapolation). For symmetric (mirrored) studies, the
`j_lo` / `j_hi` faces become symmetry walls externally.

## See also

- [Algorithm — C-grid construction](/algorithm/c-grid) — full derivation
  of TE detection, span slicing, C-curve build.
- [`shore mesh` reference](/reference/mesh) — full CLI command reference.
- [Pipeline overview](/guide/pipeline) — where the C-grid fits in the
  body-fitted workflow.
