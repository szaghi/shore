# First-tentative assembly: sphere body in a cylindrical domain

This guide walks through SHORE's first end-to-end **assembly** example
— two structured meshes, generated independently, sharing the same
world-coordinate domain. Run via `examples/05-chimera-assembly/assembly_sphere_in_cylinder.py`.

The example is *first-tentative* in a deliberately limited sense:
SHORE produces both meshes, writes them to disk, and bundles them into
a single ParaView-friendly `.vtm` for visual inspection. **Hole
cutting, fringe identification, and donor search are explicitly the
solver's responsibility** — that's the boundary between what SHORE
ships today and what a full Chimera assembly framework would do, and
where it stays for now.

## What the example produces

Two meshes covering the same domain `[-10, 10]² × [-10, 50]`:

| Mesh | Topology | Blocks | Role |
|---|---|---|---|
| Near-body | [Cubed sphere](/topologies/cubed-sphere) | 6 | Wraps a unit-radius sphere body at the origin, marched outward to total wall-normal extent ~0.2 (outer shell at r ≈ 0.7). |
| Background | [Flat-caps cylinder](/topologies/primitives) | 5 | Fills the cylindrical domain (r ≤ 10, z ∈ [-10, 50]) with no central hole — central square + 4 trapezoidal sub-blocks. No Chimera fringe at the inner boundary. |

The two meshes overlap in the near-body region (r < ~0.7 around the
origin). Resolving that overlap is the downstream solver's job; SHORE
just generates the geometry.

```
                background flat-caps cylinder (z ∈ [-10, 50], r ≤ 10)
              ┌─────────────────────────────────────────────────────┐
              │                                                     │
              │              near-body cubed sphere                 │
              │              (6 blocks, r ∈ [0.5, 0.7])             │
              │                       ╱──╲                          │
              │                      ╱    ╲                         │
              │                     ╲      ╱                        │
              │                      ╲────╱                         │
              │                                                     │
              │                                                     │
              └─────────────────────────────────────────────────────┘
                       ←  overlap region (Chimera fringe)  →
```

## Running the example

```bash
# Run with the built-in defaults; opens a PyVista viewer at the end.
python examples/05-chimera-assembly/assembly_sphere_in_cylinder.py

# Skip the viewer (CI-friendly).
python examples/05-chimera-assembly/assembly_sphere_in_cylinder.py --no-view

# Custom output directory.
python examples/05-chimera-assembly/assembly_sphere_in_cylinder.py -o my_assembly --no-view

# Drive both pieces from TOML configs.
python examples/05-chimera-assembly/assembly_sphere_in_cylinder.py \
    --near-body-config examples/05-chimera-assembly/config_assembly_near_body.toml \
    --background-config examples/05-chimera-assembly/config_assembly_background.toml \
    --no-view
```

By default the script:

1. Builds an in-memory icosphere of radius 0.5 (diameter 1.0) and
   writes it to `examples_out/sphere.stl`.
2. Builds the cubed-sphere near-body mesh through `Mesh.from_stl(...)`
   with `CubedSphereTopology()`, marches it outward, writes 6 `.geo`
   files (`wall_sub0..3.geo`, `wall_cap_north.geo`,
   `wall_cap_south.geo`).
3. Builds the flat-caps cylinder background through
   `FlatCapsCylinder.build(...)`, writes 5 `.geo` files
   (`background_center.geo`, `background_sub_e/n/w/s.geo`).
4. Bundles all 11 blocks into one `assembly.vtm` MultiBlock for
   ParaView. PyVista creates a sister `assembly/` directory holding
   the per-block `.vts` files; the `.vtm` is the manifest.

Output layout:

```
examples_out/
  sphere.stl                         # in-memory icosphere
  wall_sub0..sub3.geo            # 4 cubed-sphere equatorial blocks
  wall_cap_north/south.geo       # 2 cubed-sphere cap blocks
  background_center.geo              # flat-caps central square
  background_sub_e/n/w/s.geo         # flat-caps 4 trapezoidal sub-blocks
  assembly.vtm                       # 11-block ParaView MultiBlock manifest
  assembly/                          # 11 per-block .vts files (PyVista internals)
```

## Inspecting the result in ParaView

Open `examples_out/assembly.vtm`. ParaView pulls in all 11 blocks
automatically. Useful operations:

- **Surface with edges** representation — see the cell structure of
  every block.
- **Color by Jacobian** — every block carries a per-cell Jacobian as
  CellData, attached at write time. Healthy meshes are uniformly
  green on the default `RdYlGn` palette.
- **Block visibility toggle** — hide the background sub-blocks to see
  the near-body mesh in isolation, or vice versa.
- **Slice through z = 0** — the cleanest way to see the overlap
  region; the cubed-sphere outer shell at r ≈ 0.7 sits inside the
  background's central-square block (which extends to r ≈ 4 with the
  default `inner_frac = 0.4`).

## TOML configs

Both configs ship next to the script and are validated like any other
SHORE config:

```bash
shore config validate examples/05-chimera-assembly/config_assembly_near_body.toml
shore config validate examples/05-chimera-assembly/config_assembly_background.toml
```

Each config can also be used standalone with the corresponding CLI:

```bash
shore mesh -c examples/05-chimera-assembly/config_assembly_near_body.toml -o wall
# (will fail because input STL doesn't exist; the orchestrator overrides
#  the input key with the in-memory icosphere it generates)

shore primitive flat-caps -c examples/05-chimera-assembly/config_assembly_background.toml -o background
# 5 .geo files written; same cylinder as the assembly produces
```

## Parameter-tuning lessons baked into the defaults

Two parameter choices in this example are not obvious from first
principles. Both are documented as comments inside the TOML configs;
this section explains *why* they matter so users adapting the example
to their own bodies don't have to rediscover them.

### `ni` / `nj` and the C1 pinning window

The cubed-sphere topology pins each equator meridian's first cell to
the matching cap-edge first cell ([per-j C1 seam pin](/algorithm/k0-construction#per-j-c1-seam-pinning)).
The pin is active when the cap radial step
$\theta_{\text{cap}} / (n_c - 1)$ and the equator's natural radial
step $(\pi - 2\theta_{\text{cap}}) / (n_i - 1)$ are within a factor of 2
of each other, where $n_c = n_j / 4 + 1$.

For this example with `theta_cap = 30°`:

| `ni` | `nj` | `nc` | cap step | equator step | ratio |
|---|---|---|---|---|---|
| 30 | **60** | 16 | 0.0349 | 0.0722 | 0.48 (just outside) |
| 30 | **40** | 11 | 0.0524 | 0.0722 | 0.72 (inside) |

`30 × 60` would be a more natural choice on its own merits (more cells
on the equator and on the cap edges), but it falls outside the pinning
window — the equator falls back to unpinned spacing and emits a
`UserWarning`. `30 × 40` keeps the pinning active.

### Geometric vs. tanh wall-clustering on small bodies

Wall-clustered tanh spacing (`k_spacing = "tanh"` with `total_thickness`)
is the obvious default for boundary-layer meshes — it gives a fixed
outer shell and clusters cells near the wall. But on a *small* body
(here, radius 0.5), pushing `nk` and `beta` typical for production
wall-resolved LES (`nk = 25, beta = 4`) gives a first cell of order
$10^{-4}$, which is < 1 % of the surface cell scale. The cubed-sphere
cap-equator corner cells then invert on the very first march step.

Geometric spacing (`ds, growth`) gives a uniform cell aspect ratio
across layers and avoids this trap. The example uses
`ds = 0.01, growth = 1.1, nk = 13` — total thickness ≈ 0.21,
first cell ~30 % of the surface cell scale. The corner cells march
cleanly.

For larger bodies (wall, wing) the tanh-clustered law works as
expected — it's the geometric-scale ratio of `ds` to surface cell
size that matters, not the law per se.

## Generating Xall `cc.par` files

Once the geometry is verified in ParaView, convert each mesh to Xall
inputs with `examples/05-chimera-assembly/chimera_assembly.py` — a direct continuation of
this example that adds four steps: adjacency JSON sidecars, VTK
MultiBlock exports, binary `.grd` files, and `cc.par` generation.

```bash
# Generate everything headlessly (no viewer).
python examples/05-chimera-assembly/chimera_assembly.py

# Custom output directory.
python examples/05-chimera-assembly/chimera_assembly.py -o my_run

# Open the interactive viewer after generation.
python examples/05-chimera-assembly/chimera_assembly.py --view
```

The script uses `shore.io.cc_par.BlockMeta` to assign **per-block** BC
codes — demonstrating that `write_cc_par`'s `block_metadata` argument
(and the matching CLI `--meta` flag) give full per-block control:

```python
from shore.io.cc_par import BlockMeta, write_cc_par
from shore.io.grd import write_grd

# Pack each mesh into a binary .grd — block order must match cc.par.
write_grd("wall.grd", [
    "examples_out/wall_sub0.geo", "examples_out/wall_sub1.geo",
    "examples_out/wall_sub2.geo", "examples_out/wall_sub3.geo",
    "examples_out/wall_cap_north.geo", "examples_out/wall_cap_south.geo",
])
write_grd("background.grd", [
    "examples_out/background_center.geo",
    "examples_out/background_sub_e.geo", "examples_out/background_sub_n.geo",
    "examples_out/background_sub_w.geo", "examples_out/background_sub_s.geo",
])

# Near-body: equatorial sub-blocks and caps all use chimera BC 40.
nb_meta = {
    "sub0":      BlockMeta(free_face_bc=40),
    "sub1":      BlockMeta(free_face_bc=40),
    "sub2":      BlockMeta(free_face_bc=40),
    "sub3":      BlockMeta(free_face_bc=40),
    "cap_north": BlockMeta(free_face_bc=40),
    "cap_south": BlockMeta(free_face_bc=40),
}

# Background: all blocks get far-field extrapolation BC -9.
bg_meta = {blk.label: BlockMeta(free_face_bc=-9) for blk in bg_blocks}

write_cc_par(nb_adj, "wall.cc.par", grd_filename="wall.grd", block_metadata=nb_meta)
write_cc_par(bg_adj, "background.cc.par", grd_filename="background.grd", block_metadata=bg_meta)
```

The same result from the CLI:

```bash
# Binary grids
shore grd wall_sub0.geo wall_sub1.geo wall_sub2.geo wall_sub3.geo \
          wall_cap_north.geo wall_cap_south.geo \
          --output wall.grd

shore grd background_center.geo background_sub_e.geo background_sub_n.geo \
          background_sub_w.geo background_sub_s.geo \
          --output background.grd

# Chimera descriptors (block order must match the grd invocations above)
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd

shore cc-par background.adjacency.json background.cc.par \
  --grd background.grd --free-bc -9
```

Additional output files:

```
examples_out/
  wall.adjacency.json          # near-body adjacency sidecar
  background.adjacency.json    # background adjacency sidecar
  wall.vtm                     # near-body ParaView MultiBlock
  background.vtm               # background ParaView MultiBlock
  wall.grd                     # near-body Xall binary grid
  background.grd               # background Xall binary grid
  wall.cc.par                  # near-body Xall cc.par
  background.cc.par            # background Xall cc.par
```

See [`shore grd`](/reference/grd) for the binary format reference and
[`shore cc-par`](/reference/cc-par) for the full CLI reference including
the `--meta` flag format and BC code table.

## What's deliberately not in this example

- **Hole cutting / fringe identification / donor search.** The
  assembly is "side by side" — both meshes are written to disk in
  their world coordinates, the solver does the rest. SHORE has no
  in-tree assembly logic and there are no plans to add it.
- **Multi-body assemblies.** Single body. The pattern would extend
  naturally to N near-body meshes plus 1 background, but the
  orchestrator script doesn't loop.
- **Adaptive overlap sizing.** The background's `inner_frac = 0.4`
  was chosen so the central-square block extends well clear of the
  near-body mesh's outer shell (4.0 vs. 0.7). For tighter overlap
  budgets, reduce `inner_frac` and re-check the geometry in
  ParaView.

## See also

- [Chimera context](/guide/chimera) — what overset meshing is
  and why SHORE produces these particular outputs.
- [Cubed-sphere topology](/topologies/cubed-sphere) — the near-body
  mesh.
- [Primitives](/topologies/primitives) — the flat-caps cylinder used
  for the background.
- [Config files](/reference/config) — TOML schema for the mesh
  side; primitives use the parallel `kind = "primitive"` schema.
- [`shore grd`](/reference/grd) — pack per-block `.geo` files into
  the Xall binary `.grd`.
- [`shore export`](/reference/export) — the CLI command that bundles
  multiple `.geo` files into a `.vtm`. The orchestrator script does
  the same in-process so all 11 blocks share one manifest.
