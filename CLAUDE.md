# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Development Commands

```bash
make dev          # create .venv + install editable with [dev] extras
make test         # pytest with coverage (--cov=shore --cov-report=term-missing)
make lint         # ruff check + ruff format --check (read-only)
make fmt          # ruff check --fix + ruff format (destructive)
make clean        # remove dist/, build/, .egg-info, caches
```

Run a single test: `.venv/bin/pytest tests/test_cap_gnomonic.py -k "TestSeamStepMatching"`

## Architecture

SHORE converts triangulated STL surfaces into structured body-fitted volumes via hyperbolic normal extrusion, and packages the result for the Xall/Xnavis Chimera (overset) CFD solver chain. The package name is `shore-mesh`; the importable package and CLI are both `shore`.

Two halves:

1. **Mesh generation** — STL → structured multi-block volume (`mesh/`, `surface/`, `volume/`, `topology/`, `primitive/`). Several topologies (cubed-sphere, O-grid, C-grid, OH/CH hybrids) plus parametric primitives.
2. **Xall preprocessing pipeline** — turn the block list into the files the solver reads: ASCII `.geo` per block, binary `.grd` packs, `cc.par` chimera connectivity, `proc.input` MPI work distribution, marker boxes for hole-cutting, and a load-balance split planner (`io/`, `split.py`, `balance.py`, `boxes.py`). There are also `merge` variants (`grd-merge`, `cc-par-merge`, `proc-input-merge`) for combining independent Chimera components into one Xall run.

Geometry generators (`shapes/`, `profiles/`) round out the front end: analytic STL bodies (sphere/cylinder/torus/box/capsule), NACA airfoil sections, and section-extrusion wings — so a test case can be built end-to-end without external CAD.

### Two-Stage Pipeline (body-fitted topologies)

1. **Surface projection** (`surface/remesh.py`): Ray-cast from body centroid through a structured lat-lon sphere template onto the STL. Produces the k=0 surface layer `(ni, nj, 3)`. Requires star-shaped bodies.

2. **Hyperbolic marching** (`volume/hyperbolic.py`): Extrudes outward layer-by-layer: `P[i,j,k+1] = P[i,j,k] + step[k] * n_hat[i,j,k]`. Step sizes come from `spacing.build_k_steps` (geometric, tanh, or tanh2). Laplacian smoothing per layer. Aborts immediately on negative Jacobian.

### OOP Block Model

The mesh layer is an explicit object model, not bare arrays. `Block`
(ABC) / `HexBlock` (`mesh/block.py`) own a `(ni, nj, nk, 3)` `nodes`
array (often a numpy *view* into a topology's shared buffer so adjacent
blocks alias their seam data) plus typed `Vertex` / `Edge` / `Face`
lists. `Face` (`mesh/face.py`) carries a `FaceBC`
(`WALL | FREE | SHARED | DIRICHLET | PERIODIC | SYMMETRY`), a `partner` Face, and
an `AxisMap`. **Wiring invariant:** `HexBlock._set_face_bc` *requires*
an `AxisMap` for every SHARED/DIRICHLET/PERIODIC seam and validates it
against the partner face's along-axes at construction — a topology
builder that wires a seam wrong fails loudly then, not silently in the
downstream `cc.par`. `AxisMap` records, per along-seam axis, the
partner axis and a flip flag; `cc_par.py` turns it into Xnavis 3-digit
orientation codes.

Per-edge spacing lives on `Edge.spacing` (a `Spacing1D`). Override
detection is by **identity**: `HexBlock.__init__` assigns the same
`Spacing1D` instance to all edges of an axis, so `mesh.py`'s
`_resolve_k_steps` treats any k-edge whose `spacing is not
block._default_spacing_k` as a user override. All blocks in one `Mesh`
share a single k-step schedule (k-coupled seams must stay coincident
layer-by-layer); block 0 wins, with a `UserWarning` if other blocks
disagree.

### Topology

The `shore mesh` CLI supports five topologies, selected with
`--topology cubed_sphere|ogrid|cgrid|ch|oh` (default `cubed_sphere`).
`CGridTipCap` is available from the Python API only (and via the
dedicated `shore cap-tip` command / `--cap-tips` flag on the C-grid).
Hybrid topologies (`OHTopology`, `CHTopology`) compose a marched wrap
with non-marching far-field channel blocks; the channels are
SHARED-by-view to the wrap envelope where geometry permits and filled
by `topology.fill_channel_blocks` after marching:

**`cubed_sphere`** (6-block, `mesh/topology.py:CubedSphereTopology`, `surface/cap.py`,
`surface/equator.py`, `topology/sphere_cap.py`): 4 equatorial sub-blocks
(j-split) + 2 polar caps. Implementation modules retain "cap" as the
sub-block role name (`cap_north`, `cap_south`, `surface/cap.py`) — the
topology name is the cubed-sphere; "cap" inside it is a part.

**`ogrid`** (single-block, `mesh/topology.py:OGridTopology`,
`surface/remesh.py:sphere_latlon_directions`): single hexahedral
block, j-periodic, with polar singular axis rows at i=0 and i=ni-1.
Marched by `Mesh._march_ogrid` using `_compute_normals(periodic=(False, True))`
and the j-periodic `jacobian_field`. Writes one `<stem>.geo`.

**`cgrid`** (1- / 3- / 4-block C-grid, `mesh/topology.py:CGridTopology`,
`surface/cgrid.py`): C-grid wrap for sharp- or blunt-TE bodies
(airfoils, hydrofoils).  3-block (sharp): `upper`, `body`, `lower`
with SHARED TE-column seams and a wake-cut SHARED at `upper.k_lo ↔
lower.k_lo`.  4-block (blunt): adds `wake_base` bridging the TE-base
gap.  Optional conformal wing-tip cap blocks via `--cap-tips`
(per-tip 5 cap blocks + 1 wake-base cap; auto-rounds sharp TE).
Marched by `Mesh._march_cgrid_three` (TE-aware ghost rows).

**`ch`** (12-block sharp-TE / 14-block blunt-TE hybrid,
`mesh/topology.py:CHTopology`, `topology/ch_channel.py`): C-grid wrap
(3 or 4 blocks, sharp or blunt) + 5 k_hi-side channel blocks + 4 (or
5) outflow sub-blocks, bridging the wrap envelope to a
chord/normal/span-aligned far-field box.  The 5 k_hi-side channel
blocks `chan_upper` / `chan_body_top` / `chan_body_front` /
`chan_body_bot` / `chan_lower` are SHARED-by-view with their wrap k_hi
slices; the 4 outflow sub-blocks `chan_outflow_top_upper` /
`chan_outflow_wake_upper` / `chan_outflow_wake_lower` /
`chan_outflow_top_lower` (plus `chan_outflow_wake_base` for CH-4)
bridge the wrap's i-direction far-tip faces to the rear frame edge
with cross-axis SHARED face mapping (k=0 data copied — channel-i =
wrap-k).
Body wrap is split into 3 i-segments; the splits' outer-quad
endpoints land at the frame's top-front and bottom-front corners so
all channel outer quads are planar.  Frame is defined by
`chord_min/max`, `normal_min/max` and `freestream_world` (chord_axis
= largest non-span freestream component) + `span_axis`.  Channels
do not march; `Mesh.march` runs `_march_cgrid_three` on the wrap and
then `topology.fill_channel_blocks(blocks)`.  Optional C1 seam pin
matches each channel's first step to the wrap's last marching step
per (i, j).

**`oh`** (12-block hybrid, `mesh/topology.py:OHTopology`,
`topology/oh_channel.py`): cubed-sphere wrap (6 blocks, marched from
the STL) + 6 H-channel frustum blocks (`chan_sub0..chan_sub3`,
`chan_north`, `chan_south`) bridging the wrap k_hi envelope to an
axis-aligned far-field box.  Channels are SHARED-by-view to the wrap
k_hi slices (deeper buffer allocation) and filled by
`fill_channel_blocks` after `_march_caps_blocks`.  The wrap pins
`theta_cap = pi/4` and `phi_offset = pi/4` so each sub-block maps onto
exactly one cube face (any CLI `--theta-cap` is ignored, with a
`UserWarning`).  CLI frame is the box `--oh-box-{x,y,z}{min,max}`;
channel depth/law via the shared `--nk-channel` / `--k-law-channel` /
`--k-beta-channel` / `--c1-seam-channel` flags.  Star-shaped bodies
only (centroid ray-cast, like `cubed_sphere`).

- **Cap k=0 seed**: gnomonic flat-square projection
  (`surface/cap.py:build_cap_k0`). A flat diamond is placed in the
  tangent plane at each cap pole with corners on the cardinal axes at
  distance `tan(theta_cap)` (4 corners at azimuths
  `{0°, 90°, 180°, 270°}`, aligned with sub-block junctions). Interior
  by bilinear interpolation, then radially projected through the cached
  BVH onto the STL. Cap-quad AR ≈ 1.0 by construction; no concentric
  pole pattern.
- **Equator k=0**: `surface/equator.py:build_equator`. Seam rings are
  bit-exact copies of the cap edges (DIRICHLET-bound). Interior nodes
  along each j-meridian are great-circle slerp from south to north
  seam, ray-cast onto the STL. Each meridian's first/last interval is
  **per-j pinned** to the cap-edge first cell at the corresponding
  seam vertex (C1 spacing match at every seam node, including the
  3-edge cap-corner junctions). `_constrained_distribution` shapes the
  interior `ni − 2` cells using `i_spacing`. When the cap and equator
  radial step sizes differ by more than ~2×, the equator falls back to
  unpinned `i_spacing` and emits a `UserWarning`.
- **Marching**: `mesh/mesh.py::Mesh._march_caps` orchestrates all 6
  blocks per layer. **Seam-aware normals**: each block's
  `_compute_normals` call receives `ghost_rows` from the partner blocks
  (sub-sub j-seams + cap-equator i-seams, with index-mapping/reversal
  per `_overwrite_cap_boundaries`). The boundary tangent at every seam
  becomes a centred difference against the partner row, eliminating
  the spacing jump at multi-block edges that one-sided differences
  would accumulate over many marching layers. Cap blocks are
  Dirichlet-boundary-driven (4 edges pinned from sub-block pole rows
  at every layer, post-smoothing).
- **Spacing knobs**: `j_spacing` (`'uniform'` or `'tanh2'`) controls
  cap-edge tangent-plane distribution, inherited by the equator seam.
  `i_spacing` (`'uniform' | 'tanh' | 'tanh2'`) controls equator
  meridional distribution. `volume_k_spacing`
  (`'geometric' | 'tanh' | 'tanh2'`) controls wall-normal extrusion.
- **Deprecated CLI flags**: `--cap-pull`, `--cap-squareness`,
  `--j-cluster-lon` accept values but emit `DeprecationWarning` and
  are ignored. The gnomonic algorithm makes them unnecessary by
  construction.

### Key Conventions

- **Array shape**: grids are always `(ni, nj, nk, 3)` with last axis `(x, y, z)`
- **Periodicity**: `_compute_normals` default is `periodic=(False, True)` (legacy single-block O-grid). Cubed-sphere topology uses `(False, False)` for both equatorial sub-blocks and cap blocks, with seam-aware ghost rows replacing one-sided boundary differences.
- **Geo file ordering**: Fortran column-major (i-fastest) via `order="F"` — matching Xall_GPU
- **Intermediate format**: `.npz` files with key `"k0"` for surface layers; readable via `np.load("file.npz")["k0"]`
- **Jacobian check**: scalar triple product `(e_i x e_j) . e_k`; periodic j wraps via modular indexing; open blocks use `_jacobian_field_open` which has `nj-1` cells in j
- **Theta-cap aliases**: the CLI exposes `--theta-cap` (radians) and `--theta-cap-deg` (degrees), mutually exclusive on every command that accepts them. At the Python level `surface/remesh.py` also accepts a deprecated `theta_margin` radians kwarg (same exclusivity), but it is not surfaced as a CLI flag
- **Spacing laws**: low-level 1D distributions in `volume/spacing.py` (pure numpy — no scipy): `uniform` (equal spacing, default), `tanh` (one-sided Roberts, cluster near start), `tanh2` (two-sided Roberts, cluster at both ends or at `cluster_at`), and `pinned_distribution` (interior law with caller-given first/last interval lengths, for C1 endpoint matching at multi-block seams). `mesh/spacing.py:Spacing1D` is the OOP wrapper carried by every `Edge`: `law ∈ {uniform, tanh, tanh2, geometric, pinned}` plus parameters (`beta`, `ds`, `growth`, `total_thickness`, `cluster_at`, `first_step`/`last_step`, `interior`). Surface i/j directions and wall-normal k direction each have independent CLI controls (`--surface-i-spacing`, `--surface-j-spacing`, `--volume-k-spacing`). Tanh/tanh2 k-spacing requires `--volume-k-thickness`.
- **Xall file formats** (the solver target, all in `io/`): `.geo` ASCII per-block coordinates (i-fastest); `.grd` Fortran-unformatted binary pack of N blocks (header = block count + per-block cell-count triples, then X/Y/Z records); `cc.par` ASCII chimera connectivity + per-block metadata + marker boxes; `proc.input` MPI rank/group/body assignment. `.adjacency.json` is SHORE's own sidecar (the SHARED/DIRICHLET/PERIODIC face graph) and the bridge from the block model to `cc-par`. SHORE node counts (`ni`, `nj`, `nk`) are always **node** counts; the `.grd` header is the only place cell counts (= node − 1) appear, converted automatically. **Cubed-sphere `nj` caveat**: under `cubed_sphere` the CLI `--nj` is the full equator-*ring* node count, consumed as `nj // 4` — each of the 4 equatorial sub-blocks and both polar caps gets `nj // 4` cells (`nj // 4 + 1` nodes per sub-block, `nc = nj // 4 + 1`). So per-block j-cells = `nj // 4`, *not* `nj − 1`; e.g. `--nj 128` → 32 cells per sub-block. `_adjust_nj` rounds `nj` down to a multiple of 4. `ni`/`nk` follow the plain node = cell + 1 rule.
- **Smoothing**: `march()` applies Laplacian smoothing after each extrusion step: `smooth_iters` sweeps of `_laplacian_smooth` with blend factor `smooth_strength`. Non-periodic boundaries and the k=0 wall are always clamped. `smooth_mode='limited'` (default; `--smooth-mode`, `[march] smooth_mode`) then limits the per-node smoothing motion with `volume.hyperbolic._limit_smoothing`: pull-back ≤ 0.25·step along the marching normal, sideways ≤ 0.25·step; unbound nodes keep the plain Laplacian result bit-for-bit. Unlimited smoothing shrinks curved layers by ~strength·h²·κ and folds small first steps; `smooth_mode='legacy'` keeps that behaviour to reproduce old grids.
- **Exceptions**: `ShoreError` base with `MeshInputError(ValueError)`, `ParameterError(ValueError)`, `MarchingError(RuntimeError)`. Dual inheritance preserves backward-compat `isinstance` checks.

### Module Responsibilities

| Module | Role |
|---|---|
| `cli.py` | Typer CLI. Top-level commands: `mesh`, `project`, `collar` (full/half collars, ADR 0018), `xall` (project → merged Xall run with an options TOML, ADR 0018), `cgns` / `plot3d` (project → CGNS / Plot3D file, `-c` components, `--no-replay`, ADR 0035), `run`, `config`, `info`, `check`, `export`, `view`, `split`, `cap-tip`, `grd`, `grd-merge`, `cc-par`, `cc-par-merge`, `proc-input`, `proc-input-merge`, `balance`, `boxes-from-stl`, `boxes-export`. Sub-apps: `primitive` (`box`/`annulus`/`flat-caps`/`h-grid-rect`/`h-grid-tube`), `profile` (`naca`), `body` (`sphere`/`cylinder`/`torus`/`box`/`capsule`/`extrude`). `mesh` selects topology via `--topology cubed_sphere\|ogrid\|cgrid\|ch\|oh` (cgrid/ch need `--ni-body`/`--ni-wake`/`--n-stations`/`--wake-length`; ch adds `--ch-chord-min/max`/`--ch-normal-min/max`/`--nk-channel`, plus optional `--ch-freestream-{x,y,z}`/`--ch-i-split-body-upper/lower`/`--ch-ni-outflow`; oh needs `--oh-box-{x,y,z}{min,max}` and reuses `--nk-channel`/`--k-law-channel`/`--k-beta-channel`/`--c1-seam-channel`). `--project DIR --component C` on `mesh`/`collar`/`primitive *` loads (or creates) a project, builds component C and saves it (load-apply-save; geometries reused by SHA-256, a stem collision is an error). Parameter precedence everywhere: explicit CLI flag → `--config` TOML → built-in default (`_make_picker` + click's `ParameterSource`). Shared option aliases (`_NI`, `_NJ`, `_NK`, `_DS`, `_GROWTH`, `_SMOOTH`, `_SITERS`, spacing flags). Deprecated `mesh` flags: `--cap-pull`, `--cap-squareness`, `--j-cluster-lon` accept values, warn, are ignored. |
| `config.py` | TOML config loaders/validators + template emitters. `MeshConfig` (`[grid]`/`[march]`/`[surface]`/`[volume]`) for `shore mesh`; `PrimitiveConfig` (`[geometry]`/`[grid]`/`[spacing]`) for `shore primitive`; `SplitConfig` (`[[splits]]`) for `shore split`. `load_config` / `load_primitive_config` / `load_split_config`, `dump_template` / `dump_primitive_template`. Rejects unknown keys; resolves `--theta-cap`/`--theta-cap-deg`. |
| `shapes/` | trimesh-backed 3D STL body builders. `_primitives.py`: `sphere`, `cylinder`, `torus`, `box`, `capsule` (thin `trimesh.creation` wrappers). `_extrude.py`: `extrude_section(section_2d, span, n_stations, taper, sweep_deg, twist_deg, span_axis, end_caps)` — stack a 2D Selig section into a wing STL. Backs `shore body`. |
| `profiles/` | Pure-numpy 2D airfoil section generators. `_naca.py`: `naca` / `naca4` / `naca5` → `(N, 2)` Selig-ordered coords (upper-TE → LE → lower-TE), sharp (`closed_te=True`) or blunt TE. `io.py`: `read_section` / `write_section` for Selig `.dat`. Backs `shore profile naca`; consumed by `shapes.extrude_section`. |
| `geometry/` | `Geometry` protocol (ADR 0002): `base.py` (`Geometry`, `EntityKind`/`EntityRef`/`BODY`, `Projection`, `AnchoredCaster` — Projector-compatible anchored caster that can snap hits to exact geometry), `stl.py` (`StlGeometry(mesh, crease_deg=30)`, trimesh; points and ray casts bit-identical to `Projector`; crease-aware corner normals with planar-cut correction, ADR 0014; `Body.normals` keeps `_stl_normals`, ADR 0015), `occt.py` (`OcctGeometry`, STEP/IGES/BREP via build123d/OCP, optional `shore[cad]`: faces/edges/vertices, tier-1 tessellation projection with (u,v) + oriented normals, tier-2 exact projection with trimmed-face fallback, IGES sewn, `merged_proxy` watertight tessellation), `load_geometry(path)` suffix dispatch; `display()` on both backends returns a `base.Display` (float64 triangles with their owning entity, CAD edge polylines and vertices) for the GUI (ADR 0038). `Body.from_file` / `Mesh.from_file` / CLI `mesh`+`project` accept CAD. Tests: `test_geometry.py`, `test_geometry_occt.py`, `test_cad_pipeline.py` (CAD ones skip without build123d). |
| `mesh/association.py` | Block-to-geometry association (ADR 0003): `Association` (geometry key, `EntityKind`, `EntityKey` fingerprints; `Association.body()`), `AssociationResolver` (re-finds keys via `geometry.keys.match_keys`, projects nodes; `exact=` tier), `snap_block` (face → edge → vertex projection, vertex wins), `associate_wall`. `Vertex`/`Edge`/`Face.association`; `body_projected` is a read-only derived property. `Mesh(body=...)` associates every WALL k_lo (face + 4 edges + 4 corners) with the body; `Mesh.snap(exact=)`. `split_blocks` propagates associations (lowest-dimensional source entity); `.adjacency.json` v3 stores them (deduplicated table), v2 still read. Fingerprints: `geometry/keys.py` (`EntityKey`, bbox-normalised; `match_keys` = per-(kind,type) `linear_sum_assignment`, statuses same/changed/lost). |
| `pipeline.py` | Topology construction shared by `shore mesh` and the command layer: `k_spacing` (wall-normal `Spacing1D`), `make_topology(name, options)` (constructor args for `ch`/`oh`), `build_topology(body, name, nk=, spacing=, options=, **params)` (params checked against the named `build` parameters; their `**kwargs` would hide typos). |
| `project/` | Command layer (ADR 0004; multi-component, portable since ADR 0016: named geometries `load_geometry(path, name=)`, named components `use_component(name)` (default `"main"`; mesh commands act on the current one; `mesh`/`build` are its), `geometry=` on `build_topology`/`surface_grid`/`build_collar`, `write_xall(dir, grd=, np=, components={name: {free_face_bc, meta, boxes, body, proc_meta}}, output_basename=, multigrid_levels=, boundary_layer_thickness=, numerical_beach=)` = per-component files + merged `.grd`/`cc.par` + one `proc.input` (ADR 0017: header options run-wide, block options per component); `save` copies inputs into `inputs/` with SHA-256, drops output commands, format 2, `load` resolves inputs against the project dir and still reads format 1; `save` also writes each component's face graph `grids/<component>/adjacency.json`, and `load(replay=False)` (ADR 0033, `--no-replay`) records the journal instead of replaying it and reads each component's blocks from its grids + face graph: no march or undo of the loaded history). `commands.py`: `Command` (frozen dataclass, `apply` → `Undo`, `to_source` = journal line, `inputs` = input-file fields); state commands `LoadGeometry`, `UseComponent`, `BuildTopology`, `March` (progress + `CancelToken`, restores on failure/cancel; refuses non-marchable meshes), `Snap` (resolves every project geometry by name; `"body"` = the mesh's body, or the single geometry its component was built on after a load without replay, ADR 0039), `Associate` (`p.associate(block, entity, geometry, kind, indices)`: a face/edge/vertex `v0`..`v7` with body/surface/curve/point or `none`, own dimension or higher, fingerprinted, ADR 0039), `FillFace` (`p.surface_grid(block, face, surface=k)`: rebuild a block face on a CAD SURFACE, ADR 0007), `BuildCollar` (`p.build_collar(other, ...)`, ADR 0010), `BuildPrimitive`, `BuildTipCap`, `LoadBlocks`, `Split`, the blocking commands `CreateBlocking` / `SplitBlocking` / `Ogrid` / `MoveVertex` / `DeleteBlocks` / `AssociateBlocking` / `Premesh` (ADR 0042, 0043, 0045, see `mesh/blocking.py`); output-only (never undone) `WriteGeo` (naming: primitive rule / `per_block` / `shore mesh` rule), `WriteAdjacency`, `WriteGrd` (packs from `.geo` files, like `shore grd`), `WriteCcPar`, `WriteProcInput`, `WriteXall`, `WriteCgns` / `WritePlot3d` (ADR 0035); `UndoCache` (on-disk `.npz` snapshots). CLI `mesh`, `primitive *`, `cap-tip`, `split`, `grd`, `cc-par`, `proc-input` run over a `Project` (`--journal`); `*-merge`, `balance`, `boxes-*`, generators and inspection stay file tools. `project.py`: `Project` (state + history, one method per command, `undo`/`redo`, `journal()`/`save_journal`, `save(dir)` = `project.toml` + `journal.py` + `inputs/` + `adjacency.json` + `grids/<component>/*.npz`, `load(dir)` = input-hash check + replay + bit-identical grid check), `run_journal`. `shore mesh` runs over a `Project` (`--journal` writes it); `shore run JOURNAL [--save-project]` replays. |
| `quality.py` | Grid-quality metrics per block node array: `scaled_jacobian` (min corner, Verdict definition), `k_orthogonality_deg` (k-segment vs layer normal), `k_spacing`, `k_growth`, `aspect_ratio`, `summary(blocks)`. **Frozen**: the acceptance specs are written in these metrics. |
| `surface/io.py` | `load_stl` — loads an STL file into a `trimesh.Trimesh`; validates it is a single mesh; warns on non-watertight. |
| `surface/projector.py` | `Projector` — ray-cast helper around `trimesh.RayMeshIntersector` with cached BVH; `Projector.project(dirs)` returns hit positions + per-ray status. Used by `Body` (lazy) and by `build_cap_k0` / `build_equator` directly. |
| `surface/anchor.py` | `resolve_anchor(mesh)` — chooses an anchor inside the body for ray casting (centroid-based, with star-shaped checks). |
| `surface/cap.py` | `build_cap_k0(theta_cap, nc, pole, projector, j_spacing, j_beta)` — gnomonic flat-square cap k=0 seed. Diamond in tangent plane (corners at azimuth 0/90/180/270°, distance `tan(theta_cap)`), bilinear interior, radial projection via `Projector`. |
| `surface/equator.py` | `build_equator(cap_south, cap_north, ni, theta_cap, projector, i_spacing, i_beta)` — equator k=0 from cap seam rings. `_seam_ring(cap, pole)` extracts the `nj`-node seam in sub-block order. `_constrained_distribution(ni, start, stop, first_step, last_step, spacing, beta)` builds a 1D distribution with pinned first/last intervals. Per-j slerp-meridian + ray-cast onto STL. C1 seam pinning when cap/equator radial steps are within ~2x; warns + falls back otherwise. |
| `surface/remesh.py` | Single-block lat-lon path: `sphere_latlon_directions(ni, nj, theta_cap, i_spacing, i_beta, j_spacing, j_beta, j_cluster_lon)`, `project_sphere_to_surface(mesh, ...)`. Used by `Body.project()`, the `shore project` CLI command, and `OGridTopology` for the k=0 surface ray-cast. **Not** used by the cubed-sphere topology. |
| `surface/intersect.py` | Intersection loop of two bodies, the start curve of a collar (ADR 0009): `intersection_loop(body_a, body_b, seed, n, *, entity_a=BODY, entity_b=BODY) -> (n, 3)`. Alternating projections take the seed onto the curve; predictor steps along `n_b x n_a` with alternating-projection correctors trace it until it closes; arc-length resample to `n` nodes, polished onto both bodies. Ordered so `march_surface(A, loop)` / `march_surface(B, loop[::-1])` march the collar halves away from the other body. `Geometry.project` only (STL, CAD, analytic). `TraceParams`; grazing intersections raise. `intersection_curve(body_a, body_b, seed, n, *, plane=(point, normal))`: the open curve on the seed's side of a symmetry plane (ADR 0011), traced both ways until it reaches the plane (within tolerance: half-model bodies stop there), ends put on the plane by Newton along the tangent; works on half models. Both take `spacing=Spacing1D` (loop: uniform/tanh2; curve: also tanh) and `curvature=alpha` (ADR 0020): polyline refined where it turns > 2 deg, windowed curvature, equidistribution of `1 + alpha*kappa*L/2pi`, adjacent cells <= 1.2x (`TraceParams`); defaults bit-identical to the plain arc-length resampling. Corners (ADR 0021): a failed step (no progress, or tangent turn >= `corner_deg` 30) is redone by a recursive 1/16-step walk; a failure at `corner_scale` (1e-6 diag) is a corner, left via curve-sphere points (`_ring`, alternating projections) and located exactly by crossing branch tangent lines ever closer; corners are fixed nodes, segments get cells by length (`_allocate`) and the law per segment; smooth traces (tight bends too) bit-identical. With corners and `curvature > 0` (ADR 0022), `_graded_spacing` (g = C/m capped by cones of slope ln 1.2, C by bisection) places nodes over the whole curve, so cells stay within 1.2x across corners; smooth curves keep ADR 0020's limiter (pinned). |
| `surface/march.py` | Hyperbolic surface marching, the collar-grid kernel (ADR 0008): `march_surface(geometry, curve, steps, *, periodic=False, entity=BODY) -> (n, m+1, 3)`. Per layer one sparse solve of the linearised surface system (orthogonality, cell area, tangent plane) about pure `n x t` marching, area smoothing + implicit/explicit dissipation + contraction sensor (`MarchParams`), then projection onto the surface (which also gives the next layer's normals: one geometry query per node and layer). Works on STL and CAD bodies. Converging fronts overshoot the schedule by ~3.5 % (sphere). `plane=(point, normal)` holds both ends of an open curve on a symmetry plane (ADR 0012): mirror ghost nodes, ends put on surface and plane by alternating projections; `plane=None` is bit-identical to the free-end march. At creases (a front corner on a surface edge, `|n . t| > MarchParams.crease_tolerance` 0.1) the frame's normal is made orthogonal to the front tangent (ADR 0023); smooth marches bit-identical. At a sharp crease (front folded back, `|t|` < `MarchParams.sharp_crease` 0.5 x the arc spacing) the tangent is lengthened to the arc spacing (ADR 0024). At clustered corners (ADR 0027): a layer whose step exceeds `max_step_ratio` 24 x the smallest spacing of the nodes not closer than on the starting curve is marched in sub-steps of <= `substep_ratio` 2 x that spacing (`_substeps`); a layer that folds against the previous one (`_inverted`: a cell reversed or a front segment pointing backwards) is smoothed along the front round the fold and re-projected until it does not (`_untangle`); marches where neither fires are bit-identical. |
| `collar.py` | The collar block of two intersecting bodies (ADR 0010): `build_collar(body_a, body_b, seed, *, n_loop, steps_a, steps_b, steps_k, label="collar") -> HexBlock`. Wall = `intersection_loop` + `march_surface` on both bodies (junction row `i = len(steps_b)`, j periodic around the loop); hyperbolic layers while thinner than half the junction spacing; per loop station, an outer boundary offset in the station plane with its swallowtail trimmed and the crease mapped to the junction row; TFI + Winslow/Thomas-Middlecoff interior (Picard + SuperLU, all sheets at once). `CollarParams`. `extrude_collar(wall, steps_k, junction)` is the public volume stage. `plane=(point, normal)`: the half-collar of a half model (ADR 0013): `intersection_curve` + `march_surface(plane=)`, the volume extruded from the wall mirrored into a closed loop and halved, j-faces on the plane as `FaceBC.SYMMETRY`; warns when the curve meets the plane obliquely. `Project.build_collar(other, ..., plane=None)` (`BuildCollar`): body B loaded from a file path, the collar replaces the mesh (no body association). `loop_spacing` / `loop_curvature` forward ADR 0020's node distribution along the junction curve (`Project.build_collar(loop_law=, loop_beta=, loop_cluster_at=, loop_curvature=)`, `shore collar --loop-*`). `split_collar(...) -> list[HexBlock]` (ADR 0024): the loop cut at sharp corners (tracer corners turning > `CollarParams.sharp_deg` 105, via the private `intersect._loop`), each piece extended by a wake line on A (bisector, cells from the corner spacing growing <= `wake_ratio` 1.2, `wake_length` default `sum(steps_a)`) and a virtual sheet (the corner's B-line translated; B marched with its ends on the two sheets, `march._march_between`/`_Ends`), extruded open in j, split into six blocks `{label}_{p}_{wake_lo|body|wake_hi}_{b|a}`; the sheet blocks' k_lo SHARED across the cut (AxisMap i same, j reversed); no sharp corner: `build_collar`'s block; `Project.build_collar(split=True, wake_length=)` / `shore collar --split --wake-length`. `CollarParams.volume` (ADR 0025, 0026): `"smoothed"` (default since ADR 0026; `_smoothed_outer`: each station's outer line on its 3-D offset polyline (`_offset_line`, `_trim_origins`), station offsets smoothed along the loop by a Gaussian of `offset_smoothing` 0.5 x thickness with one kernel for all i, only near corners at the collar's scale (`_near_corners`: the wall loop turns >= 60 deg within +- the thickness), local 1-2-1 smoothing where the boundary steps back (`backstep_ratio` 0.5), each line trimmed again) or `"stations"` (ADR 0010's per-station outer boundary, bit-identical to before); `BuildCollar.volume` / `--collar-volume` default `None` = the default, except while a project saved before format 3 replays (`legacy_collar_volume`: `"stations"`, so old projects still load). Projects are format 3. Valid thick collars at convex corners; concave corners >= 4.7 x and the blunt base still fold a few cells. |
| `surface/tfi.py` | Surface grids on one parametric face (ADR 0007): `surface_grid(geometry, surface, boundary, smoothing="elliptic" or "none") -> SurfaceGrid(points, params)`. Boundary rows projected once to (u, v), TFI in (u, v), Winslow in (u, v) with Christoffel terms + Thomas-Middlecoff control functions (clustering kept), metric from a once-evaluated parameter lattice (`EllipticParams`), Picard + SuperLU; one evaluation at the end. Needs a parametric SURFACE (`OcctGeometry`); seam-crossing / singular faces out of scope and rejected (boundary parameter-step ratio > 50x median, or a collapsed boundary edge). |
| `surface/cgrid.py` | C-grid k=0 construction from an STL: `detect_te_edges` (sharp/blunt TE by dihedral), `slice_by_span` (per-station cross-sections via `mesh_plane`), `build_c_curve` / `build_cgrid_k0` (wrap + wake-branch polylines, arc-length resample), `round_body_te_in_surface` (synthetic TE rounding). Feeds `CGridTopology` / `CHTopology`. |
| `volume/spacing.py` | Low-level 1D spacing laws: `uniform`, `tanh_one_sided`, `tanh_two_sided`; domain builders `build_i_coords` (`uniform`/`tanh`/`tanh2`), `build_j_coords` (`uniform`/`tanh2`, periodic), `build_k_steps` (`geometric`/`tanh`/`tanh2`), `pinned_distribution` (interior law + pinned first/last intervals). Pure numpy — no scipy. |
| `mesh/spacing.py` | `Spacing1D` dataclass — the OOP spacing law carried by every `Edge`. `build(n, start, stop)` dispatches to `volume/spacing.py`; `build_steps(nk)` returns wall-normal step sizes. See the **Spacing laws** convention above. |
| `mesh/block.py` | `Block` ABC + `HexBlock` (8 vertices / 12 edges / 6 faces over a `(ni, nj, nk, 3)` node array, often a view). `_set_face_bc(name, bc, partner, axis_map)` enforces the AxisMap-required-for-connection-seams invariant. Default face BCs: k_lo=WALL, all others FREE; topology rewires after construction. |
| `mesh/face.py` | `FaceBC` enum (`WALL`/`FREE`/`SHARED`/`DIRICHLET`/`PERIODIC`/`SYMMETRY`), `Face` dataclass (name/bc/partner/axis_map/body_projected), and `AxisMap` (per-along-axis partner-axis + flip; `identity`/`from_pairs`/`inverse`/`validate_against_partner`). `FACE_ALONG_AXES` maps each face name to its two in-plane axes. |
| `mesh/edge.py` / `mesh/vertex.py` | `Edge` (axis, endpoint vertices, node count, `Spacing1D`, `body_projected`) and `Vertex` (corner index, ijk sign tuple, `body_projected`) — the rest of the block element model. |
| `mesh/body.py` | `Body` dataclass — STL wrapper; an unresolvable anchor is kept in `anchor_error` and raised by `projector` (open half-model bodies serve collars, ADR 0018). `Body.normals` (the reference of ADR 0005's k-line angle metric; the marcher's `blend_normals_k` that read them is removed, ADR 0034: 0 warns, anything else raises `ParameterError`): `_stl_normals` vertex-normal blend for STL, exact normals for CAD (ADR 0015). `Body.from_stl(path, skip_hygiene=...)` loads + (optionally) repairs; `.projector` (lazy cached `Projector`), `.project(...)` (lat-lon k=0), `.project_caps(...)`, `.directions(...)`, `save_surface`/`load_surface` (`.npz` with key `"k0"`). |
| `mesh/inspect.py` | File-level inspection/export backing the CLI: `info_geo`/`info_surface`, `check_geo`, `export_geo`/`export_surface`, `view_geo_file`/`view_surface_file`. |
| `topology/cgrid_cap.py` | Wing-tip cap builder for the C-grid: `build_cap_blocks` (3/4 H-style blocks closing the airfoil section at a span extreme), `find_symmetric_le_cuts`, `detect_te_thickness`, `auto_round_sharp_te`. Backs `CGridTipCap` and the `shore cap-tip` CLI. |
| `primitive/` | Parametric structured-grid primitives (no STL). `base.py` (`PrimitiveMesh` ABC, `resolve_axis`/`axis_frame`/`coords_1d`), `box.py` (`BoxMesh`), `annulus.py` (`AnnularCylinderMesh`, optional periodic-j sector), `flat_caps.py` (`FlatCapsCylinder`, 5-block Chimera-free background), `h_grid_rect.py` (`HGridRectMesh`, rectangular duct), `h_grid_tube.py` (`HGridTubeMesh`, circular bore, j-periodic, k=0 polar singular axis). All expose `.build(...)` + `.blocks` + `write_geo`. |
| `exceptions.py` | `ShoreError` base + `MeshInputError(ShoreError, ValueError)`, `ParameterError(ShoreError, ValueError)`, `MarchingError(ShoreError, RuntimeError)`. |
| `volume/hyperbolic.py` | `_compute_normals(layer, periodic, ghost_rows=None)`: outward unit normals via cross of i- and j-tangents. When `ghost_rows` is supplied (dict with optional keys `'i_lo'` / `'i_hi'` / `'j_lo'` / `'j_hi'`), the corresponding boundary tangent uses a centred difference against the ghost row instead of one-sided. `_laplacian_smooth(layer, strength, periodic)`: per-iteration `strength * (laplacian - node)` blend; non-periodic boundary rows clamped, k=0 wall always clamped. `_stl_normals` (vertex-normal blend at the closest point; `Body.normals` on STL). |
| `volume/jacobian.py` | `jacobian_field` (periodic-j, shape `(ni-1, nj, nk-1)`), `_jacobian_field_open` (non-periodic, shape `(ni-1, nj-1, nk-1)`), `min_jacobian_two_layers` (used during marching), `jacobian_per_layer` (used by `check`). |
| `topology/sphere_cap.py` | Low-level march helpers for the 6-block cubed-sphere topology: `_adjust_nj` (rounds nj to a multiple of 4), `_overwrite_cap_boundaries` (Dirichlet pin from equator pole rows), `_march_one_layer` / `_march_cap_layer` (hyperbolic step + Laplacian, with optional `ghost_rows` forwarded to `_compute_normals`), `_check_jacobian`. Cap **seed** face lives in `surface/cap.py`; equator seed in `surface/equator.py`. |
| `topology/oh_channel.py` | OH (and shared by CH) channel-block frustum builder: `build_channel_block(inner_face, outer_quad, nk_channel, k_law, k_beta, c1_first_step)` — bilinear sweep of `outer_quad` to produce the channel's k=nk-1 face, ruled-surface interior, optional per-(i, j) C1 first-step pin via geometric tail (`_solve_geometric_ratio` Newton solver). Pure numpy; no marching. |
| `topology/ch_channel.py` | CH-specific frame and outer-quad helpers: `CHFrame` dataclass (chord/normal/span axis indices + box extents), `resolve_frame(...)` (validates inputs, picks chord_axis from freestream), `outer_quads_for_cgrid(...)` (returns 7 or 8 planar `(4, 3)` quads keyed by channel label — `upper`, `body_top`, `body_front`, `body_bot`, `lower`, `outflow_upper`, `outflow_lower`, +`outflow_wake_base` for CH-4). Body-side splits land at frame corners so all outer quads are planar and adjacent quads agree at shared columns. |
| `mesh/topology.py` | `Topology` ABC and `CubedSphereTopology.build(surface, nk, ni, nj, spacing_k, projector, theta_cap, i_spacing, i_beta, j_spacing, j_beta)` — calls `build_cap_k0` for both poles, `build_equator`, then wires the 6 `HexBlock`s with SHARED views (sub-sub j-seams) and DIRICHLET partners (cap-equator i-seams). The `surface` arg is ignored (kept for backwards-compat); the gnomonic algorithm projects through `projector` directly.  Also defines `OGridTopology` (single-block, j-periodic), `CGridTopology` (sharp/blunt 1- / 3- / 4-block C-grid; conformal cap-tips path), `CGridTipCap` (C-grid wing-tip cap, fully volumetric — no marching), `OHTopology` (cubed-sphere wrap + 6-face H-channel; 12 blocks, channel views aliased to wrap.k_hi via deeper buffer allocation, `fill_channel_blocks(...)` after marching), and `CHTopology` (C-grid wrap + H-channel; 12 blocks for sharp TE, 14 for blunt — 5 k_hi-side channels SHARED-by-view via single `nodes_all` buffer of depth `nk_wrap+nk_chan-1`; 4 (or 5) outflow sub-blocks per side `top_upper / wake_upper / wake_lower / top_lower` (+ `wake_base` for blunt) with conformal SHARED to wrap.upper.i_lo / lower.i_hi (nk = nk_wrap) and to chan_upper.i_lo / chan_lower.i_hi (nk = nk_channel)).  All non-cubed-sphere topologies expose `is_blunt` and `fill_channel_blocks` where applicable. |
| `mesh/mesh.py` | `Mesh.from_stl(...)` orchestrates body load + topology build for `(ni, nj, nk)`-signature topologies; `Mesh.from_array(surface, topology, nk, ...)` builds from a pre-made surface array (e.g. a cross-section polyline) with `body=None` (used by `CGridTipCap`). `Mesh.march()` dispatches by topology: `_march_caps()` for `CubedSphereTopology` (per layer: snapshots all 6 blocks, builds ghost rows via `_sub_ghost_rows` / `_cap_ghost_rows`, advances each block with seam-aware normals + Laplacian + cap-block DIRICHLET pin, restitches sub-sub seams, runs Jacobian guards, commits), `_march_ogrid()` for `OGridTopology` (single-block, `periodic=(False, True)`, j-periodic Jacobian guard), `_march_cgrid()` (1/3/4-block C-grid, dispatching to `_march_cgrid_three`), the OH/CH hybrids (run the underlying wrap marcher via `_march_caps_blocks` / `_march_cgrid_three`, then `topology.fill_channel_blocks(blocks)`), and `CGridTipCap` (no-op — fully volumetric at build). `_resolve_k_steps` picks the one shared k-schedule (see OOP Block Model). In-memory inspection: `.info()`, `.check()`, `.export()`/`.write_vtk()`, `.view()`. `Mesh.from_stl` rejects `CGridTopology` and `CHTopology` because their build params (`ni_body`/`ni_wake`/`n_stations`/`wake_length`) don't fit the `(ni, nj, nk)` signature — build the block list explicitly. |
| `io/geo.py` | ASCII geo format (Xall_GPU compatible): `read_geo`, `write_geo`. |
| `io/cgns.py` | CGNS/HDF5 export (ADR 0035, `h5py`, extra `[cgns]`): `write_cgns({base: blocks}, path)`: a base per component, a structured zone per block, `GridConnectivity1to1_t` in both zones for SHARED/DIRICHLET/PERIODIC (`transform_vector` from the `AxisMap`; `_links` pairs faces once and `_locate` places the smaller face inside the larger by its nodes, so CH's multi-partner `body.k_hi` works; unmatched nodes refused), `BC_t` for WALL (`BCWall`), SYMMETRY (`BCSymmetryPlane`), FREE (`BCTypeUserDefined` + `shore:FREE`); library version 3.4. |
| `mesh/seams.py` | Point-matched seams (ADR 0035, 0036): `seams(blocks, tol, only=)` pairs connected faces once and locates the smaller face inside the larger by its nodes (SIDS ranges + `Transform`); `shared_faces` (faces with several partners); `transform_vector`. Used by `io/cgns.py` and by `write_adjacency_json` (the `patches` of shared faces). `Aliases(blocks, tol)` / `sync(changed)` (ADR 0040): groups of nodes that are one point across SHARED/DIRICHLET seams; after an in-place edit each group takes its last changed member's value (as views would); `Snap` and `FillFace` edit inside `_seams_written_through`, so one-sided associations and loaded projects keep seams closed; wall-only snaps of built meshes stay bit-identical. |
| `mesh/redistribute.py` | Spacing editor core (ADR 0041): `redistribute(blocks, block, edge, Spacing1D)`: the law over the edge's arc length → fractions `t`; every line of the block along the axis resampled at `t` on its old polyline (`resample`, batched), ends kept; `plan` propagates through SHARED/DIRICHLET faces containing the axis (reversed across flipped axis maps), refusing multi-partner and periodic faces and conflicting distributions. Command `Redistribute` / `Project.redistribute(block, edge, law, *, beta, cluster_at, first_step, last_step, interior)`, inside `_seams_written_through`. |
| `mesh/blocking.py` | Blocking (ADR 0042): `Blocking` (skeleton `points`, `blocks` of 8 vertex ids in `HexBlock.vertices()` order, `nodes` per edge, `associations` keyed `("v", id)` / `("e", a, b)` / `("f", a, b, c, d)`; `labels`/`block`/`shape`/`entity`/`edges`/`faces`/`polylines`), `box`, `split` (at a node, carried through every block sharing a cut edge; high parts `b<n>`), `move` (ADR 0045: one vertex, projected onto its own or inherited association), `ogrid` (ADR 0043: a layer on a block set's faces, `faces` or its boundary less `exclude`; vertex copies moved by `(1 - factor)` x half the sum of their across edges; open faces passed through; refuses faces touching the layer shared with blocks outside the set), `delete`, `associate`, `effective_associations` (inherited lowest dimension first when they agree), `premesh` (edges and faces once, projected, copied into every block: seams bit for bit; `tfi_face` / `tfi_volume`; SHARED with axis maps from corner ids, WALL if associated, else FREE). Commands `CreateBlocking`, `SplitBlocking`, `Ogrid`, `MoveVertex`, `DeleteBlocks`, `AssociateBlocking`, `Premesh`; `Project.blocking` / `blockings` (per component; rebuilt by replay only). GUI: `Session.skeleton()`, `PygfxViewport.set_skeleton`, `SkeletonPick`, Blocking dock (`split_selection`); ADR 0044: `Blocking.centres`, `Session.skeleton_markers`, `PygfxViewport.set_markers` / `MarkerPick`, one `blocking_selection` of `SkeletonKey`s (edge mode / face and block markers; shift adds or toggles), `ogrid_selection`, `delete_selection`, `associate_blocking_selection`; vertex mode, `move_selection`, `move_selection_to_pick` (ADR 0045). |
| `io/plot3d.py` | Plot3D export (ADR 0035): `write_plot3d(blocks, path)`, multi-block whole 3-D grid, Fortran unformatted, little endian, double precision. |
| `io/grd.py` | Fortran-unformatted binary `.grd` pack of N `.geo` blocks: `write_grd`, `read_grd_metadata` (per-block cell-count triples — also the weight source for `balance`/`proc-input`), `merge_grd` (byte-passthrough concat for multi-component Xall runs). Backs `shore grd` / `grd-merge`. |
| `io/vtk.py` | VTK XML export (optional `pyvista` dep): `write_vts`, `write_surface_vts`, `write_caps_vtm` (MultiBlock for multi-block topologies), `write_boxes_vtm` (marker boxes), `geo_to_pyvista` / `surface_to_pyvista` / `open_to_pyvista`. |
| `io/adjacency.py` | JSON sidecar for block lists: `write_adjacency_json` / `read_adjacency_json` serialise the per-face BC + partner graph (FaceBC name, `{label, face}` partner reference, AxisMap). `summarise_adjacency` for a human summary. The bridge from the block model to `cc-par`. |
| `io/cc_par.py` | Xnavis `cc.par` reader-aware writer/merger. `write_cc_par(adjacency_payload, path, *, grd_filename, ..., block_metadata, boxes)` translates `.adjacency.json` → Fortran ASCII: maps `WALL/FREE/SHARED/DIRICHLET/SYMMETRY` to Xnavis BC codes (SYMMETRY → natural-BC index 2), emits paired connection patches with 3-digit orientation codes derived from the `AxisMap` (`_orientation_code_from_axis_map`; the normal digit is odd for faces on opposite sides of their blocks and even for the same side, as Xall's overset pre-processor reads it, ADR 0036 — before it lo faces were reversed) and one patch per partner on a face shared with several (`patches` in the sidecar, CH's `body.k_hi`), and writes a marker-boxes section. `BlockMeta` = per-block level/group/priority/comment + `free_face_bc` (40 chimera vs negative natural BC). `merge_cc_par` shifts block/patch indices across components (parses inputs back in). |
| `io/proc_input.py` | `proc.input` (MPI work distribution) writer/merger: `write_proc_input(grd, out, np, ...)` greedily weight-balances blocks across ranks (weights from the `.grd` cell counts), with per-block `proc`/`group`/`body` overrides; `merge_proc_input` concatenates components (no re-balancing). `BalanceSummary` reports per-rank load + imbalance. |
| `io/boxes.py` | `Box` dataclass (8-vertex hexahedron, canonical i-fastest order) + `box_to_json` / `box_from_json`. The cc.par marker-box value type. |
| `split.py` | `BlockSplit` dataclass + `split_blocks(blocks, splits)`: subdivides a list of `HexBlock`s along ijk vertex indices, wires internal SHARED seams between consecutive chunks (numpy-view aliased; bit-exact), rewires external SHARED partner pointers. Raises `ParameterError` on cross-block split mismatches (e.g. splitting `sub0` along k without paired splits on `sub1`/`sub2`/`sub3`). Used by `shore split` and by `balance`-generated plans. |
| `balance.py` | Iterative split *planner* for MPI load balancing. `plan_balance(grd, np, tolerance, max_iterations, adjacency_path)` reads `.grd` cell counts (+ optional adjacency seam graph), greedily schedules topology-safe cuts on the heaviest block (k-axis preferred since k-faces are WALL/FREE, never SHARED; falls back to j/i with paired co-splits), returns a `BalancePlan` that serialises to a `shore split` `split.toml`. Plans only — never mutates files. Backs `shore balance`. |
| `boxes.py` | `boxes_from_stl(stl, block, max_outward_gap, method)` — generate chimera hole-cutting marker boxes from a body STL. `voxel-fill` (default, watertight interior cover; resolution from the gap budget), `aabb` / `obb` (single best-effort box with an outward-gap post-check + warning). Backs `shore boxes-from-stl`. |
| `gui/` | Desktop GUI (ADR 0006, 0037; extra `[gui]`): `session.py` (Qt-free `Session` over one `Project`: `run(method, ...)` = one command = one journal line, `undo`/`redo`/`open`/`save`/`new`, `blocks()` → `(BlockKey(component, label), block)`, `tree()`, listeners), `forms.py` (`Form`/`Field`/`FORMS`: a form per command, `arguments(values)`), `viewport/__init__.py` (`Viewport` protocol, `PickResult`, structured arithmetic `quad_cell`/`quad_centre`/`boundary_values`, `colormap`), `viewport/pygfx_backend.py` (the only pygfx/wgpu/rendercanvas import: `PygfxViewport` with float64→float32 rebase over blocks, slices and geometries; quads with per-face colours; pick → `PickResult` (block, face, cell, xyz; slices `SliceKey`) or `GeometryPick` (geometry, `EntityRef`, xyz), ADR 0038), `window.py` (`MainWindow`: model tree (with an association column) + Geometry tab, viewport, command forms, Display dock (visibility, slices: view state, not journaled), Association dock (selection mode face/edge/vertex via `nearest_edge`/`nearest_vertex`, geometry picks, shift adds, Associate / Clear → `p.associate`, ADR 0039), Spacing dock (`edge_info`, `redistribute_selection` → `p.redistribute`, ADR 0041), journal/log console; long forms on a one-thread `QThreadPool` with progress/cancel; `main`). `shore gui` / `shore-gui`. Tests in `tests_gui/` (outside `testpaths`; `test-gui` CI job, offscreen Qt + lavapipe). |
| `vis.py` | PyVista interactive visualisation: `view_geo`, `view_surface`, `save_screenshot`. Optional dep (`shore[vis]`); Jupyter backend via `jupyter=True`. |

### Frozen acceptance specs (phase 4 discipline)

`tests/spec/` holds acceptance specs written and validated **before** an
algorithm is implemented (ADR 0005; the `ai-research-verification`
discipline; ADR 0007 surface grids, ADR 0008 surface marching, ADR 0009 intersection loops, ADR 0010 collar block, ADR 0011 open intersection curves, ADR 0012 marcher ends on a plane, ADR 0013 half-collar block, ADR 0014 STL normals, ADR 0016 multi-component projects, ADR 0017 Xall options, ADR 0018 CLI projects, ADR 0019 blend deprecation (spec replaced by ADR 0034's), ADR 0020 junction-curve spacing, ADR 0021 junction curves with corners, ADR 0022 corner grading, ADR 0023 marching at creases, ADR 0024 split collar, ADR 0025 collar volume at corners, ADR 0026 smoothed volume by default, ADR 0027 marching at clustered corners, ADR 0033 load without replay, ADR 0034 blend removal, ADR 0035 CGNS and Plot3D exports, ADR 0036 cc.par as Xall reads it, ADR 0037 GUI slice 1, ADR 0038 GUI slice 2, ADR 0039 association editing, ADR 0040 seams kept closed, ADR 0041 redistribute, ADR 0042 blocking, ADR 0043 O-grid, ADR 0044 blocking selection, ADR 0045 vertex moves, ADR 0048 orthogonal wall lines). While implementing, do **not** edit `tests/spec/**`,
`tests/conftest.py`, the pytest config or `src/shore/quality.py`. The only
permitted edit is deleting a case's `@PENDING` (strict xfail) line once that
case passes. `python scripts/check_spec.py` verifies the hashes in
`tests/spec/SHA256SUMS` (spec modules hashed without `@PENDING` lines; CI
lint job). Changing a threshold or metric is a re-freeze: a separate,
user-approved commit (ADR 0005 records the first one). Report the verdict as
counts (e.g. "13 passed, 0 skipped, 0 xfailed"), not "it passes".
Verdicts are per module and per CI job: `test_surface_tfi.py`,
`test_surface_march.py`, `test_intersection.py`, `test_collar.py`, `test_open_curve.py`, `test_march_ends.py`, `test_half_collar.py`, `test_stl_normals.py`, `test_project_components.py`, `test_xall_options.py`, `test_cli_projects.py`, `test_blend_removal.py` (with its data `data/blend_removal_default.npz`, `data/blended_project.zip`), `test_loop_spacing.py` (with its golden `data/loop_default.npz`) `test_loop_corners.py` (golden `data/loop_bends_default.npz`) `test_corner_grading.py` and `test_march_creases.py` (golden `data/march_default.npz`) `test_split_collar.py` (golden `data/collar_default.npz`) `test_collar_volume.py`, `test_collar_default.py`, `test_march_corners.py`, `test_fast_load.py`, `test_exports.py` (needs `pycgns`, in `dev`), `test_cc_par_xall.py` (data `data/xall_before.zip`) , `test_seam_sync.py`, `test_redistribute.py`, `test_blocking.py`, `test_ogrid.py`, `test_move_vertex.py` and `test_orthogonalize.py` run in the main job, `tests_gui/spec/test_gui_shell.py`, `test_gui_inspect.py`, `test_gui_associate.py`, `test_gui_spacing.py`, `test_gui_blocking.py`, `test_gui_ogrid.py`, `test_gui_blocking_select.py`, `test_gui_move_vertex.py` and `test_gui_orthogonalize.py` (need `shore[gui]`, `test_gui_associate.py` and `test_gui_spacing.py` also `shore[cad]`; outside `testpaths`) in the GUI job, `tests/spec/test_associate.py` (needs `shore[cad]`) in the CAD job, `test_surface_tfi_cad.py` and `test_loop_corners_cad.py` (need
`shore[cad]`) in the CAD job.

### Test Suite

Tests in `tests/` use `conftest.py` session fixtures/helpers:
- `sphere_stl` — exports a unit icosphere STL once per session
- `naca0012_stl` / `naca0012_blunt_stl` — extruded NACA 0012 wing STLs (sharp / blunt TE) for the C-grid and hybrid tests
- `sphere_layer(ni, nj, ...)` — plain function returning a `(ni, nj, 3)` lat-lon sphere layer (legacy single-block tests)

Test files (grouped by subsystem):
- **CLI / config**: `test_cli.py`, `test_primitive_cli.py`, `test_config.py`
- **Mesh + topologies**: `test_mesh.py`, `test_cap_gnomonic.py` (gnomonic cap corners, aspect ratio, seam C0/C1, volume seam continuity, marching robustness), `test_ogrid.py`, `test_cgrid.py`, `test_cgrid_tip_cap.py`, `test_cgrid_cap.py`, `test_oh_topology.py`, `test_oh_channel.py`, `test_ch_topology.py`, `test_ch_channel.py`, `test_axis_map.py`, `test_edge_spacing.py`
- **Primitives**: `test_primitive.py`, `test_h_grid_rect.py`, `test_h_grid_tube.py`
- **Geometry front end**: `test_profiles.py`, `test_shapes.py`
- **Surface / numerics**: `test_anchor.py`, `test_load_stl_hygiene.py`, `test_projector.py`, `test_projector_validation.py`, `test_spacing.py`, `test_jacobian.py`
- **Xall pipeline**: `test_geo.py`, `test_grd.py`, `test_cc_par.py`, `test_proc_input.py`, `test_balance.py`, `test_boxes.py`, `test_boxes_from_stl.py`, `test_split.py`, `test_vtk.py`, `test_vis.py`
- **Worked examples** (run the `examples/` scripts end-to-end): `test_assembly_example.py`, `test_split_example.py`
- `test_exceptions.py` — exception hierarchy

Full suite: 999 tests, ~1 min on a modern laptop (the marching and example tests dominate; coverage ~85%).
