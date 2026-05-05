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

Run a single test: `.venv/bin/pytest tests/test_hyperbolic.py -k "test_march_sphere"`

## Architecture

SHORE converts triangulated STL surfaces into structured body-fitted O-grid volumes via hyperbolic normal extrusion. It is the foundational meshing layer for Chimera (overset) CFD mesh generation.

### Two-Stage Pipeline

1. **Surface projection** (`surface/remesh.py`): Ray-cast from body centroid through a structured lat-lon sphere template onto the STL. Produces the k=0 surface layer `(ni, nj, 3)`. Requires star-shaped bodies.

2. **Hyperbolic marching** (`volume/hyperbolic.py`): Extrudes outward layer-by-layer: `P[i,j,k+1] = P[i,j,k] + step[k] * n_hat[i,j,k]`. Step sizes come from `spacing.build_k_steps` (geometric, tanh, or tanh2). Laplacian smoothing per layer. Aborts immediately on negative Jacobian.

### Topology

The `shore mesh` CLI supports four topologies, selected with
`--topology cubed_sphere|ogrid|cgrid|ch` (default `cubed_sphere`).
Two further topologies (`OHTopology`, `CGridTipCap`) are available
from the Python API only.  Hybrid topologies (`OHTopology`,
`CHTopology`) compose a marched wrap with non-marching far-field
channel blocks; the channels are SHARED-by-view to the wrap envelope
where geometry permits and filled by `topology.fill_channel_blocks`
after marching:

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

**`ch`** (10- / 12-block hybrid, `mesh/topology.py:CHTopology`,
`topology/ch_channel.py`): C-grid wrap (3 or 4 blocks, sharp or
blunt) + 7 or 8 H-channel blocks bridging the wrap envelope to a
chord/normal/span-aligned far-field box.  Channel blocks: `chan_upper`
/ `chan_body_top` / `chan_body_front` / `chan_body_bot` / `chan_lower`
SHARED-by-view with their wrap k_hi slices; `chan_outflow_upper` /
`chan_outflow_lower` (and `chan_outflow_wake_base` for CH-4) bridge
the wrap's i-direction far-tip faces to the rear frame edge with
cross-axis SHARED face mapping (k=0 data copied — channel-i = wrap-k).
Body wrap is split into 3 i-segments; the splits' outer-quad
endpoints land at the frame's top-front and bottom-front corners so
all channel outer quads are planar.  Frame is defined by
`chord_min/max`, `normal_min/max` and `freestream_world` (chord_axis
= largest non-span freestream component) + `span_axis`.  Channels
do not march; `Mesh.march` runs `_march_cgrid_three` on the wrap and
then `topology.fill_channel_blocks(blocks)`.  Optional C1 seam pin
matches each channel's first step to the wrap's last marching step
per (i, j).

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
- **Theta-cap aliases**: `--theta-cap` (radians), `--theta-cap-deg` (degrees), `--theta-margin` (deprecated radians alias) are mutually exclusive on every command that accepts them
- **Spacing laws**: Three 1D distributions in `volume/spacing.py`: `uniform` (equal spacing, default), `tanh` (one-sided Roberts, cluster near start), `tanh2` (two-sided Roberts, cluster at both ends or at `cluster_at`). Surface i/j directions and wall-normal k direction each have independent spacing controls via CLI flags (`--surface-i-spacing`, `--surface-j-spacing`, `--volume-k-spacing`). Tanh/tanh2 k-spacing requires `--volume-k-thickness`.
- **Smoothing**: `march()` applies Laplacian smoothing after each extrusion step: `smooth_iters` sweeps of `_laplacian_smooth` with blend factor `smooth_strength`. Non-periodic boundaries and the k=0 wall are always clamped.
- **Exceptions**: `ShoreError` base with `MeshInputError(ValueError)`, `ParameterError(ValueError)`, `MarchingError(RuntimeError)`. Dual inheritance preserves backward-compat `isinstance` checks.

### Module Responsibilities

| Module | Role |
|---|---|
| `cli.py` | Typer CLI — 10 commands: `mesh`, `primitive` (sub-app: `box`/`annulus`/`flat-caps`), `split`, `cc-par`, `project`, `config`, `info`, `check`, `export`, `view`. `mesh` selects topology via `--topology cubed_sphere|ogrid` (default `cubed_sphere`); `project` runs the single-block lat-lon path. Shared option aliases (`_NI`, `_NJ`, `_NK`, `_DS`, `_GROWTH`, `_SMOOTH`, `_SITERS`, spacing flags) enforce consistent help text. Deprecated flags: `--cap-pull`, `--cap-squareness`, `--j-cluster-lon` accept values, warn, and are ignored. |
| `surface/io.py` | `load_stl` — loads an STL file into a `trimesh.Trimesh`; validates it is a single mesh. |
| `surface/projector.py` | `Projector` — ray-cast helper around `trimesh.RayMeshIntersector` with cached BVH; `Projector.project(dirs)` returns hit positions + per-ray status. Used by `Body` (lazy) and by `build_cap_k0` / `build_equator` directly. |
| `surface/anchor.py` | `resolve_anchor(mesh)` — chooses an anchor inside the body for ray casting (centroid-based, with star-shaped checks). |
| `surface/cap.py` | `build_cap_k0(theta_cap, nc, pole, projector, j_spacing, j_beta)` — gnomonic flat-square cap k=0 seed. Diamond in tangent plane (corners at azimuth 0/90/180/270°, distance `tan(theta_cap)`), bilinear interior, radial projection via `Projector`. |
| `surface/equator.py` | `build_equator(cap_south, cap_north, ni, theta_cap, projector, i_spacing, i_beta)` — equator k=0 from cap seam rings. `_seam_ring(cap, pole)` extracts the `nj`-node seam in sub-block order. `_constrained_distribution(ni, start, stop, first_step, last_step, spacing, beta)` builds a 1D distribution with pinned first/last intervals. Per-j slerp-meridian + ray-cast onto STL. C1 seam pinning when cap/equator radial steps are within ~2x; warns + falls back otherwise. |
| `surface/remesh.py` | Single-block lat-lon path: `sphere_latlon_directions(ni, nj, theta_cap, i_spacing, i_beta, j_spacing, j_beta, j_cluster_lon)`, `project_sphere_to_surface(mesh, ...)`. Used by `Body.project()`, the `shore project` CLI command, and `OGridTopology` for the k=0 surface ray-cast. **Not** used by the cubed-sphere topology. |
| `volume/spacing.py` | 1D spacing laws: low-level `uniform`, `tanh_one_sided`, `tanh_two_sided`; domain builders `build_i_coords` (`uniform`/`tanh`/`tanh2`), `build_j_coords` (`uniform`/`tanh2`, periodic), `build_k_steps` (`geometric`/`tanh`/`tanh2`). Pure numpy — no scipy. |
| `exceptions.py` | `ShoreError` base + `MeshInputError(ShoreError, ValueError)`, `ParameterError(ShoreError, ValueError)`, `MarchingError(ShoreError, RuntimeError)`. |
| `volume/hyperbolic.py` | `_compute_normals(layer, periodic, ghost_rows=None)`: outward unit normals via cross of i- and j-tangents. When `ghost_rows` is supplied (dict with optional keys `'i_lo'` / `'i_hi'` / `'j_lo'` / `'j_hi'`), the corresponding boundary tangent uses a centred difference against the ghost row instead of one-sided. `_laplacian_smooth(layer, strength, periodic)`: per-iteration `strength * (laplacian - node)` blend; non-periodic boundary rows clamped, k=0 wall always clamped. `_stl_normals` for blending. |
| `volume/jacobian.py` | `jacobian_field` (periodic-j, shape `(ni-1, nj, nk-1)`), `_jacobian_field_open` (non-periodic, shape `(ni-1, nj-1, nk-1)`), `min_jacobian_two_layers` (used during marching), `jacobian_per_layer` (used by `check`). |
| `topology/sphere_cap.py` | Low-level march helpers for the 6-block cubed-sphere topology: `_adjust_nj` (rounds nj to a multiple of 4), `_overwrite_cap_boundaries` (Dirichlet pin from equator pole rows), `_march_one_layer` / `_march_cap_layer` (hyperbolic step + Laplacian, with optional `ghost_rows` forwarded to `_compute_normals`), `_check_jacobian`. Cap **seed** face lives in `surface/cap.py`; equator seed in `surface/equator.py`. |
| `topology/oh_channel.py` | OH (and shared by CH) channel-block frustum builder: `build_channel_block(inner_face, outer_quad, nk_channel, k_law, k_beta, c1_first_step)` — bilinear sweep of `outer_quad` to produce the channel's k=nk-1 face, ruled-surface interior, optional per-(i, j) C1 first-step pin via geometric tail (`_solve_geometric_ratio` Newton solver). Pure numpy; no marching. |
| `topology/ch_channel.py` | CH-specific frame and outer-quad helpers: `CHFrame` dataclass (chord/normal/span axis indices + box extents), `resolve_frame(...)` (validates inputs, picks chord_axis from freestream), `outer_quads_for_cgrid(...)` (returns 7 or 8 planar `(4, 3)` quads keyed by channel label — `upper`, `body_top`, `body_front`, `body_bot`, `lower`, `outflow_upper`, `outflow_lower`, +`outflow_wake_base` for CH-4). Body-side splits land at frame corners so all outer quads are planar and adjacent quads agree at shared columns. |
| `mesh/topology.py` | `Topology` ABC and `CubedSphereTopology.build(surface, nk, ni, nj, spacing_k, projector, theta_cap, i_spacing, i_beta, j_spacing, j_beta)` — calls `build_cap_k0` for both poles, `build_equator`, then wires the 6 `HexBlock`s with SHARED views (sub-sub j-seams) and DIRICHLET partners (cap-equator i-seams). The `surface` arg is ignored (kept for backwards-compat); the gnomonic algorithm projects through `projector` directly.  Also defines `OGridTopology` (single-block, j-periodic), `CGridTopology` (sharp/blunt 1- / 3- / 4-block C-grid; conformal cap-tips path), `CGridTipCap` (C-grid wing-tip cap, fully volumetric — no marching), `OHTopology` (cubed-sphere wrap + 6-face H-channel; 12 blocks, channel views aliased to wrap.k_hi via deeper buffer allocation, `fill_channel_blocks(...)` after marching), and `CHTopology` (C-grid wrap + H-channel; 12 blocks for sharp TE, 14 for blunt — 5 k_hi-side channels SHARED-by-view via single `nodes_all` buffer of depth `nk_wrap+nk_chan-1`; 4 (or 5) outflow sub-blocks per side `top_upper / wake_upper / wake_lower / top_lower` (+ `wake_base` for blunt) with conformal SHARED to wrap.upper.i_lo / lower.i_hi (nk = nk_wrap) and to chan_upper.i_lo / chan_lower.i_hi (nk = nk_channel)).  All non-cubed-sphere topologies expose `is_blunt` and `fill_channel_blocks` where applicable. |
| `mesh/mesh.py` | `Mesh.from_stl(...)` orchestrates body load + topology build. `Mesh.march()` dispatches by topology: `_march_caps()` for `CubedSphereTopology` (per layer: snapshots all 6 blocks, builds ghost rows via `_sub_ghost_rows` / `_cap_ghost_rows`, advances each block with seam-aware normals + Laplacian + cap-block DIRICHLET pin, restitches sub-sub seams, runs Jacobian guards, commits), `_march_ogrid()` for `OGridTopology` (single-block hyperbolic march with `periodic=(False, True)` and j-periodic Jacobian guard), `_march_cgrid()` (1/3/4-block C-grid, dispatches to `_march_cgrid_three` for the multi-block path), and the OH/CH hybrids (run the underlying wrap marcher, then call `topology.fill_channel_blocks(blocks)` to fill the channel interiors). `Mesh.from_stl` rejects `CGridTopology` and `CHTopology` because their build params (`ni_body`/`ni_wake`/`n_stations`/`wake_length`) don't fit the `(ni, nj, nk)` signature — build the block list explicitly. |
| `io/geo.py` | ASCII geo format (Xall_GPU compatible): `read_geo`, `write_geo`. |
| `io/vtk.py` | VTK XML StructuredGrid export: `write_vts`, `write_surface_vts`, `write_caps_vtm` (MultiBlock for the 6-block cubed-sphere topology), `geo_to_pyvista`, `surface_to_pyvista`. |
| `io/adjacency.py` | JSON sidecar for split block lists: `write_adjacency_json` / `read_adjacency_json` serialise the per-face BC + partner graph (FaceBC name + `{label, face}` partner reference). Read by downstream solvers to learn which faces are SHARED between which split chunks. |
| `io/cc_par.py` | Xnavis `cc.par` writer (the chimera pre-processor input format used by Xall_GPU): `write_cc_par(adjacency_payload, path, *, grd_filename, ..., block_metadata)` translates a SHORE `.adjacency.json` into the Fortran ASCII format. Maps `WALL/FREE/SHARED/DIRICHLET` to Xnavis BC codes; emits paired connection patches with the 3-digit orientation codes (135 for same-axis seams, lookup table from `_overwrite_cap_boundaries` for cubed-sphere cap-equator cross-axis pairs). `BlockMeta` dataclass exposes per-block level/group/priority/comment + `free_face_bc` (40 chimera vs negative natural BC for far-field). One-way bridge — no reader. |
| `split.py` | `BlockSplit` dataclass + `split_blocks(blocks, splits)`: post-processing pass that subdivides any list of `HexBlock`s along ijk vertex indices, wires internal SHARED seams between consecutive chunks (numpy-view aliased; bit-exact node coincidence), and rewires external SHARED partner pointers. v1 supports single-axis-per-block per call; raises `ParameterError` on cross-block split mismatches (e.g. splitting `sub0` along k without paired splits on `sub1`/`sub2`/`sub3`). Used by the `shore split` CLI and reusable from Python. |
| `vis.py` | PyVista interactive visualisation: `view_geo`, `view_surface`, `save_screenshot`. Optional dep (`shore[vis]`); Jupyter backend via `jupyter=True`. |

### Test Suite

Tests in `tests/` use `conftest.py` helpers:
- `sphere_stl` — `pytest.fixture(scope="session")`, exports a unit icosphere STL once per session
- `sphere_layer(ni, nj, ...)` — plain function returning a `(ni, nj, 3)` lat-lon sphere layer (used by legacy single-block tests)

Test files:
- `test_cli.py` — all CLI commands via `typer.testing.CliRunner`
- `test_mesh.py` — OOP `Mesh` API + `CubedSphereTopology.build` block wiring
- `test_cap_gnomonic.py` — gnomonic cap k=0 corners on a unit sphere; cap aspect ratio; seam C0 conformity (north + south); equator meridian colatitude monotonicity; per-j seam C1 step matching (`TestSeamStepMatching`); volume seam continuity across marching layers (`TestVolumeSeamContinuity`); marching robustness for `nk` ∈ {5, 10, 20}
- `test_anchor.py` / `test_load_stl_hygiene.py` — anchor resolution, STL preconditioning
- `test_projector.py` / `test_projector_validation.py` — `Projector` ray-cast contract + status enum
- `test_spacing.py` — 1D spacing laws, domain builders, error cases
- `test_jacobian.py` — Jacobian field calculations, inversion detection
- `test_geo.py` / `test_vtk.py` / `test_vis.py` — file I/O round-trips, VTK export, PyVista helpers
- `test_split.py` — `BlockSplit` validation, single- and multi-cut splits, internal-seam wiring, external-face BC preservation, cross-block adjacency rewiring, JSON sidecar round-trip, `shore split` CLI
- `test_exceptions.py` — exception hierarchy

Full suite: 194 tests, runs in ~2.5 s on a modern laptop.
