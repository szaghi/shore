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

The `shore mesh` CLI builds the **caps topology** only (6 blocks). Other
topologies (single-block lat-lon O-grid, Yin-Yang) still exist as Python
APIs (`Body.project`, `surface/remesh.py`, `surface/yinyang.py`) but are
not wired to the CLI.

**`caps`** (6-block, `mesh/topology.py:CapsTopology`, `surface/cap.py`,
`surface/equator.py`, `topology/sphere_cap.py`): 4 equatorial sub-blocks
(j-split) + 2 polar caps.

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
- **Periodicity**: `_compute_normals` default is `periodic=(False, True)` (legacy single-block O-grid). Caps topology uses `(False, False)` for both equatorial sub-blocks and cap blocks, with seam-aware ghost rows replacing one-sided boundary differences.
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
| `cli.py` | Typer CLI — 6 commands: `mesh`, `project`, `info`, `check`, `export`, `view`. `mesh` always builds the caps topology (no `--topology` flag); `project` runs the legacy single-block lat-lon path. Shared option aliases (`_NI`, `_NJ`, `_NK`, `_DS`, `_GROWTH`, `_SMOOTH`, `_SITERS`, spacing flags) enforce consistent help text. Deprecated flags: `--cap-pull`, `--cap-squareness`, `--j-cluster-lon` accept values, warn, and are ignored. |
| `surface/io.py` | `load_stl` — loads an STL file into a `trimesh.Trimesh`; validates it is a single mesh. |
| `surface/projector.py` | `Projector` — ray-cast helper around `trimesh.RayMeshIntersector` with cached BVH; `Projector.project(dirs)` returns hit positions + per-ray status. Used by `Body` (lazy) and by `build_cap_k0` / `build_equator` directly. |
| `surface/anchor.py` | `resolve_anchor(mesh)` — chooses an anchor inside the body for ray casting (centroid-based, with star-shaped checks). |
| `surface/cap.py` | `build_cap_k0(theta_cap, nc, pole, projector, j_spacing, j_beta)` — gnomonic flat-square cap k=0 seed. Diamond in tangent plane (corners at azimuth 0/90/180/270°, distance `tan(theta_cap)`), bilinear interior, radial projection via `Projector`. |
| `surface/equator.py` | `build_equator(cap_south, cap_north, ni, theta_cap, projector, i_spacing, i_beta)` — equator k=0 from cap seam rings. `_seam_ring(cap, pole)` extracts the `nj`-node seam in sub-block order. `_constrained_distribution(ni, start, stop, first_step, last_step, spacing, beta)` builds a 1D distribution with pinned first/last intervals. Per-j slerp-meridian + ray-cast onto STL. C1 seam pinning when cap/equator radial steps are within ~2x; warns + falls back otherwise. |
| `surface/remesh.py` | Legacy single-block lat-lon path: `sphere_latlon_directions(ni, nj, theta_cap, i_spacing, i_beta, j_spacing, j_beta, j_cluster_lon)`, `project_sphere_to_surface(mesh, ...)`. Used by `Body.project()` and the `shore project` CLI command. **Not** used by the caps topology. |
| `surface/yinyang.py` | `YANG_ROTATION` matrix; `project_yang_to_surface`. Yin-Yang topology, Python-API only (not wired to CLI). |
| `volume/spacing.py` | 1D spacing laws: low-level `uniform`, `tanh_one_sided`, `tanh_two_sided`; domain builders `build_i_coords` (`uniform`/`tanh`/`tanh2`), `build_j_coords` (`uniform`/`tanh2`, periodic), `build_k_steps` (`geometric`/`tanh`/`tanh2`). Pure numpy — no scipy. |
| `exceptions.py` | `ShoreError` base + `MeshInputError(ShoreError, ValueError)`, `ParameterError(ShoreError, ValueError)`, `MarchingError(ShoreError, RuntimeError)`. |
| `volume/hyperbolic.py` | `_compute_normals(layer, periodic, ghost_rows=None)`: outward unit normals via cross of i- and j-tangents. When `ghost_rows` is supplied (dict with optional keys `'i_lo'` / `'i_hi'` / `'j_lo'` / `'j_hi'`), the corresponding boundary tangent uses a centred difference against the ghost row instead of one-sided. `_laplacian_smooth(layer, strength, periodic)`: per-iteration `strength * (laplacian - node)` blend; non-periodic boundary rows clamped, k=0 wall always clamped. `_stl_normals` for blending. |
| `volume/jacobian.py` | `jacobian_field` (periodic-j, shape `(ni-1, nj, nk-1)`), `_jacobian_field_open` (non-periodic, shape `(ni-1, nj-1, nk-1)`), `min_jacobian_two_layers` (used during marching), `jacobian_per_layer` (used by `check`). |
| `topology/sphere_cap.py` | Low-level march helpers for the 6-block caps topology: `_adjust_nj` (rounds nj to a multiple of 4), `_overwrite_cap_boundaries` (Dirichlet pin from equator pole rows), `_march_one_layer` / `_march_cap_layer` (hyperbolic step + Laplacian, with optional `ghost_rows` forwarded to `_compute_normals`), `_check_jacobian`. Cap **seed** face lives in `surface/cap.py`; equator seed in `surface/equator.py`. |
| `mesh/topology.py` | `Topology` ABC and `CapsTopology.build(surface, nk, ni, nj, spacing_k, projector, theta_cap, i_spacing, i_beta, j_spacing, j_beta)` — calls `build_cap_k0` for both poles, `build_equator`, then wires the 6 `HexBlock`s with SHARED views (sub-sub j-seams) and DIRICHLET partners (cap-equator i-seams). The `surface` arg is ignored (kept for backwards-compat); the gnomonic algorithm projects through `projector` directly. |
| `mesh/mesh.py` | `Mesh.from_stl(...)` orchestrates body load + topology build. `Mesh.march()` calls `_march_caps()`, which per layer: snapshots all 6 blocks, builds ghost rows via `_sub_ghost_rows` / `_cap_ghost_rows`, advances each block with seam-aware normals + Laplacian + (caps) DIRICHLET pin, restitches sub-sub seams, runs Jacobian guards, commits. |
| `io/geo.py` | ASCII geo format (Xall_GPU compatible): `read_geo`, `write_geo`. |
| `io/vtk.py` | VTK XML StructuredGrid export: `write_vts`, `write_surface_vts`, `write_caps_vtm` (MultiBlock for the 6-block caps topology), `geo_to_pyvista`, `surface_to_pyvista`. |
| `vis.py` | PyVista interactive visualisation: `view_geo`, `view_surface`, `save_screenshot`. Optional dep (`shore[vis]`); Jupyter backend via `jupyter=True`. |

### Test Suite

Tests in `tests/` use `conftest.py` helpers:
- `sphere_stl` — `pytest.fixture(scope="session")`, exports a unit icosphere STL once per session
- `sphere_layer(ni, nj, ...)` — plain function returning a `(ni, nj, 3)` lat-lon sphere layer (used by legacy single-block tests)

Test files:
- `test_cli.py` — all CLI commands via `typer.testing.CliRunner`
- `test_mesh.py` — OOP `Mesh` API + `CapsTopology.build` block wiring
- `test_cap_gnomonic.py` — gnomonic cap k=0 corners on a unit sphere; cap aspect ratio; seam C0 conformity (north + south); equator meridian colatitude monotonicity; per-j seam C1 step matching (`TestSeamStepMatching`); volume seam continuity across marching layers (`TestVolumeSeamContinuity`); marching robustness for `nk` ∈ {5, 10, 20}
- `test_anchor.py` / `test_load_stl_hygiene.py` — anchor resolution, STL preconditioning
- `test_projector.py` / `test_projector_validation.py` — `Projector` ray-cast contract + status enum
- `test_spacing.py` — 1D spacing laws, domain builders, error cases
- `test_jacobian.py` — Jacobian field calculations, inversion detection
- `test_geo.py` / `test_vtk.py` / `test_vis.py` — file I/O round-trips, VTK export, PyVista helpers
- `test_exceptions.py` — exception hierarchy

Full suite: 194 tests, runs in ~2.5 s on a modern laptop.
