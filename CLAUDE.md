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

### Topology Variants

- **`single`** (`cli.py:_run_single_topology`): One O-grid block, j-periodic.
- **`double`** (Yin-Yang, `cli.py:_run_double_topology`, `surface/yinyang.py`): Two rotated O-grids covering the sphere without pole singularities. Each block is independent after projection. Default `theta_cap = pi/6`.
- **`caps`** (6-block, `cli.py:_run_caps_topology`, `topology/sphere_cap.py`): 4 equatorial sub-blocks (j-split) + 2 polar caps. Caps are Dirichlet-boundary-driven: all 4 edges pinned from equatorial sub-blocks at every layer, interior solved by Laplacian convergence. `split_ogrid_with_polar_caps` orchestrates simultaneous marching of all 6 blocks. Supports `--cap-pull` and `--cap-squareness` quality controls.

### Key Conventions

- **Array shape**: grids are always `(ni, nj, nk, 3)` with last axis `(x, y, z)`
- **Periodicity**: `march()` default is `periodic=(False, True)` (i open, j periodic) for single/Yin-Yang; `(False, False)` for caps and equatorial sub-blocks
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
| `cli.py` | Typer CLI — 7 commands: `mesh`, `project`, `march`, `info`, `check`, `export`, `view`. Routes to `_run_single_topology`, `_run_double_topology`, `_run_caps_topology`. `_resolve_spacing()` validates spacing flags and builds `surface_kwargs` dict + `k_steps` array. Shared option aliases (`_NI`, `_NJ`, `_NK`, `_DS`, `_GROWTH`, `_SMOOTH`, `_SITERS`, spacing flags) enforce consistent help text across commands. |
| `surface/io.py` | `load_stl` — loads an STL file into a `trimesh.Trimesh`; validates it is a single mesh. |
| `surface/remesh.py` | `sphere_latlon_directions` builds unit direction vectors (accepts `i_spacing`, `i_beta`, `j_spacing`, `j_beta`, `j_cluster_lon`); `project_sphere_to_surface` does the ray-casting projection (same new kwargs forwarded) |
| `surface/yinyang.py` | `YANG_ROTATION` matrix (90-degree rotation); `project_yang_to_surface` applies the rotation then same projection logic (same 5 spacing kwargs forwarded) |
| `volume/spacing.py` | 1D spacing laws: low-level `uniform`, `tanh_one_sided`, `tanh_two_sided`; domain builders `build_i_coords` (latitude, accepts `uniform`/`tanh`/`tanh2`), `build_j_coords` (longitude, periodic, `uniform`/`tanh2`), `build_k_steps` (wall-normal, `geometric`/`tanh`/`tanh2`). Pure numpy — no scipy. |
| `exceptions.py` | `ShoreError` base + `MeshInputError(ShoreError, ValueError)`, `ParameterError(ShoreError, ValueError)`, `MarchingError(ShoreError, RuntimeError)`. |
| `volume/hyperbolic.py` | `march()` is the core extruder. Accepts `steps: np.ndarray \| None` — when supplied overrides `ds`/`growth`. `_compute_normals` and `_laplacian_smooth` accept `periodic=(periodic_i, periodic_j)` tuple; default `(False, True)`. `_laplacian_smooth` blends `strength * (laplacian - node)` per iteration; boundary rows clamped in non-periodic directions and at wall (k=0). |
| `volume/jacobian.py` | `jacobian_field` (periodic-j, shape `(ni-1, nj, nk-1)`), `_jacobian_field_open` (non-periodic, shape `(ni-1, nj-1, nk-1)`), `min_jacobian_two_layers` (two-layer Jacobian used during marching), `jacobian_per_layer` (min per layer, used by `check` command) |
| `topology/sphere_cap.py` | Splits O-grid into 6 blocks, builds cap seed faces via Laplacian convergence with Dirichlet BCs, marches all blocks simultaneously with shared-face enforcement. Accepts `steps: np.ndarray \| None`. |
| `topology/sphere_cap.py` | Splits O-grid into 6 blocks, builds cap seed faces via Laplacian convergence with Dirichlet BCs, marches all blocks simultaneously with shared-face enforcement. Accepts `steps: np.ndarray \| None`. |
| `io/geo.py` | ASCII geo format (Xall_GPU compatible): `read_geo`, `write_geo` |
| `io/vtk.py` | VTK XML StructuredGrid export: `write_vts`, `write_surface_vts`, `geo_to_pyvista`, `surface_to_pyvista` |
| `vis.py` | PyVista interactive visualisation: `view_geo`, `view_surface`, `save_screenshot`. Optional dep (`shore[vis]`); Jupyter backend via `jupyter=True`. |

### Test Suite

Tests in `tests/` use `conftest.py` helpers:
- `sphere_stl` — `pytest.fixture(scope="session")`, exports a unit icosphere STL once per session
- `sphere_layer(ni, nj, ...)` — plain function returning a `(ni, nj, 3)` lat-lon sphere layer
- `yang_layer(ni, nj, ...)` — plain function returning a `(ni, nj, 3)` Yang-rotated layer

Test files: `test_cli.py` (CLI via `typer.testing.CliRunner`), `test_hyperbolic.py`, `test_jacobian.py`, `test_sphere_cap.py`, `test_yinyang.py`, `test_geo.py`, `test_vtk.py`, `test_vis.py`, `test_spacing.py` (1D spacing laws, domain builders, error cases). Coverage: geometric correctness (outer radius > inner, within 5%), Jacobian inversion detection, VTK round-trips, all topology variants, spacing monotonicity and endpoint constraints.
