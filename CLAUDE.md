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

2. **Hyperbolic marching** (`volume/hyperbolic.py`): Extrudes outward layer-by-layer: `P[i,j,k+1] = P[i,j,k] + ds * r^k * n_hat[i,j,k]`. Laplacian smoothing per layer. Aborts immediately on negative Jacobian.

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

### Module Responsibilities

| Module | Role |
|---|---|
| `cli.py` | Typer CLI — 7 commands: `mesh`, `project`, `march`, `info`, `check`, `export`, `view`. Routes to `_run_single_topology`, `_run_double_topology`, `_run_caps_topology`. |
| `surface/io.py` | `load_stl` — loads an STL file into a `trimesh.Trimesh`; validates it is a single mesh. |
| `surface/remesh.py` | `sphere_latlon_directions` builds unit direction vectors; `project_sphere_to_surface` does the ray-casting projection |
| `surface/yinyang.py` | `YANG_ROTATION` matrix (90-degree rotation); `project_yang_to_surface` applies the rotation then same projection logic |
| `volume/hyperbolic.py` | `march()` is the core extruder. `_compute_normals` and `_laplacian_smooth` accept `periodic=(periodic_i, periodic_j)` tuple; default `(False, True)`. |
| `volume/jacobian.py` | `jacobian_field` (periodic-j, shape `(ni-1, nj, nk-1)`), `_jacobian_field_open` (non-periodic, shape `(ni-1, nj-1, nk-1)`), `min_jacobian_two_layers` (two-layer Jacobian used during marching), `jacobian_per_layer` (min per layer, used by `check` command) |
| `topology/sphere_cap.py` | Splits O-grid into 6 blocks, builds cap seed faces via Laplacian convergence with Dirichlet BCs, marches all blocks simultaneously with shared-face enforcement |
| `io/geo.py` | ASCII geo format (Xall_GPU compatible): `read_geo`, `write_geo` |
| `io/vtk.py` | VTK XML StructuredGrid export: `write_vts`, `write_surface_vts`, `geo_to_pyvista`, `surface_to_pyvista` |
| `vis.py` | PyVista interactive visualisation: `view_geo`, `view_surface`, `save_screenshot`. Optional dep (`shore[vis]`); Jupyter backend via `jupyter=True`. |

### Test Suite

Tests in `tests/` use `conftest.py` helpers:
- `sphere_stl` — `pytest.fixture(scope="session")`, exports a unit icosphere STL once per session
- `sphere_layer(ni, nj, ...)` — plain function returning a `(ni, nj, 3)` lat-lon sphere layer
- `yang_layer(ni, nj, ...)` — plain function returning a `(ni, nj, 3)` Yang-rotated layer

Test files: `test_cli.py` (CLI via `typer.testing.CliRunner`), `test_hyperbolic.py`, `test_jacobian.py`, `test_sphere_cap.py`, `test_yinyang.py`, `test_geo.py`, `test_vtk.py`, `test_vis.py`. Coverage: geometric correctness (outer radius > inner, within 5%), Jacobian inversion detection, VTK round-trips, all topology variants.
