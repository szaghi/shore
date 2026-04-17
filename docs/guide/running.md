# Running Examples and Tests

This guide covers how to run SHORE's example scripts and the full test suite.

## Prerequisites

Ensure SHORE is installed in development mode with all extras:

```bash
git clone https://github.com/szaghi/shore.git
cd shore
python -m venv .venv && source .venv/bin/activate
make dev    # pip install -e ".[dev,vis]"
```

## Running Examples

SHORE ships two standalone example scripts that serve as both documentation and sanity checks:

### `examples/visualize_geo.py` — Single O-grid

The primary example generates a unit-sphere test case and optionally visualises it:

```bash
# Generate sphere.stl + sphere.geo and open the interactive viewer
python examples/visualize_geo.py --generate

# Generate and save a screenshot instead of opening a window
python examples/visualize_geo.py --generate --screenshot sphere_preview.png

# Generate only (skip the viewer)
python examples/visualize_geo.py --generate --no-view
```

Key parameters:

| Flag | Default | Description |
|---|---|---|
| `--ni` | 40 | Latitude points |
| `--nj` | 60 | Longitude points |
| `--nk` | 20 | Wall-normal layers |
| `--ds` | 0.02 | First-layer thickness (geometric) |
| `--growth` | 1.10 | Layer growth ratio (geometric) |
| `--smooth` | 0.2 | Laplacian blend strength [0–1] |
| `--smooth-iters` | 2 | Laplacian sweeps per layer |
| `--surface-i-spacing` | `uniform` | Latitude spacing: `uniform`, `tanh`, `tanh2` |
| `--surface-i-beta` | 3.0 | Clustering strength for latitude tanh/tanh2 |
| `--surface-j-spacing` | `uniform` | Longitude spacing: `uniform`, `tanh2` |
| `--surface-j-beta` | 3.0 | Clustering strength for longitude tanh2 |
| `--j-cluster-lon` | — | Longitude to cluster toward (radians) |
| `--volume-k-spacing` | `geometric` | Wall-normal spacing: `geometric`, `tanh`, `tanh2` |
| `--volume-k-beta` | 3.0 | Clustering strength for volume tanh/tanh2 |
| `--volume-k-thickness` | — | Total wall-normal thickness (required for tanh/tanh2) |

Spacing examples:

```bash
# Tanh volume spacing — fixed far-field boundary
python examples/visualize_geo.py --generate \
  --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 2.0

# Longitude clustering at the leading edge (phi=0)
python examples/visualize_geo.py --generate \
  --surface-j-spacing tanh2 --surface-j-beta 5.0 --j-cluster-lon 0.0

# Latitude clustering + tanh volume spacing
python examples/visualize_geo.py --generate \
  --surface-i-spacing tanh --surface-i-beta 4.0 \
  --volume-k-spacing tanh --volume-k-beta 3.0 --volume-k-thickness 1.5
```

See [Spacing laws](/concepts/spacing) for an explanation of each distribution and what physical region it clusters toward.

Visualise an existing `.geo` or `.npz` file directly:

```bash
python examples/visualize_geo.py sphere.geo
python examples/visualize_geo.py sphere_k0.npz
python examples/visualize_geo.py sphere.geo --screenshot out.png --no-edges
```

### `examples/visualize_geo.py` — Double O-grid (Yin-Yang)

Generate both Yin and Yang blocks for a complete sphere without polar singularities:

```bash
# Generate sphere_yin.geo and sphere_yang.geo
python examples/visualize_geo.py --generate --topology double

# With custom grid parameters
python examples/visualize_geo.py --generate --topology double --ni 40 --nj 60 --nk 20

# Validate both blocks
shore check sphere_yin.geo
shore check sphere_yang.geo
```

For the double topology, two files are produced:

- `sphere_yin.geo` — standard O-grid block
- `sphere_yang.geo` — rotated O-grid block covering the complementary polar region

Both blocks share the same `(ni, nj, nk)` dimensions and can be exported to ParaView individually or as a combined MultiBlock (see VTK export).

### `examples/visualize.ipynb` — Jupyter Notebook

A Jupyter notebook covering the full pipeline, Jacobian analysis, and spacing-law comparison:

- Sections 0–8: generation, inspection, Jacobian profiles, VTK export
- Sections 9–12: surface and volume spacing distributions, tanh vs. geometric comparison

```bash
# Launch Jupyter with the notebook
jupyter notebook examples/visualize.ipynb
```

Or run headlessly for CI environments:

```python
import pyvista as pv
pv.set_jupyter_backend("static")   # no display required

from shore.vis import view_geo
from shore.io.geo import read_geo

grid = read_geo("sphere.geo")
pl = view_geo(grid, jupyter=True)  # returns a plotter, does not call .show()
```

## Running the Test Suite

### Full suite

```bash
make test
```

Or directly via pytest:

```bash
source .venv/bin/activate
python -m pytest tests/ -v
```

Expected output:

```
84 passed
```

All 84 tests must pass before any PR is merged.

### Running a specific test file

```bash
# Unit tests for specific modules
python -m pytest tests/test_hyperbolic.py -v      # marching kernel
python -m pytest tests/test_jacobian.py -v        # Jacobian calculations
python -m pytest tests/test_yinyang.py -v         # Yin-Yang double O-grid
python -m pytest tests/test_cli.py -v            # CLI commands
python -m pytest tests/test_vtk.py -v            # VTK export

# Run from the project root (pytest needs the src/ path)
cd /path/to/shore
python -m pytest tests/test_geo.py tests/test_vis.py -v
```

### Running a single test

```bash
python -m pytest tests/test_cli.py::TestMesh::test_exit_code_ok -v
```

### Test file descriptions

| File | What it tests |
|---|---|
| `tests/conftest.py` | Shared fixtures (`sphere_stl`, `sphere_layer`, `yang_layer`) |
| `tests/test_geo.py` | Geo file read/write round-trip fidelity |
| `tests/test_vtk.py` | VTK XML export correctness |
| `tests/test_vis.py` | PyVista conversion helpers |
| `tests/test_jacobian.py` | Jacobian field calculations and inversion detection |
| `tests/test_hyperbolic.py` | Normal estimation, Laplacian smoothing, marching |
| `tests/test_yinyang.py` | Yin-Yang rotation, coverage masks, combined coverage |
| `tests/test_cli.py` | All CLI commands end-to-end |
| `tests/test_spacing.py` | 1D spacing laws, domain builders, error cases |

### Filtering tests by keyword

```bash
# Run only tests containing "sphere" in the name
python -m pytest tests/ -k sphere -v

# Run only tests for the Yin-Yang module
python -m pytest tests/ -k yinyang -v

# Skip slow tests (if marked)
python -m pytest tests/ --ignore=tests/test_performance.py -v
```

### Test parameters

SHORE's test suite uses small default grids (`ni=10, nj=12, nk=6`) to keep execution fast. The full suite runs in under 1 second on a modern laptop.

For coverage reporting:

```bash
python -m pytest tests/ --cov=shore --cov-report=html -v
```

HTML coverage report is written to `htmlcov/index.html`.

## Makefile Targets

All development tasks are available via `make`:

```bash
make dev    # install package in editable mode with dev + vis extras
make test   # run test suite with pytest
make lint   # ruff check + format check (read-only)
make fmt    # ruff fix + format (destructive — auto-fixes issues)
make clean  # remove build artefacts, .ruff_cache, .pytest_cache
```

## Troubleshooting

### pytest collects no tests

If you see `No tests collected`, ensure:

1. You are running from the project root (not inside the `tests/` directory)
2. The `.venv` is activated: `source .venv/bin/activate`
3. pytest is installed: `pip show pytest`

```bash
# Verify setup
source .venv/bin/activate
python -c "import sys; print(sys.executable)"   # should point to .venv/bin/python
python -m pytest --version                       # should show pytest version
python -c "import tests.conftest; print('conftest OK')"  # fixture import check
```

### Virtual environment issues (WSL2)

On WSL2 the system Python is externally managed. Always use `.venv`:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

If you see `error: externally-managed-environment`, the system Python is being used instead of `.venv`.

### Visualisation not available

The `shore view` command and `examples/visualize_geo.py` interactive mode require pyvista:

```bash
pip install shore[vis]
# or
pip install pyvista
```

Without pyvista, the `--screenshot` and `--no-view` flags still work for offline generation.

## CI / Automated Environments

For headless CI (no display required):

```bash
# Install
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

# Lint
make lint

# Test (headless)
export DISPLAY=  # suppress pyvista display attempts
make test

# Generate example files only (no viewer)
python examples/visualize_geo.py --generate --no-view
shore check sphere.geo
shore info sphere.geo
```