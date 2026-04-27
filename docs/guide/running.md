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

SHORE ships a set of standalone example scripts and a Jupyter
notebook that serve as both documentation and sanity checks.

### `examples/visualize_cubed_sphere.py` — cubed-sphere topology around a sphere

The script generates a unit-sphere body-fitted test case (always
6-block cubed-sphere topology) and optionally visualises it.  When called
*without* `--generate`, it acts as a generic `.geo` / `.npz` viewer.

```bash
# Generate sphere.stl + 6 sphere_<label>.geo files and open the viewer
python examples/visualize_cubed_sphere.py --generate

# Generate and save a screenshot instead of opening a window
python examples/visualize_cubed_sphere.py --generate --screenshot sphere_preview.png

# Generate only (skip the viewer)
python examples/visualize_cubed_sphere.py --generate --no-view

# Export to a single VTK MultiBlock .vtm for ParaView
python examples/visualize_cubed_sphere.py --generate --vtk
```

Key parameters:

| Flag | Default | Description |
|---|---|---|
| `--ni` | 40 | Equator meridional resolution |
| `--nj` | 60 | Longitude resolution (multiple of 4) |
| `--nk` | 20 | Wall-normal layers |
| `--ds` | 0.02 | First-layer thickness (geometric) |
| `--growth` | 1.10 | Layer growth ratio (geometric) |
| `--theta-cap` | π/6 | Polar exclusion angle in radians |
| `--theta-cap-deg` | — | Polar exclusion angle in degrees |
| `--smooth` | 0.2 | Laplacian blend strength [0–1] |
| `--smooth-iters` | 2 | Laplacian sweeps per layer |
| `--surface-i-spacing` | `uniform` | Equator meridional spacing law |
| `--surface-i-beta` | 3.0 | Clustering strength for tanh/tanh2 i-spacing |
| `--surface-j-spacing` | `uniform` | Cap-edge tangent-plane spacing law |
| `--surface-j-beta` | 3.0 | Clustering strength for tanh2 j-spacing |
| `--volume-k-spacing` | `geometric` | Wall-normal spacing law |
| `--volume-k-beta` | 3.0 | Clustering strength for tanh/tanh2 k-spacing |
| `--volume-k-thickness` | — | Total wall-normal thickness (required for tanh/tanh2) |
| `--blend-normals-k` | 0 | Layers over which to blend STL normals (0 = off) |

Spacing examples:

```bash
# Tanh wall-normal spacing — fixed far-field boundary
python examples/visualize_cubed_sphere.py --generate \
  --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 2.0

# Equator meridional clustering at the cap seams
python examples/visualize_cubed_sphere.py --generate \
  --surface-i-spacing tanh --surface-i-beta 4.0

# Cap-corner clustering (tanh2 j-spacing — inherited by equator seam)
python examples/visualize_cubed_sphere.py --generate \
  --surface-j-spacing tanh2 --surface-j-beta 5.0
```

See [Spacing laws](/algorithm/spacing-laws) for an explanation of each
distribution and what physical region it clusters toward.

Visualise an existing `.geo` or `.npz` file directly:

```bash
python examples/visualize_cubed_sphere.py sphere_sub0.geo
python examples/visualize_cubed_sphere.py sphere_k0.npz
python examples/visualize_cubed_sphere.py sphere_sub0.geo --screenshot out.png --no-edges
```

Output: 6 `.geo` files per run

```
sphere_sub0.geo .. sphere_sub3.geo   # equatorial sub-blocks
sphere_cap_north.geo                 # north polar cap
sphere_cap_south.geo                 # south polar cap
```

### `examples/config_fuselage.toml` — config-driven workflow

When the parameter set gets large enough that the CLI invocation
becomes a wall of flags, switch to a TOML config file. Every CLI flag
on `shore mesh` has an equivalent TOML key; the file is parsed,
validated, and merged with any explicit CLI flags before the run
starts (CLI flag > config value > built-in default).

```bash
# bootstrap a fresh config from defaults
shore config template > shore.toml

# or start from the annotated example
cp examples/config_fuselage.toml my_case.toml

# parse-check before running (catches typos, version mismatches, etc.)
shore config validate my_case.toml

# run the full pipeline; INPUT is taken from the config's `input` key
shore mesh -c my_case.toml -v

# CLI overrides still apply: this picks up everything from the file
# but uses --nk 80 instead of the configured value.
shore mesh -c my_case.toml --nk 80
```

The example config (`examples/config_fuselage.toml`) demonstrates
every section with a non-default choice — tanh wall-normal spacing
with fixed thickness, equator clustering at the cap seams, blended
STL normals in the near-wall region — so you can see how the knobs
combine in a single recipe.

See the full schema and precedence rules in
[Config files](/reference/config).

### Primitive examples

For the parametric primitives (box, annulus, flat-caps) there is one
standalone script per primitive plus a matching annotated TOML config.
Each script generates the primitive in-memory, writes a `.geo`
(or several, for multi-block primitives) **and a `.vtm` MultiBlock**
for ParaView, and optionally opens the PyVista viewer.

| Script | Primitive | Config | Output |
|---|---|---|---|
| `examples/visualize_box.py` | Cartesian box | `config_primitive_box.toml` | `box.geo`, `box.vtm` |
| `examples/visualize_annulus.py` | annular cylinder | `config_primitive_annulus.toml` | `annulus.geo`, `annulus.vtm` |
| `examples/visualize_flat_caps.py` | 5-block flat-caps cylinder | `config_primitive_flat_caps.toml` | `flat_caps_<label>.geo` ×5, `flat_caps.vtm` |

Each script accepts the same shape of options:

```bash
# Defaults — open the viewer
python examples/visualize_box.py

# Skip the viewer (CI-friendly)
python examples/visualize_box.py --no-view

# Save a screenshot instead of opening a window
python examples/visualize_box.py --screenshot box.png

# Drive the parameters from a TOML config (CLI flags still override)
python examples/visualize_annulus.py --config examples/config_primitive_annulus.toml --no-view
python examples/visualize_flat_caps.py --config examples/config_primitive_flat_caps.toml --no-view

# Skip the .vtm export (only .geo is written)
python examples/visualize_box.py --no-vtk --no-view
```

The `.vtm` file is a small XML manifest pointing at one `.vts` per
block stored under a sister directory (`box/`, `annulus/`, or
`flat_caps/`).  Open the `.vtm` in ParaView; it pulls in every block
automatically.

### `examples/visualize.ipynb` — Jupyter Notebook

A complete notebook covering mesh generation, Jacobian analysis,
spacing-law comparison, and the C1 seam-matching demo:

- Sections 0–7: generation, inspection, Jacobian profiles, VTK export
- Sections 8–9: surface and volume spacing distributions
- Section 10: cap/equator C1 seam matching, with side-by-side comparison
  of compatible and incompatible parameter regimes

```bash
jupyter notebook examples/visualize.ipynb
```

For headless CI environments:

```python
import pyvista as pv
pv.set_jupyter_backend("static")   # no display required
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

### Running a specific test file

```bash
python -m pytest tests/test_cap_gnomonic.py -v   # gnomonic cap construction + seam C1
python -m pytest tests/test_jacobian.py -v       # Jacobian calculations
python -m pytest tests/test_cli.py -v            # CLI commands
python -m pytest tests/test_vtk.py -v            # VTK export
```

### Running a single test

```bash
python -m pytest tests/test_cli.py::TestMesh::test_exit_code_ok -v
```

### Test file descriptions

| File | What it tests |
|---|---|
| `tests/conftest.py` | Shared fixtures (`sphere_stl`, `sphere_layer`) |
| `tests/test_geo.py` | Geo file read/write round-trip fidelity |
| `tests/test_vtk.py` | VTK XML export correctness |
| `tests/test_vis.py` | PyVista conversion helpers |
| `tests/test_jacobian.py` | Jacobian field calculations and inversion detection |
| `tests/test_cap_gnomonic.py` | Gnomonic cap, seam C0/C1, volume seam continuity |
| `tests/test_mesh.py` | OOP `Mesh` API, `CubedSphereTopology.build` |
| `tests/test_anchor.py` / `test_load_stl_hygiene.py` | Mesh preconditioning |
| `tests/test_projector.py` / `test_projector_validation.py` | Ray-cast Projector |
| `tests/test_cli.py` | All CLI commands end-to-end |
| `tests/test_spacing.py` | 1D spacing laws, domain builders, error cases |

### Filtering tests by keyword

```bash
# Run only seam-related tests
python -m pytest tests/ -k seam -v

# Run only the gnomonic cap tests
python -m pytest tests/test_cap_gnomonic.py -v
```

### Test parameters

The test suite uses small default grids to keep execution fast. The
full suite runs in under 5 seconds on a modern laptop.

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

1. You are running from the project root (not inside `tests/`).
2. The `.venv` is activated: `source .venv/bin/activate`.
3. pytest is installed: `pip show pytest`.

### Virtual environment issues (WSL2)

On WSL2 the system Python is externally managed. Always use `.venv`:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

### Visualisation not available

The `shore view` command and the example script's interactive mode
require pyvista:

```bash
pip install shore[vis]
```

Without pyvista, `--screenshot` and `--no-view` still work for offline
generation.

## CI / Automated Environments

For headless CI (no display required):

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

# Lint
make lint

# Test (headless)
export DISPLAY=
make test

# Smoke-generate the example mesh
python examples/visualize_cubed_sphere.py --generate --no-view
shore check sphere_sub0.geo
shore info sphere_sub0.geo
```
