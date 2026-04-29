# Installation

## Requirements

- **Python** 3.11 or later
- **pip** or any PEP 517-compatible installer

All Python dependencies are installed automatically:

| Package | Role |
|---|---|
| `numpy >= 1.26` | Core array numerics |
| `trimesh >= 4.0` | STL I/O, BVH ray casting |
| `rtree >= 1.0` | Spatial index for trimesh ray backend |
| `scipy >= 1.12` | Spatial queries (fallback) |
| `typer >= 0.12` | CLI framework |
| `rich >= 13.0` | Progress bar and formatted output |

## From PyPI

```bash
pip install shore
```

Verify the installation:

```bash
shore --version
shore --help
```

## From source (development)

```bash
git clone https://github.com/szaghi/shore.git
cd shore
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

The `[dev]` extra adds `pytest`, `pytest-cov`, `ruff`, and `pyvista`.

Run the test suite to confirm everything works:

```bash
make test
```

Expected output (count grows with each release; the run completes in
under a minute on a modern laptop):

```
588 passed in 53s
```

## Optional: pyvista for visualisation

```bash
pip install shore[vis]
# or
pip install pyvista
```

The `shore view` command (and the example scripts' interactive mode) require pyvista. Without it, the rest of SHORE — generation, validation, VTK / `.geo` / `.grd` / `cc.par` export — works unchanged.

## Makefile targets

All development tasks are available via `make`:

```bash
make dev    # install in editable mode with dev extras
make test   # run test suite with coverage
make lint   # ruff check + format check (read-only)
make fmt    # ruff fix + format (destructive)
make clean  # remove build artefacts
```

::: tip WSL2 / externally managed Python
On WSL2 (and most modern Linux distributions) the system Python is externally managed — direct `pip install` without a virtual environment will fail. Always activate `.venv` before running any SHORE command:

```bash
source .venv/bin/activate
shore mesh ...
```
:::
