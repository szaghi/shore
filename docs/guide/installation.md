# Installation

## Requirements

- **Python** 3.11 or later
- A POSIX shell for the one-liner (Linux / macOS / WSL); Windows users
  should run the `pipx install` command directly.

All Python dependencies are installed automatically:

| Package | Role |
|---|---|
| `numpy >= 1.26` | Core array numerics |
| `trimesh >= 4.0` | STL I/O, BVH ray casting |
| `rtree >= 1.0` | Spatial index for trimesh ray backend |
| `networkx >= 3.0` | Normal repair on STL load (trimesh) |
| `scipy >= 1.12` | Sparse solvers (surface marching, elliptic smoothing), spatial queries |
| `typer >= 0.12` | CLI framework |
| `rich >= 13.0` | Progress bar and formatted output |

::: tip Distribution name vs import name
The PyPI distribution name is `shore-mesh` (the unprefixed `shore` is
squatted on PyPI). The Python import name and CLI entry point both
stay `shore`:

```python
import shore
```

```bash
shore --help
```
:::

## One-liner (recommended for CLI use)

```bash
curl -fsSL https://raw.githubusercontent.com/szaghi/shore/main/install.sh | bash
```

The installer:

- Verifies Python ≥ 3.11 is available.
- Installs [pipx](https://pipx.pypa.io/) if missing — via `apt` on
  Debian/Ubuntu (PEP 668-compliant), via `pip --user` elsewhere.
- Runs `pipx install shore-mesh`, which puts SHORE in its own venv at
  `~/.local/share/pipx/` and symlinks the `shore` CLI to
  `~/.local/bin/shore`.
- Verifies the install and prints next steps.

Optional flags:

```bash
# install with the [vis] extras (pyvista for interactive viewers):
curl -fsSL https://raw.githubusercontent.com/szaghi/shore/main/install.sh | bash -s -- --vis

# pin a specific version:
curl -fsSL https://raw.githubusercontent.com/szaghi/shore/main/install.sh | bash -s -- --version 0.1.0
```

To upgrade later: `pipx upgrade shore-mesh`. To uninstall: `pipx uninstall shore-mesh`.

## Optional extras

The core install meshes STL bodies and writes every Xall file. Extras add
the rest:

| Extra | Adds | Needed for |
|---|---|---|
| `cad` | `build123d` (OpenCascade) | STEP / IGES / BREP geometry, exact projection, CAD faces and edges |
| `gui` | PySide6, pygfx, wgpu, rendercanvas | `shore gui` |
| `cgns` | `h5py` | `shore cgns`, `Project.write_cgns` |
| `vis` | PyVista | `shore view`, `shore export`, `shore boxes-export`, `shore.vis` |

Combine them as usual: `pipx install 'shore-mesh[cad,gui,cgns,vis]'`.

## With pipx (manual)

If you already have pipx and don't want a wrapper script:

```bash
pipx install shore-mesh
# or with the [vis] extras
pipx install 'shore-mesh[vis]'
# CGNS export (shore cgns, Project.write_cgns) needs the [cgns] extra (h5py)
pipx install 'shore-mesh[cgns]'
# the desktop GUI (shore gui) needs the [gui] extra; on Windows use native Python, not WSL2
pipx install 'shore-mesh[gui]'
```

This is what the one-liner ends up running. Same isolation
guarantees, same upgrade / uninstall commands.

## With pip (library use)

If you want to **import shore** from your own scripts or notebooks —
not just run the CLI — install into a virtual environment with pip:

```bash
python -m venv .venv && source .venv/bin/activate
pip install shore-mesh
```

::: warning Modern Linux distributions (PEP 668)
On WSL2, recent Ubuntu/Debian, and similar PEP 668 systems, the
system Python is **externally managed** — running `pip install`
without a virtual environment will fail. Always create a `.venv`
first, or use `pipx` (which manages its own venv for you).
:::

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

The `[dev]` extra adds `pytest`, `pytest-cov`, `ruff`, `pyvista`, `h5py` and
`pycgns` (the independent CGNS reader of the export specs, ADR 0035). For
the CAD and GUI tests add those extras: `pip install -e ".[dev,cad,gui]"`.

Run the test suite to confirm everything works:

```bash
make test                                         # about 2000 tests, ~9 min
QT_QPA_PLATFORM=offscreen pytest tests_gui --no-cov   # the GUI tests, needs [gui]
```

[Examples and tests](/guide/running) describes the suites and the frozen
acceptance specs.

## Optional: pyvista for visualisation

The `shore view` command (and the example scripts' interactive mode)
require pyvista. Pick whichever install path matches how you installed
SHORE:

```bash
# one-liner
curl -fsSL https://raw.githubusercontent.com/szaghi/shore/main/install.sh | bash -s -- --vis

# pipx
pipx install 'shore-mesh[vis]'

# pip (inside a venv)
pip install 'shore-mesh[vis]'
```

Without pyvista, the rest of SHORE — generation, validation,
VTK / `.geo` / `.grd` / `cc.par` export — works unchanged.

## Makefile targets

All development tasks are available via `make`:

```bash
make dev    # install in editable mode with dev extras
make test   # run test suite with coverage
make lint   # ruff check + format check (read-only)
make fmt    # ruff fix + format (destructive)
make clean  # remove build artefacts
make docs-api      # regenerate the Python API reference from the docstrings
make docs-figures  # redraw the docs figures (light and dark)
make docs-build    # build the docs site (fails on dead links)
```
