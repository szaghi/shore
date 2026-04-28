# Contributing

Contributions are welcome. This page covers the development workflow, code conventions, and what is most useful to work on.

## Development setup

```bash
git clone https://github.com/szaghi/shore.git
cd shore
python -m venv .venv && source .venv/bin/activate
make dev    # pip install -e ".[dev]"
make test   # should print "35 passed"
```

## Code conventions

- **Formatter / linter**: Ruff — `make fmt` to fix, `make lint` to check only
- **Type annotations**: all functions fully annotated; `from __future__ import annotations` in every file
- **Docstrings**: NumPy-style (`Parameters`, `Returns`, `Raises` sections)
- **Exceptions**: never raise bare `Exception`; use typed project exceptions or `RuntimeError` with a precise message
- **No print statements** in library code; CLI output goes through `rich.console.Console`

## Commit messages

Use [Conventional Commits](https://www.conventionalcommits.org/):

```
feat(surface): add STEP input via pythonOCC
fix(hyperbolic): correct Jacobian sign for non-lat-lon surface grids
docs(guide): add multi-block roadmap section
test(cli): add shore view command test
```

## Test conventions

- All tests live in `tests/`
- Shared fixtures and helpers in `tests/conftest.py`
- Mock nothing except unavoidable external I/O; the test suite runs against real geometry
- Use small grids (`ni=10, nj=12, nk=6`) to keep the suite under 1 s
- New features must include at least one geometric correctness test (not just exit-code checks)

## What to work on

The most impactful contributions, in rough priority order:

1. **Box outer boundary morphing** — the outermost layer currently matches the sphere template shape. Xall background grids are Cartesian, so the outer boundary should be morphed to a box. This is a post-processing step on `grid[:,:,-1,:]`.

2. **Steger–Sorenson orthogonality source term** — the current marching kernel uses pure normal extrusion. At high curvature the normals drift; the Steger–Sorenson correction keeps cells orthogonal. See: Steger & Sorenson (1979), *Automatic mesh-point clustering near a boundary in grid generation with elliptic partial differential equations*.

3. **`shore view` command** — pyvista-based visualisation of geo files. Should display the structured grid as a wireframe with a colourmap of the Jacobian.

4. **Multi-block surface decomposition** — decompose non-star-shaped bodies into structured patches. This is the hardest item and requires block topology reasoning.

5. **STEP / parametric surface input** — via pythonOCC or cadquery. The `surface/io.py` module should grow a `load_step()` function that returns a trimesh-compatible representation.

## Pull requests

- Open a draft PR early to discuss the approach before investing significant time
- Every PR must pass `make lint` and `make test`
- Add or update documentation in `docs/` for any user-facing change
