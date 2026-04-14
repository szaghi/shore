# shore mesh

Full pipeline: STL → structured O-grid volume → geo.

Equivalent to running `shore project` followed by `shore march` in sequence.

## Synopsis

```
shore mesh [OPTIONS] INPUT
```

## Arguments

| Argument | Description |
|---|---|
| `INPUT` | Path to the input `.stl` file (binary or ASCII) |

## Options

### Grid resolution

| Option | Default | Description |
|---|---|---|
| `--ni INTEGER` | 40 | Latitude resolution (pole-to-pole, exclusive of poles) |
| `--nj INTEGER` | 60 | Longitude resolution (full circle, periodic) |
| `--nk INTEGER` | 30 | Wall-normal layers, including k=0 (wall) |

### Marching parameters

| Option | Default | Description |
|---|---|---|
| `--ds FLOAT` | 0.001 | First-layer thickness in body length units |
| `--growth FLOAT` | 1.15 | Geometric growth ratio between consecutive layers |
| `--smooth FLOAT` | 0.2 | Laplacian blend strength per layer [0–1] |
| `--smooth-iters INTEGER` | 2 | Laplacian sweeps applied after each layer |

### I/O and display

| Option | Default | Description |
|---|---|---|
| `-o` / `--output PATH` | `<INPUT stem>.geo` | Output geo file path |
| `-v` / `--verbose` | off | Show Rich progress bar with per-layer stats |
| `--help` | — | Show this help and exit |

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Input file not found, projection failed, or negative Jacobian during march |

## Examples

Minimal invocation — defaults for all grid parameters:

```bash
shore mesh sphere.stl
# → sphere.geo  (40×60×30, ds=1e-3, growth=1.15)
```

High-resolution boundary layer grid with verbose output:

```bash
shore mesh fuselage.stl \
  --ni 80 --nj 120 --nk 50 \
  --ds 1e-5 --growth 1.12 --smooth 0.3 \
  -o fuselage_bl.geo -v
```

Custom output name:

```bash
shore mesh wing.stl --ni 60 --nj 90 --nk 40 -o grids/wing_ogrid.geo
```

## Notes

- The STL must be a single connected, watertight mesh. Multi-body STL files (multiple shells) will fail during projection — split them first with a CAD tool or trimesh.
- `--nj` should be even to avoid degenerate cells at the poles of the sphere template.
- The outer boundary (k=nk-1) is currently sphere-shaped. Use `shore check` to verify the grid before passing it to Xall_GPU.

## See also

- [`shore project`](project) — run stage 1 only
- [`shore march`](march) — run stage 2 only
- [`shore check`](check) — validate Jacobians
- [Quick Start](/guide/quickstart) — worked examples
- [Hyperbolic marching](/concepts/hyperbolic-marching) — algorithm details
