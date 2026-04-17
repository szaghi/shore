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

### Topology

| Option | Default | Description |
|---|---|---|
| `--topology TEXT` | `single` | `'single'`, `'double'` (Yin-Yang), or `'caps'` (6-block) |
| `--theta-cap FLOAT` | — | Polar exclusion angle in radians |
| `--theta-cap-deg FLOAT` | — | Same as `--theta-cap` but in degrees |
| `--cap-pull FLOAT` | 0.0 | [caps] Dome cap interior toward pole [0, 1) |
| `--cap-squareness FLOAT` | 0.0 | [caps] Straighten cap boundary edges [0, 1] |
| `--correct-spacing` | off | [caps] Equalise cell widths at all block interfaces after marching |

### Surface spacing

Surface spacing redistributes the angular positions of the projection rays. It is set at projection time (stage 1) and baked into the `.npz`.

| Option | Default | Description |
|---|---|---|
| `--surface-i-spacing TEXT` | `uniform` | Latitude spacing: `'uniform'`, `'tanh'` (cluster at cap boundaries), `'tanh2'` (both ends) |
| `--surface-i-beta FLOAT` | 3.0 | Clustering strength for latitude tanh/tanh2 |
| `--surface-j-spacing TEXT` | `uniform` | Longitude spacing: `'uniform'`, `'tanh2'` (cluster at a target longitude) |
| `--surface-j-beta FLOAT` | 3.0 | Clustering strength for longitude tanh2 |
| `--j-cluster-lon FLOAT` | — | Target longitude $\phi_0$ in radians where cells are finest (tanh2 only). `0.0` = +x, `pi/2` = +y, `pi` = -x. Default `pi` when omitted. |

### Volume spacing

| Option | Default | Description |
|---|---|---|
| `--volume-k-spacing TEXT` | `geometric` | Wall-normal spacing: `'geometric'`, `'tanh'`, `'tanh2'` |
| `--volume-k-beta FLOAT` | 3.0 | Clustering strength for tanh/tanh2 |
| `--volume-k-thickness FLOAT` | — | Total wall-normal thickness (required for tanh/tanh2) |
| `--ds FLOAT` | 0.001 | First-layer thickness (geometric spacing) |
| `--growth FLOAT` | 1.15 | Geometric growth ratio |

### Smoothing

| Option | Default | Description |
|---|---|---|
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
# → sphere.geo  (40x60x30, ds=1e-3, growth=1.15)
```

High-resolution boundary layer grid with verbose output:

```bash
shore mesh fuselage.stl \
  --ni 80 --nj 120 --nk 50 \
  --ds 1e-5 --growth 1.12 --smooth 0.3 \
  -o fuselage_bl.geo -v
```

Tanh-spaced boundary layer with fixed total thickness:

```bash
shore mesh wing.stl \
  --ni 60 --nj 90 --nk 40 \
  --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 2.0 \
  -o wing_tanh.geo
```

Latitude clustering at the cap boundaries (caps topology — finer cells where equatorial and polar blocks meet):

```bash
shore mesh body.stl \
  --topology caps --ni 60 --nj 80 --nk 30 \
  --surface-i-spacing tanh --surface-i-beta 4.0 \
  --cap-pull 0.4 \
  -o body_caps.geo
```

Caps topology with interface spacing equalisation:

```bash
shore mesh body.stl \
  --topology caps --ni 60 --nj 80 --nk 30 \
  --correct-spacing --cap-pull 0.4 --cap-squareness 0.5 \
  -o body_caps_corrected.geo
```

Longitude clustering toward the leading edge. `--j-cluster-lon 0.0` places the dense band at $\phi = 0$ (positive x-axis); cells are coarsest at the opposite side ($\phi = \pi$):

```bash
shore mesh airfoil.stl \
  --surface-j-spacing tanh2 --surface-j-beta 5.0 --j-cluster-lon 0.0 \
  --ni 40 --nj 80 --nk 30 \
  -o airfoil_le.geo
```

## Notes

- The STL must be a single connected, watertight mesh. Multi-body STL files (multiple shells) will fail during projection — split them first with a CAD tool or trimesh.
- `--nj` should be even to avoid degenerate cells at the poles of the sphere template.
- The outer boundary (k=nk-1) is currently sphere-shaped. Use `shore check` to verify the grid before passing it to Xall_GPU.
- Surface spacing is set at projection time (stage 1) and baked into the `.npz`; volume spacing is applied at march time (stage 2). Use `shore project` + `shore march` separately to experiment with volume spacing without re-projecting.
- See [Spacing laws](/concepts/spacing) for a detailed explanation of each distribution and when to use it.

## See also

- [`shore project`](project) — run stage 1 only
- [`shore march`](march) — run stage 2 only
- [`shore check`](check) — validate Jacobians
- [Spacing laws](/concepts/spacing) — clustering distributions explained
- [Quick Start](/guide/quickstart) — worked examples
- [Hyperbolic marching](/concepts/hyperbolic-marching) — algorithm details
