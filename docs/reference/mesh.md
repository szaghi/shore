# shore mesh

Full pipeline: STL → 6-block O-grid volume → geo files.

The `mesh` command always builds the [caps topology](/topologies/caps)
(4 equatorial sub-blocks + 2 polar caps, 6 `.geo` files written). The cap
k=0 face is a gnomonic flat-square projection; the equator is built from
the cap seam rings with C1 spacing match at every seam vertex; marching
uses seam-aware normals so multi-block edge spacing stays continuous
across all layers.

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
| `--ni INTEGER` | 40 | Equator meridional resolution (south seam → north seam) |
| `--nj INTEGER` | 60 | Longitude resolution (multiple of 4; rounded down with a warning if not) |
| `--nk INTEGER` | 30 | Wall-normal layers, including k=0 (wall) |

### Cap geometry

| Option | Default | Description |
|---|---|---|
| `--theta-cap FLOAT` | $\pi/6$ | Polar exclusion angle in radians (cap corners sit at this colatitude) |
| `--theta-cap-deg FLOAT` | — | Same as `--theta-cap` but in degrees (mutually exclusive) |

::: tip Cap / equator compatibility
Choose `ni`, `nj`, `theta_cap` so the cap radial step
`theta_cap / (nc − 1)` (with `nc = nj // 4 + 1`) is within ~2× of the
equator's natural radial step `(π − 2*theta_cap) / (ni − 1)`. This keeps
the per-j C1 seam pinning active. When out of range, the equator falls
back to its unpinned `i_spacing` law and emits a single `UserWarning`.
:::

### Surface spacing

| Option | Default | Description |
|---|---|---|
| `--surface-i-spacing TEXT` | `uniform` | Equator meridional spacing: `'uniform'`, `'tanh'`, `'tanh2'` |
| `--surface-i-beta FLOAT` | 3.0 | Clustering strength for tanh/tanh2 i-spacing |
| `--surface-j-spacing TEXT` | `uniform` | Cap-edge tangent-plane spacing: `'uniform'` or `'tanh2'`; inherited by the equator seam |
| `--surface-j-beta FLOAT` | 3.0 | Clustering strength for tanh2 j-spacing |

### Volume spacing

| Option | Default | Description |
|---|---|---|
| `--volume-k-spacing TEXT` | `geometric` | Wall-normal spacing: `'geometric'`, `'tanh'`, `'tanh2'` |
| `--volume-k-beta FLOAT` | 3.0 | Clustering strength for tanh/tanh2 |
| `--volume-k-thickness FLOAT` | — | Total wall-normal thickness (required for tanh/tanh2) |
| `--ds FLOAT` | 0.001 | First-layer thickness (geometric spacing) |
| `--growth FLOAT` | 1.15 | Geometric growth ratio |

### Smoothing & normals

| Option | Default | Description |
|---|---|---|
| `--smooth FLOAT` | 0.2 | Laplacian blend strength per layer [0–1] |
| `--smooth-iters INTEGER` | 2 | Laplacian sweeps applied after each layer |
| `--blend-normals-k INTEGER` | 0 | Number of layers over which to blend STL surface normals into the grid-computed normals (0 = off; useful for sharp-feature bodies) |

### I/O and display

| Option | Default | Description |
|---|---|---|
| `-o` / `--output PATH` | `<INPUT stem>` | Output stem; writes `<stem>_sub0.geo` ... `<stem>_cap_south.geo` (6 files) |
| `-v` / `--verbose` | off | Show Rich progress bar with per-layer stats |
| `--help` | — | Show this help and exit |

### Deprecated

| Option | Status |
|---|---|
| `--cap-pull FLOAT` | Ignored. Cap interior is now built by gnomonic flat-square projection — square cells by construction. |
| `--cap-squareness FLOAT` | Ignored. Same reason. |
| `--j-cluster-lon FLOAT` | Ignored. The gnomonic cap has 4-fold symmetry; use `--surface-j-spacing tanh2` to cluster nodes near the cap-corner longitudes. |

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Input file not found, projection failed, or negative Jacobian during march |

## Output

The command writes 6 `.geo` files (one per block):

```
<stem>_sub0.geo        sub-block 0 (phi 0°→90°)         (ni × nc × nk)
<stem>_sub1.geo        sub-block 1 (phi 90°→180°)       (ni × nc × nk)
<stem>_sub2.geo        sub-block 2 (phi 180°→270°)      (ni × nc × nk)
<stem>_sub3.geo        sub-block 3 (phi 270°→360°)      (ni × nc × nk)
<stem>_cap_north.geo   north polar cap                  (nc × nc × nk)
<stem>_cap_south.geo   south polar cap                  (nc × nc × nk)
```

with `nc = nj // 4 + 1`.

## Examples

Minimal invocation — defaults for all grid parameters:

```bash
shore mesh sphere.stl
# → sphere_sub*.geo, sphere_cap_*.geo  (40 x 60 x 30, theta_cap = pi/6)
```

High-resolution boundary layer with verbose output:

```bash
shore mesh fuselage.stl \
  --ni 80 --nj 120 --nk 50 \
  --ds 1e-5 --growth 1.12 --smooth 0.3 \
  -o fuselage_bl -v
```

Tanh-spaced boundary layer with fixed total thickness:

```bash
shore mesh wing.stl \
  --ni 60 --nj 80 --nk 40 \
  --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 2.0 \
  -o wing_tanh
```

Equator meridional clustering at the cap seams:

```bash
shore mesh body.stl \
  --ni 60 --nj 80 --nk 30 \
  --surface-i-spacing tanh --surface-i-beta 4.0 \
  --theta-cap-deg 30 \
  -o body
```

Cap-corner refinement (j-spacing tanh2):

```bash
shore mesh body.stl \
  --ni 40 --nj 60 --nk 30 \
  --surface-j-spacing tanh2 --surface-j-beta 4.0 \
  --theta-cap-deg 30 \
  -o body
```

STL-normal blending for sharp features (first 5 layers blend toward the
true surface normal of the underlying triangulation):

```bash
shore mesh airfoil.stl \
  --ni 40 --nj 60 --nk 30 \
  --blend-normals-k 5 \
  -o airfoil
```

## Notes

- The STL must be a single connected, watertight mesh. Multi-body STL files
  (multiple shells) will fail during projection — split them first with a
  CAD tool or trimesh.
- `--nj` must be a multiple of 4 (otherwise rounded down with a warning,
  so the equatorial band can be split into 4 equal sectors).
- The outer boundary (k=nk-1) is currently sphere-shaped after marching.
  Use `shore check` to verify the per-block Jacobians.
- See [Spacing laws](/algorithm/spacing-laws) for a detailed explanation of each
  distribution and when to use it.

## See also

- [`shore project`](project) — generate the k=0 surface only (single-block,
  legacy `body.project` path; not the caps topology)
- [`shore check`](check) — validate Jacobians on existing `.geo` files
- [`shore export`](export) — convert `.geo` to VTK MultiBlock for ParaView
- [Caps topology](/topologies/caps) — block layout, cap construction,
  seam matching
- [Spacing laws](/algorithm/spacing-laws) — clustering distributions
- [Hyperbolic marching](/algorithm/hyperbolic-marching) — algorithm details
