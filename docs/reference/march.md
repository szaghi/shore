# shore march

Extrude a structured surface layer (`.npz` from `shore project`) into a structured volume (`.geo`).

## Synopsis

```
shore march [OPTIONS] INPUT
```

## Arguments

| Argument | Description |
|---|---|
| `INPUT` | Path to a `.npz` file produced by `shore project` |

## Options

### Grid resolution

| Option | Default | Description |
|---|---|---|
| `--nk INTEGER` | 30 | Wall-normal layers, including k=0 |

### Volume spacing

| Option | Default | Description |
|---|---|---|
| `--volume-k-spacing TEXT` | `geometric` | Wall-normal spacing law: `'geometric'` (cluster at wall, unbounded outer), `'tanh'` (cluster at wall, fixed outer), `'tanh2'` (cluster at wall and far-field, fixed outer) |
| `--volume-k-beta FLOAT` | 3.0 | Clustering strength for tanh/tanh2 |
| `--volume-k-thickness FLOAT` | — | Total wall-normal thickness (required for tanh/tanh2 — sets the outer boundary distance) |
| `--ds FLOAT` | 0.001 | First-layer thickness (geometric only) |
| `--growth FLOAT` | 1.15 | Geometric growth ratio (geometric only) |

### Smoothing

| Option | Default | Description |
|---|---|---|
| `--smooth FLOAT` | 0.2 | Laplacian blend strength [0–1] |
| `--smooth-iters INTEGER` | 2 | Laplacian sweeps per layer |

### I/O and display

| Option | Default | Description |
|---|---|---|
| `-o` / `--output PATH` | `<INPUT stem>.geo` | Output `.geo` file |
| `-v` / `--verbose` | off | Show Rich progress bar |
| `--help` | — | Show this help and exit |

## Spacing laws

The `--volume-k-spacing` flag selects how wall-normal step sizes are distributed. The "start" of the interval is always the body wall (k=0); the "end" is the far-field boundary (k=nk-1).

| Law | Dense region | Outer boundary | Requires |
|---|---|---|---|
| `geometric` (default) | Wall (k=0) | Unbounded — depends on `ds`, `growth`, `nk` | `--ds`, `--growth` |
| `tanh` | Wall (k=0) | Fixed at `--volume-k-thickness` | `--volume-k-thickness` |
| `tanh2` | Wall (k=0) **and** far-field (k=nk-1) | Fixed at `--volume-k-thickness` | `--volume-k-thickness` |

When using `tanh` or `tanh2`, `--ds` and `--growth` are ignored; the total extent is fixed by `--volume-k-thickness` and the distribution is controlled by `--volume-k-beta`. Use `geometric` when you need to control first-layer thickness for $y^+$ resolution; use `tanh`/`tanh2` when you need a known far-field distance.

See [Spacing laws](/concepts/spacing) for the mathematical details and visual comparison.

## Progress bar

With `-v`, a live bar is displayed during marching:

```
Marching       12/29  jmin=3.41e-04  ds=4.05e-03  0:00:01
```

Fields:

| Field | Description |
|---|---|
| `12/29` | Completed layers / total layers to march |
| `jmin=3.41e-04` | Minimum cell Jacobian in the last completed layer |
| `ds=4.05e-03` | Physical thickness of the last completed layer |

A `jmin` value dropping toward zero is a warning sign — increase `--smooth` or reduce the layer thickness.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Input file not found, invalid `.npz` (missing `k0` key), missing `--volume-k-thickness`, or negative Jacobian detected |

## Examples

Basic march with geometric spacing:

```bash
shore march sphere_k0.npz --nk 30 -o sphere.geo
```

Fine boundary-layer grid:

```bash
shore march fuselage_k0.npz --nk 50 --ds 1e-5 --growth 1.10 --smooth 0.3 -v
```

Tanh-spaced volume with fixed total thickness (good for CFD far-field control):

```bash
shore march k0.npz \
  --nk 60 --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 5.0 \
  -o grid_tanh.geo
```

Two-sided tanh for a bounded domain (wall-to-far-field):

```bash
shore march k0.npz \
  --nk 40 --volume-k-spacing tanh2 --volume-k-beta 3.0 --volume-k-thickness 3.0 \
  -o grid_tanh2.geo
```

Experiment with different growth ratios without re-projecting:

```bash
shore march k0.npz --nk 40 --ds 1e-3 --growth 1.10 -o grid_r110.geo
shore march k0.npz --nk 40 --ds 1e-3 --growth 1.20 -o grid_r120.geo
shore check grid_r110.geo
shore check grid_r120.geo
```

## Error handling

If a negative Jacobian is detected, marching aborts:

```
Error: Negative Jacobian (-1.23e-04) detected at layer k=7.
Reduce ds, growth ratio, or increase smooth_strength.
```

Typical remedies:

1. Reduce `--ds` (smaller first layer)
2. Reduce `--growth` (slower layer expansion)
3. Increase `--smooth` (more Laplacian damping)
4. Increase `--smooth-iters` (more sweeps per layer)
5. Increase `--volume-k-beta` (tighter clustering at wall for tanh/tanh2)

## See also

- [`shore project`](project) — produce the `.npz` input
- [`shore check`](check) — validate the output
- [Spacing laws](/concepts/spacing) — clustering distributions explained
- [Hyperbolic marching](/concepts/hyperbolic-marching) — algorithm and failure modes
