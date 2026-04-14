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

| Option | Default | Description |
|---|---|---|
| `--nk INTEGER` | 30 | Wall-normal layers, including k=0 |
| `--ds FLOAT` | 0.001 | First-layer thickness |
| `--growth FLOAT` | 1.15 | Geometric growth ratio |
| `--smooth FLOAT` | 0.2 | Laplacian blend strength [0–1] |
| `--smooth-iters INTEGER` | 2 | Laplacian sweeps per layer |
| `-o` / `--output PATH` | `<INPUT stem>.geo` | Output `.geo` file |
| `-v` / `--verbose` | off | Show Rich progress bar |
| `--help` | — | Show this help and exit |

## Progress bar

With `-v`, a live bar is displayed during marching:

```
⠸ Marching       ████████████░░░░░░░░  12/29  jmin=3.41e-04  Δs=4.05e-03  0:00:01
```

Fields:

| Field | Description |
|---|---|
| `12/29` | Completed layers / total layers to march |
| `jmin=3.41e-04` | Minimum cell Jacobian in the last completed layer |
| `Δs=4.05e-03` | Physical thickness of the last completed layer |

A `jmin` value dropping toward zero is a warning sign — increase `--smooth` or reduce `--ds` or `--growth`.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Input file not found, invalid `.npz` (missing `k0` key), or negative Jacobian detected |

## Examples

Basic march with defaults:

```bash
shore march sphere_k0.npz --nk 30 -o sphere.geo
```

Fine boundary-layer grid:

```bash
shore march fuselage_k0.npz --nk 50 --ds 1e-5 --growth 1.10 --smooth 0.3 -v
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

1. Reduce `--ds` (smaller first layer → less aggressive initial step)
2. Reduce `--growth` (slower layer expansion)
3. Increase `--smooth` (more Laplacian damping)
4. Increase `--smooth-iters` (more sweeps per layer)

## See also

- [`shore project`](project) — produce the `.npz` input
- [`shore check`](check) — validate the output
- [Hyperbolic marching](/concepts/hyperbolic-marching) — algorithm and failure modes
