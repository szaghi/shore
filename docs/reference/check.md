# shore check

Validate the cell Jacobians of a structured `.geo` file layer by layer.

Reports the minimum Jacobian per layer in a table and exits with code 1 if any layer contains an inverted cell.

## Synopsis

```
shore check [OPTIONS] INPUT
```

## Arguments

| Argument | Description |
|---|---|
| `INPUT` | Path to a `.geo` file |

## Options

| Option | Default | Description |
|---|---|---|
| `--warn-threshold FLOAT` | 1e-6 | Layers with jmin below this value are flagged as warnings even if positive |
| `--help` | — | Show this help and exit |

## Output

### Clean grid

```
┌────────────────────────────────────────────────┐
│ Jacobian check — sphere.geo                    │
├───────┬────────────┬────────┤
│ Layer │ jmin       │ Status │
├───────┼────────────┼────────┤
│     1 │ 3.4112e-04 │ ok     │
│     2 │ 3.6889e-04 │ ok     │
│   ... │        ... │ ok     │
│    29 │ 8.2341e-03 │ ok     │
└───────┴────────────┴────────┘
✓ All Jacobians positive.
```

### Grid with near-inversion warning

```
│    14 │ 3.2100e-07 │ warn   │
```

### Grid with inverted cells

```
│     7 │ -1.2300e-04 │ INVERTED │
✗ Inverted cells detected.
```

## Exit codes

| Code | Meaning |
|---|---|
| 0 | All Jacobians positive (warnings are non-fatal) |
| 1 | File not found or at least one layer has jmin ≤ 0 |

The exit code makes `shore check` composable in shell pipelines:

```bash
shore mesh body.stl && shore check body.geo && xall_gpu body.geo
```

## Jacobian definition

The Jacobian for each cell is the scalar triple product:

$$J_{i,j} = (\mathbf{e}_2 \times \mathbf{e}_1) \cdot \mathbf{e}_3$$

where $\mathbf{e}_1$, $\mathbf{e}_2$, $\mathbf{e}_3$ are the three edge vectors of the cell (j-tangent, i-tangent, wall-normal). The reported `jmin` per layer is the minimum over all $(i, j)$ cells in that layer interval.

See [Hyperbolic marching — Jacobian check](/algorithm/hyperbolic-marching#jacobian-check) for the full derivation.

## Examples

Validate after generation:

```bash
shore mesh sphere.stl -o sphere.geo
shore check sphere.geo
```

Use a stricter warning threshold:

```bash
shore check sphere.geo --warn-threshold 1e-4
```

In a CI/CD script:

```bash
#!/bin/bash
set -e
shore mesh "$INPUT" -o "$OUTPUT"
shore check "$OUTPUT"
echo "Grid $OUTPUT validated."
```

## See also

- [`shore mesh`](mesh) — generate the grid
- [Hyperbolic marching](/algorithm/hyperbolic-marching) — why Jacobians matter
