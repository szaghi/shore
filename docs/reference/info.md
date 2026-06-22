# shore info

Print metadata of a structured `.geo` block or an `.npz` surface file: dimensions, point and cell counts, and axis-aligned bounding box.

## Synopsis

```
shore info [OPTIONS] INPUT
```

## Arguments

| Argument | Description |
|---|---|
| `INPUT` | Path to a `.geo` block file or an `.npz` surface file |

## Options

| Option | Description |
|---|---|
| `--help` | Show this help and exit |

## Output

```
┌─────────────────────────────────────────────┐
│ sphere.geo                                  │
├─────────────┬──────────────────────────────┤
│ Dimensions  │ ni=40  nj=60  nk=30          │
│ Points      │ 72,000                        │
│ Cells       │ 47,151                        │
│ x bounds    │ [-1.042000,  1.042000]        │
│ y bounds    │ [-1.042000,  1.042000]        │
│ z bounds    │ [-1.042000,  1.042000]        │
└─────────────┴──────────────────────────────┘
```

**Cells** = `(ni-1) × (nj-1) × (nk-1)` — the number of hexahedral cells in the grid.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | File not found |

## Examples

```bash
shore info sphere.geo
shore info grids/fuselage.geo
```

## See also

- [`shore check`](check) — validate Jacobians
- [geo format](/reference/api-geo) — file format details
