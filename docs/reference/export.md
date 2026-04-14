# shore export

Export a `.geo` volume or `.npz` surface layer to VTK StructuredGrid format (`.vts`) for ParaView.

## Synopsis

```
shore export [OPTIONS] INPUT
```

## Arguments

| Argument | Description |
|---|---|
| `INPUT` | Path to a `.geo` volume file or a `.npz` surface layer file |

## Options

| Option | Default | Description |
|---|---|---|
| `-o` / `--output PATH` | `<stem>.vts` | Output `.vts` file path |
| `--no-jacobian` | off | Omit Jacobian cell data (only applies to `.geo` input) |
| `--help` | — | Show this help and exit |

## Output

A VTK StructuredGrid file (`.vts`) readable by ParaView, VisIt, and PyVista. The `.vts` extension is appended automatically if the output path has no extension.

**From a `.geo` file:**

- Grid dimensions: `(ni, nj+1, nk)` — the j-direction is closed (periodic wrap repeated).
- Point data: x, y, z coordinates.
- Cell data: `Jacobian` — scalar triple product $J = (\mathbf{e}_2 \times \mathbf{e}_1) \cdot \mathbf{e}_3$ for each cell (unless `--no-jacobian`).

**From an `.npz` file:**

- Grid dimensions: `(ni, nj+1, 1)` — a single surface layer, j-closed.
- No cell data (a single layer has no volume cells).

## Exit codes

| Code | Meaning |
|---|---|
| 0 | File written successfully |
| 1 | Input not found, unrecognised extension, or missing `k0` key in `.npz` |

## Examples

Export a volume grid and open in ParaView:

```bash
shore export sphere.geo
# → sphere.vts
paraview sphere.vts
```

Export with a custom output path:

```bash
shore export sphere.geo -o results/sphere_grid.vts
```

Export a projected surface layer:

```bash
shore export sphere_k0.npz -o sphere_surface.vts
paraview sphere_surface.vts
```

Export without Jacobian cell data (smaller file):

```bash
shore export sphere.geo --no-jacobian
```

## Viewing in ParaView

After loading a `.vts` file in ParaView:

1. Apply **Surface with Edges** representation to inspect cell structure.
2. Color by **Jacobian** to identify low-quality cells (blue = low, red = high with the default `RdYlGn` palette).
3. Use **Clip** or **Slice** to inspect interior layers.

## See also

- [`shore check`](check) — validate Jacobians without exporting
- [`shore view`](view) — interactive PyVista viewer (no ParaView required)
- [`shore.io.vtk`](api-vtk) — Python API for VTK export
