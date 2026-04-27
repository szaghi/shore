# shore export

Export `.geo` volumes or a `.npz` surface layer to VTK XML for ParaView.

- **One** `.geo` or `.npz` input → one `.vts` StructuredGrid.
- **Two or more** `.geo` inputs → one `.vtm` MultiBlock, with one
  PyVista block per input. Block labels are derived from each input's
  filename stem.

## Synopsis

```
shore export [OPTIONS] INPUTS...
```

## Arguments

| Argument | Description |
|---|---|
| `INPUTS...` | One or more input files. A single `.geo` / `.npz` produces a `.vts`; multiple `.geo` inputs produce a `.vtm` MultiBlock. |

## Options

| Option | Default | Description |
|---|---|---|
| `-o` / `--output PATH` | `<first-input-stem>.vts` (single) or `<first-input-stem>.vtm` (multi) | Output file path. The right extension (`.vts` / `.vtm`) is forced if missing. |
| `--no-jacobian` | off | Omit Jacobian cell data (applies to `.geo` inputs in both single and multi modes) |
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

## Multi-block export

Pass multiple `.geo` files in one invocation; SHORE bundles them into
a single ParaView-friendly `.vtm` MultiBlock. Block keys inside the
`.vtm` come from each file's stem.

```bash
# 6-block caps topology — sphere_<label>.geo  →  sphere.vtm
shore export sphere_*.geo -o sphere.vtm
paraview sphere.vtm

# 5-block flat-caps cylinder primitive
shore export bg_*.geo -o background.vtm
```

The shell glob is the natural way to drive this; the CLI itself
takes a list of paths, so `shore export a.geo b.geo c.geo -o out.vtm`
works equally well.

If you pass a `.npz` together with `.geo` files, the command rejects
the mix — multi-input mode is `.geo`-only by design (a single `.vtm`
cannot represent both volume blocks and a surface layer
unambiguously). Convert the surface separately if you need both.

## See also

- [`shore check`](check) — validate Jacobians without exporting
- [`shore view`](view) — interactive PyVista viewer (no ParaView required)
- [`shore.io.vtk`](api-vtk) — Python API for VTK export
- [`shore.topology.sphere_cap`](api-sphere-cap) — 6-block polar cap topology
