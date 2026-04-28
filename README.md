<div align="center">

# SHORE — Structured Hyperbolic Overlapping Remeshing Engine

[![Latest Version](https://img.shields.io/pypi/v/shore.svg)](https://pypi.org/project/shore/)
[![Build Status](https://github.com/szaghi/shore/actions/workflows/python-package.yml/badge.svg)](https://github.com/szaghi/shore/actions)
[![Coverage](https://img.shields.io/endpoint?url=https://szaghi.github.io/shore/coverage-badge.json)](https://szaghi.github.io/shore/)
[![Supported Python versions](https://img.shields.io/badge/Py-%203.11,%203.12-blue.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)]()

> Generate structured body-fitted O-grid volumes from raw triangulated surfaces via hyperbolic normal extrusion.
> Building block for Chimera (overset) mesh generation in computational fluid dynamics.

<div>
<table>
<tr>
<td><b>📐 STL → structured volume in one command</b><br><sub>Feed any raw triangulated STL — SHORE projects it onto a structured sphere template, marches layers outward along smooth vertex normals, and writes a Xall_GPU-compatible geo file with a single <code>shore mesh</code> call. <a href="https://szaghi.github.io/shore/guide/quickstart">Quick start</a></sub></td>
<td><b>🌊 Hyperbolic marching kernel</b><br><sub>Each wall-normal layer is advanced by a geometric step along area-weighted vertex normals, with per-layer Laplacian smoothing and a Jacobian check that aborts immediately on cell inversion — no silent bad meshes. <a href="https://szaghi.github.io/shore/algorithm/hyperbolic-marching">Algorithm</a></sub></td>
</tr>
<tr>
<td><b>🔍 Pipeline-aware CLI</b><br><sub>Eight focused commands — <code>mesh</code>, <code>primitive</code>, <code>project</code>, <code>config</code>, <code>info</code>, <code>check</code>, <code>export</code>, <code>view</code> — let you inspect and re-run any individual stage, drive runs from TOML, and emit ParaView-ready output. A Rich progress bar shows layer count, minimum Jacobian, and step thickness in real time. <a href="https://szaghi.github.io/shore/reference/mesh">CLI reference</a></sub></td>
<td><b>🏗️ Multi-block topologies</b><br><sub>Two body-fitted topologies — single-block O-grid and the 6-block cubed sphere (4 equatorial sub-blocks + 2 hexahedral cap blocks) — plus three parametric primitives (box, annular cylinder, flat-caps cylinder) for backgrounds and buffers. All share conforming node-for-node faces and are exported as ParaView MultiBlock (.vtm). <a href="https://szaghi.github.io/shore/topologies/cubed-sphere">Cubed-sphere topology</a></sub></td>
</tr>
<tr>
<td><b>🎨 ParaView and Jupyter visualization</b><br><sub>Export to VTK StructuredGrid (<code>shore export</code>) for ParaView with per-cell Jacobian colour maps. Or use <code>shore view</code> for an immediate PyVista window, or <code>shore.vis</code> directly in Jupyter notebooks with the HTML backend. <a href="https://szaghi.github.io/shore/reference/api-vis">API reference</a></sub></td>
<td><b>✅ Fully tested pipeline</b><br><sub>Pytest suite covers every CLI command, geometric correctness (outer radius > inner, surface projection within 5 % of target), Jacobian inversion detection, VTK export round-trip, and geo format fidelity. <a href="https://szaghi.github.io/shore/guide/contributing">Contributing</a></sub></td>
</tr>
<tr>
<td><b>🆓 Free and open source</b><br><sub>Released under the MIT license. Pure Python — depends only on NumPy, trimesh, SciPy, Typer, and Rich. PyVista is an optional dependency for visualization. No compiled extensions, no C++ meshing backends.</sub></td>
<td></td>
</tr>
</table>
</div>

**[Full documentation](https://szaghi.github.io/shore/)**

</div>

## What is SHORE?

SHORE is the foundational meshing layer for a Chimera (overset) CFD mesh generator. Given a body surface — a raw triangulated STL file for now, parametric and STEP inputs planned — it produces a structured body-fitted volume grid by hyperbolic extrusion:

1. **Project** a structured sphere template onto the body surface → k=0 layer
2. **March** outward along smooth vertex normals, layer by layer, with geometric growth
3. **Write** the structured `(ni, nj, nk)` volume to a geo file ready for Xall_GPU

The result is a body-fitted structured volume that wraps tightly around the body, with high near-wall resolution and smoothly growing cells toward the far field. Two topologies are available: the 6-block cubed sphere (default — no polar singularity, conforming SHARED seams) and the single-block, j-periodic O-grid (`--topology ogrid`).

```bash
# Full pipeline — one command
shore mesh fuselage.stl --ni 60 --nj 80 --nk 40 --ds 1e-4 --growth 1.15 -o fuselage.geo

# Inspect the result
shore info fuselage.geo
shore check fuselage.geo
```

## Quick start

### Installation

```bash
pip install shore
```

Or from source for development:

```bash
git clone https://github.com/szaghi/shore.git
cd shore
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

### Your first O-grid

Given a body surface `sphere.stl`, generate a 40×60×30 O-grid with a first-layer thickness of 1 mm and growth ratio 1.15:

```bash
shore mesh sphere.stl --ni 40 --nj 60 --nk 30 --ds 1e-3 --growth 1.15 -o sphere.geo -v
```

The `-v` flag shows a live progress bar:

```
⠸ Loading STL     ████████████████████████  1/1   320 triangles        0:00:00
⠸ Projecting      ████████████████████████  1/1   (40×60) surface layer 0:00:00
⠸ Marching        ████████████░░░░░░░░░░░░  14/29  jmin=3.41e-04  Δs=4.05e-03  0:00:01
```

Validate the result:

```bash
shore check sphere.geo
```

```
┌──────────────────────────────────────────────┐
│ Jacobian check — sphere.geo                  │
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

## Pipeline stages

The `shore mesh` command builds the
[6-block cubed-sphere topology](docs/topologies/cubed-sphere.md)
end-to-end by default (gnomonic flat-square cap projection +
great-circle meridian equator + seam-aware hyperbolic marching) and
writes 6 `.geo` files. Pass `--topology ogrid` to build a
single-block, j-periodic O-grid instead — written as one `.geo`. For
debugging or downstream tools that consume a single-block lat-lon
surface, `shore project` writes the `(ni, nj, 3)` surface layer as
`.npz` without marching.

```bash
# Project the lat-lon surface only (single-block, .npz)
shore project sphere.stl --ni 40 --nj 60 -o sphere_k0.npz -v

# Inspect any geo file
shore info sphere_sub0.geo
```

The `.npz` intermediate format is a plain `numpy.savez_compressed` file — inspect it with:

```python
import numpy as np
k0 = np.load("sphere_k0.npz")["k0"]   # shape (ni, nj, 3)
```

## CLI reference

| Command | Description |
|---|---|
| `shore mesh INPUT.stl` | Full pipeline: STL → cubed-sphere (6 blocks, default) or O-grid (1 block, `--topology ogrid`) volume → `.geo` file(s) |
| `shore project INPUT.stl` | Project a single-block lat-lon surface onto the STL (`.npz`) |
| `shore info INPUT.geo` | Print dimensions, point count, bounding box |
| `shore check INPUT.geo` | Validate per-layer Jacobians, exit 1 on inversion |
| `shore export INPUT` | Export `.geo` or `.npz` to VTK `.vts` for ParaView |
| `shore view INPUT` | Open interactive PyVista viewer (requires `shore[vis]`) |

All commands accept `--help` for full option lists.

### Key options for `shore mesh`

| Option | Default | Description |
|---|---|---|
| `--ni` | 40 | Latitude resolution (pole-to-pole) |
| `--nj` | 60 | Longitude resolution (periodic) |
| `--nk` | 30 | Wall-normal layers (includes k=0) |
| `--ds` | 1e-3 | First-layer thickness (body units) |
| `--growth` | 1.15 | Geometric growth ratio |
| `--smooth` | 0.2 | Laplacian smooth strength [0–1] |
| `--smooth-iters` | 2 | Laplacian sweeps per layer |
| `-o` / `--output` | `<stem>.geo` | Output path |
| `-v` / `--verbose` | off | Rich progress bar |

## Algorithm overview

### Surface projection (stage 1)

A structured `(ni, nj)` lat-lon sphere is centred at the STL body centroid. For each grid point on the sphere, a ray is cast from the centroid outward; the first intersection with the STL triangulation is taken as the k=0 surface point. This requires the body to be **star-shaped** with respect to its centroid — a valid assumption for most CFD bodies (fuselages, wings, cylinders, nacelles).

Ray casting uses trimesh's BVH-accelerated intersector backed by `rtree`.

### Hyperbolic marching (stage 2)

Given the k=0 layer, each subsequent layer k+1 is computed as:

$$P_{i,j}^{k+1} = P_{i,j}^{k} + \Delta s \cdot r^k \cdot \hat{n}_{i,j}^{k}$$

where $\Delta s$ is the first-layer thickness, $r$ is the growth ratio, and $\hat{n}$ is the area-weighted unit normal estimated from the four in-plane neighbours (central differences, periodic in j). Laplacian smoothing is applied after each layer, followed by a Jacobian check; if any cell inverts, marching aborts with a clear error message.

## Roadmap

- [ ] Box-shaped outer boundary morphing (Xall_GPU background grid compatibility)
- [ ] Multi-block surface decomposition for non-star-shaped bodies
- [ ] Steger–Sorenson orthogonality source term for high-curvature bodies
- [ ] STEP and parametric surface input (via pythonOCC or cadquery)
- [x] PyVista visualization (`shore view`, `shore export`, `shore.vis`)
- [x] 6-block cubed-sphere topology (production default for `shore mesh`)
- [x] Gnomonic flat-square cap k=0 projection with C1 seam-step pinning
- [x] Seam-aware hyperbolic marching (centred-difference normals across block seams)
- [x] Single-block O-grid topology (`shore mesh --topology ogrid`)
- [x] Parametric primitives (box, annular cylinder, flat-caps cylinder) for backgrounds and buffers
- [ ] GitHub Actions CI

## Author

**Stefano Zaghi** — [GitHub](https://github.com/szaghi)

## License

SHORE is released under the [MIT License](LICENSE).
