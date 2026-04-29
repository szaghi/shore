<div align="center">

# SHORE — Structured Hyperbolic Overlapping Remeshing Engine

[![Latest Version](https://img.shields.io/pypi/v/shore.svg)](https://pypi.org/project/shore/)
[![Build Status](https://github.com/szaghi/shore/actions/workflows/python-package.yml/badge.svg)](https://github.com/szaghi/shore/actions)
[![Coverage](https://img.shields.io/endpoint?url=https://szaghi.github.io/shore/coverage-badge.json)](https://szaghi.github.io/shore/)
[![Supported Python versions](https://img.shields.io/badge/Py-%203.11,%203.12-blue.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)]()

> Hyperbolic mesh generation for Chimera CFD — body-fitted near-body grids and parametric backgrounds, in one tool.

<div>
<table>
<tr>
<td><b>🏗️ Full Chimera assembly in one tool</b><br><sub>Generate every grid an overset assembly needs — body-fitted near-body meshes and parametric backgrounds — in the same world frame. Pack them into the Xall solver bundle (.grd + cc.par + marker boxes) without leaving the SHORE CLI. <a href="https://szaghi.github.io/shore/guide/cc-par-pipeline">Pipeline walkthrough</a></sub></td>
<td><b>📐 Body-fitted near-body meshes</b><br><sub>Two topologies for surface inputs (STL today; CAD on the roadmap) — single-block O-grid and 6-block cubed sphere with gnomonic cap k=0 seed and C1 cap-equator seam pinning. One <code>shore mesh</code> call from a triangulated body to Xall-compatible .geo files. <a href="https://szaghi.github.io/shore/topologies/cubed-sphere">Cubed sphere</a></sub></td>
</tr>
<tr>
<td><b>🧱 Parametric backgrounds</b><br><sub>Box, annular cylinder, and 5-block flat-caps cylinder primitives with per-axis spacing laws — fill the far field without an STL. Per-edge spacing on Box and Annulus, transfinite-interpolation blends; flat-caps cylinder is Chimera-fringe-free for GPU-friendly backgrounds. <a href="https://szaghi.github.io/shore/topologies/primitives">Primitives</a></sub></td>
<td><b>🌊 Seam-aware hyperbolic marching</b><br><sub>Each wall-normal layer advances along outward normals computed with centred differences across multi-block seams (ghost rows from partner blocks). Spacing continuity across block edges holds layer by layer; per-layer Laplacian smoothing + Jacobian guard catch any cell inversion immediately. <a href="https://szaghi.github.io/shore/algorithm/hyperbolic-marching">Algorithm</a></sub></td>
</tr>
<tr>
<td><b>🔍 Focused CLI</b><br><sub>Nineteen commands cover the full chimera pipeline — generate (<code>mesh</code>, <code>primitive</code>), subdivide and balance for parallel solvers (<code>split</code>, <code>balance</code>), pack for the solver (<code>grd</code>, <code>cc-par</code>, <code>proc-input</code>, <code>boxes-from-stl</code>), merge multi-component runs (<code>grd-merge</code>, <code>cc-par-merge</code>, <code>proc-input-merge</code>), inspect (<code>info</code>, <code>check</code>), and visualise (<code>export</code>, <code>view</code>, <code>boxes-export</code>). A Rich progress bar shows per-layer Jacobians during marching. <a href="https://szaghi.github.io/shore/reference/mesh">CLI reference</a></sub></td>
<td><b>🎨 ParaView and Jupyter visualization</b><br><sub>Export to VTK StructuredGrid (<code>shore export</code>) and MultiBlock (<code>.vtm</code>) for ParaView with per-cell Jacobian colour maps and overlaid marker boxes. Or use <code>shore view</code> for an immediate PyVista window, or <code>shore.vis</code> directly in Jupyter notebooks with the HTML backend. <a href="https://szaghi.github.io/shore/reference/api-vis">API reference</a></sub></td>
</tr>
<tr>
<td><b>✅ Fully tested pipeline</b><br><sub>Pytest suite covers every CLI command, geometric correctness, Jacobian inversion detection, multi-block adjacency wiring, axis-map orientation, VTK export round-trip, and geo / grd / cc.par format fidelity. <a href="https://szaghi.github.io/shore/guide/contributing">Contributing</a></sub></td>
<td><b>🆓 Free and open source</b><br><sub>Released under the MIT license. Pure Python — depends only on NumPy, trimesh, Typer, and Rich. PyVista is an optional dependency for visualization. No compiled extensions, no C++ meshing backends.</sub></td>
</tr>
</table>
</div>

**[Full documentation](https://szaghi.github.io/shore/)**

</div>

## What is SHORE?

SHORE produces every structured grid a Chimera (overset) CFD assembly needs and packages them for the solver. Two product lines, one tool:

- **Body-fitted near-body meshes** by hyperbolic normal extrusion of a body surface. Surface input is a triangulated STL today; CAD (STEP / IGES / BREP) is on the roadmap.
- **Parametric backgrounds** that fill the far field — no STL needed, just analytic geometry. Box, annular cylinder, and flat-caps cylinder primitives.

Hole cutting, fringe identification, and donor search are the **solver's** responsibility — SHORE produces every grid and every descriptor the solver needs, then steps out of the way.

```bash
# Body-fitted near-body mesh (default: 6-block cubed sphere)
shore mesh fuselage.stl --ni 60 --nj 80 --nk 40 --ds 1e-4 --growth 1.15 -o fuselage -v

# Single-block O-grid alternative (one .geo, j-periodic)
shore mesh fuselage.stl --topology ogrid --ni 60 --nj 80 --nk 40 -o fuselage

# Parametric background — no STL needed
shore primitive flat-caps --r-out 10 --z0 -10 --z1 50 --ni 16 --nc 24 --nk 60 -o background

# Inspect
shore info fuselage_sub0.geo
shore check fuselage_sub0.geo
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
make dev    # pip install -e ".[dev]"
make test   # currently ~588 passing
```

### Your first body-fitted mesh

Given a body surface `sphere.stl`, generate a 6-block cubed-sphere mesh with 40 equator latitude points, 60 longitude points, and 30 wall-normal layers:

```bash
shore mesh sphere.stl --ni 40 --nj 60 --nk 30 --ds 1e-3 --growth 1.15 -o sphere -v
```

The `-v` flag shows a live progress bar:

```
⠸ Loading STL     ████████████████  1/1   5,120 triangles  0:00:00
⠸ Marching        ████████████████ 29/29  jmin=2.30e-05    0:00:00
✓ Wrote 6 .geo files  (sub: 40x16x30, caps: 16x16x30)
```

Validate the result:

```bash
shore check sphere_sub0.geo
```

```
✓ All Jacobians positive.
```

### Your first parametric background

No STL required — primitives build their nodes from analytic geometry:

```bash
# 5-block flat-caps cylinder (Chimera-fringe-free; central square + 4 trapezoidal sub-blocks)
shore primitive flat-caps --r-out 10 --z0 -10 --z1 50 --ni 16 --nc 24 --nk 60 -o background
# → background_center.geo, background_sub_e.geo, background_sub_n.geo,
#   background_sub_w.geo, background_sub_s.geo
```

Other primitives:

```bash
shore primitive box --min -10 -10 -10 --max 10 10 50 --ni 30 --nj 30 --nk 60 -o domain
shore primitive annulus --r-in 1.5 --r-out 5 --z0 -2 --z1 2 --ni 30 --nj 64 --nk 40 -o buffer
```

## Xall pipeline

Once each component (near-body and / or background) is generated, package it for the Xall solver:

```bash
# Pack each component into binary .grd
shore grd wall_*.geo -o wall.grd
shore grd background_*.geo -o background.grd

# Generate the chimera descriptor (consumed by overset preprocessor)
shore boxes-from-stl fuselage.stl --block sub0 --max-outward-gap 0.10 -o boxes.json
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd --boxes boxes.json
shore cc-par background.adjacency.json background.cc.par --grd background.grd --free-bc -9

# Generate the runtime work-distribution descriptor (consumed by Xall)
shore proc-input wall.grd       wall.proc.input       --np 16
shore proc-input background.grd background.proc.input --np 16

# Merge multi-component runs into a single Xall job
shore grd-merge        background.grd          wall.grd          -o assembly.grd
shore cc-par-merge     background.cc.par       wall.cc.par       -o assembly.cc.par --grd assembly.grd
shore proc-input-merge background.proc.input   wall.proc.input   -o assembly.proc.input
```

The full walkthrough — including the `examples/05-chimera-assembly/cc_par_assembly.py` reference script — is at the [cc.par + boxes pipeline guide](https://szaghi.github.io/shore/guide/cc-par-pipeline).

## CLI reference

Nineteen commands, grouped by purpose:

| Group | Commands |
|---|---|
| **Mesh generation** | `shore mesh`, `shore primitive`, `shore project` |
| **Adjacency / parallelism** | `shore split`, `shore balance` |
| **Xall pipeline** | `shore grd`, `shore grd-merge`, `shore cc-par`, `shore cc-par-merge`, `shore proc-input`, `shore proc-input-merge`, `shore boxes-from-stl` |
| **Inspection** | `shore info`, `shore check` |
| **Visualisation** | `shore export`, `shore view`, `shore boxes-export` |
| **Configuration** | `shore config` |

All commands accept `--help` for full option lists. See the [CLI reference](https://szaghi.github.io/shore/reference/mesh) for full documentation.

### Key options for `shore mesh`

| Option | Default | Description |
|---|---|---|
| `--topology` | `cubed_sphere` | `cubed_sphere` (6 blocks) or `ogrid` (1 block, j-periodic) |
| `--ni` | 40 | Latitude / equator meridional **node count** |
| `--nj` | 60 | Longitude **node count** (multiple of 4 for cubed sphere) |
| `--nk` | 30 | Wall-normal **node count** (includes the wall row at k=0; the marcher extrudes `nk - 1` cells) |
| `--ds` | 1e-3 | First-layer thickness (geometric spacing) |
| `--growth` | 1.15 | Geometric growth ratio |
| `--theta-cap-deg` | 30° | Polar exclusion angle (cubed sphere) |
| `--volume-k-spacing` | `geometric` | `geometric`, `tanh`, or `tanh2` |
| `--smooth` | 0.2 | Laplacian smooth strength [0–1] |
| `--smooth-iters` | 2 | Laplacian sweeps per layer |
| `-o` / `--output` | `<stem>` | Output stem (multi-block writes `<stem>_<label>.geo`) |
| `-v` / `--verbose` | off | Rich progress bar |

## Algorithm overview

### Body-fitted: surface projection + hyperbolic marching

For body-fitted meshes, SHORE runs two stages:

1. **Surface projection.** The body's anchor (centroid by default) is the origin of a structured ray-cast onto the STL. The 6-block cubed sphere uses gnomonic flat-square cap projection at each pole + great-circle slerp meridians for the equator, with C1 step matching at every cap-equator seam vertex. The single-block O-grid uses a uniform lat-lon ray template. Both require the body to be **star-shaped** from the anchor.
2. **Hyperbolic marching.** Each wall-normal layer is computed as

   $$P_{i,j}^{k+1} = P_{i,j}^{k} + \Delta s_k \cdot \hat{n}_{i,j}^{k}$$

   where $\Delta s_k$ comes from the chosen spacing law (geometric, tanh, or tanh2) and $\hat{n}$ is the outward normal. At multi-block seams, the boundary tangent uses a centred difference against ghost rows from partner blocks instead of a one-sided fallback — this is the **seam-aware kernel** that keeps spacing continuous across block edges layer by layer. Per-layer Laplacian smoothing and a Jacobian guard catch any cell inversion immediately.

### Parametric backgrounds

Primitives build their nodes from analytic geometry — no projection, no marching. Each axis has its own spacing law (uniform / tanh / tanh2); per-edge overrides on Box and Annulus produce a transfinite-interpolation blend between the four parallel edges. The 5-block flat-caps cylinder pins its radial and tangential distributions geometrically (central square + 4 trapezoidal sub-blocks with conforming SHARED seams).

### Marker boxes for hole-cutting

For chimera assemblies, SHORE generates marker boxes from the body STL via voxel-fill + greedy AABB merge: voxelise the body's AABB, include every voxel that overlaps the body, then greedy-merge face-adjacent grid blocks into slab-shaped AABBs whenever the merged corners stay within the user's `max_outward_gap` budget. The output is a small slab cover that contains the body and fits inside the near-mesh — both invariants enforced by construction. See [Algorithm — marker boxes](https://szaghi.github.io/shore/algorithm/boxes).

## Roadmap

- [x] Body-fitted topologies: 6-block cubed sphere (production default) + single-block O-grid
- [x] Parametric primitives: box, annular cylinder, 5-block flat-caps cylinder
- [x] Gnomonic flat-square cap k=0 projection with C1 seam-step pinning
- [x] Seam-aware hyperbolic marching (centred-difference normals across block seams)
- [x] Per-edge spacing data model (`Spacing1D`, `Edge.first_step`/`last_step`, `pinned` law)
- [x] Adjacency JSON sidecar with `axis_map` orientation
- [x] Xall pipeline: `.geo` → `.grd` (`shore grd`), chimera descriptor (`shore cc-par`), marker boxes (`shore boxes-from-stl`)
- [x] Multi-component merge (`shore grd-merge`, `shore cc-par-merge`) for Xall assemblies
- [x] PyVista visualization (`shore view`, `shore export`, `shore boxes-export`, `shore.vis`)
- [x] GitHub Actions CI + VitePress documentation site
- [ ] CAD surface input (STEP / IGES / BREP via pythonOCC or cadquery)
- [ ] Multi-block surface decomposition for non-star-shaped bodies
- [ ] Steger–Sorenson orthogonality source term for high-curvature bodies
- [ ] Box-shaped outer-boundary morphing (alternative to background overlap)
- [ ] Other CFD formats (Plot3D, CGNS)

## Author

**Stefano Zaghi** — [GitHub](https://github.com/szaghi)

## License

SHORE is released under the [MIT License](LICENSE).
