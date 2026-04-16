---
layout: home

hero:
  name: SHORE
  text: Structured Hyperbolic Overlapping Remeshing Engine
  tagline: Generate structured body-fitted O-grid volumes from raw triangulated surfaces — one command, zero manual topology.
  actions:
    - theme: brand
      text: Quick Start
      link: /guide/quickstart
    - theme: alt
      text: Guide
      link: /guide/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/shore

features:
  - icon: 📐
    title: STL → structured volume in one command
    details: Feed any raw triangulated STL — SHORE projects it onto a structured sphere template, marches layers outward along smooth vertex normals, and writes a Xall_GPU-compatible geo file with a single `shore mesh` call.
    link: /guide/quickstart
    linkText: Quick start
  - icon: 🌊
    title: Hyperbolic marching kernel
    details: Each wall-normal layer is advanced by a geometric step along area-weighted vertex normals, with per-layer Laplacian smoothing and a Jacobian guard that aborts immediately on cell inversion — no silent bad meshes.
    link: /concepts/hyperbolic-marching
    linkText: Algorithm
  - icon: 🔍
    title: Pipeline-aware CLI
    details: Five focused commands — mesh, project, march, info, check — let you inspect and re-run any individual stage. A Rich progress bar shows layer count, minimum Jacobian, and step thickness in real time.
    link: /reference/mesh
    linkText: CLI reference
  - icon: 🏗️
    title: Chimera-ready output
    details: Writes single-block structured grids in the Xall_GPU geo format (ni, nj, nk header + Fortran-ordered x, y, z arrays). Multi-block output and block decomposition are on the roadmap.
    link: /concepts/chimera
    linkText: Chimera meshing
  - icon: ✅
    title: Fully tested pipeline
    details: 35 pytest tests cover every CLI command, geometric correctness, Jacobian inversion detection, and geo round-trip fidelity. The synthetic sphere fixture runs the full pipeline in under 300 ms.
    link: /guide/contributing
    linkText: Contributing
  - icon: 🆓
    title: Free and open source
    details: Released under the MIT license. Pure Python — depends only on NumPy, trimesh, SciPy, Typer, and Rich. No compiled extensions, no C++ meshing backends.
---

## What is SHORE?

SHORE is the foundational meshing layer for a Chimera (overset) CFD mesh generator. Given a body surface — a raw triangulated STL file for now, parametric and STEP inputs planned — it produces a structured body-fitted volume grid by hyperbolic extrusion.

```bash
# Full pipeline — one command
shore mesh fuselage.stl --ni 60 --nj 80 --nk 40 --ds 1e-4 --growth 1.15 -o fuselage.geo -v
```

The result is an O-grid: a single structured block that wraps tightly around the body, with high near-wall resolution and smoothly growing cells toward the far field — exactly what a Chimera solver needs for its near-body block.

## Pipeline at a glance

```
STL (raw triangles)
      │
      ▼  shore project
  k=0 layer (ni × nj × 3)   ← structured surface on the body
      │
      ▼  shore march
  volume (ni × nj × nk × 3) ← O-grid block
      │
      ▼  shore mesh (= project + march)
  fuselage.geo               ← Xall_GPU structured grid
```


## Topology

SHORE supports two O-grid topologies for closing polar singularities and handling body-surface topology:

### Single O-grid

The default topology projects a single lat-lon sphere grid onto the body surface. The poles are excluded to avoid singular convergence points:

```bash
# Single O-grid (default)
shore mesh sphere.stl --ni 40 --nj 60 --nk 30
```

The single O-grid leaves two polar caps uncovered. For most aerodynamic bodies (wing tips, nacelle ends) this is acceptable since poles are handled separately. For closed bodies (spheres, fuselages), use the double O-grid.

### Double O-grid (Yin-Yang)

The **Yin-Yang double O-grid** resolves polar singularities without placing any grid node at a singular pole. Two standard O-grid blocks are rotated 90° relative to each other; each block's excluded caps fall inside the other block's coverage zone.

```bash
# Double O-grid (Yin-Yang) with default theta_m = pi/6
shore mesh sphere.stl --ni 40 --nj 60 --nk 30 --topology double

# Double O-grid with custom theta_m
shore mesh sphere.stl --ni 40 --nj 60 --nk 30 --topology double --theta-margin 0.6

# Validate both blocks
shore check sphere_yin.geo
shore check sphere_yang.geo

# Export to VTK for ParaView
shore export sphere_yin.geo
shore export sphere_yang.geo
```

The overlap band (default ~73% at theta_m = pi/6) is the natural chimera interpolation zone for SHORE's overset context.

#### Mathematical foundation

- **Yin piece**: standard lat-lon directions with colatitude theta in (theta_m, pi - theta_m)
- **Yang piece**: same grid rotated by R_y(+pi/2); its pseudo-poles land on the equatorial x-axis
- **Full coverage** is guaranteed if and only if theta_m < pi/4

Reference: Kageyama & Sato (2004), *Geochem. Geophys. Geosyst.* 5(9).

### Choosing theta_m

| theta_m (rad) | Approx. degrees | Overlap fraction | Notes |
|---|---|---|---|
| pi/6 | 30° | ~73% | Default. Wide overlap, comfortable margin. |
| pi/5 | 36° | ~62% | Good balance. |
| 0.65 | ~37° | ~59% | Minimum tested in examples. |
| pi/4 - epsilon | →45° | →41% | Minimum for full coverage; narrow margin. |

### Limitations

The two blocks are currently independent grids (no automatic stitching of the overlap). The overlap zone must be handled by the downstream chimera solver.

## Author

**Stefano Zaghi** — [GitHub](https://github.com/szaghi)
