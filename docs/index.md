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

## Author

**Stefano Zaghi** — [GitHub](https://github.com/szaghi)
