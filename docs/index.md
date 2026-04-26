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
    title: STL → 6-block structured volume in one command
    details: Feed any star-shaped triangulated STL — SHORE builds the 6-block caps topology with a single `shore mesh` call. Gnomonic flat-square cap k=0 seed, great-circle equator, C1 spacing match at every cap-equator seam vertex, six Xall_GPU-compatible .geo files out.
    link: /guide/quickstart
    linkText: Quick start
  - icon: 🌊
    title: Seam-aware hyperbolic marching
    details: Each wall-normal layer advances along outward normals computed with centred differences across multi-block seams (ghost rows from partner blocks). Spacing continuity across block edges holds layer by layer; per-layer Laplacian smoothing + Jacobian guard catch any cell inversion immediately.
    link: /algorithm/hyperbolic-marching
    linkText: Algorithm
  - icon: 🧩
    title: Topology catalogue
    details: The production caps topology coexists with the conceptual O-grid primitive and the legacy Yin-Yang option. Each topology is one page; new topologies (C-grid, H-grid) plug into the same algorithm pipeline.
    link: /topologies/
    linkText: Topologies
  - icon: 🔍
    title: Focused CLI
    details: Six commands — mesh, project, info, check, export, view — let you generate, inspect, validate, and visualise. A Rich progress bar shows per-layer Jacobians during marching.
    link: /reference/mesh
    linkText: CLI reference
  - icon: 🏗️
    title: Chimera-ready output
    details: Writes a 6-block structured volume per body in the Xall_GPU .geo format. The conforming SHARED + DIRICHLET seams between sub-blocks and caps slot directly into a downstream Chimera (overset) workflow.
    link: /concepts/chimera
    linkText: Chimera meshing
  - icon: 🆓
    title: Free and open source
    details: Released under the MIT license. Pure Python — depends only on NumPy, trimesh, Typer, and Rich. No compiled extensions, no C++ meshing backends.
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

`shore mesh` builds a **6-block caps topology**: 4 equatorial sub-blocks
plus 2 polar cap blocks, all conforming through SHARED-memory seams
(sub-sub) and DIRICHLET pinning (cap-equator). Cap k=0 seeds use
**gnomonic flat-square projection**: a flat diamond in the tangent plane
at each cap pole, radially projected onto the body. Cap cells are square
by construction, with no concentric pole pattern.

```bash
shore mesh sphere.stl --ni 40 --nj 60 --nk 30
# → sphere_sub0.geo .. sphere_sub3.geo, sphere_cap_north.geo, sphere_cap_south.geo
```

Spacing continuity is enforced at two levels:

- **k=0 seam C1**: the equator's first/last meridional cell is per-j
  pinned to the cap-edge first cell at every seam vertex (so the 3
  edges meeting at each cap-corner share a single first-step length).
- **Volume seam continuity**: each block computes its outward normal
  using centred differences across multi-block seams (ghost rows from
  partner blocks), so the spacing stays continuous on the outer surface
  layer-by-layer.

See [Caps topology](topologies/caps.md) for the full algorithm.

## Author

**Stefano Zaghi** — [GitHub](https://github.com/szaghi)
