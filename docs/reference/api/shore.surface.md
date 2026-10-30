---
title: shore.surface
outline: [2, 3]
---

# `shore.surface`

Surface I/O and unstructured-to-structured remeshing.

**Modules**

- [`shore.surface.anchor`](./shore.surface.anchor) — Anchor-point resolution strategies for ray-cast projections
- [`shore.surface.cap`](./shore.surface.cap) — Polar cap k=0 seed face via gnomonic projection of a flat tangent square
- [`shore.surface.cgrid`](./shore.surface.cgrid) — C-grid surface construction
- [`shore.surface.equator`](./shore.surface.equator) — Equator band construction from cap seam boundaries
- [`shore.surface.intersect`](./shore.surface.intersect) — Intersection loop of two bodies, the start curve of a collar grid
- [`shore.surface.io`](./shore.surface.io) — STL reader
- [`shore.surface.march`](./shore.surface.march) — Hyperbolic surface marching
- [`shore.surface.projector`](./shore.surface.projector) — Ray-mesh projector with cached BVH and per-ray status tracking
- [`shore.surface.remesh`](./shore.surface.remesh) — Unstructured STL surface → structured (ni, nj) quad grid via sphere projection
- [`shore.surface.tfi`](./shore.surface.tfi) — Surface grids on one parametric face
