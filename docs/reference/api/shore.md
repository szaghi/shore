---
title: shore
outline: [2, 3]
---

# `shore`

SHORE — Structured Hyperbolic Overlapping Remeshing Engine.

**Modules**

- [`shore.balance`](./shore.balance) — Multi-way equal-partition split planner for Xall MPI load balancing
- [`shore.boxes`](./shore.boxes) — Generate :class:`shore.io.boxes.Box` instances from an STL surface
- [`shore.collar`](./shore.collar) — The collar block of two intersecting bodies
- [`shore.config`](./shore.config) — TOML-driven configuration for the `shore mesh` pipeline
- [`shore.exceptions`](./shore.exceptions) — SHORE-specific exception hierarchy
- [`shore.geometry`](./shore.geometry) — Geometry backends behind one protocol
- [`shore.gui`](./shore.gui) — SHORE's desktop GUI
- [`shore.io`](./shore.io) — I/O
- [`shore.mesh`](./shore.mesh) — OOP mesh layer for SHORE
- [`shore.pipeline`](./shore.pipeline) — Meshing pipeline steps shared by the CLI and the command layer
- [`shore.primitive`](./shore.primitive) — Parametric structured-grid primitives for Chimera/overset assembly
- [`shore.profiles`](./shore.profiles) — 2D airfoil profile generators (NACA 4-digit, 5-digit)
- [`shore.project`](./shore.project) — Command layer
- [`shore.quality`](./shore.quality) — Structured-grid quality metrics (roadmap phase 4)
- [`shore.shapes`](./shore.shapes) — 3D body STL builders backed by trimesh
- [`shore.split`](./shore.split) — Split structured blocks along ijk vertex indices, preserving adjacency
- [`shore.surface`](./shore.surface) — Surface I/O and unstructured-to-structured remeshing
- [`shore.topology`](./shore.topology) — Topology decomposition utilities for multi-block structured grids
- [`shore.vis`](./shore.vis) — Visualisation helpers for SHORE grids
- [`shore.volume`](./shore.volume) — Hyperbolic normal-extrusion

**Re-exported**

- [`MarchingError`](./shore.exceptions#MarchingError)
- [`MeshInputError`](./shore.exceptions#MeshInputError)
- [`ParameterError`](./shore.exceptions#ParameterError)
- [`ShoreError`](./shore.exceptions#ShoreError)
