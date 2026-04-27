---
title: Changelog
---

# Changelog

## [Unreleased]

### Added
- Full pipeline: `shore mesh`, `shore project`, `shore march`, `shore info`, `shore check`
- Hyperbolic marching kernel with geometric growth, Laplacian smoothing, and per-layer Jacobian check
- Sphere-projection surface remeshing via BVH-accelerated ray casting (trimesh + rtree)
- Xall_GPU geo format reader/writer (ASCII, Fortran column-major ordering)
- Rich progress bar with live jmin and step display
- 35 pytest tests covering all commands and geometric correctness
- VitePress documentation
- Per-edge spacing data model: `Edge.first_step` / `Edge.last_step` fields, `Spacing1D(law='pinned')`, public `pinned_distribution` (lifted from `surface/equator.py`).
- Per-block k-spacing override through `Mesh.march()` — replace the four k-edges' `Spacing1D` instance to override `(ds, growth)`. Identity-check resolver enforces a single schedule across all blocks (k-coupled SHARED / DIRICHLET seams). Works on cubed-sphere and O-grid topologies.
- Per-edge i/j/k on `BoxMesh` and i/j on `AnnularCylinderMesh` via `edge_pins` kwarg; constructor builds a transfinite-interpolation (TFI) blend of the four parallel edges' distributions when they disagree, collapses to legacy tensor-product when they don't.
- Examples: `per_edge_ij_box.py`, `per_edge_k_spacing_cubed_sphere.py` (renamed from `per_edge_k_spacing.py`), `per_edge_k_spacing_ogrid.py`.

### Documented
- Per-edge spacing section in `algorithm/spacing-laws.md` with the full support table by topology / primitive.
- `topologies/cubed-sphere.md`, `topologies/o-grid.md`, `topologies/primitives.md`, `reference/primitive.md` updated with per-edge guidance and the structural reasons each axis is or isn't wired.
