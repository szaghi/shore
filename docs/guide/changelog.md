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
