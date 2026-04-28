---
title: Changelog
---

# Changelog

## [Unreleased]

### Added
- Full pipeline: `shore mesh`, `shore project`, `shore march`, `shore info`, `shore check`
- Hyperbolic marching kernel with geometric growth, Laplacian smoothing, and per-layer Jacobian check
- Sphere-projection surface remeshing via BVH-accelerated ray casting (trimesh + rtree)
- Xall geo format reader/writer (ASCII, Fortran column-major ordering)
- Rich progress bar with live jmin and step display
- 35 pytest tests covering all commands and geometric correctness
- VitePress documentation
- Per-edge spacing data model: `Edge.first_step` / `Edge.last_step` fields, `Spacing1D(law='pinned')`, public `pinned_distribution` (lifted from `surface/equator.py`).
- Per-block k-spacing override through `Mesh.march()` — replace the four k-edges' `Spacing1D` instance to override `(ds, growth)`. Identity-check resolver enforces a single schedule across all blocks (k-coupled SHARED / DIRICHLET seams). Works on cubed-sphere and O-grid topologies.
- Per-edge i/j/k on `BoxMesh` and i/j on `AnnularCylinderMesh` via `edge_pins` kwarg; constructor builds a transfinite-interpolation (TFI) blend of the four parallel edges' distributions when they disagree, collapses to legacy tensor-product when they don't.
- Examples: `per_edge_ij_box.py`, `per_edge_k_spacing_cubed_sphere.py` (renamed from `per_edge_k_spacing.py`), `per_edge_k_spacing_ogrid.py`.
- `shore cc-par` command + `shore.io.cc_par` module: translate the SHORE `.adjacency.json` sidecar into the Xall `cc.par` chimera pre-processor input format. WALL/FREE/SHARED/DIRICHLET map to Xall BC codes; SHARED/DIRICHLET orientation codes are derived from each face's `axis_map` (no per-topology lookup tables). `BlockMeta` dataclass exposes per-block level / group / priority / comment + `free_face_bc` for body-fitted (40) vs background (-9, -3, -4) meshes.

### Changed
- **Adjacency JSON bumped from v1 to v2.** Every connection face (`SHARED` / `DIRICHLET` / `PERIODIC`) now carries an `axis_map` field describing how its along-seam axes map onto the partner face's along-seam axes. v1 sidecars are rejected with a regenerate-the-sidecar message. Sidecars are derived artefacts — re-run the topology / split that produced them.
- **`cc.par` orientation codes derived from `axis_map` (no per-topology tables).** The previously hard-coded `_CAP_EQUATOR_PARTNER` table in `shore.io.cc_par` is gone. Orientation codes are now derived purely from the face's `axis_map` plus the partner face name; any topology that stamps valid `axis_map` values on its seams produces correct `cc.par` output without writer changes. As a side effect, same-axis seams between two blocks now emit asymmetric codes (`135` on the `*_lo` side, `145` on the `*_hi` side) — both correctly describe the same seam from each block's own frame.
- **Cubed-sphere cap-equator partners point at cap *lateral* faces.** The runtime `Face.partner` for a sub-block's `i_lo` / `i_hi` is now the cap's matching lateral face (`i_lo` / `i_hi` / `j_lo` / `j_hi`), not the cap's `k_lo`. This matches the geometric face-to-face seam and lets a single rule cover both sub-side and cap-side patches in the cc.par writer.

### Fixed
- **`FlatCapsCylinder` adjacency was empty.** The 5-block flat-caps cylinder primitive built its blocks but never wired the 8 SHARED seams (4 center↔sub + 4 sub↔sub corners), so `background.adjacency.json` had every face as `FREE` / `WALL` and `background.cc.par` emitted zero connection patches. All 8 seams are now wired with the correct partners and `axis_map` values.

### Added
- `shore.mesh.face.AxisMap`: frozen dataclass describing how a face's along-seam axes map onto its partner's along-seam axes (axis correspondence + per-axis flip flags). Required on every connection face; `_set_face_bc` validates it against the partner face name at wiring time so a topology that omits the field fails loudly instead of silently producing wrong cc.par downstream.
- `shore.mesh.face.FACE_ALONG_AXES`: canonical along-axis order per face name, used for `AxisMap` interpretation.

### Documented
- Per-edge spacing section in `algorithm/spacing-laws.md` with the full support table by topology / primitive.
- `topologies/cubed-sphere.md`, `topologies/o-grid.md`, `topologies/primitives.md`, `reference/primitive.md` updated with per-edge guidance and the structural reasons each axis is or isn't wired.
- New `reference/cc-par.md` covering the `shore cc-par` command, BC-code mapping, orientation-digit decoder, and the cubed-sphere cap-equator orientation table.
- `reference/split.md` and `reference/cc-par.md` updated for the v2 adjacency JSON schema and axis_map-driven orientation codes.
