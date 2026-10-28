# Changelog

All notable changes to SHORE are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.4.0] — 2026-10-28
### Added
- **collar**: Cut the collar at sharp junction corners with split_collar

- **cli**: Build split collars from projects and shore collar --split

- **collar**: Smooth the collar's outer boundary along the loop at corners

- **collar**: Build collars with the smoothed outer boundary by default


## [0.3.0] — 2026-10-27
### Added
- **collar**: Distribute nodes along the junction curve by law and curvature

- **surface**: Trace junction curves through corners of either body

- **surface**: Grade the junction-curve spacing across corners

- **surface**: March surfaces cleanly from fronts with corners on creases


### Fixed
- **geometry**: Route OCCT reader console output to the debug log

- **mesh**: Zero-initialise topology node buffers before marching

- **surface**: Trace junction corners on thin wedges and sliver STLs

- **shapes**: Write extruded wings with outward, consistent normals

- **io**: Declare networkx, needed by load_stl's normal repair


## [0.2.0] — 2026-10-27
### Added
- **geometry**: Add Geometry protocol with STL and OCCT backends

- **mesh**: Accept STEP/IGES/BREP input across the meshing pipeline

- **mesh**: Associate block entities with geometry ⚠ BREAKING CHANGE

- **project**: Add command layer with replayable journal and shore run

- **project**: Run the Xall pipeline commands over the command layer

- **quality**: Add grid metrics and freeze hyperbolic-marching spec

- **volume**: Add Steger-Chaussee hyperbolic marching scheme

- **surface**: Add TFI and elliptic surface-grid kernel in (u, v)

- **project**: Rebuild a block face on a CAD surface with surface_grid

- **surface**: Add hyperbolic surface marching for collar grids

- **surface**: Trace the intersection loop of two bodies for collars

- **collar**: Build the collar block of two intersecting bodies

- **surface**: Trace open intersection curves ending on a symmetry plane

- **surface**: Hold surface-marching ends on a symmetry plane

- **collar**: Build half-collar blocks on a symmetry plane

- **project**: Hold several geometries and components in a portable project

- **project**: Pass per-component Xall options through write_xall

- **cli**: Build collars and multi-component projects from the CLI

- **mesh**: Deprecate blend_normals_k and refuse it on C-grids


### Documentation
- **adr**: Add architecture decisions and GUI spike for the mesher roadmap


### Fixed
- **volume**: Limit Laplacian smoothing so thin first layers cannot fold ⚠ BREAKING CHANGE

- **volume**: Keep hyperbolic explicit dissipation from acting across seams

- **volume**: March the full step on convex walls in hyperbolic scheme

- **surface**: Reject seam-crossing and singular faces in surface_grid

- **geometry**: Give STL bodies crease-aware surface normals

- **mesh**: Restore vertex normals for the marcher's normal blending


## [0.1.2] — 2026-06-22
### Added
- **balance**: Multi-way equal-partition planner and multi-axis split


## [0.1.1] — 2026-06-18
### Added
- **cli**: Add curl|bash installer + shore --version flag

- **body**: Add shore.profiles + shore.shapes for body generation

- **cli**: Emit .adjacency.json sidecar from mesh + primitive

- **examples**: Add bash sibling for the C-grid example

- **mesh**: Split C-grid into 3-block (sharp) / 4-block (blunt) topology ⚠ BREAKING CHANGE

- **primitive**: Add H-grid baseline for ducted internal flow

- **mesh**: Add OHTopology — cubed-sphere wrap + H-channel hybrid

- **mesh**: Add C-grid wing-tip caps with conformal + chimera paths

- **mesh**: Add CHTopology — C-grid wrap + H-channel far-field hybrid

- **skill**: Ship Claude Code skill for operating the shore CLI

- **cli**: Expose --topology oh and remaining CHTopology knobs


### Documentation
- **claude**: Rewrite CLAUDE.md for two-half architecture + new topologies

- **skill**: Sync shore skill with --topology oh and new CH knobs


### Fixed
- **mesh**: Make OHTopology channel routing body-centre relative

- **deps**: Declare click as a direct dependency

- **cli**: Stop depending on click for CLI-flag precedence

- **spacing**: Correct tanh_two_sided interior clustering


## [0.1.0] — 2026-04-30
### Added
- **cli**: Add yinyang surface generation and visualization support ⚠ BREAKING CHANGE

- **topology**: Add 6-block sphere polar cap topology with quality controls

- **scaffold**: Add exception handling and tests for surface generation ⚠ BREAKING CHANGE

- Add author information to pyproject.toml

- **spacing**: Add tanh/tanh2 spacing laws for surface and volume distributions

- Remove bugged correct-spacing algorithm ⚠ BREAKING CHANGE

- **mesh**: Implement comprehensive mesh topology and geometry components ⚠ BREAKING CHANGE

- **surface**: Implement comprehensive cap and anchor surface generation with validation

- **caps**: Rebuild caps topology with gnomonic cap and seam-aware marching ⚠ BREAKING CHANGE

- **mesh**: Smooth STL normals via barycentric vertex-normal blend

- **cli**: Add TOML config-file driven invocation for `shore mesh`

- **primitive**: Add box, annulus, and flat-caps cylinder primitives

- **mesh**: Promote single-block O-grid to a first-class topology

- **examples**: Add first-tentative chimera assembly example

- **split**: Add BlockSplit + adjacency sidecar for parallel decomp

- **spacing**: Add per-edge spacing data model with TFI on primitives

- **io**: Add cc.par writer to bridge SHORE adjacency to Xnavis

- **cli**: Expose --meta flag on shore cc-par for per-block BlockMeta

- Add shore grd command and rename Xnavis references to Xall

- **mesh**: Track seam orientation via AxisMap on every connection face ⚠ BREAKING CHANGE

- **io**: Add merge_grd and merge_cc_par for multi-mesh Xall jobs

- **io**: Add marker boxes for cc.par chimera hole-cutting

- **io**: Add proc.input writer for Xall MPI work distribution

- **balance**: Add iterative split planner for MPI load balancing

- **mesh**: Add C-grid topology for sharp-TE bodies


### Changed
- **mesh**: Rename CapsTopology to CubedSphereTopology ⚠ BREAKING CHANGE

- **examples**: Organize into 5 numbered subdirs by feature

- **examples**: Rename cc_par_assembly.py to chimera_assembly.py


### Documentation
- Restructure into algorithm + topologies sections, add figures

- **topologies**: Drop the yin-yang topology page and references

- Rename caps-sphere example, sweep up residual cap-topology refs

- Add assembly walkthrough; refresh chimera page

- Refactor for first release — accuracy, scope, structure

- **cli**: Clarify ni/nj/nk are node counts, not cell counts

- Expand authors section to credit Andrea Di Mascio and Claude


### Fixed
- **tests**: Force PyVista off-screen rendering on headless CI

- **cli**: Make help text ASCII-safe on Windows + sort config imports

- **io**: Make voxel-fill cover the whole body and merge into slabs


### build
- Drop Python 3.10 support, require >=3.11 ⚠ BREAKING CHANGE



