# Changelog

All notable changes to SHORE are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

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



