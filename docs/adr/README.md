# Architecture decision records

These records fix the decisions that take SHORE from a semi-automatic
STL → Chimera pipeline to a CAD-associative structured/overset mesher with an
interactive GUI. Each record states the decision, the forces behind it, the
alternatives rejected and the consequences accepted. A record is changed only
by a later record that supersedes it.

| # | Title | Status |
|---|---|---|
| [0001](./0001-scope.md) | Scope: a CAD-associative structured/overset mesher, not a general mesher | Proposed |
| [0002](./0002-geometry-backend.md) | Geometry protocol with STL and OCCT (build123d/OCP) backends | Proposed |
| [0003](./0003-association.md) | Block-to-geometry association replaces `body_projected` | Proposed |
| [0004](./0004-command-layer.md) | Command layer and replayable journal under CLI, API and GUI | Proposed |
| [0005](./0005-marching-schemes.md) | Name the current marcher "algebraic"; add a true hyperbolic scheme | Proposed |
| [0006](./0006-gui-stack.md) | GUI: PySide6 + pygfx/wgpu, desktop, Python + numpy only | Accepted |
| [0007](./0007-surface-grids.md) | Surface grids: TFI + elliptic smoothing in (u, v) per face | Proposed |
| [0008](./0008-surface-marching.md) | Hyperbolic surface marching: the collar-grid kernel | Proposed |
| [0009](./0009-intersection-curves.md) | Intersection loops of two bodies, for collar grids | Proposed |
| [0010](./0010-collar-block.md) | The collar block of two intersecting bodies | Accepted |
| [0011](./0011-open-intersection-curves.md) | Open intersection curves ending on a symmetry plane | Accepted |
| [0012](./0012-march-ends-on-a-plane.md) | Surface marching with ends on a symmetry plane | Accepted |
| [0013](./0013-half-collar-block.md) | The half-collar block on a symmetry plane | Accepted |
| [0014](./0014-crease-aware-stl-normals.md) | Crease-aware STL surface normals | Accepted |
| [0015](./0015-stl-normals-for-marching-blend.md) | STL normals for the marcher's normal blending | Accepted |
| [0016](./0016-multi-component-projects.md) | Multi-component, portable projects | Accepted |
| [0017](./0017-xall-options-per-component.md) | Per-component Xall options in `write_xall` | Accepted |
| [0018](./0018-cli-collars-and-projects.md) | CLI for collars and multi-component projects | Accepted |
| [0019](./0019-deprecate-normal-blending.md) | Deprecate the marcher's normal blending (`blend_normals_k`) | Accepted (removal: 0034) |
| [0020](./0020-junction-curve-spacing.md) | Spacing laws along the junction curve | Accepted |
| [0021](./0021-junction-curves-with-corners.md) | Junction curves with corners | Accepted |
| [0022](./0022-collars-at-blunt-corners.md) | Graded spacing across corners; collars at blunt corners | Accepted |
| [0023](./0023-surface-marching-at-creases.md) | Surface marching at creases | Accepted |
| [0024](./0024-split-collar-at-sharp-corners.md) | The split collar at sharp corners | Accepted |
| [0025](./0025-collar-volume-at-corners.md) | The collar volume at loop corners | Accepted |
| [0026](./0026-smoothed-collar-volume-by-default.md) | The smoothed collar volume by default | Accepted |
| [0027](./0027-surface-marching-at-clustered-corners.md) | Surface marching at clustered corners | Accepted |
| [0028](./0028-collar-volume-at-extreme-clustering.md) | The collar volume at extreme clustering (investigation) | Proposed |
| [0029](./0029-spreading-converging-fronts.md) | Spreading converging fronts (not adopted); a four-corner surface metric | Rejected / Accepted |
| [0030](./0030-concave-notches.md) | Collars in concave notches (investigation) | Rejected |
| [0031](./0031-polishing-exact-cycles.md) | Polishing: stop at exact cycles | Accepted |
| [0032](./0032-stl-projection-cache.md) | STL projection: each query point computed once | Accepted |
| [0033](./0033-load-without-replay.md) | Loading a saved project without replay | Accepted |
| [0034](./0034-remove-normal-blending.md) | Remove the marcher's normal blending (`blend_normals_k`) | Accepted |
| [0035](./0035-cgns-plot3d-export.md) | CGNS and Plot3D exports | Accepted |
| [0036](./0036-cc-par-as-xall-reads-it.md) | cc.par connections as Xall reads them | Accepted |
| [0037](./0037-gui-viewer-command-shell.md) | GUI slice 1: viewer and command shell | Accepted |
| [0038](./0038-gui-geometry-and-slices.md) | GUI slice 2: geometry display and slice planes | Accepted |
| [0039](./0039-association-editing.md) | Association editing (GUI slice 3) | Accepted |
| [0040](./0040-seams-kept-closed.md) | Seams kept closed by in-place edits | Accepted |
| [0041](./0041-redistribute-edge.md) | Redistributing the nodes along a block edge (GUI spacing editor) | Accepted |
| [0042](./0042-blocking-skeleton-premesh.md) | Blocking: a block skeleton, split and premeshed (phase 7, slice 1) | Accepted |
| [0043](./0043-ogrid-insertion.md) | O-grid insertion in a blocking (phase 7, slice 2) | Accepted |
| [0044](./0044-blocking-selection.md) | Selecting blocking faces and blocks in the GUI (phase 7, slice 3) | Accepted |
| [0045](./0045-move-blocking-vertex.md) | Moving blocking vertices (phase 7, slice 4) | Accepted |
| [0046](./0046-wing-blocking-investigation.md) | A wing on a plane with the blocking tools | Investigation |
| [0047](./0047-walls-from-outside-rejected.md) | Wall entities projected from across the block | Rejected |
| [0048](./0048-orthogonal-wall-lines.md) | Grid lines leaving the walls along their normal | Accepted |
| [0049](./0049-stl-feature-curves.md) | Feature curves and points of STL geometry | Accepted |
| [0050](./0050-wall-layer.md) | Wall layers with an absolute first cell (algebraic, after premesh) | Proposed |

The phased plan that implements them is in [roadmap.md](./roadmap.md).

## Template

```markdown
# NNNN — Title

- **Status:** Proposed | Accepted | Superseded by NNNN
- **Date:** YYYY-MM-DD

## Context
## Decision
## Alternatives considered
## Consequences
## Verification
```
