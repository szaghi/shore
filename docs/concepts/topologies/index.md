# Topologies

A **topology** is a fixed block layout: how the grid around a class of
bodies is split into hexahedral blocks, and how the blocks are connected
at their faces. The [near-body pipeline](/concepts/near-body-pipeline)
pages describe the stages (projection, k=0 construction, marching); this
section describes the layouts built from them. All of them produce
`HexBlock`s with typed faces and axis maps ([Block model](/concepts/block-model)).

## Catalogue

| Topology | Blocks | Built by | For |
|---|---|---|---|
| [Cubed sphere](cubed-sphere) | 6: 4 equatorial + 2 caps | `shore mesh` (default, `--topology cubed_sphere`) | closed star-shaped bodies; no singular axis |
| [O-grid](o-grid) | 1, j-periodic | `shore mesh --topology ogrid` | a single periodic block; leaves the two polar regions uncovered |
| [C-grid](c-grid) | 3 (sharp TE, default), 4 (blunt TE), 1 (legacy) | `shore mesh --topology cgrid` | wings and foils, with wake branches |
| [C-grid tip caps](c-grid#tip-caps) | 6 per tip (conformal) or 5 per tip (overset) | `shore mesh --cap-tips`, `shore cap-tip` | closing a finite wing at its tips |
| [OH hybrid](oh-hybrid) | 12: cubed-sphere wrap + 6 channels | `shore mesh --topology oh` | one grid from the wall to an axis-aligned box |
| [CH hybrid](ch-hybrid) | 12 (sharp TE) / 14 (blunt TE) | `shore mesh --topology ch` | one grid from a wing's wall to a chord/normal frame |
| [H-grids](h-grid) | 1 | `shore primitive h-grid-rect` / `h-grid-tube` | straight ducts without a body |
| [Primitives](primitives) | 1 or 5 | `shore primitive box` / `annulus` / `flat-caps` | backgrounds and buffers |

The topologies cover bodies of a known shape class. Junctions between
components are closed with [collar grids](/concepts/collars), and
arbitrary layouts are built by hand as a [blocking](/concepts/blocking)
on CAD or STL geometry.

## Choosing

1. **Closed body, star-shaped from an interior point** (spheres,
   fuselages, pods): the cubed sphere. Add channel blocks to a box with
   the OH hybrid when the grid must end on a box instead of overlapping a
   background.
2. **Wing or foil** with a sharp or blunt trailing edge: the C-grid, with
   tip caps for a finite wing; the CH hybrid for a single grid out to a
   frame.
3. **Duct** with no body: an H-grid.
4. **Background or buffer**: a primitive.

The O-grid is the simplest single block around a body, but it does not
cover the poles; prefer the cubed sphere when the wall grid must be
closed.

For the cubed sphere and the OH wrap, check the cap–equator spacing match
when choosing `(ni, nj, theta_cap)`:
[k=0 construction](/concepts/k0-construction#compatibility-window-and-fallback).
The CLI flags are in [`shore mesh`](/reference/mesh) and
[`shore primitive`](/reference/primitive).

## See also

- [Near-body pipeline](/concepts/near-body-pipeline), [Pipeline limits](/concepts/pipeline-limits)
- [Chimera mesh generation](/concepts/chimera): how the grids are combined
- [`shore.mesh.topology`](/reference/api/shore.mesh.topology)
