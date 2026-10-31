# Pipeline limits

The [near-body pipeline](/concepts/near-body-pipeline) is fast and
automatic because it assumes a lot about the body. This page lists those
assumptions, what happens when a body breaks them, and which part of SHORE
handles the case instead.

## Star-shaped bodies

The cubed-sphere, O-grid and OH topologies cast every k=0 node along a ray
from one anchor point inside the body. Each ray must hit the body once, so
the body must be **star-shaped from the anchor**. The anchor is chosen by
`resolve_anchor` (default strategy `auto`: centre of the minimum enclosing
sphere, then the bounding-box centre, each checked to lie inside the body).
Re-entrant cavities, through-holes, appendages that hide part of the
surface from the anchor, and multi-component STLs break this. The
projector then raises `MeshInputError` with the per-ray failure counts
([Surface projection](/concepts/projection#failure-modes)).

What to do instead:

- **Wings and fins**: the C-grid family slices the body at span stations
  and needs no anchor ([C-grid](/concepts/topologies/c-grid)).
- **Bodies meeting bodies** (wing–fuselage, appendages): grid each
  component separately and close the junction with a
  [collar grid](/concepts/collars), in a multi-component
  [project](/concepts/projects).
- **Anything else**: build a [blocking](/concepts/blocking) by hand (split,
  O-grid insertion, association to CAD entities) and premesh it.

## Thin first layers and concave walls

The default algebraic marcher extrudes along the grid normal and smooths.
Converging fronts at concave edges fold the grid: the marcher stops with
`MarchingError` at the first layer with a negative Jacobian, it never writes
a folded grid. The hyperbolic scheme (`--march-scheme hyperbolic`) spreads
converging fronts and completes the concave cases of its acceptance spec
where the algebraic scheme fails ([Marching](/concepts/hyperbolic-marching),
[ADR 0005](/adr/0005-marching-schemes)). Its acceptance excludes blunt-TE
C-grids, the OH/CH channels, non-star-shaped bodies and fine grids; those
remain unverified.

## Cap–equator spacing match

The cubed sphere pins the first and last meridional interval of the
equator to the cap edge (C1 spacing across the seam) only when the two
steps are within a factor 2 of each other. Outside that window the equator
falls back to its unpinned law and emits one `UserWarning`. The default
cap, 30° since [ADR 0054](/adr/0054-default-cap-30-degrees), with the
default `--ni 40 --nj 60` is inside the window; the old 9° default was
outside it. See
[k=0 construction](/concepts/k0-construction#compatibility-window-and-fallback)
for the rule.

## Outer boundary

The outer layer of a marched grid follows the body shape, offset by the
step schedule. For an overset run that is what the background grid
expects. When the volume must end on a box instead, the hybrid topologies
add non-marching channel blocks out to an axis-aligned box
([OH](/concepts/topologies/oh-hybrid)) or a chord/normal frame
([CH](/concepts/topologies/ch-hybrid)).

## O-grid poles

The single-block O-grid covers colatitudes $[\theta_{cap},\, \pi -
\theta_{cap}]$ only: its first and last i-rows are rings, and two polar
regions of half-angle `theta_cap` are left without cells
([O-grid](/concepts/topologies/o-grid)). Use the cubed sphere for a closed
wall grid.

## See also

- [Near-body pipeline](/concepts/near-body-pipeline)
- [ADR 0005](/adr/0005-marching-schemes) (marching schemes and their
  stated scope), [ADR 0010](/adr/0010-collar-block) (collars),
  [ADR 0042](/adr/0042-blocking-skeleton-premesh) (blocking)
- [Roadmap](/adr/roadmap): open items of each phase
