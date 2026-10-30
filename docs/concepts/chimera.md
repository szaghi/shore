# Chimera mesh generation

Chimera (overset) grid generation covers a complex domain with several
structured grids that **overlap** instead of one grid that conforms
everywhere. Each grid is built for one job; the flow solver connects them
at run time by interpolating between overlapping cells. SHORE targets the
Xall/Xnavis chain, and writes the files its overset pre-processor reads.

## The model

An overset assembly has:

- **Near-body grids**: body-fitted blocks around each component, fine at
  the wall.
- **A background grid**: simple blocks that fill the far field.
- **Buffer grids** (optional): intermediate grids that bridge the
  resolution jump between the two.
- **Overlap**: regions covered by more than one grid. The solver cuts
  holes where a grid lies inside a body (driven by marker boxes), marks
  fringe cells at hole and outer boundaries, and finds donor cells in the
  other grids to interpolate them from.

<Figure name="chimera-assembly" caption="A cubed-sphere near-body grid inside a cylindrical background, cut away to show the overlap." />

The gain is decoupling: each grid is generated independently and kept
simple, so a multi-component body becomes a set of near-body grids plus a
background rather than one block topology for the whole configuration.

## What SHORE produces

| Role | Tool |
|---|---|
| Near-body grids | `shore mesh` ([topologies](/concepts/topologies/)), [collar grids](/concepts/collars) at component junctions, [blockings](/concepts/blocking) for bodies the topologies cannot wrap |
| Background and buffers | `shore primitive`: [box, annulus, flat-caps, H-grids](/concepts/topologies/primitives) |
| Hole cutting input | [marker boxes](/concepts/marker-boxes) from the body (`shore boxes-from-stl`) |
| Solver input | `.grd` (all blocks), `cc.par` (BCs, connections, marker boxes), `proc.input` (rank assignment) |
| Load balance | [split planning](/concepts/load-balancing) (`shore balance`, `shore split`) |

SHORE does **not** cut holes, identify fringes or search donors: that is
the overset pre-processor's job. Its boundary is to deliver grids in world
coordinates with a correct face graph, so the pre-processor reads every
block connection and boundary condition as intended
([ADR 0036](/adr/0036-cc-par-as-xall-reads-it)).

## From blocks to `cc.par`

Each block face carries a `FaceBC`. `shore cc-par` translates the
`.adjacency.json` face graph:

| `FaceBC` | `cc.par` |
|---|---|
| `WALL` | 1 (viscous wall) |
| `SYMMETRY` | natural BC 2 |
| `SHARED` / `DIRICHLET` | a connection patch with a 3-digit orientation code from the `AxisMap` |
| `FREE` | the block's `free_face_bc`: 40 (chimera fringe) by default; a positive natural code in `[1, 19]` (e.g. 9 extrapolation, 3/4 inflow/outflow) for outer boundaries |

Several components (a near-body grid, a background, a collar) are merged
with `shore grd-merge` / `shore cc-par-merge` / `shore proc-input-merge`,
or written in one step from a multi-component
[project](/concepts/projects) with `shore xall`
([ADR 0016](/adr/0016-multi-component-projects),
[ADR 0017](/adr/0017-xall-options-per-component)). The end-to-end runs are
in the [chimera assembly](/tutorials/chimera-assembly) and
[Xall bundle](/tutorials/xall-bundle) tutorials.

## Requirements on near-body grids

1. **On the wall.** The k=0 layer lies on the body: ray hits on an STL,
   snapped to the exact surface for CAD bodies.
2. **Enough overlap.** The outer layer must reach far enough into the
   background for the solver's interpolation stencils. How far depends on
   the solver's stencil and on the multigrid levels in use: every block
   must keep enough coarse-level cells (the `--ngr` / `--min-coarse-cells`
   floor of [load balancing](/concepts/load-balancing)).
3. **No folded cells.** The marcher stops at the first negative Jacobian;
   [`shore check`](/reference/check) and the
   [quality metrics](/concepts/quality) are the independent check before
   the grids go to the solver.

## Requirements on the background

1. **Resolution matched at the overlap.** Background cells comparable to
   the outer cells of the near-body grids give the donor search
   well-shaped donors.
2. **Far boundary far away.** Far-field conditions imposed too close to
   the body change the solution.

The `flat-caps` primitive fills a cylinder without an inner hole, so the
background itself needs no internal fringe; an `annulus` plus a box needs
one between them but imposes no alignment.

## Bodies the pipeline cannot wrap

The ray-cast topologies need a body star-shaped from one point. For
assemblies of bodies (wing–fuselage, appendages) SHORE grids each
component on its own and closes each junction with a collar grid built
from the bodies' intersection curve ([ADR 0010](/adr/0010-collar-block));
half models use half-collars on the symmetry plane
([ADR 0013](/adr/0013-half-collar-block)). Where no topology fits, a
blocking is built by hand, associated to the CAD entities and premeshed
([ADR 0042](/adr/0042-blocking-skeleton-premesh)); the phase-7 acceptance
case is a half wing–body built that way
([ADR 0051](/adr/0051-wing-body-acceptance)). Open items are on the
[roadmap](/adr/roadmap).

## See also

- [Chimera assembly tutorial](/tutorials/chimera-assembly), [Xall bundle tutorial](/tutorials/xall-bundle)
- [Marker boxes](/concepts/marker-boxes), [Load balancing](/concepts/load-balancing)
- [`shore cc-par`](/reference/cc-par), [`shore xall`](/reference/xall), [File formats](/reference/formats)
