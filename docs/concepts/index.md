---
title: How SHORE works
---

# How SHORE works

SHORE makes structured multi-block grids for overset (Chimera) CFD and
packages them for the Xall/Xnavis solver chain. It has two halves:

1. **Grid generation**: geometry (STL or CAD) in, structured hexahedral
   blocks out. Near-body grids are marched from the wall; junctions get
   collar grids; anything else is blocked by hand, as in ICEM CFD Hexa.
2. **Solver packaging**: the block list out to the files the solver reads:
   `.geo` per block, a `.grd` pack, `cc.par` connectivity with marker boxes
   for hole cutting, and a `proc.input` MPI distribution, with a
   load-balancing split planner. CGNS, Plot3D and VTK cover other solvers
   and viewers.

## The data model

```mermaid
flowchart LR
    G[Geometry<br/>STL / STEP / IGES] --> A[Association<br/>fingerprints]
    G --> T[Topology or primitive<br/>or collar]
    G --> B[Blocking<br/>skeleton]
    T --> M[Blocks<br/>HexBlock + face graph]
    B -- premesh --> M
    A --> M
    M --> P[Project<br/>commands + journal]
    P --> X[Xall: .geo .grd<br/>cc.par proc.input]
    P --> C[CGNS / Plot3D]
    M --> V[VTK]
```

- **Geometry** is reached through one protocol with an STL and a CAD
  backend; entities (surfaces, curves, points) carry fingerprints that
  survive a re-read of an edited file ([Geometry](/concepts/geometry)).
- **Blocks** are `HexBlock`s: a `(ni, nj, nk, 3)` node array with typed
  vertices, edges and faces; faces carry boundary conditions and seams to
  other blocks ([Block model](/concepts/block-model)).
- **Associations** tie block faces, edges and vertices to geometry entities;
  `snap` projects the nodes onto them ([Association](/concepts/association)).
- **The project** applies every change as a command, can undo it, and
  records a journal that replays the grid bit for bit
  ([Projects and journals](/concepts/projects)).
- **Quality** is measured by a fixed set of metrics that the acceptance
  specs are written in ([Grid quality](/concepts/quality)).

## Four ways to make a grid

| Way | What it does | Use it for | Entry points |
|---|---|---|---|
| Marched topologies | project a template onto the body, then extrude layer by layer along the wall normal | near-body grids of star-shaped bodies (cubed sphere, O-grid) and wings (C-grid, CH/OH hybrids with far-field channels) | `shore mesh`, `p.build_topology` + `p.march` |
| Primitives | analytic structured blocks, no geometry | backgrounds, ducts, test grids (box, annulus, flat-caps cylinder, H-grids) | `shore primitive`, `p.build_primitive` |
| Collars | intersection curve of two bodies, a surface grid marched on each, a volume block over the junction | wing–body, appendage–hull junctions in an overset assembly | `shore collar`, `p.build_collar` |
| Interactive blocking | a skeleton of blocks split, O-gridded, associated and premeshed onto the geometry | bodies the templates cannot handle, single-grid junction regions, CAD-driven blocking | `p.create_blocking` … `p.premesh`, the GUI |

The four write the same block model, so the outputs, the quality checks,
the project commands and the GUI work alike on all of them. A Chimera run
usually combines several, each as a project **component**: a body-fitted
grid, a collar, a background box ([Choosing a workflow](/guide/workflows)).

## Concepts

<div class="shore-cards">
<a class="shore-card" href="/shore/concepts/block-model"><strong>Block model</strong><span>HexBlock, faces, boundary conditions, axis maps, seams</span></a>
<a class="shore-card" href="/shore/concepts/geometry"><strong>Geometry</strong><span>STL and CAD behind one protocol; fingerprints</span></a>
<a class="shore-card" href="/shore/concepts/association"><strong>Association</strong><span>block entities linked to geometry, and snapping</span></a>
<a class="shore-card" href="/shore/concepts/projects"><strong>Projects and journals</strong><span>commands, undo, replay, save and load</span></a>
<a class="shore-card" href="/shore/concepts/quality"><strong>Grid quality</strong><span>the frozen metrics and what the specs require</span></a>
<a class="shore-card" href="/shore/concepts/near-body-pipeline"><strong>Near-body pipeline</strong><span>from STL to a marched volume, stage by stage</span></a>
<a class="shore-card" href="/shore/concepts/projection"><strong>Surface projection</strong><span>ray casting a template onto the body</span></a>
<a class="shore-card" href="/shore/concepts/k0-construction"><strong>k=0 construction</strong><span>the cubed-sphere wall grid: gnomonic caps and seams</span></a>
<a class="shore-card" href="/shore/concepts/hyperbolic-marching"><strong>Marching schemes</strong><span>algebraic and hyperbolic layer extrusion</span></a>
<a class="shore-card" href="/shore/concepts/seam-aware-normals"><strong>Seam-aware normals</strong><span>centred tangents across block seams</span></a>
<a class="shore-card" href="/shore/concepts/spacing-laws"><strong>Spacing laws</strong><span>uniform, tanh, tanh2, geometric, pinned</span></a>
<a class="shore-card" href="/shore/concepts/pipeline-limits"><strong>Limits</strong><span>what the marched pipeline does not handle</span></a>
<a class="shore-card" href="/shore/concepts/topologies/"><strong>Topologies</strong><span>O-grid, cubed sphere, C-grid, H-grid, hybrids, primitives</span></a>
<a class="shore-card" href="/shore/concepts/c-grid-construction"><strong>C-grid construction</strong><span>trailing-edge detection, sections, wake</span></a>
<a class="shore-card" href="/shore/concepts/airfoil-profiles"><strong>Airfoil profiles</strong><span>NACA sections and extruded wings</span></a>
<a class="shore-card" href="/shore/concepts/surface-grids"><strong>Surface grids</strong><span>TFI and elliptic smoothing on a CAD face</span></a>
<a class="shore-card" href="/shore/concepts/intersection-curves"><strong>Intersection curves</strong><span>tracing where two bodies meet</span></a>
<a class="shore-card" href="/shore/concepts/surface-marching"><strong>Surface marching</strong><span>hyperbolic grids grown on a surface</span></a>
<a class="shore-card" href="/shore/concepts/collars"><strong>Collar blocks</strong><span>volume grids over a junction</span></a>
<a class="shore-card" href="/shore/concepts/blocking"><strong>Blocking</strong><span>skeleton, split, O-grid, premesh, wall layers</span></a>
<a class="shore-card" href="/shore/concepts/chimera"><strong>Chimera grids</strong><span>overset components and how they connect</span></a>
<a class="shore-card" href="/shore/concepts/marker-boxes"><strong>Marker boxes</strong><span>hole-cutting boxes from a body STL</span></a>
<a class="shore-card" href="/shore/concepts/load-balancing"><strong>Load balancing</strong><span>splitting blocks across MPI ranks</span></a>
</div>

## See also

- [ADR index](/adr/) and [roadmap](/adr/roadmap): the decisions behind each
  concept, with their measured limits
- [Tutorials](/tutorials/), [Cookbook](/cookbook/), [CLI reference](/reference/),
  [Python API](/reference/api/)
