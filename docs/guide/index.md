# What is SHORE?

SHORE (**S**tructured **H**yperbolic **O**verlapping **R**emeshing
**E**ngine) makes the structured multi-block grids of an overset (Chimera)
CFD case and writes the files the solver reads. It is a Python package with
a command-line interface (`shore`), a scriptable `Project` API and a desktop
GUI; all three issue the same commands and record the same replayable
journal.

<Figure name="hero-wing-body" caption="A half wing–body blocked, premeshed and given a wall layer with SHORE's blocking commands (tutorial 7)." />

## What it makes

| Grid | How | When to use it |
|---|---|---|
| Near-body grids | a body surface marched outwards, layer by layer, in one of six topologies: O-grid, cubed sphere, C-grid (with wing-tip caps), OH and CH hybrids | compact bodies and wings; minutes from an STL to a grid |
| Background grids | parametric primitives: box, annular cylinder, flat-caps cylinder, H-grid duct and tube | the far field of an overset assembly |
| Collars | the intersection curve of two bodies, a surface grid marched on each, and the volume around the junction | wing–fuselage and other junctions in an overset assembly |
| Blocked grids | a block skeleton split, holed and O-gridded interactively, associated with the geometry and premeshed | bodies no single topology fits; point-matched multi-block grids |

Geometry comes as STL, or as STEP, IGES or BREP through OpenCascade
(`shore[cad]`). CAD entities are projected exactly and fingerprinted, so a
journal written on one revision of a model replays on the next.

## What it writes

- The **Xall/Xnavis** bundle: `.geo` blocks, the `.grd` pack, `cc.par`
  (boundary conditions, connection patches with orientation codes, marker
  boxes) and `proc.input` (MPI ranks), merged across components for a
  multi-component overset run.
- **CGNS** (structured zones, 1-to-1 connectivity, wall and symmetry BCs),
  **Plot3D** and **VTK** for other solvers and for ParaView.

Hole cutting, fringe identification and donor search belong to the solver:
SHORE writes every grid and descriptor it needs, and stops there.

## How it is built

- **Explicit stages.** Geometry, association, topology or blocking,
  marching or premesh, outputs: each is a command you can inspect, undo and
  rerun.
- **Loud failures.** Marching stops at the first inverted cell; seams are
  validated when blocks are wired; commands that would invert a block are
  refused. SHORE does not write a grid it knows is invalid.
- **Frozen acceptance specs.** Each algorithm was accepted against a spec
  frozen before it was written (quality metrics, seam closure, replay bit
  for bit); the [design records](/adr/) say what was decided, measured and
  left out.

## Limits

- The marched topologies need a body that the topology fits: star-shaped
  from its centroid (O-grid, cubed sphere, OH) or wing-like (C-grid, CH).
  Other shapes are [blocked](/concepts/blocking) or built from components
  and collars.
- Blocking is interactive: SHORE does not decompose an arbitrary body into
  blocks by itself.
- Volume grids are marched or algebraic (TFI, Hermite wall lines,
  geometric wall layers); there is no 3-D elliptic smoothing of block
  interiors yet. The [roadmap](/adr/roadmap) lists what is open.

## Next

- [Install SHORE](/guide/installation), then the [quick start](/guide/quickstart).
- [Choose a workflow](/guide/workflows) for your body.
- Learn end to end with the [tutorials](/tutorials/).
