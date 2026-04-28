# Chimera mesh generation

Chimera (overset) mesh generation is a strategy for discretising
complex geometries in CFD by assembling multiple overlapping
structured grids rather than building a single monolithic mesh.

## The Chimera concept

A Chimera mesh consists of:

- **Near-body grids** — body-fitted structured blocks around each
  surface component (one per body in a multi-body case).
- **Background grid** — a simple structured block that fills the far
  field.
- **Buffer grids** — *optional* intermediate meshes that bridge a
  resolution jump between near-body and background.
- **Overlap region** — cells from different grids that spatially
  overlap; the CFD solver interpolates between them at run time.

```
         Background grid (cylindrical or Cartesian)
    ┌──────────────────────────────────────┐
    │                                      │
    │    ╔══════════════╗                  │
    │    ║  Near-body   ║  ← cubed-sphere  │
    │    ║  cubed-sphere║    or O-grid     │
    │    ╚══════════════╝                  │
    │                                      │
    └──────────────────────────────────────┘
              overlap zone
```

The key advantage: each grid is generated independently and optimised
for its purpose (body-fitted near-body, simple far-field). Complex
multi-component geometries become a collection of independent
near-body grids rather than one intractable topology problem.

## SHORE's role

SHORE produces **all the structured grids** an overset assembly
needs:

- **Near-body** via `shore mesh` — choose
  [cubed sphere](/topologies/cubed-sphere) (default, 6 blocks, no
  polar singularity) or [O-grid](/topologies/o-grid) (single block,
  j-periodic, with polar singular axes).
- **Background and buffer** via `shore primitive` — pick from
  [box](/topologies/primitives) (Cartesian), [annulus](/topologies/primitives)
  (annular cylinder, optionally a partial sector), or
  [flat-caps cylinder](/topologies/primitives) (5-block conforming
  cylinder with no central hole — the Chimera-fringe-free
  alternative when you want to avoid the cost of overset
  interpolation in the radial direction).

What SHORE **deliberately does not do**: hole cutting, fringe
identification, donor search. Those are the solver's job. SHORE's
boundary is "produce the grids in their world coordinates, write
them to disk, optionally bundle them into one ParaView-friendly
`.vtm`." The downstream Chimera-capable solver (Xall,
OVERFLOW, etc.) does the assembly at run time.

## A worked end-to-end assembly

The first-tentative assembly example
(`examples/assembly_sphere_in_cylinder.py`) produces the simplest
non-trivial overset stack: a 6-block cubed-sphere mesh wrapping a
unit-radius sphere, plus a 5-block flat-caps cylinder background
spanning `[-10, 10]² × [-10, 50]`. See
[Assembly walkthrough](/guide/assembly) for the parameter choices,
output layout, and ParaView inspection workflow.

```bash
python examples/assembly_sphere_in_cylinder.py --no-view
# → examples_out/{fuselage_*,background_*}.geo  +  examples_out/assembly.vtm
```

The `.vtm` is the load-bearing piece: opening it in ParaView shows
all 11 blocks in their actual world positions, lets you slice
through the overlap region, and verifies that the near-body
outer shell sits well inside the background's central block.

## Workflow with a full Chimera solver

```
1.  shore mesh fuselage.stl --ni 80 --nj 120 --nk 50 -o fuselage
        → fuselage_sub0..3.geo, fuselage_cap_north/south.geo  (6 files)
2.  shore primitive flat-caps -c background.toml -o background
        → background_center.geo, background_sub_e/n/w/s.geo    (5 files)
3.  shore export fuselage_*.geo background_*.geo -o assembly.vtm
        → ParaView-friendly bundle for visual sanity check
4.  Hand to Chimera-capable solver for the assembly + CFD run.
```

## Near-body grid requirements

For a Chimera solver to accept a near-body grid:

1. **Body conforming** — the k=0 layer must lie on the body surface.
   SHORE's projection guarantees this within ray-casting precision.
2. **Adequate overlap** — the outer boundary (k=nk-1) must extend
   far enough into the background grid to establish a valid
   interpolation stencil. Typical Chimera solvers require 3–5 cells
   of overlap.
3. **No inverted cells** — `shore check` must exit 0 before passing
   the grid to the solver. The marcher's per-layer Jacobian guard
   already catches inversions during generation; `check` is the
   independent cross-check.

## Background grid requirements

Background grids have looser constraints than near-body grids since
they don't conform to a body, but two things still matter:

1. **Cells coarse enough to be cheap, fine enough to overlap** with
   the near-body grids' outer shells. The simplest robust choice is
   a background cell size ~2–4× the near-body outer-shell cell size.
2. **Domain large enough to hold the physics** — far-field BCs
   imposed too close to the body affect the solution. A rule of
   thumb is 10–20 body-lengths in each direction for external
   aerodynamics.

The `shore primitive flat-caps` primitive is particularly suited to
GPU-friendly backgrounds: the 5-block conforming layout has no inner
fringe to interpolate across (the solver only needs Chimera fringes
between background and near-body).

## Multi-block extension (roadmap)

For bodies that are not star-shaped w.r.t. their centroid (e.g. an
aircraft configuration with a wing-fuselage junction), a single
near-body grid cannot wrap the full surface. The planned multi-block
extension will:

1. Decompose the surface into structured patches by tracing feature
   lines.
2. Generate one near-body block per patch.
3. Stitch adjacent blocks at their shared faces (C0 continuity).
4. Output multiple geo files for the solver's multi-block assembly.

This is a significantly harder problem — contributions and
discussion are welcome.
