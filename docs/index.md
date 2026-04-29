---
layout: home

hero:
  name: SHORE
  text: Structured Hyperbolic Overlapping Remeshing Engine
  tagline: Hyperbolic mesh generation for Chimera CFD — body-fitted near-body grids and parametric backgrounds, in one tool.
  actions:
    - theme: brand
      text: Quick Start
      link: /guide/quickstart
    - theme: alt
      text: Guide
      link: /guide/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/shore

features:
  - icon: 🏗️
    title: Full Chimera assembly in one tool
    details: Generate every grid an overset assembly needs — body-fitted near-body meshes and parametric backgrounds — in the same world frame. Pack them into the Xall solver bundle (.grd + cc.par + marker boxes) without leaving the SHORE CLI.
    link: /guide/cc-par-pipeline
    linkText: Pipeline walkthrough
  - icon: 📐
    title: Body-fitted near-body meshes
    details: Three topologies for surface inputs (STL today; CAD on the roadmap) — single-block O-grid, 6-block cubed sphere with gnomonic cap k=0 seed and C1 cap-equator seam pinning, and a single-block C-grid for sharp-TE bodies (airfoils, hydrofoils) with auto TE detection and downstream wake branches. One `shore mesh` call from a triangulated body to Xall-compatible .geo files.
    link: /topologies/cubed-sphere
    linkText: Cubed sphere
  - icon: 🧱
    title: Parametric backgrounds
    details: Box, annular cylinder, and 5-block flat-caps cylinder primitives with per-axis spacing laws — fill the far field without an STL. Per-edge spacing on Box and Annulus, transfinite-interpolation blends; flat-caps cylinder is Chimera-fringe-free for GPU-friendly backgrounds.
    link: /topologies/primitives
    linkText: Primitives
  - icon: 🌊
    title: Seam-aware hyperbolic marching
    details: Each wall-normal layer advances along outward normals computed with centred differences across multi-block seams (ghost rows from partner blocks). Spacing continuity across block edges holds layer by layer; per-layer Laplacian smoothing + Jacobian guard catch any cell inversion immediately.
    link: /algorithm/hyperbolic-marching
    linkText: Algorithm
  - icon: 🔍
    title: Focused CLI
    details: Nineteen commands cover the full chimera pipeline — generate (`mesh`, `primitive`), subdivide and balance for parallel solvers (`split`, `balance`), pack for the solver (`grd`, `cc-par`, `proc-input`, `boxes-from-stl`), merge multi-component runs (`grd-merge`, `cc-par-merge`, `proc-input-merge`), inspect (`info`, `check`), and visualise (`export`, `view`, `boxes-export`). A Rich progress bar shows per-layer Jacobians during marching.
    link: /reference/mesh
    linkText: CLI reference
  - icon: 🆓
    title: Free and open source
    details: Released under the MIT license. Pure Python — depends only on NumPy, trimesh, Typer, and Rich. No compiled extensions, no C++ meshing backends.
---

## What is SHORE?

SHORE produces every structured grid a Chimera (overset) CFD assembly needs and packages them for the solver. Two product lines, one tool:

**Body-fitted near-body meshes** by hyperbolic normal extrusion of a body surface. Surface input is a triangulated STL today; CAD (STEP / IGES / BREP) is on the roadmap. Two topologies:

```bash
# 6-block cubed sphere (default) — no polar singularity, conforming seams
shore mesh fuselage.stl --ni 60 --nj 80 --nk 40 --ds 1e-4 --growth 1.15 -o fuselage -v

# Single-block O-grid — j-periodic, one .geo file
shore mesh fuselage.stl --topology ogrid --ni 60 --nj 80 --nk 40 -o fuselage
```

**Parametric backgrounds** that fill the far field — no STL needed, just analytic geometry:

```bash
# Cartesian box
shore primitive box --min -10 -10 -10 --max 10 10 50 --ni 30 --nj 30 --nk 60 -o background

# 5-block flat-caps cylinder (Chimera-fringe-free; central square + 4 trapezoidal sub-blocks)
shore primitive flat-caps --r-out 10 --z0 -10 --z1 50 --ni 16 --nc 24 --nk 60 -o background
```

Then assemble with the Xall pipeline:

```bash
# Pack each component into binary .grd
shore grd wall_*.geo -o wall.grd
shore grd background_*.geo -o background.grd

# Generate the chimera descriptor with marker boxes for hole-cutting
shore boxes-from-stl fuselage.stl --block sub0 --max-outward-gap 0.10 -o boxes.json
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd --boxes boxes.json
shore cc-par background.adjacency.json background.cc.par --grd background.grd --free-bc -9

# Merge multi-component runs into a single Xall job
shore grd-merge    background.grd    wall.grd    -o assembly.grd
shore cc-par-merge background.cc.par wall.cc.par -o assembly.cc.par --grd assembly.grd
```

Hole cutting, fringe identification, and donor search are the **solver's** responsibility — SHORE produces every grid and every descriptor the solver needs, then steps out of the way.

## Pipeline at a glance

```
Body-fitted (STL)                 Parametric backgrounds
  STL (raw triangles)               (no surface input)
        │                                  │
        ▼  shore mesh                      ▼  shore primitive
  cubed sphere (6 blocks)            box | annulus | flat-caps
  or O-grid (1 block)                (1 or 5 blocks)
        │                                  │
        ▼                                  ▼
  wall_*.geo                          background_*.geo
        │                                  │
        ▼  shore grd                       ▼  shore grd
  wall.grd                            background.grd
        │                                  │
        ▼  shore cc-par                    ▼  shore cc-par
  + shore boxes-from-stl              (no boxes for far-field)
  wall.cc.par                         background.cc.par
        │                                  │
        ▼  shore proc-input                ▼  shore proc-input
  wall.proc.input                     background.proc.input
        └─────────────┬────────────────────┘
                      ▼  shore grd-merge / cc-par-merge / proc-input-merge
            assembly.grd  +  assembly.cc.par  +  assembly.proc.input
                       (Xall multi-component job)
```


## Body-fitted topologies

The cubed-sphere topology is SHORE's flagship and the default for `shore mesh`: 4 equatorial sub-blocks plus 2 polar cap blocks, all conforming through SHARED-memory seams (sub-sub) and DIRICHLET pinning (cap-equator). Cap k=0 seeds use **gnomonic flat-square projection** — a flat diamond in the tangent plane at each cap pole, radially projected onto the body — so cap cells are square by construction with no concentric pole pattern.

```bash
shore mesh sphere.stl --ni 40 --nj 60 --nk 30
# → sphere_sub0.geo .. sphere_sub3.geo, sphere_cap_north.geo, sphere_cap_south.geo
```

Spacing continuity is enforced at two levels:

- **k=0 seam C1**: the equator's first/last meridional cell is per-j pinned to the cap-edge first cell at every seam vertex (so the 3 edges meeting at each cap-corner share a single first-step length).
- **Volume seam continuity**: each block computes its outward normal using centred differences across multi-block seams (ghost rows from partner blocks), so the spacing stays continuous on the outer surface layer-by-layer.

The single-block O-grid is the simpler alternative — one j-periodic block, polar singular axis rows, one `.geo` file. Suitable for downstream tools that need a single periodic-j block.

See [Cubed-sphere topology](topologies/cubed-sphere.md) and [O-grid topology](topologies/o-grid.md) for the full algorithm.

## Authors

**[Stefano Zaghi](https://github.com/szaghi)** · stefano.zaghi@cnr.it
> *Chief Yak Shaver, Accidental Research Scientist, and HPC Farmer* — CFD researcher who decided that one more day debugging Fortran build systems was one day too many, opened a Python REPL "just to prototype," and now finds himself maintaining a meshing library, a chimera assembler, an MPI load balancer, and the seven blog tabs he keeps meaning to read.

**Andrea Di Mascio** · andrea.dimascio@univaq.it
> *Sommo Vate & Resident Gandalf* — responsible for the deep CFD intuitions SHORE pretends to know.

**[Claude](https://claude.ai)** (Anthropic)
> *Omniscient Code Oracle & Tireless Rubber Duck* — AI pair programmer, responsible for writing the boring parts so humans don't have to.
