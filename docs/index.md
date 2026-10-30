---
layout: home

hero:
  name: SHORE
  text: Structured grids for overset CFD
  tagline: Hyperbolic near-body grids, collars and interactive multi-block blocking on STL and CAD bodies, packaged for the Xall/Xnavis Chimera solver, with CGNS and Plot3D exports.
  image:
    light: /figures/light/hero-wing-body.webp
    dark: /figures/dark/hero-wing-body.webp
    alt: Multi-block grid of a half wing–body, blocked and premeshed with SHORE
  actions:
    - theme: brand
      text: Get started
      link: /guide/
    - theme: alt
      text: Tutorials
      link: /tutorials/
    - theme: alt
      text: Cookbook
      link: /cookbook/

features:
  - title: Marched near-body grids
    details: Six topologies from a body surface — O-grid, cubed sphere, C-grid with tip caps, OH and CH hybrids out to a far-field box — extruded layer by layer with seam-aware normals and a Jacobian guard, by an algebraic or a Steger–Chaussee hyperbolic scheme.
    link: /concepts/near-body-pipeline
    linkText: How marching works
  - title: Interactive blocking
    details: Split a box into blocks, cut a hole, insert O-grids, associate faces, edges and vertices with the geometry, premesh with seams matched bit for bit, then wall-normal lines and a wall layer with an exact first cell.
    link: /concepts/blocking
    linkText: Blocking
  - title: STL and CAD geometry
    details: Triangulated STL with feature curves, or STEP / IGES / BREP through OpenCascade with exact projection. Entities are fingerprinted, so a journal replays on an edited model.
    link: /concepts/geometry
    linkText: Geometry
  - title: Collars at junctions
    details: Intersection curves of two bodies, hyperbolic surface marching on both, and the collar volume around the junction — closed loops, half models on a symmetry plane, and split collars at sharp corners.
    link: /concepts/collars
    linkText: Collars
  - title: Projects you can replay
    details: Every change is a command; the journal is a Python script. Undo, redo, save, reload and replay bit for bit — from the CLI, the Python API or the GUI, which records the same journal.
    link: /concepts/projects
    linkText: Projects
  - title: Solver-ready output
    details: The Xall bundle (.grd, cc.par with orientation codes, proc.input, marker boxes) for multi-component overset runs, plus CGNS with 1-to-1 connectivity, Plot3D and VTK.
    link: /tutorials/xall-bundle
    linkText: The Xall bundle
---

<div class="vp-doc" style="max-width: 1152px; margin: 64px auto 0; padding: 0 24px;">

## A first grid in three commands

::: code-group

```bash [CLI]
pipx install shore-mesh
shore body sphere -o sphere.stl
shore mesh sphere.stl --ni 33 --nj 64 --nk 25 --ds 0.01 --growth 1.15 \
    --theta-cap-deg 35 -o sphere --journal sphere.py
```

```python [Python]
import numpy as np
from shore.project import Project
from shore.shapes import sphere

sphere().export("sphere.stl")
p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=25, ni=33, nj=64, ds=0.01, growth=1.15,
                 theta_cap=np.deg2rad(35))
p.march()
p.write_geo("sphere")
p.save_journal("sphere.py")
```

:::

Six `.geo` blocks of a cubed-sphere O-grid, and a journal that rebuilds them
bit for bit with `shore run sphere.py`. The [quick start](/guide/quickstart)
goes on from here; the [tutorials](/tutorials/) take you to collars,
projects and a blocked wing–body.

## Where to go

<div class="shore-cards">
  <a class="shore-card" href="/shore/guide/"><strong>Get started</strong><span>What SHORE is, installation, the quick start and which workflow fits your body.</span></a>
  <a class="shore-card" href="/shore/tutorials/"><strong>Tutorials</strong><span>Eight end-to-end lessons, from a sphere to a blocked CAD wing–body in the GUI.</span></a>
  <a class="shore-card" href="/shore/cookbook/"><strong>Cookbook</strong><span>Short recipes for one task each: spacing, associations, exports, splitting.</span></a>
  <a class="shore-card" href="/shore/concepts/"><strong>Concepts</strong><span>The block model, geometry, marching, collars and blocking explained.</span></a>
  <a class="shore-card" href="/shore/reference/"><strong>Reference</strong><span>Every CLI command, the Python API, project scripts, file formats.</span></a>
  <a class="shore-card" href="/shore/adr/"><strong>Design</strong><span>51 architecture decision records and the roadmap: what was decided and why.</span></a>
</div>

</div>
