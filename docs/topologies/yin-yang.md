# Yin-Yang topology

The **Yin-Yang topology** resolves the polar singularity of a single
[O-grid](o-grid) without introducing cap blocks: instead of closing
the polar holes, it covers them with a *second* O-grid rotated 90°,
so each block's polar exclusion zone falls inside the other block's
covered region.

::: warning Status: Python-API only
The Yin-Yang topology is **not wired to the `shore mesh` CLI** in the
current release — the CLI builds the [caps topology](caps)
unconditionally. The Yin-Yang building blocks remain in the Python
API (`shore.surface.yinyang`) for users who want to construct the
topology manually or experiment with downstream Chimera-style
post-processing.
:::

## Concept

Two single-block O-grids:

- **Yin** — standard lat-lon directions with colatitude
  $\theta \in (\theta_m, \pi - \theta_m)$.
- **Yang** — same grid, rotated by `R_y(+π/2)`. Its pseudo-poles land
  on the equatorial x-axis of the Yin grid.

Each block leaves a polar exclusion of half-angle $\theta_m$ from its
own pole. With $\theta_m < \pi/4$ the two polar exclusions of one
block fall strictly inside the other block's covered region, so the
two-block union is a complete covering of the sphere.

## Overlap zone

Because the two blocks are independent O-grids, their interiors
overlap in a band determined by $\theta_m$:

| $\theta_m$ (rad) | Approx. degrees | Overlap fraction |
|---|---|---|
| π/6 | 30° | ~73% |
| π/5 | 36° | ~62% |
| 0.65 | ~37° | ~59% |
| π/4 − ε | →45° | →41% |

Larger overlap means more redundancy (more cells) but a more
comfortable Chimera interpolation margin. $\theta_m \to \pi/4$ is the
minimum for full coverage.

## Limitations vs. caps

| | Yin-Yang | Caps |
|---|---|---|
| Polar singularity | Avoided | Avoided |
| Block count | 2 | 6 |
| Block conformity | Independent grids (no automatic stitching) | Conforming SHARED + DIRICHLET seams |
| Suited for | Chimera-capable solvers that interpolate the overlap | Single-grid solvers; downstream Chimera between body and far-field |
| Cell count | Higher (overlap is duplicated) | Lower (no duplication) |

For SHORE's overset CFD context the **caps topology** is the better
default — the conforming 6-block layout matches the downstream
Chimera-between-bodies pattern without an extra Yin-Yang
interpolation step inside each body's near-field.

## Mathematical foundation

Reference: Kageyama, A. & Sato, T. (2004). "Yin-Yang grid: An
overset grid in spherical geometry." *Geochem. Geophys. Geosyst.*
**5**(9), Q09005. https://doi.org/10.1029/2004GC000734

## Python-API access

```python
import numpy as np
from shore.surface.io import load_stl
from shore.surface.remesh import project_sphere_to_surface
from shore.surface.yinyang import YANG_ROTATION, project_yang_to_surface

mesh = load_stl("sphere.stl")
yin  = project_sphere_to_surface(mesh, ni=40, nj=60, theta_cap_deg=30.0)
yang = project_yang_to_surface(mesh,  ni=40, nj=60, theta_cap_deg=30.0)
```

`YANG_ROTATION` is the 3×3 matrix that takes a Yin direction to its
Yang counterpart. `project_yang_to_surface` is a convenience wrapper
that applies the rotation and casts the rays.

## See also

- [O-grid topology](o-grid) — the single-block primitive Yin-Yang
  duplicates
- [Caps topology](caps) — SHORE's production alternative
- [Chimera context](/concepts/chimera) — overset mesh generation
