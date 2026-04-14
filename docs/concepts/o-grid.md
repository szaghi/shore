# O-grid topology

An O-grid (also called *O-type grid* or *O-mesh*) is a structured grid topology that wraps completely around a body, forming a closed loop in one parametric direction. The name comes from the shape of the grid lines around a cross-section, which resemble the letter O.

## Topology

In SHORE's implementation:

- The **i-direction** spans from one pole to the other (latitude, open, not periodic)
- The **j-direction** wraps completely around the body (longitude, periodic: `j=0` is adjacent to `j=nj-1`)
- The **k-direction** goes wall-normal from body surface (k=0) to far field (k=nk-1)

```
          k=nk-1 (far field)
         ╔══════════════════╗
         ║  O-grid block    ║
         ╠══════════════════╣  k=2
         ╠══════════════════╣  k=1
         ╠══════════════════╣  k=0 (wall)
         ╚══════════════════╝
              body surface
```

In 3D, the periodic j-direction forms a closed shell around the body:

```
i=0 (north pole)
     │
     │  ← j loops around (0 … nj-1, periodic)
     │
i=ni-1 (south pole)
```

## Grid dimensions

| Index | Direction | Periodicity | Typical range |
|---|---|---|---|
| `i` | Latitude (pole to pole) | Not periodic | 20 – 80 |
| `j` | Longitude (around body) | Periodic: `j=0 ↔ j=nj-1` | 30 – 120 |
| `k` | Wall normal | Not periodic | 20 – 60 |

The total number of grid points is `ni × nj × nk`. For typical values (40 × 60 × 30) this is 72,000 points — trivial to hold in memory.

## Near-wall resolution

The most important grid parameter for CFD is the first-layer thickness `ds`. For wall-resolved Large Eddy Simulation (LES) or Direct Numerical Simulation (DNS), the requirement is:

$$\Delta y^+ = \frac{\Delta y \, u_\tau}{\nu} \approx 1$$

where $u_\tau$ is the friction velocity and $\nu$ is the kinematic viscosity. For RANS simulations with wall functions, $\Delta y^+ \approx 30$ – 100 is acceptable.

Given a target $\Delta y^+$ and estimated flow conditions, the required first-layer thickness is:

$$\Delta s = \frac{\Delta y^+ \, \nu}{u_\tau}$$

## O-grid vs. C-grid vs. H-grid

| Topology | Shape | When to use |
|---|---|---|
| **O-grid** | Closed loop around body | Bodies with no sharp trailing edges (cylinders, fuselages, nacelles) |
| **C-grid** | Open C-shape, wake cut | Airfoils and wings with sharp trailing edges |
| **H-grid** | Rectangular with body cut-out | Bluff bodies, internal flows |

SHORE currently implements the O-grid only. C-grid and H-grid topologies require different surface parametrisation strategies and are on the roadmap.

## Periodicity in j

Because the j-direction is periodic, normal estimation and Laplacian smoothing use circular indexing:

```python
jp1 = (np.arange(nj) + 1) % nj   # j+1, wrapping nj-1 → 0
jm1 = (np.arange(nj) - 1) % nj   # j-1, wrapping 0 → nj-1
```

This is applied consistently in `_compute_normals`, `_laplacian_smooth`, and `_min_jacobian`. The i-direction uses one-sided differences at boundaries (no wrapping).

## Visualisation

The structured grid can be visualised with pyvista (optional dependency):

```python
import pyvista as pv
import numpy as np
from shore.io.geo import read_geo

grid = read_geo("sphere.geo")
ni, nj, nk, _ = grid.shape

# Convert to pyvista StructuredGrid
points = grid.reshape(-1, 3)
sg = pv.StructuredGrid()
sg.dimensions = (ni, nj, nk)
sg.points = points

sg.plot(show_edges=True, opacity=0.5)
```
