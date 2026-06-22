# O-grid topology

An O-grid (also called *O-type grid* or *O-mesh*) is a structured grid
topology that wraps completely around a body, forming a closed loop in
one parametric direction. The name comes from the shape of the grid
lines around a cross-section, which resemble the letter O.

The single-block O-grid is a **first-class topology** in SHORE: it
can be built and marched directly through the CLI with
`shore mesh --topology ogrid` (or `topology = "ogrid"` in a TOML
config), producing a single ``.geo`` file. It coexists with the
default [cubed-sphere topology](cubed-sphere) as a peer — the two
topologies are structurally independent and pick different
trade-offs:

| | O-grid | Cubed sphere |
|---|---|---|
| Blocks | 1 | 6 (4 sub + 2 caps) |
| Polar coverage | singular axis at i=0, i=ni-1 | non-singular gnomonic caps |
| j direction | periodic | conforming SHARED seams |
| Output | one `<stem>.geo` | six `<stem>_<label>.geo` |
| Best for | Chimera background, downstream tools that expect a periodic-j block | production near-body meshing |

## Building from the CLI

```bash
# Default invocation: uniform spacing, theta_cap=9 deg.
shore mesh sphere.stl --topology ogrid -o sphere

# With a TOML config (`[grid] topology = "ogrid"`):
shore mesh -c ogrid.toml
```

The output is one `<stem>.geo` file with shape `(ni, nj, nk, 3)` and
j-periodic. All ``shore mesh`` flags apply equally — `--ni`, `--nj`,
`--nk`, `--ds`, `--growth`, `--smooth`, `--smooth-iters`, the spacing
laws, `--blend-normals-k`, `--theta-cap` / `--theta-cap-deg`. The
legacy cap-shaping flags (`--cap-pull`, `--cap-squareness`) are
deprecated for both topologies and silently ignored (the cap k=0 is now
built by gnomonic flat-square projection).

The Python API is the symmetric of the cubed sphere's:

```python
from shore.mesh import Mesh, OGridTopology, Spacing1D

m = Mesh.from_stl(
    "sphere.stl",
    topology=OGridTopology(),
    ni=40, nj=60, nk=30,
    theta_cap_deg=9.0,
    spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
)
m.march(ds=0.02, growth=1.1)
m.write_geo("sphere")  # → sphere_ogrid.geo (with the auto-suffix writer)
```

## Topology

In SHORE's single-block model:

- The **i-direction** spans from one pole to the other (latitude, open, not periodic)
- The **j-direction** wraps completely around the body (longitude, periodic: `j=0` is adjacent to `j=nj-1`)
- The **k-direction** goes wall-normal from body surface (k=0) to far field (k=nk-1)

In the cubed-sphere topology the j-direction is split into 4 sub-blocks (each
spanning 90°) plus 2 cap blocks; the periodic j-wrap of the single-block
view is realised as a SHARED-memory sub3↔sub0 seam.

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

SHORE implements the O-grid (this page) and the cubed sphere as
peers, plus the [C-grid](c-grid) (`--topology cgrid`) for sharp-/blunt-TE
bodies and the H-grid analytic primitives (`shore primitive h-grid-rect`
/ `h-grid-tube`) for body-less ducts.

## Per-edge spacing

Per-edge spacing on the O-grid topology has one viable axis: **k**. The
i and j directions are foreclosed by the geometry:

- **i**: the j direction is periodic, so the four parallel i-edges
  (`j_lo_k_lo`, `j_hi_k_lo`, `j_lo_k_hi`, `j_hi_k_hi`) lie at adjacent
  azimuths under the j-wrap rather than at distinct transverse
  corners. There is no degree of freedom over the per-axis
  `i_spacing` argument.
- **j**: periodic — the "first" and "last" cells of any j-edge are
  the same cell, so endpoint pinning is geometrically meaningless.
- **k**: the four parallel k-edges (`i_lo_j_lo`, `i_hi_j_lo`,
  `i_lo_j_hi`, `i_hi_j_hi`) of the single block can be overridden via
  the same identity-check resolver as the cubed sphere. Replace the
  `spacing` field on all four with a shared `Spacing1D` instance and
  `Mesh.march()` reads it instead of `(ds, growth)`:

```python
from shore.mesh.spacing import Spacing1D

override = Spacing1D(law="geometric", ds=5e-4, growth=1.2)
for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
    m.blocks[0].edge(name).spacing = override

m.march(ds=99, growth=99)   # ds, growth ignored — override wins
```

See [Per-edge spacing](/algorithm/spacing-laws#per-edge-spacing) for
the full mechanism, and `examples/03-per-edge-spacing/per_edge_k_spacing_ogrid.py` for a
side-by-side comparison.

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

## See also

- [Cubed-sphere topology](cubed-sphere) — the production 6-block
  alternative without polar singularities
- [Surface projection](/algorithm/projection) — how the body STL is
  ray-cast onto the structured i × j sphere template
- [Hyperbolic marching](/algorithm/hyperbolic-marching) — how the
  k-layers are extruded
