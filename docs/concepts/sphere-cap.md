# Sphere polar cap topology

A standard lat-lon O-grid wraps the equatorial band of a sphere but leaves
open holes at both poles. The **sphere-cap topology** closes those holes with
two hexahedral cap blocks, one per pole, producing a watertight 6-block
structured grid with no singular axes.

## Motivation

In a conventional lat-lon O-grid the latitude bounds are chosen to avoid the
poles (typically `theta ∈ [0.05π, 0.95π]`). This leaves two open edges at
`i = 0` (south) and `i = ni - 1` (north) that cannot be extruded into a
conforming volume — the normal field degenerates as the ring of grid points
collapses to a single pole node.

Two standard workarounds exist:

- **Pole singularity** — the two pole edges are collapsed to a single node.
  This produces degenerate cells adjacent to the pole and is numerically
  problematic for most CFD solvers.
- **Yin-Yang** — two overlapping curvilinear O-grids shifted 90° avoid the
  poles entirely. Chimera interpolation handles the overlap. This requires a
  solver with Chimera capability and adds complexity.

The **sphere-cap topology** takes a third route: the equatorial O-grid is
split into four 90° sub-blocks along the azimuthal direction, and each polar
hole is closed by a dedicated cap block whose seed face is determined by the
four shared edges that meet at the pole.

## Block decomposition

The result is six structured blocks:

```
              cap_north (nc × nc × nk)
              ┌─────────────┐
              │             │
  sub3  ┌─────┤     N       ├─────┐  sub1
(ni×nc) │     │             │     │ (ni×nc)
        │     └─────────────┘     │
        │       sub0 (ni×nc)      │
        └──────────────────────── ┘
              ┌─────────────┐
              │     S       │
              │             │
              └─────────────┘
              cap_south (nc × nc × nk)
```

where `nc = nj // 4 + 1`. Each sub-block spans one 90° azimuthal sector.
Adjacent sub-blocks share their boundary j-column exactly (node-for-node
conformity). Each cap shares its four boundary edges with the pole row of the
four adjacent sub-blocks.

## Seed face construction

The k=0 face of each cap block is built by discrete Laplacian interpolation:

1. The four boundary edges are set from the equatorial sub-block pole rows
   (`i = ni-1` for the north cap, `i = 0` for the south cap).
2. The interior nodes are initialised to zero and solved to convergence with
   Dirichlet boundary conditions using the discrete Laplacian
   (strength = 0.5, tolerance = $10^{-12}$).

This produces the minimum-energy surface that interpolates the four boundary
edges while keeping all boundary nodes exactly coincident with the equatorial
blocks.

::: info Why not gnomonic projection?
A gnomonic (great-circle) projection of the four corner points onto the sphere
is the natural choice for a spherical cap, but the winding correction required
to make normals point outward changes the interior orientation relative to
what the boundary edge mapping expects. The Laplacian approach bypasses this
by building the seed face directly from the boundary constraints.
:::

## Simultaneous layer-by-layer marching

All six blocks are marched simultaneously, one layer at a time:

```
for k in 0 .. nk-2:
    1. Advance 4 equatorial sub-blocks (normal march + Laplacian smooth)
    2. Re-stitch shared edges between adjacent sub-blocks (average)
    3. Jacobian check on sub-blocks
    4. Extract north/south pole rows from sub-blocks (the cap boundaries at k+1)
    5. Advance 2 cap blocks:
       a. Normal-based displacement (same step and normals as equatorial blocks)
       b. Laplacian smoothing (same strength and iterations)
       c. Overwrite boundary edges with equatorial pole rows (Dirichlet pinning)
    6. Jacobian check on cap interior cells
    7. Commit all 6 blocks to layer k+1
```

The Dirichlet boundary pinning in step 5c is what enforces conformity at every
layer: the cap boundary is forced to match the equatorial sub-block pole rows
exactly, regardless of what the cap's interior normal march would have produced.

## Jacobian convention

Equatorial sub-blocks and cap blocks use `_jacobian_field_open` (non-periodic
j) because all six blocks have open boundaries. The four **corner cells** at
each cap's `(i=0,j=0)`, `(i=0,j=nc-2)`, `(i=nc-2,j=0)`, `(i=nc-2,j=nc-2)`
are excluded from the Jacobian check: these cells border the junction of four
sub-blocks at a single shared node, and their signed volume is near-zero by
construction — not a mesh-quality defect, but a topological consequence of
four blocks meeting at one point.

## Controlling cap size with `theta_cap`

`theta_cap` (in radians, or `theta_cap_deg` in degrees) is the polar exclusion
angle passed to `project_sphere_to_surface` when building the k=0 surface
layer. It determines how far the equatorial band extends toward the poles, and
therefore the physical size of the cap faces.

The equatorial band spans colatitudes `[theta_cap, pi - theta_cap]`. The cap
boundary parallel has radius `sin(theta_cap)` on the unit sphere:

| `theta_cap` | degrees | Cap boundary radius / sphere radius |
|---|---|---|
| 0.087 | 5° | 0.087 |
| 0.157 | 9° (default) | 0.156 |
| 0.262 | 15° | 0.259 |
| 0.349 | 20° | 0.342 |
| 0.524 | 30° | 0.500 |

Larger `theta_cap` gives a physically larger cap face and a shorter equatorial
band. The grid-space cap size `nc = nj // 4 + 1` is independent of `theta_cap`
— it is set by `nj`. `theta_cap` controls only the physical arc length of the
cap boundary.

```python
# Default ~9° cap
surface = sphere_layer(40, 60)                              # theta_cap = 0.05*pi

# 20° cap — larger polar region
surface = sphere_layer(40, 60, theta_cap_deg=20.0)

# Explicit radians
surface = sphere_layer(40, 60, theta_cap=0.349)
```

Or via `project_sphere_to_surface` for a real STL body:

```python
k0 = project_sphere_to_surface(mesh, ni=40, nj=60, theta_cap_deg=20.0)
```

## Cap quality parameters {#quality}

Two parameters control the mesh quality of the polar cap faces.  They are
independent and can be combined freely.

### `cap_squareness` — straighten boundary edges {#cap_squareness}

Without any shaping the 4 cap boundary edges are arcs of the polar parallel
circle, each spanning 90° of longitude.  Two adjacent arcs meet at a corner
with an angle of roughly 157°.  For CFD solvers this wide angle produces
poor aspect-ratio cells near the corners.

The `cap_squareness` parameter blends each boundary edge from the circular arc
toward the **geodesic chord** (great-circle arc) connecting its two corner
nodes — the edge of the square inscribed in the parallel circle.  Corner nodes
are fixed (they are shared between adjacent equatorial sub-blocks and cannot
move); only the intermediate nodes along each edge are displaced.

The blend uses spherical linear interpolation (slerp) between the two corner
unit vectors, ensuring the displaced nodes remain on the sphere surface at
every layer.  The squareness is re-applied at every layer so the cap boundary
stays in the inscribed-square geometry as the grid grows outward.

**Effect on corner angle:**

| `cap_squareness` | Corner angle |
|---|---|
| 0.0 (default) | ~157° circular arc |
| 0.25 | ~140° |
| 0.50 | ~120° |
| 0.75 | ~100° |
| 1.0 | ~90° inscribed square |

The full range `[0, 1]` is numerically stable — slerp is geometrically
monotone and does not create concavities.

### `cap_pull` — dome the cap interior {#cap_pull}

Without `cap_pull` the cap k=0 seed face is the discrete Laplacian solution
interpolating the four boundary arcs.  This produces wide, fan-shaped cells
near the cap centre.

The `cap_pull` parameter domes the interior toward the geographic pole by
blending each interior node from its Laplacian position toward the pole point
and re-projecting onto the sphere.  The 4 boundary edges are kept exactly
unchanged, so the cap remains fully conforming with the equatorial sub-blocks.

**Weight function** — a product of two sine profiles, symmetric around the
face centre:

```
w(i, j) = cap_pull × sin(π × i/(nc−1)) × sin(π × j/(nc−1))
```

Zero on all 4 boundary edges, peak at the cap centre.

**Effect on interior cell quality:**

| `cap_pull` | Interior cell quality |
|---|---|
| 0.0 (default) | flat Laplacian face |
| 0.3 | noticeably improved |
| 0.5 | significantly improved |
| ≥ 0.7 | may fail with aggressive `ds`/`growth` |

**Practical recommendation:** use `cap_pull = 0.3`–`0.5`.  Values above ~0.6
may trigger Jacobian failures at early layers as the dome transition becomes
too steep for the normal march to follow.

### Combined usage

```python
blocks = split_ogrid_with_polar_caps(
    surface, nk=10, ds=0.02, growth=1.10,
    cap_squareness=0.5,   # straighten boundary edges toward inscribed square
    cap_pull=0.4,         # dome interior toward pole
)

# Equalise cell widths at all block interfaces after marching
from shore.topology.sphere_cap import correct_interface_spacing
correct_interface_spacing(blocks)
```

Or via the example script:

```bash
.venv/bin/python examples/visualize_geo.py \
    --generate --topology caps --vtk --no-view --nk 10 \
    --cap-squareness 0.5 --cap-pull 0.4 --correct-spacing
```

## Grid parameters

The cap block size `nc = nj // 4 + 1` is determined by `nj`. If `nj` is not
divisible by 4, it is rounded down to the nearest multiple of 4 with a
`warnings.warn`. All six blocks use the same `nk`, `ds`, and `growth` as the
underlying single-block O-grid.

**Practical limit:** because the equatorial blocks have open i-boundaries and
the cap blocks have open ij-boundaries, the cell quality degrades faster than
in a periodic O-grid as the grid expands outward. For a unit sphere, `nk ≤ 10`
with `ds = 0.02` and `growth = 1.10` gives a total thickness of ~0.32, which
is well-behaved. Pushing beyond ~0.6× the sphere radius will invert equatorial
cells near the poles.

## API

The topology is implemented in `shore.topology.sphere_cap`. The main entry
point is:

```python
from shore.topology.sphere_cap import split_ogrid_with_polar_caps

blocks = split_ogrid_with_polar_caps(
    surface, nk=10, ds=0.02, growth=1.10, cap_pull=0.4  # optional: dome interior
)
# blocks[0..3]: equatorial sub-blocks, shape (ni, nc, nk, 3)
# blocks[4]:   cap_north, shape (nc, nc, nk, 3)
# blocks[5]:   cap_south, shape (nc, nc, nk, 3)
```

See [`shore.topology.sphere_cap`](/reference/api-sphere-cap) for the full API
reference.

## Example

```python
import numpy as np
from tests.conftest import sphere_layer
from shore.topology.sphere_cap import split_ogrid_with_polar_caps
from shore.io.vtk import write_caps_vtm

# Build a lat-lon surface layer with default cap size (~9°)
surface = sphere_layer(40, 60)

# Or with an explicit 20° cap
surface = sphere_layer(40, 60, theta_cap_deg=20.0)

# Generate 6-block topology (with optional interior dome for quality)
blocks = split_ogrid_with_polar_caps(surface, nk=10, ds=0.02, growth=1.10, cap_pull=0.4)

# Export as a ParaView MultiBlock (.vtm) with Jacobian cell data
write_caps_vtm("sphere_caps.vtm", blocks, with_jacobian=True)
```

Or via the example script:

```bash
# Default cap (~9°), no dome
.venv/bin/python examples/visualize_geo.py \
    --generate --topology caps --vtk --no-view --nk 10

# Larger cap (20°) with interior dome and spacing correction
.venv/bin/python examples/visualize_geo.py \
    --generate --topology caps --vtk --no-view --nk 10 \
    --theta-cap-deg 20 --cap-pull 0.4 --correct-spacing

paraview sphere_caps.vtm
```

## See also

- [`shore.topology.sphere_cap`](/reference/api-sphere-cap) — Python API
- [Hyperbolic marching](/concepts/hyperbolic-marching) — underlying algorithm
- [Chimera mesh generation](/concepts/chimera) — where this topology fits
