# Cubed-sphere topology (6-block)

The **cubed-sphere topology** is SHORE's production topology and the
one the `shore mesh` CLI builds. It is the standard cubed-sphere
construction from atmospheric and planetary CFD (Sadourny 1972;
Ronchi, Iacono & Paolucci 1996): the body surface is decomposed into
6 conforming hexahedral blocks — 4 equatorial sub-blocks plus 2 polar
caps — corresponding to the 6 faces of a cube projected outward onto
the body. There is no singular axis, all seams are conforming, and
each cap is a regular structured block in its own right.

It is **not** a wrapped O-grid: the cubed sphere does not extend the
single-block [O-grid](o-grid) by augmenting it with caps. They are
two distinct constructions that share only the underlying spacing
primitives.

![Cubed-sphere topology — perspective view](/figures/cubed-sphere-topology-blocks.png)

This page describes the topology — the block layout, the multi-block
connectivity, and the boundary conditions. The algorithmic details
(how each block is built and marched) live in the
[Algorithm](/algorithm/) section.

## Block decomposition

Six structured blocks:

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

with `nc = nj // 4 + 1`. Each sub-block spans one 90° azimuthal
sector; `nj` must be a multiple of 4 (rounded down with a warning if
not).

| Block | Shape | Role |
|---|---|---|
| `sub0` | `(ni, nc, nk, 3)` | Equatorial sector, $\phi \in [0°, 90°]$ |
| `sub1` | `(ni, nc, nk, 3)` | Equatorial sector, $\phi \in [90°, 180°]$ |
| `sub2` | `(ni, nc, nk, 3)` | Equatorial sector, $\phi \in [180°, 270°]$ |
| `sub3` | `(ni, nc, nk, 3)` | Equatorial sector, $\phi \in [270°, 360°]$ |
| `cap_north` | `(nc, nc, nk, 3)` | North polar cap (`z < 0` of body) |
| `cap_south` | `(nc, nc, nk, 3)` | South polar cap (`z > 0` of body) |

## Top-down view

The cap is a square-celled patch that meets the four equatorial sub-
blocks at the colatitude-`theta_cap` parallel. The four cap corners
sit on that parallel; cap edges dip slightly poleward at midpoints
(intrinsic to the gnomonic-square geometry):

![Top-down view](/figures/cubed-sphere-topology-topdown.png)

## Multi-block connectivity

### Sub-sub seams (j-direction)

The four equatorial sub-blocks share their `j_lo` / `j_hi` columns
with their azimuthal neighbours. The four sub-blocks live in a single
contiguous **ring** allocation `(ni, 4·(nc-1)+1, nk, 3)`; each
sub-block is a numpy view into a 90° slice. Adjacent sub-blocks share
their boundary j-column by **memory aliasing** (no reconciliation
needed):

| Seam | Mechanism | Boundary condition |
|---|---|---|
| sub0.j_hi ↔ sub1.j_lo | Memory alias (numpy view) | `FaceBC.SHARED` |
| sub1.j_hi ↔ sub2.j_lo | Memory alias | `FaceBC.SHARED` |
| sub2.j_hi ↔ sub3.j_lo | Memory alias | `FaceBC.SHARED` |
| sub3.j_hi ↔ sub0.j_lo | Memory alias (periodic wrap) | `FaceBC.SHARED` |

### Cap-equator seams (i-direction)

The cap's four boundary edges are **DIRICHLET-bound** to the four
equatorial sub-block pole rows:

**North cap (z<0 of body):**

| Cap edge | ← Equator pole row |
|---|---|
| `cap_north[0, :]` | `sub0.nodes[ni-1, :]` |
| `cap_north[:, -1]` | `sub1.nodes[ni-1, :]` |
| `cap_north[-1, :]` | `sub2.nodes[ni-1, ::-1]` |
| `cap_north[:, 0]` | `sub3.nodes[ni-1, ::-1]` |

**South cap (z>0 of body):**

| Cap edge | ← Equator pole row |
|---|---|
| `cap_south[:, 0]` | `sub0.nodes[0, :]` |
| `cap_south[-1, :]` | `sub1.nodes[0, :]` |
| `cap_south[:, -1]` | `sub2.nodes[0, ::-1]` |
| `cap_south[0, :]` | `sub3.nodes[0, ::-1]` |

The `_overwrite_cap_boundaries` helper writes the sub-block pole rows
into the cap edges at every layer of marching, keeping cap-equator
conformity bit-exact regardless of what the cap's interior march would
have produced.

## Construction order

The `CubedSphereTopology.build` orchestrator:

1. Builds **both cap k=0 layers** by gnomonic flat-square projection
   (see [k=0 mesh construction](/algorithm/k0-construction)). The cap
   edges are now known.
2. Builds the **equator k=0** by extracting the cap seam rings and
   slerp-interpolating great-circle meridians from south seam to north
   seam. Each meridian's first/last interval is per-j pinned to the
   cap-edge first cell at the corresponding seam vertex (C1 spacing
   match).
3. Allocates the contiguous equator ring + cap arrays and wires the
   `HexBlock`s with `FaceBC.SHARED` (sub-sub) and `FaceBC.DIRICHLET`
   (cap-equator) face BCs.

`Mesh.march` then advances all 6 blocks simultaneously layer by layer
with [seam-aware normals](/algorithm/seam-aware-normals) so spacing
continuity holds across every multi-block seam at every layer.

## Multi-block edge spacing rule

At each cap-equator seam vertex, three multi-block edges meet:
two cap-seam edges and one equator meridional edge. The
[per-j C1 seam pin](/algorithm/k0-construction#per-j-c1-seam-pinning)
forces all three first-cells to share one length:

![3-edge cap-corner junction](/figures/junction-3edge.png)

The pin is active when the cap radial step and the equator's natural
radial step are within ~2× of each other. Beyond that, the equator
falls back to its unpinned `i_spacing` law and emits one
`UserWarning` per build. See
[k=0 mesh construction](/algorithm/k0-construction) for the rule of
thumb on choosing `(ni, nj, theta_cap)`.

## Volume seam continuity

Without seam-aware normals, each block computes its boundary normal
from a one-sided difference internal to itself. The two sides of a
seam disagree, and the disagreement accumulates over many layers as a
visible spacing jump on the outer surface. With seam-aware normals
(centred differences against ghost rows from the partner block) both
sides agree on the normal at every shared node and the spacing stays
continuous all the way to the outer layer:

![Volume seam continuity (k=0 vs outer layer)](/figures/volume-seam-continuity.png)

The blue (inner) and orange (outer) layers both show smooth grids
across the block seams; cell-size variation on the outer surface
follows the natural growth law, not the block decomposition. See
[Seam-aware normals](/algorithm/seam-aware-normals) for the
implementation and quantitative verification.

## Cap geometry knobs

| Parameter | Effect |
|---|---|
| `theta_cap` | Polar exclusion angle in radians. The 4 cap corners sit at colatitude `theta_cap` on a unit sphere; cap edges dip slightly poleward at midpoints. Default $\pi/6$ (30°). |
| `j_spacing` | Cap-edge tangent-plane node distribution (`'uniform'` or `'tanh2'`). The equator inherits the same distribution at the seam by DIRICHLET binding. |
| `j_beta` | Clustering strength for `tanh2` j-spacing. |

For the equator's meridional distribution (`i_spacing`, `i_beta`) and
the wall-normal stepping (`spacing_k`), see
[Spacing laws](/algorithm/spacing-laws).

## Public API

```python
import numpy as np
from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

m = Mesh.from_stl(
    "sphere.stl",
    topology=CubedSphereTopology(),
    ni=40, nj=60, nk=20,
    theta_cap=float(np.pi / 6),
    spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
)
m.march(ds=0.02, growth=1.1)
m.write_geo("sphere")  # 6 sphere_<label>.geo files
```

The 6 `.geo` files can be bundled into one ParaView MultiBlock with
`shore.io.vtk.write_caps_vtm`:

```python
from shore.io.geo import read_geo
from shore.io.vtk import write_caps_vtm

labels = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]
blocks = [read_geo(f"sphere_{lbl}.geo") for lbl in labels]
write_caps_vtm("sphere.vtm", blocks, with_jacobian=True)
```

## See also

- [Surface projection](/algorithm/projection) — the ray-cast kernel
- [k=0 mesh construction](/algorithm/k0-construction) — gnomonic cap
  + great-circle equator + per-j C1 seam pin
- [Hyperbolic marching](/algorithm/hyperbolic-marching) — outward
  extrusion
- [Seam-aware normals](/algorithm/seam-aware-normals) — what keeps
  spacing continuous across multi-block seams during marching
- [O-grid topology](o-grid) — the single-block sibling topology
  (structurally independent of the cubed sphere)
- [`shore.surface.cap`](/reference/api-sphere-cap),
  [`shore.surface.equator`](/reference/api-sphere-cap),
  [`shore.topology.sphere_cap`](/reference/api-sphere-cap) — Python API
