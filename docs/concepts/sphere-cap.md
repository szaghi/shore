# Sphere polar cap topology

A standard lat-lon O-grid wraps the equatorial band of a sphere but leaves
open holes at both poles. The **caps topology** closes those holes with
two hexahedral cap blocks, one per pole, producing a watertight 6-block
structured grid with no singular axes.

## Motivation

In a conventional lat-lon O-grid the latitude bounds are chosen to avoid the
poles (typically `theta ∈ [theta_cap, π − theta_cap]`). This leaves two open
edges at `i = 0` (south) and `i = ni − 1` (north) that cannot be extruded
into a conforming volume — the normal field degenerates as the ring of grid
points collapses to a single pole node.

The **caps topology** splits the equatorial band into four 90° sub-blocks
along the azimuthal direction and closes each polar hole with a dedicated
cap block whose seed face is built by **gnomonic projection of a flat
tangent square**. This is the same construction used by classic cubed-sphere
schemes (e.g. for atmospheric-flow simulations on the sphere) and produces
square-celled cap blocks with no concentric pole singularity.

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

where `nc = nj // 4 + 1`. Each sub-block spans one 90° azimuthal sector;
`nj` must be a multiple of 4 (rounded down with a warning if not). Adjacent
sub-blocks share their boundary j-column exactly (node-for-node conformity
via numpy views into a contiguous ring allocation). Each cap shares its
four boundary edges with the pole row of the four adjacent sub-blocks.

## Cap k=0 seed: gnomonic flat-square projection

The cap is constructed by placing a flat **diamond** (a square rotated 45°)
in 3D, tangent to the sphere at the cap pole, with corners on the cardinal
axes at distance `tan(theta_cap)`. The 4 corners are the cap-equator
junction points at azimuths `{0°, 90°, 180°, 270°}` — aligned with the 4
sub-block junctions by construction. Interior nodes of the diamond are
filled by bilinear interpolation, then every node is radially projected
onto the body surface through the cached BVH ray-cast.

```
         tangent plane at pole
       (z = 1 in pole-local frame)

              C1 (0, h, 1)
                ╱╲
               ╱  ╲
              ╱    ╲
   C2 ◄─────╳──────► C0   (corners on cardinal axes,
   (-h,0,1) ╲      ╱  (h,0,1)   distance h = tan(theta_cap))
              ╲    ╱
               ╲  ╱
                ╲╱
              C3 (0, -h, 1)
```

After radial projection (`p / |p|` then scale to body radius), the 4 corners
land exactly at colatitude `theta_cap` on a unit sphere. Cap edges are
great-circle arcs between adjacent corners. Cap cells are square-shaped in
the bulk (modest gnomonic-square anisotropy at corners, ~√2). There is no
concentric ring pattern around the pole — the cap is just a smooth
square-topology patch on the body.

### Why a diamond, not an axis-aligned square

The 4 sub-block junctions on the seam parallel sit at azimuths
`{0°, 90°, 180°, 270°}`. An axis-aligned square placed at the pole would
have corners at azimuths `{±45°, ±135°}` — a 45° offset from the
junctions. Rotating the square 45° (turning it into a diamond) places its
corners exactly on the sub-block junctions, and the bilinear interpolation
of the 4 corners
$$
P(a,b) = (1-a)(1-b)\,C_0 + a(1-b)\,C_1 + ab\,C_2 + (1-a)b\,C_3
$$
yields the closed-form expression
$$
x_{\text{local}} = h(1 - a - b),\qquad y_{\text{local}} = h(a - b)
$$
with `h = tan(theta_cap)` and `(a, b) ∈ [0, 1]²`.

### `j_spacing` inheritance

The tangent-plane parametrisation uses 1D laws (`uniform` or `tanh2`) on
`[−h, +h]`, so the cap-edge node distribution along the seam matches the
equator's `j_spacing`. The equator's pole rows are filled bit-exactly from
the cap edges (DIRICHLET binding), so the equator inherits the cap's
distribution at the seam without any reconciliation step.

## Equator: great-circle meridians

The equator's k=0 layer is built **after** the cap (so the cap edges are
known). Each j-column is a great-circle arc from `seam_south[j]` to
`seam_north[j]`, parametrised in colatitude by the user's `i_spacing` law
(`uniform`, `tanh`, `tanh2`). Interior nodes are slerp-interpolated and
ray-cast onto the body to handle non-spherical surfaces.

The first/last interval on each meridian is **pinned** to the cap-edge
first cell at the corresponding seam vertex, giving C1 spacing continuity
across the seam (see next section).

## Multi-block edge spacing match (C1 at seam vertices)

At each cap-equator seam vertex three multi-block edges meet:

1. The cap-seam edge running along one side of the cap (e.g. `cap[:, 0]`).
2. The cap-seam edge running along the perpendicular side (e.g. `cap[0, :]`).
3. The equatorial-meridian edge running radially inward (e.g. sub_q's
   `[:, j_local]` column).

The contract: **all three edges share a single first-cell length at the
shared vertex**. This is the standard rule for clean multi-block
edge-spacing in a structured grid.

The cap is built first, fixing edges 1 and 2 by gnomonic-square geometry
(both equal to the great-circle arc from seam[j] to seam[j+1]). The
equator's per-meridian first step is then pinned to that arc length:

```python
first_step[j] = arccos(u_seam[j] · u_seam[j+1])  *  interval_j / omega_j
```

where `u_seam` are unit directions from anchor through the seam nodes,
`omega_j` is the meridional arc angle, and `interval_j` is the colatitude
span used by the meridian's spacing law. The factor
`interval_j / omega_j` converts between the two angular parametrisations.

The interior `ni − 2` cells of each meridian are then graded by `i_spacing`
with these endpoint steps fixed.

### Compatibility window and fallback

Pinning the first/last step is only stable when the cap's seam step and
the equator's natural first step (on the chosen `i_spacing` law) are
within a factor of ~2 of each other. Beyond that, pinning produces a
sharp discontinuity in cell size between the first and second cells of
each meridian that destabilises hyperbolic marching.

When the steps differ by more than 2×, the equator falls back to its
unpinned `i_spacing` law and a single `UserWarning` is emitted with
guidance on choosing compatible parameters.

**Rule of thumb** (uniform i_spacing): for the per-j pinning to be active
across the whole seam, choose `ni`, `nj`, `theta_cap` so that

$$
\frac{\theta_{\text{cap}}}{n_c - 1} \approx \frac{\pi - 2\theta_{\text{cap}}}{n_i - 1}
$$

with `n_c = n_j // 4 + 1`. For example, `theta_cap = π/6 (30°)` with
`n_i = 40, n_j = 60` gives cap step ~2° and equator natural step ~3.1°
(ratio 1.55× — within the 2× window).

## Hyperbolic marching with seam-aware normals

All six blocks are marched simultaneously, one layer at a time. Without
seam-aware normal computation, each block would compute its boundary
tangent from a one-sided difference internal to itself — sub_q's tangent
at its `j_hi` seam uses sub_q's `layer[:, -1] - layer[:, -2]`, while
sub_{q+1}'s tangent at its `j_lo` seam (the same seam) uses sub_{q+1}'s
`layer[:, 1] - layer[:, 0]`. These disagree, the two blocks compute
different normals at the same shared node, and the disagreement
accumulates over many marching layers as a visible spacing jump on the
outer surface.

**Fix.** Each block receives a `ghost_rows` dict containing the partner
block's "row one inward" data at the current layer, mapped through the
inverse of `_overwrite_cap_boundaries` for cap-equator seams (with
reversal for sub2/sub3). The boundary tangent then uses a centred
difference against the ghost row, exactly matching what the partner
block would compute on its side.

The result: at layer k=0 the j-step ratio across each sub-sub seam is
within ~1% of 1.0; over 20 marching layers it **converges** toward
1.0000 (centred differences damp the small initial disagreement).

```
for k in 0 .. nk−2:
    1. Snapshot the k-layer of all 6 blocks (consistent state for ghosts).
    2. For each sub_q:
       - gather ghost rows from sub_{q-1}, sub_{q+1}, cap_north, cap_south
       - march with seam-aware normals + Laplacian smoothing
    3. Re-stitch shared sub-sub edges (sub-FP-noise residual averaging).
    4. Jacobian check on sub-blocks (corners excluded).
    5. Extract pole rows from the new sub-block layers.
    6. For each cap:
       - gather ghost rows from the 4 sub-blocks
       - march with seam-aware normals + Laplacian smoothing
       - DIRICHLET-pin the 4 boundary edges from sub-block pole rows
    7. Jacobian check on cap blocks (corners excluded).
    8. Commit all 6 blocks to layer k+1.
```

## Jacobian convention

All six blocks use `_jacobian_field_open` (non-periodic j) because all
boundaries are open. The four **corner cells** at each cap's
`(0,0), (0,nc-2), (nc-2,0), (nc-2,nc-2)` are excluded from the check:
these cells border the junction of four sub-blocks at a single shared
node, and their signed volume is near-zero by construction — not a
mesh-quality defect, but a topological consequence of four blocks
meeting at one point. The tolerance for the remaining cells is `−1e−6`
(floating-point noise floor for the seam-aware normal computation).

## Controlling cap size with `theta_cap`

`theta_cap` (radians) is the polar exclusion angle. The 4 cap corners
sit at colatitude `theta_cap` on a unit sphere; cap-edge midpoints sit
at colatitude `atan(tan(theta_cap)/√2)` (smaller than `theta_cap` —
edges dip toward the pole, the natural shape of a great-circle arc
between two parallel-circle corners).

| `theta_cap` | degrees | Corner colatitude | Edge-midpoint colatitude |
|---|---|---|---|
| π/12 | 15° | 15.00° | 10.71° |
| π/6 | 30° (default) | 30.00° | 22.21° |
| π/4 | 45° | 45.00° | 35.26° |
| π/3 | 60° | 60.00° | 50.77° |

Larger `theta_cap` gives a physically larger cap face (more of the body
covered by caps, less by the equatorial band). The grid-space cap size
`nc = nj // 4 + 1` is set by `nj` and is independent of `theta_cap`.

## API

The topology is implemented across three modules:

- `shore.surface.cap.build_cap_k0(theta_cap, nc, pole, projector, ...)`
  builds the gnomonic-square cap k=0 seed.
- `shore.surface.equator.build_equator(cap_south, cap_north, ni, theta_cap,
  projector, ...)` builds the equator from the cap seam rings.
- `shore.mesh.topology.CapsTopology` orchestrates them and wires the 6
  HexBlocks (SHARED memory views for sub-sub seams; DIRICHLET pinning for
  cap-equator seams).

The high-level entry point is `Mesh.from_stl(...)` with
`topology=CapsTopology()`:

```python
from shore.mesh import CapsTopology, Mesh, Spacing1D
import numpy as np

m = Mesh.from_stl(
    "sphere.stl",
    topology=CapsTopology(),
    ni=40, nj=60, nk=20,
    theta_cap=float(np.pi / 6),
    spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
)
m.march(ds=0.02, growth=1.1)
m.write_geo("sphere")  # writes 6 sphere_<label>.geo files
```

See [`shore.topology.sphere_cap`](/reference/api-sphere-cap) for the
marching helpers (`_march_one_layer`, `_march_cap_layer`,
`_overwrite_cap_boundaries`, `_check_jacobian`).

## See also

- [Hyperbolic marching](/concepts/hyperbolic-marching) — underlying
  layer-by-layer algorithm
- [Spacing laws](/concepts/spacing) — `i_spacing` / `j_spacing` /
  `volume_k_spacing` knobs
- [Chimera mesh generation](/concepts/chimera) — where this topology fits
