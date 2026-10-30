---
title: shore.surface.cap
outline: [2, 3]
---

# `shore.surface.cap`

Polar cap k=0 seed face via gnomonic projection of a flat tangent square.

The cap is constructed by placing a flat planar square (a diamond, with its 4
corners on the tangent-plane cardinal axes at distance `tan(theta_cap)`) at
the cap pole and then projecting every node radially onto the body surface
through the cached `Projector`.  The 4 corners land at colatitude
`theta_cap` exactly (on a sphere); cap edges are great-circle arcs between
adjacent corners; the interior is filled by bilinear interpolation in the
tangent plane (which produces equispaced corner-to-corner straight lines and,
after radial projection, equispaced great-circle arcs on the sphere).

Why a diamond and not an axis-aligned square: the 4 cap corners must coincide
with the 4 sub-block-junction longitudes `phi ∈ {0, π/2, π, 3π/2}` to
satisfy the existing `_overwrite_cap_boundaries` contract.  A diamond with
corners on the cardinal axes places the cap-corners exactly at those
longitudes.  An axis-aligned square would put corners at ``phi ∈ {±π/4,
±3π/4}`` — a 45° offset from the sub-block junctions.

Spacing inheritance
-------------------
The tangent-plane parametric grid uses `build_j_coords`-style 1D laws on
`[-tan(theta_cap), +tan(theta_cap)]`, so the cap edges' node distribution
matches the equator's `j_spacing` law along the seam.  Cap edges are then
the equator pole rows (bit-exact) and the equator inherits this distribution
on its first/last latitude row.

Pole convention
---------------
- `pole='south'` (shore: `cap_south`, sphere top, z &gt; 0): no rotation.
- `pole='north'` (shore: `cap_north`, sphere bottom, z &lt; 0): cap-local
  `+z'` axis maps to world `-z`; the (u, v) tangent-plane coords are
  reflected so the cap's index orientation matches the existing override
  convention (`cap[0, :]` → sub0 pole row, etc.).

## `build_cap_k0` <Badge type="info" text="function" /> {#build_cap_k0}

```python
build_cap_k0(
    theta_cap: float,
    nc: int,
    pole: str,
    projector: object,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    phi_offset: float = 0.0,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/cap.py#L64">source</a>

Build the `(nc, nc, 3)` cap k=0 seed face via gnomonic projection.

**Parameters**

- **`theta_cap`** — Polar exclusion angle in radians.  The 4 cap corners sit at colatitude `theta_cap` on a sphere; cap-edge midpoints sit at colatitude `atan(tan(theta_cap)/sqrt(2))` (smaller than `theta_cap` — edges dip toward the pole, which is the natural shape of a great-circle arc between two parallel-circle corners).
- **`nc`** — Cap resolution.  Must equal `nj // 4 + 1` so the cap edges have the same node count as the sub-block pole rows.
- **`pole`** — `'north'` (cap_north, `z<0` of the body) or `'south'` (cap_south, `z>0`).
- **`projector`** — `Projector` with cached BVH.
- **`j_spacing`** — `'uniform'` (default) or `'tanh2'`.  The tangent-plane node distribution along the cap edges; inherited by the equator seam.
- **`j_beta`** — Clustering strength for `'tanh2'` (ignored for `'uniform'`).
- **`phi_offset`** — Azimuthal rotation (radians) applied to the gnomonic diamond about the polar axis.  Default `0.0` places cap corners at world azimuths `{0, pi/2, pi, 3*pi/2}` (the CubedSphereTopology convention). :class:`shore.mesh.OHTopology` passes `pi/4` so the four sub-block *centres* (rather than corners) face the world cardinal axes — which is what makes a 1:1 sub-block-to-box-face mapping work for the wrap-channel hybrid.  The cap array's parametric (a, b) layout is unchanged; only the world-space orientation rotates.

**Returns**

- `numpy.ndarray, shape ``(nc, nc, 3)``` — Cap k=0 layer with nodes lying on the body surface.
