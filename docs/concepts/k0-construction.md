# k=0 construction (cubed sphere)

The k=0 layer is the wall grid from which the volume is marched. For the
[cubed sphere](/concepts/topologies/cubed-sphere) (and the wrap of the
[OH hybrid](/concepts/topologies/oh-hybrid)) it is built in two stages:

1. **Caps**: a flat diamond in the tangent plane at each pole, cast onto
   the body (gnomonic projection).
2. **Equator**: great-circle meridians from the south cap seam to the
   north cap seam, cast onto the body, with the first and last interval of
   each meridian pinned to the cap edge (C1 spacing at the seam).

Both stages cast rays through the [projector](/concepts/projection). The
result is six arrays: four equatorial sub-blocks of `(ni, nc)` nodes and
two caps of `(nc, nc)` nodes, with `nc = nj // 4 + 1`.

## Caps: gnomonic flat-square projection

At each pole a flat **diamond** is placed in the tangent plane $z = 1$,
corners on the cardinal axes at distance $h = \tan\theta_{cap}$.
Parametrised by $(a, b) \in [0, 1]^2$, bilinear interpolation gives

$$
x(a, b) = h\,(1 - a - b), \qquad y(a, b) = h\,(a - b), \qquad z = 1 .
$$

Each tangent-plane point is normalised to a direction and cast through
the anchor. On a sphere centred on the anchor, the four corners land at
colatitude $\theta_{cap}$, the cap edges are great-circle arcs, and the
edge midpoints sit at colatitude $\arctan(\tan\theta_{cap}/\sqrt{2})$,
closer to the pole than the corners. The cells come out nearly square
by construction (aspect ratio about 1), with no concentric pattern
around the pole.

<Figure name="gnomonic-cap" ext="svg" caption="The cap diamond in the tangent plane (left) and its radial projection onto the body (right): the index structure stays that of the flat square." />

**Why a diamond.** The four sub-block junctions sit at azimuths 0°, 90°,
180°, 270°. An axis-aligned square would put its corners at ±45° and
±135°, between the junctions; turned by 45° its corners land exactly on
them, which is the contract `_overwrite_cap_boundaries` expects.

**Cap-edge spacing.** The 1D coordinates on $[-h, +h]$ follow
`j_spacing`: `uniform`, or `tanh2` to cluster nodes toward the cap
corners. Because the cap edges become the equator's pole rows, the
equator inherits the same distribution along its seams. (`tanh` is not
offered: the law must be symmetric in the four quadrants.)

The OH hybrid passes `phi_offset = pi/4` so the four sub-block centres,
not their corners, face the world axes.

## Equator: great-circle meridians

The seam ring of each cap (its four edges in sub-block order, shared
corners dropped) gives `nj` seam nodes per pole. For every column $j$:

1. Slerp from `seam_south[j]` to `seam_north[j]` along the great circle
   through the two anchor directions, at the colatitudes of a 1D
   distribution (below).
2. Cast the interior directions onto the body.
3. Copy the two seam rows from the caps bit for bit (no re-cast), so the
   cap–equator seams are conforming to round-off.

## Per-j C1 seam pinning

At each seam vertex three block edges meet: two cap-edge edges and one
equator meridian. For the spacing to be continuous across the seam all
three must start with the same first cell.

<Figure name="junction-3edge" ext="svg" caption="A cap corner: two cap-edge edges and one equator meridian meet, and all three take the same first cell." />

The cap-edge first cell at column $j$ is the angle, seen from the anchor,
between adjacent seam nodes:

$$
\Delta s(j) = \arccos\!\left(\hat u_{seam}(j)\cdot \hat u_{seam}(j+1)\right).
$$

The meridian is parametrised in colatitude, not in anchor angle, so the
step is converted per column,

$$
\Delta\theta_{first}(j) = \Delta s(j)\,\frac{\text{interval}_j}{\omega_j},
$$

with $\text{interval}_j$ the colatitude span of the meridian and
$\omega_j$ the angle it subtends at the anchor (equal on a sphere centred
on the anchor; the ratio keeps it exact elsewhere). The same is done at
the north seam for the last interval. `pinned_distribution` then places
the first and last interval exactly and shapes the interior `ni − 2` cells
with `i_spacing` (`uniform`, `tanh` toward the south seam, `tanh2` toward
both seams).

### Compatibility window and fallback

Pinning is used on a column only when each pinned step is within a factor
2 of the step the unpinned `i_spacing` law would give there
($0.5 \le \Delta\theta_{pin}/\Delta\theta_{natural} \le 2$), and the two
pins take less than 90 % of the meridian. Otherwise that column uses the
unpinned law, and the build emits one `UserWarning` ("Cap radial step ...
is incompatible ..."). A forced pin outside the window would put a large
cell-size jump just inside the seam, which the marcher then carries
outward.

<Figure name="per-j-pinning" ext="svg" caption="Compatible (left): the cap step is within a factor 2 of the natural equator step and is pinned; incompatible (right): the equator falls back to its unpinned law." />

For `uniform` spacing on a sphere, the cap step is the corner-to-corner
arc divided by the cap cells and the natural equator step is the meridian
divided by its cells, so the rule is

$$
\frac{1}{2} \;\le\; \frac{\arccos(\cos^2\theta_{cap})/(n_c - 1)}{(\pi - 2\theta_{cap})/(n_i - 1)} \;\le\; 2,
\qquad n_c = \frac{n_j}{4} + 1 .
$$

The defaults (`--ni 40 --nj 60`, $\theta_{cap} = 9°$) give a ratio of
about 0.2: they are **outside** the window and the CLI warns. With
`--theta-cap-deg 30` the ratio is about 0.9 and every column is pinned.

## Final layout

`CubedSphereTopology.build` allocates one contiguous ring of shape
`(ni, 4·(nc − 1) + 1, nk, 3)` for the four sub-blocks; each sub-block is a
numpy view of a 90° slice, so adjacent sub-blocks share their `j_hi` /
`j_lo` columns in memory (`SHARED`). The caps are separate `(nc, nc, nk,
3)` arrays bound to the sub-block pole rows (`DIRICHLET`);
`_overwrite_cap_boundaries` rewrites the cap edges from the pole rows at
every marched layer.

<Figure name="sphere-blocks" caption="The six blocks of the cubed sphere at k=0: four equatorial sub-blocks and a polar cap." />

<Figure name="sphere-topdown" caption="Top-down view of a cap: square cells, corners on the sub-block junctions, edges dipping toward the pole." />

## See also

- [Surface projection](/concepts/projection): the ray-casting kernel
- [Marching](/concepts/hyperbolic-marching) and [Seam-aware normals](/concepts/seam-aware-normals): how the pinned spacing is carried outward
- [Cubed-sphere topology](/concepts/topologies/cubed-sphere): block layout and connectivity
- [`build_cap_k0`](/reference/api/shore.surface.cap#build_cap_k0), [`build_equator`](/reference/api/shore.surface.equator#build_equator), [`pinned_distribution`](/reference/api/shore.volume.spacing#pinned_distribution)
