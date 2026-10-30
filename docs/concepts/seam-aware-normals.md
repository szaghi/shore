# Seam-aware normals

Every block computes the marching normals from its own grid. At a block
seam (a row or column of nodes shared by two blocks) a one-sided
difference on each side gives two different tangents, hence two different
normals for the same node. The difference repeats at every layer and shows
on the outer layer as a step in cell size along every seam.

Seam-aware normals remove it: on an edge with a partner block, the
boundary tangent is a **centred difference against a ghost row** taken
from the partner. Both blocks then compute the same normal at the shared
node and move it the same way.

<Figure name="seam-aware-normals" ext="svg" caption="At a seam, one-sided tangents differ on the two sides (left); a centred difference against the partner's ghost row gives both sides the same tangent (right)." />

## One-sided differences

The `j_hi` column of sub-block $q$ is the `j_lo` column of sub-block
$q+1$. Without ghost rows:

- sub $q$ uses the backward difference `layer[:, -1] - layer[:, -2]`;
- sub $q+1$ uses the forward difference `layer[:, 1] - layer[:, 0]`.

The two tangents differ in length and, on a curved wall, in direction.
Averaging the shared column after each step makes the *positions* agree,
but the *directions* each block marches in still disagree, and the
disagreement accumulates.

## Centred differences against ghost rows

`_compute_normals(layer, periodic, ghost_rows)` takes an optional dict
with one row per face:

```python
ghost_rows = {
    "j_lo": ...,  # (ni, 3): the partner's row next to the shared column
    "j_hi": ...,
    "i_lo": ...,  # (nj, 3)
    "i_hi": ...,
}
```

For a face with a ghost row, the tangent at the boundary node is
$\tfrac12(\mathbf{P}_{1} - \mathbf{P}_{ghost})$: the same centred formula
the interior uses, and the same one the partner computes from its side.
The hyperbolic scheme uses the same ghost rows for its tangents and its
implicit terms ([Marching](/concepts/hyperbolic-marching)).

## Where the ghost rows come from

- **Sub–sub seams** (j): the `j_hi` ghost row of sub $q$ is column 1 of
  sub $q+1$, one in from the shared column, and symmetrically.
- **Cap–equator seams** (i): the cap and the sub-blocks are rotated and
  partly reversed relative to each other (the `_overwrite_cap_boundaries`
  mapping). `_sub_ghost_rows(q, ...)` returns, for each sub-block, the cap
  row that lines up with its pole row after that mapping;
  `_cap_ghost_rows(pole, ...)` does the inverse for the caps (the
  sub-blocks' second row from the pole).

All ghost rows of a layer are taken from one snapshot of the six blocks
at layer k, so no block reads a neighbour that has already been advanced.

The other topologies use the same mechanism: the O-grid wraps its
periodic j-direction, the C-grid passes ghost rows across the upper-TE
and lower-TE columns (`_march_cgrid_three`), and the hybrid wraps reuse the
cubed-sphere and C-grid marchers.

## Measured effect

On a unit icosphere (subdivision 4), cubed sphere `ni=60, nj=60, nk=20`,
$\theta_{cap} = \pi/6$, geometric steps `ds=0.02`, `growth=1.1`, the ratio
of the first j-cells on the two sides of each sub–sub seam,

$$
r_k(i) = \frac{\lVert \mathbf{P}^{(k)}_{q+1}[i, 1] - \mathbf{P}^{(k)}_{q+1}[i, 0]\rVert}
              {\lVert \mathbf{P}^{(k)}_{q}[i, -1] - \mathbf{P}^{(k)}_{q}[i, -2]\rVert},
$$

stays close to 1 and the deviation decreases outward:

| Layer k | max over seams and i of $\lvert r_k - 1\rvert$ |
|---:|---:|
| 0 | 0.0164 |
| 5 | 0.0067 |
| 10 | 0.0033 |
| 15 | 0.0015 |
| 19 | 0.0008 |

The k = 0 value is the wall grid's own mismatch; marching with centred
seam tangents damps it instead of growing it.

<Figure name="seam-continuity" caption="The outer layer over the wall: cell size varies smoothly across the block seams, following the growth law rather than the block decomposition." />

## Cap corners

At a cap corner three block edges meet (two cap edges and one equator
meridian). The [per-j C1 pin](/concepts/k0-construction#per-j-c1-seam-pinning)
gives them one first cell at k = 0; the seam-aware normals on both sides
of the cap–equator seam carry that agreement to the outer layers, and the
Dirichlet overwrite of the cap edges keeps the seam nodes identical at
every layer.

## See also

- [Marching](/concepts/hyperbolic-marching): where the normals are used
- [k=0 construction](/concepts/k0-construction): the C1 seam pin at the wall
- [Seams kept closed by edits](/adr/0040-seams-kept-closed): the related rule for in-place grid edits
- [`shore.topology.sphere_cap`](/reference/api/shore.topology.sphere_cap): `_march_one_layer` / `_march_cap_layer`
