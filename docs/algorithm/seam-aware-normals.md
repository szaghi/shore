# Seam-aware normals

When SHORE marches the 6 blocks of the cubed-sphere topology outward, each
block computes its own outward normals from its own grid. At a
**multi-block seam** — a column or row shared between two blocks —
the two blocks' boundary tangents disagree if each uses a one-sided
difference internal to itself. Marching with mismatched normals
accumulates that disagreement layer by layer as a visible spacing
jump on the outer surface.

**Seam-aware normals** fix this by computing the boundary tangent as a
**centred difference against a ghost row** drawn from the partner
block's data. Both blocks then arrive at the same normal at the shared
node and march in step.

![Seam-aware normals: one-sided vs ghost-row centred](/figures/seam-aware-normals.png)

## The problem (one-sided differences)

Sub_q's `j_hi` boundary is the same column of nodes as sub_{q+1}'s
`j_lo` boundary. Without ghost rows:

- Sub_q computes the j-tangent at its `j_hi` boundary as
  `layer[:, -1] - layer[:, -2]` — backward, internal to sub_q.
- Sub_{q+1} computes the j-tangent at its `j_lo` boundary as
  `layer[:, 1] - layer[:, 0]` — forward, internal to sub_{q+1}.

These two tangents have different magnitudes and (generally) different
directions, so the cross-product normals that drive the march also
differ. Both blocks then advance the same shared node in two
inconsistent directions, the post-march averaging of the SHARED column
absorbs the discrepancy in *position*, but the underlying *velocity*
mismatch persists. After many layers the outer surface shows a stripe
across each block seam.

## The fix (centred difference via ghost rows)

Each block's `_compute_normals` accepts an optional `ghost_rows` dict
with one entry per face:

```python
ghost_rows = {
    "j_lo": <(ni, 3) row from the partner block>,
    "j_hi": <(ni, 3) row from the next partner>,
    "i_lo": <(nj, 3) row from cap_south>,
    "i_hi": <(nj, 3) row from cap_north>,
}
```

When a face has a ghost row, that boundary's tangent becomes a
**centred difference** between the next interior row and the ghost row
— the same formula both blocks would compute on their side. Both
blocks now agree on the normal at the shared node, and march along the
same velocity by construction.

## Where the ghost rows come from

Sub-sub j-seams are easy: sub_q's `j_hi` ghost row is sub_{q+1}'s
column-1 (one in from sub_{q+1}'s `j_lo` boundary), and symmetrically
on the other side.

Cap-equator i-seams are harder because the cap and sub-block grids are
**rotated and possibly reversed** relative to each other, per the
existing `_overwrite_cap_boundaries` mapping. The
`_sub_ghost_rows(q, sub_k, cap_n_k, cap_s_k)` helper bakes in this
mapping, so for each `q` it returns the cap row that, after the
appropriate reversal, lines up with sub_q's i-boundary. Same for
`_cap_ghost_rows(pole, sub_k, ni)` from the cap's perspective.

The ghost rows are gathered **once per layer**, from a snapshot of all
6 blocks at the current k. This gives a consistent set of partner data
across all 6 blocks (no race condition where one block has already
been advanced when its neighbour reads the ghost).

## Numerical evidence

On a unit-icosphere mesh (`ni=60, nj=60, nk=20, theta_cap=π/6`), the
j-step ratio across the sub_q/sub_{q+1} seam at marching layer `k`,
defined as

$$
r_k(i) = \frac{|\mathbf{P}^{(k)}_{q+1}[i, 1] - \mathbf{P}^{(k)}_{q+1}[i, 0]|}
              {|\mathbf{P}^{(k)}_{q}[i, -1] - \mathbf{P}^{(k)}_{q}[i, -2]|},
$$

stays close to 1 at every layer:

| Layer `k` | Max $|r_k - 1|$ |
|---:|---:|
| 0 | 0.0163 |
| 5 | 0.0067 |
| 10 | 0.0033 |
| 15 | 0.0015 |
| 19 | 0.0008 |

The deviation **shrinks** with `k` — centred differences damp the
small initial disagreement out of the system over the march.

## Cap interior at seam vertices

Where 3 multi-block edges meet (cap corner: 2 cap-seam edges + 1
equator meridian), the seam-aware normals propagate the
[per-j C1 seam pinning](k0-construction#per-j-c1-seam-pinning) outward
with bit-exact agreement. The `cap` blocks' i-boundary normals draw
their ghost row from the appropriate sub-block's i-2 row (mapped per
`_overwrite_cap_boundaries`); the sub-block's i-boundary normals draw
theirs from the cap's row 1 (the same mapping inverted).

The result: at every layer, cap-edge first cells equal the equator
meridional first cells equal each other to FP precision, even after 20
layers of marching.

## See also

- [Hyperbolic marching](hyperbolic-marching) — where the normals are
  used
- [k=0 mesh construction](k0-construction) — the C1 seam pin that
  this propagates outward
- [`shore.volume.hyperbolic`](/reference/api-hyperbolic) — `_compute_normals`
  with the `ghost_rows` kwarg
- [`shore.topology.sphere_cap`](/reference/api-sphere-cap) —
  `_march_one_layer` / `_march_cap_layer` with `ghost_rows`
- `shore.mesh.mesh._sub_ghost_rows` and `_cap_ghost_rows` — the helpers
  that build the ghost-row dicts per layer
