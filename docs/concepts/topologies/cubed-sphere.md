# Cubed-sphere topology

The cubed sphere is the default topology of `shore mesh`. It splits the
grid around a closed body into six blocks, the faces of a cube projected
onto the body: four equatorial sub-blocks and two polar caps (Sadourny
1972; Ronchi, Iacono & Paolucci 1996). No block has a singular axis, every
seam is point-matched, and each cap is a regular structured block.

It is not an O-grid with caps added: the two topologies are built
independently and share only the spacing kernels.

<Figure name="sphere-blocks" caption="The six cubed-sphere blocks at the wall: four equatorial sub-blocks and a polar cap." />

## Blocks

With `nc = nj // 4 + 1` (`nj` the equator-ring node count, rounded down
to a multiple of 4 with a warning):

| Block | Shape | Covers |
|---|---|---|
| `sub0` | `(ni, nc, nk)` | azimuth 0°–90° |
| `sub1` | `(ni, nc, nk)` | 90°–180° |
| `sub2` | `(ni, nc, nk)` | 180°–270° |
| `sub3` | `(ni, nc, nk)` | 270°–360° |
| `cap_north` | `(nc, nc, nk)` | the cap at $z < 0$ of the body |
| `cap_south` | `(nc, nc, nk)` | the cap at $z > 0$ of the body |

(The cap names are historical: `cap_north` sits at negative $z$.) In
each sub-block i runs along the meridian from the `cap_south` seam
(`i = 0`) to the `cap_north` seam, j along the azimuth, k away from the
wall. `--ni 65 --nj 128 --nk 65` gives 64×32×64-cell sub-blocks and
32×32×64-cell caps.

<Figure name="sphere-topdown" caption="Top-down view of a cap: four corners on the sub-block junctions at colatitude theta_cap, edges dipping toward the pole between them." />

## Connectivity

**Sub–sub seams (j).** The four sub-blocks are views into one contiguous
ring of shape `(ni, 4·(nc − 1) + 1, nk, 3)`, so the shared j-columns are
the same memory. Faces `sub_q.j_hi ↔ sub_{q+1}.j_lo` (and `sub3.j_hi ↔
sub0.j_lo`) are `SHARED`.

**Cap–equator seams (i).** Each cap edge is bound to one sub-block pole
row by a `DIRICHLET` face, and `_overwrite_cap_boundaries` copies the
rows into the cap at every marched layer:

| `cap_north` edge | ← | `cap_south` edge | ← |
|---|---|---|---|
| `[0, :]` | `sub0[ni-1, :]` | `[:, 0]` | `sub0[0, :]` |
| `[:, -1]` | `sub1[ni-1, :]` | `[-1, :]` | `sub1[0, :]` |
| `[-1, :]` | `sub2[ni-1, ::-1]` | `[:, -1]` | `sub2[0, ::-1]` |
| `[:, 0]` | `sub3[ni-1, ::-1]` | `[0, :]` | `sub3[0, ::-1]` |

The two caps are wound oppositely (their outward normals point opposite
ways), hence the two mappings. Every block's `k_lo` is `WALL`, every
`k_hi` is `FREE`.

## Construction

`CubedSphereTopology.build`:

1. builds both caps' k=0 by gnomonic projection of a tangent-plane
   diamond ([k=0 construction](/concepts/k0-construction));
2. builds the equator k=0 from great-circle meridians between the two cap
   seam rings, pinning each meridian's end intervals to the cap edge when
   the steps are compatible (C1 spacing across the seam);
3. allocates the ring and the caps and wires the faces.

`Mesh.march` then advances the six blocks together, layer by layer, with
[seam-aware normals](/concepts/seam-aware-normals) and the cap edges
re-pinned at every layer ([Marching](/concepts/hyperbolic-marching)).

## Spacing at the seams

At each cap corner three block edges meet: two cap edges and one
meridian. The per-j pin makes their first cells equal; it applies when
the cap step and the equator's natural step are within a factor 2.
Outside that window the equator falls back to its unpinned law with one
`UserWarning`. The defaults (`--ni 40 --nj 60`, `theta_cap` 9°) are
outside it; `--theta-cap-deg 30` is inside. The rule is in
[k=0 construction](/concepts/k0-construction#compatibility-window-and-fallback).

<Figure name="junction-3edge" ext="svg" caption="A cap corner: the two cap edges and the meridian share one first cell." />

<Figure name="seam-continuity" caption="Outer layer over the wall: the cell size follows the growth law across the block seams, not the block decomposition." />

## Parameters

| Parameter | CLI | Default | Effect |
|---|---|---|---|
| `theta_cap` | `--theta-cap` / `--theta-cap-deg` | $0.05\,\pi$ (9°) | colatitude of the cap corners; sets the cap size |
| `j_spacing`, `j_beta` | `--surface-j-spacing`, `--surface-j-beta` | `uniform`, 3.0 | cap-edge distribution (`uniform` or `tanh2` toward the corners), inherited by the equator seams |
| `i_spacing`, `i_beta` | `--surface-i-spacing`, `--surface-i-beta` | `uniform`, 3.0 | meridian interior (`tanh` toward the south seam, `tanh2` toward both) |
| k law | `--volume-k-spacing`, `--ds`, `--growth`, ... | geometric, 0.001, 1.15 | the step schedule shared by the six blocks |

Per-edge spacing: the k-edges of all six blocks can carry one shared
`Spacing1D` override, read by `Mesh.march`; per-edge i/j is not available
here, because the seam pin is geometric
([Spacing laws](/concepts/spacing-laws#per-edge-spacing)).

## Limits

- The body must be star-shaped from the anchor
  ([Pipeline limits](/concepts/pipeline-limits#star-shaped-bodies)).
- The outer layer follows the body; for a grid that ends on a box use the
  [OH hybrid](oh-hybrid).

## See also

- [Sphere tutorial](/tutorials/sphere): a cubed sphere end to end
- [`shore mesh`](/reference/mesh): every flag
- [k=0 construction](/concepts/k0-construction), [Seam-aware normals](/concepts/seam-aware-normals), [Marching](/concepts/hyperbolic-marching)
- [O-grid](o-grid): the single-block alternative
- [`CubedSphereTopology`](/reference/api/shore.mesh.topology#CubedSphereTopology), [`build_cap_k0`](/reference/api/shore.surface.cap#build_cap_k0), [`build_equator`](/reference/api/shore.surface.equator#build_equator)
