# Marker boxes

The overset pre-processor cuts holes in the grids that pass through a
body using **marker boxes**: hexahedra listed in the `cc.par` boxes
section. They must **cover the body interior** and **stay inside the
body's near-body grid**. A single bounding box never satisfies both. SHORE
builds a cover of axis-aligned boxes in two stages, a voxel fill and a
greedy merge (`shore boxes-from-stl`, method `voxel-fill`, the default).

<Figure name="marker-boxes" caption="A body and its voxel-fill marker boxes after the greedy merge (translucent): slabs that follow the surface within the gap budget." />

For the data model and CLI / Python API, see the
[Marker boxes reference](/reference/boxes); for the underlying
generators, [`shore.boxes`](/reference/api/shore.boxes) and
[`shore.io.boxes`](/reference/api/shore.io.boxes).

## Two constraints

1. **Coverage of body volume.** Every body-interior cell-centre of
   the background mesh must lie inside at least one box, so the
   chimera marks it as donors-only. Equivalently, the union of the
   boxes must be a superset of the body.
2. **Confinement to the near-mesh.** Every box vertex must lie
   inside the body's near-mesh (the wall mesh). Box vertices that
   sit past the near-mesh outer extent spuriously hole-cut
   legitimate flow cells.

Constraint (2) is parameterised by `max_outward_gap` — the maximum
allowed distance from any box vertex to the body surface, measured
outward only (points on the box surface that lie inside the body
don't count).

A single AABB or OBB satisfies (1) trivially but fails (2): the
bounding cube's corners stick out by `r * (sqrt(3) - 1) ≈ 0.732 r`
past a circumscribed sphere of radius `r` — much further than any
practical near-mesh thickness.

## Stage 1 — voxel-fill seed

Voxelise the body's axis-aligned bounding box at side

$$
s = \frac{\text{max\_outward\_gap}}{\sqrt{3}}
$$

so the voxel half-diagonal is exactly `max_outward_gap / 2`. Include
a voxel iff its centre satisfies

$$
\text{signed\_distance}(\text{centre}) > -\tfrac{1}{2}\,\text{diag}
$$

with `signed > 0` inside the body and `< 0` outside (per
`trimesh.proximity.signed_distance`). This admits two classes of
voxels:

1. **Centre inside the body** (`signed > 0`).
2. **Body-shell voxels**: centre outside but within one half-diagonal
   of the surface (`-half_diag < signed <= 0`). These voxels
   *intersect* the body even though their centre sits just past the
   surface.

### Why both classes are needed

Class (1) alone leaves an uncovered shell of width up to `s/2`
along the body surface — body protrusions (sphere poles, equator
extremes) sit outside the box union. Adding class (2) closes that
shell: any point on the body surface lies within `half_diag` of
some grid cell centre, so it lands in at least one included voxel.

### Why the budget still holds

The maximum outward gap from any included voxel's corner to the
body surface is bounded by `2 * half_diag = s * sqrt(3) =
max_outward_gap`: a body-shell voxel whose centre is exactly
`half_diag` outside contributes a corner up to one further
half-diagonal beyond. So every box vertex satisfies constraint (2)
by construction.

## Stage 2 — greedy AABB merge

The voxel-fill seed is correct but wasteful. A cubic grid forces
all three axes to the same fine resolution even when the local body
geometry would tolerate a coarser slab — the seed cover of a
unit icosphere at `max_outward_gap = 0.3` is 1280 voxels.

Two adjacent voxels that both touch the same nearly-flat patch of
the body can be merged into one larger box without violating either
constraint. The cubic grid cannot represent that — but a greedy
post-process can.

### The merge

Each merge candidate is a contiguous **grid block** `[i0:i1,
j0:j1, k0:k1]`. Two grid blocks are merge-eligible iff:

1. They are face-adjacent along one axis `a` (one's high coord on
   `a` equals the other's low coord on `a`).
2. Their grid ranges on the other two axes match exactly.
3. The merged AABB's 8 world-space corners *all* have outward gap
   to the body ≤ `max_outward_gap`.

Condition (3) preserves constraint (2) at the enlarged corners.
Condition (1) plus the body-overlap cover from stage 1 preserves
constraint (1) — a merge is a strict superset of its parts, so
the union doesn't shrink.

### Pass order

For each axis in turn, in **descending order of the seed grid's
extent** on that axis: group active blocks by their (other-two-axes)
range, sort each group by the merge-axis low coord, and walk a
sweep — merge consecutive face-adjacent blocks whenever the merged
corners pass the gap check, restart on rejection or gap. Repeat
the sweep until a pass makes no merges.

Sweeping the longest axis first biases the output toward
**slab-shaped boxes** when the body is elongated. For a sphere all
three axes are equivalent and the result is roughly cubic per
local curvature.

### Result

Measured on a unit icosphere (subdivision 3):

| `max_outward_gap` | Seed voxels | Boxes after merge | Largest outward gap of a box corner |
|---|---|---|---|
| 0.30 | 1280 | 49 | 0.283 |
| 0.15 | 8624 | 219 | 0.138 |

The merged count tracks the **body surface area**, not its volume: the
merge fuses interior cells freely and stops near the surface at the local
curvature scale. Halving the budget multiplied the seed by 6.7 but the
merged count by 4.5.

## Limits

- **Non-watertight STL.** The voxel-fill inclusion test relies on
  `trimesh.contains` / `signed_distance`, which assume a closed
  surface. The function rejects non-watertight STLs with a clear
  `MeshInputError`. Repair the STL or fall back to single-box
  methods (`aabb`, `obb`) which don't need inside / outside
  classification.
- **Tiny `max_outward_gap`.** The seed grid scales as
  `body-volume * (sqrt(3) / gap)^3`. At sub-percent gap budgets
  the seed dominates memory; the merger is linear in seed size, so
  the limiting factor is the per-merge `signed_distance` query
  budget. Practically: pick `max_outward_gap` smaller than the
  near-mesh thickness, but no smaller than necessary.

## Single-box fallbacks

Two best-effort methods are kept for quick previews:

- **`aabb`** — one world-axis-aligned box bounding the body.
- **`obb`** — one oriented bounding box from area-weighted PCA.
  Better than `aabb` for elongated bodies because the OBB aligns
  with the body's principal axis.

Both measure the actual outward gap on a sampled grid of box-surface
points and emit a `UserWarning` when it exceeds the budget. Useful
for very loose budgets on near-cubic bodies, or to eyeball whether
a body needs the full voxel-fill machinery.

## See also

- [Marker boxes reference](/reference/boxes) — data model, CLI,
  ParaView visualisation.
- [`shore.boxes`](/reference/api/shore.boxes) — `boxes_from_stl` Python
  API.
- [`shore.io.boxes`](/reference/api/shore.io.boxes) — the `Box`
  dataclass.
- [Chimera mesh generation](/concepts/chimera), [Xall bundle tutorial](/tutorials/xall-bundle)
- [`shore cc-par --boxes`](/reference/cc-par) — bake the boxes
  into `cc.par`.
