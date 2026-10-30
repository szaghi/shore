# Build and edit a blocking

Use a blocking when no built-in topology fits the body and you want to lay
out the blocks yourself: a skeleton of hex blocks, edited by commands and
premeshed into ordinary blocks.

The recipes below build one case step by step: a unit sphere in a box from
-3 to 3, cut into a 3 x 3 x 3 lattice, the centre block removed, an O-grid
on the hole, wall-normal lines and a wall layer. Each recipe continues from
the previous one. The last two recipes are independent small cases. The
whole script runs in about 2 s; it writes `sphere.stl` and the journal
`sphere_blocking.journal.py`, which `shore run` replays in about 5 s.

Full script: `examples/cookbook/blocking.py`. Two helpers in it are used
below: `block_at(blocking, point)` returns the block whose skeleton box
holds a point, and `min_jacobian(blocks)` the smallest
[scaled Jacobian](/concepts/quality) over every cell.

## Split a block at a node

Start from one block around the region and cut it into a lattice.

<<< @/../examples/cookbook/blocking.py#split{python}

```text
(13, 13, 13)
27 blocks
```

- `split_blocking(block, axis, fraction)` cuts at node
  `round(fraction * (n - 1))` of the block's `n` nodes along `axis`. Here
  13 nodes over 6 units put a node every 0.5, so `x = -1` is node 4 exactly.
- The cut runs through every block sharing a cut edge (the sheet), so the
  blocking stays conforming: one `split_blocking` along i cuts the whole
  slab of blocks.
- The low part keeps its label; the high part of each cut block is `b<n>`,
  numbered up over the project. Track blocks by position (as `block_at`
  does) rather than by guessing labels.
- `create_blocking(ni, nj, nk, margin=1.0)` without `vmin`/`vmax` takes the
  geometry's bounding box, grown by `margin` times its extent on each side.
- Refused (`ParameterError`, project unchanged): a fraction outside (0, 1)
  or one that does not fall on an interior node, a bad axis or block.

## Delete a block to make a hole round the body

Associate the centre block's faces with the sphere, then delete the block:
its faces become the hole's boundary, which premesh projects onto the body
and marks WALL.

<<< @/../examples/cookbook/blocking.py#hole{python}

```text
26 blocks, 6 WALL faces
min scaled Jacobian 0.699, wall angle 40.8 deg
```

- `associate_blocking(block, entity)` defaults to geometry `"body"` and
  kind `"body"`. Edges and vertices of the associated faces inherit the
  association when their faces agree, so the hole's corners land on the
  sphere too.
- Premesh spaces edges equally, fills faces and blocks by transfinite
  interpolation, and wires shared faces as SHARED seams (bit for bit),
  associated boundary faces as WALL, the rest FREE.
- The plain lattice is valid but poor at the wall: the grid lines meet the
  sphere up to 41 degrees off its normal (the angle printed is the largest
  between a wall node's first segment and the radius). The next recipe
  fixes that.
- Deleting every block is refused.

## Insert an O-grid around the body

Layer the faces of the hole. Name them as `(block, face)` pairs: here, the
faces of the remaining blocks that coincide with the deleted centre's.

<<< @/../examples/cookbook/blocking.py#ogrid{python}

```text
32 blocks
min scaled Jacobian 0.439, wall angle 5.7 deg
```

<Figure name="blocking-ogrid" caption="O-grid insertion: a block set before and after a layer of blocks on its faces." />

- Each layered face gets a block between it and an inner copy of the face,
  with `nodes` nodes across (17 here, so the wall layer below has room);
  the six layer blocks bring the count from 26 to 32.
- `factor` sets how far the copies move: every vertex of a layered face
  moves into the set by `(1 - factor)` times half the sum of its edges that
  leave the layered faces. Inside one block, the core is the block scaled
  by `factor` about its centre.
- A layered face keeps its association on the layer's outer face, so the
  wall stays on the sphere.
- Instead of layering the hole, you can O-grid the centre block first,
  associate its core's faces and delete the core: the same 32 blocks.

To O-grid a single block, or leave one face open (a half-O, for a body on a
wall or a symmetry plane):

<<< @/../examples/cookbook/blocking.py#ogrid-box{python}

```text
['b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6'] (4, 9, 9)
['b0', 'b1', 'b2', 'b3', 'b4', 'b5']
```

- With `blocks` given and no `faces`, every face on the set's boundary is
  layered, less `exclude`. An excluded (open) face stays in its plane, and
  its association passes to the layer faces lying on it.
- Refused: `factor` outside (0, 1), `nodes` below 2, a face not on the
  set's boundary, an open face shared with a block outside the set (add
  that block to the set, or layer the face), nothing left to layer.
- Not supported: O-grids through faces shared with blocks outside the set,
  smoothing of the layer ([ADR 0043](/adr/0043-ogrid-insertion)).

## Make grid lines leave the wall along its normal

After premesh, the lines across the wall follow the transfinite
interpolation and leave a curved wall obliquely. `orthogonalize` reshapes
them.

<<< @/../examples/cookbook/blocking.py#orthogonalize{python}

```text
min scaled Jacobian 0.439, wall angle 1.9 deg
```

- In every block with exactly one WALL face, each line across the wall
  becomes a cubic Hermite curve leaving the wall along its normal (tangent
  `strength` times its chord). The wall and the opposite face do not move;
  the nodes keep their arc-length fractions.
- The wall normal is common to every block holding a node, so seams stay
  closed; neighbours whose faces move are re-filled.
- On the plain 26-block lattice it also helps (41 to below 25 degrees in
  the acceptance spec), but the O-grid comes first for a body-fitted wall.
- Refused: `strength <= 0`, or one that would invert a valid block; the
  mesh and journal are then unchanged.
- Not handled: sharp wall edges (no single normal) and blocks with two or
  more walls ([ADR 0048](/adr/0048-orthogonal-wall-lines)).

## Cluster a wall layer with an absolute first cell

Prescribe the first cell height at the wall; the cells grow geometrically
along each wall-normal line.

<<< @/../examples/cookbook/blocking.py#wall-layer{python}

```text
b27: first cell 1.00e-03, ratio 1.391
min scaled Jacobian 0.439
```

- On every line across the wall of a block with one WALL face, the first
  cell is exactly `first_step` long and the ratio is solved per line from
  its length. Here the lines are 0.5 to 1.6 long; with 17 nodes the
  printed (face-centre) line grows by 1.39 per cell: give the layer more
  nodes for a gentler growth. There is no
  growth-ratio limit.
- Run it after `orthogonalize`: the order is
  `premesh → orthogonalize → wall_layer(ds)`. The wall node and the far
  node stay put; seams stay closed.
- It needs a layer topology (an O-grid round the body). On the plain
  lattice it would invert the corner blocks and is refused, as are a first
  step `<= 0` and one not below a line's uniform spacing.
- The layer is algebraic, not marched ([ADR 0050](/adr/0050-wall-layer)).
  For other spacing laws on an edge, use
  [redistribute](/cookbook/redistribute) after premesh.

## Move a vertex, freely or onto the geometry

Move a skeleton vertex to shape the blocks. A free vertex goes where you
say; an associated one slides onto its geometry.

<<< @/../examples/cookbook/blocking.py#move{python}

```text
[1.2 1.3 1.1]
[0.22562802 0.68167172 0.69153652] |x| = 0.9969
```

- Vertices are `v0`..`v7` in `HexBlock.vertices()` order (index
  `4 i + 2 j + k`, so `v7` is the `(i_hi, j_hi, k_hi)` corner). Every block
  sharing the vertex follows: `b1`'s `v3` is `b0`'s `v7`.
- An associated vertex (its own association, or inherited from its edges
  or faces) moves to the projection of the point onto the geometry: here
  onto the facetted sphere STL, hence `|x|` slightly below 1. After a
  geometry edit, replaying the journal re-projects it.
- A move is not checked for validity: it can invert a block, and premesh
  reports nothing. Check with [`shore.quality`](/cookbook/quality).
- Refused: a bad block or vertex name, a face or edge name, `xyz` that is
  not three finite numbers ([ADR 0045](/adr/0045-move-blocking-vertex)).

## Notes

- Every command is journaled and undoable (`p.undo()`); the journal replays
  with `shore run` ([Projects](/cookbook/projects)). A project loaded
  without replay has its mesh but no blocking.
- The GUI edits the same blocking with picks on the skeleton
  ([GUI reference](/reference/gui)).
- Not done: merge and collapse of blocks, elliptic smoothing, spacing laws
  on blocking edges ([ADR 0042](/adr/0042-blocking-skeleton-premesh)).

## See also

- [Blocking](/concepts/blocking): the model, premesh and inheritance.
- [Association](/concepts/association) and
  [Associate with geometry](/cookbook/associate).
- [`shore run` → Blocking](/reference/run#blocking): every blocking command.
- API: [`Project.create_blocking`](/reference/api/shore.project.project#Project-create_blocking),
  [`Project.ogrid`](/reference/api/shore.project.project#Project-ogrid),
  [`Project.wall_layer`](/reference/api/shore.project.project#Project-wall_layer),
  [`Blocking`](/reference/api/shore.mesh.blocking#Blocking).
- ADRs: [0042](/adr/0042-blocking-skeleton-premesh),
  [0043](/adr/0043-ogrid-insertion), [0045](/adr/0045-move-blocking-vertex),
  [0048](/adr/0048-orthogonal-wall-lines), [0050](/adr/0050-wall-layer).
