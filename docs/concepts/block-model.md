---
title: Block model
---

# Block model

Every grid SHORE makes, whichever way it is made, ends as a list of
structured hexahedral blocks. A block is an object, not a bare array: it
owns its node array and a typed description of its 8 corners, 12 edges and
6 faces. The face descriptions carry the boundary conditions and the seams
to other blocks; the writers (`.geo`, `.grd`, `cc.par`, CGNS, Plot3D) read
everything they need from them.

The classes live in `shore.mesh`: `Block` (abstract) and `HexBlock`
([`shore.mesh.block`](/reference/api/shore.mesh.block#HexBlock)), `Vertex`,
`Edge`, `Face`, `FaceBC` and `AxisMap`
([`shore.mesh.face`](/reference/api/shore.mesh.face#AxisMap)), and
`Spacing1D` ([`shore.mesh.spacing`](/reference/api/shore.mesh.spacing#Spacing1D)).

## Nodes

A `HexBlock` holds `nodes`, an array of shape `(ni, nj, nk, 3)`: node
`(i, j, k)` has coordinates `nodes[i, j, k] = (x, y, z)`.

- `ni`, `nj`, `nk` are **node** counts. A block has `(ni - 1)(nj - 1)(nk - 1)`
  cells. The `.grd` header is the only place where cell counts appear, and
  the writer converts.
- `k` is the wall-normal direction of body-fitted blocks: `k = 0` is the
  wall, `k = nk - 1` the outer boundary.
- A valid cell is right-handed: $(\mathbf{e}_i \times \mathbf{e}_j) \cdot
  \mathbf{e}_k > 0$, with $k$ pointing away from the wall. The quality
  metrics and the marchers' Jacobian guards use this orientation (see
  [Grid quality](/concepts/quality)).
- `.geo` files store the nodes i-fastest (Fortran order), as Xall reads them.

`nodes` is often a numpy *view* into a buffer the topology builder allocated
for several blocks, so that neighbouring blocks share their seam nodes in
memory (see [Seams](#seams-views-and-write-through) below).

## Vertices, edges and faces

<Figure name="block-model" ext="svg" caption="A HexBlock: vertex v = 4i + 2j + k over the corner bits, 12 edges named by their two fixed sides, 6 faces i_lo … k_hi; k_lo is the wall by default." />

| Entity | Count | Naming | Example |
|---|---|---|---|
| Vertex | 8 | `v0` … `v7`, index $4i + 2j + k$ with $i, j, k \in \{0, 1\}$ for low / high corner | `v5` is `nodes[-1, 0, -1]` |
| Edge | 12 | the two **fixed** sides, e.g. `j_lo_k_lo` runs along `i` at `j = 0`, `k = 0` | `i_lo_j_lo` is the k-edge `nodes[0, 0, :]` |
| Face | 6 | `i_lo`, `i_hi`, `j_lo`, `j_hi`, `k_lo`, `k_hi` | `k_lo` is `nodes[:, :, 0]` |

Each entity has a `nodes(block)` accessor returning a view of its nodes, and
an `association` slot that links it to geometry (see
[Association](/concepts/association)). The same vertex order (index
$4i + 2j + k$) is used by the blocking skeleton
([Blocking](/concepts/blocking)), so `v0` … `v7` mean the same corners in
journal lines like `p.associate('box', 'v5', ...)`.

Each face has two **along-axes**, the block axes that lie in it, in a fixed
order (`FACE_ALONG_AXES`): `i_*` faces have `(j, k)`, `j_*` faces `(i, k)`,
`k_*` faces `(i, j)`.

## Face boundary conditions

A face carries one `FaceBC`. A new `HexBlock` has `k_lo = WALL` and every
other face `FREE`; the topology builder rewires the rest.

| `FaceBC` | Meaning | `cc.par` (Xall) | CGNS |
|---|---|---|---|
| `WALL` | body surface (`k = 0` of a body-fitted block) | BC `1` (viscous wall) | `BCWall` |
| `FREE` | outer boundary, decided downstream: overset fringe or natural BC | per block `free_face_bc`, default `40` (chimera); a background uses a positive natural code in `[1, 19]` (e.g. `9` far field); a face's `bc_code` overrides | `BCTypeUserDefined`, descriptor `shore:FREE` |
| `SHARED` | the same nodes as a partner face, usually shared in memory | connection patch, 3-digit orientation code from the `AxisMap` | `GridConnectivity1to1_t` |
| `DIRICHLET` | the same nodes as a partner face, held as a copy (the cubed sphere's cap–equator seams, whose axes are permuted) | as `SHARED` | `GridConnectivity1to1_t` |
| `PERIODIC` | wraps to a face of the same block | no native code: `40` with a `UserWarning` | `GridConnectivity1to1_t` |
| `SYMMETRY` | lies on a symmetry plane (half models, [ADR 0013](/adr/0013-half-collar-block)) | natural BC `2` | `BCSymmetryPlane` |

The `cc.par` orientation digits follow the reading of Xall's overset
pre-processor ([ADR 0036](/adr/0036-cc-par-as-xall-reads-it)); the CGNS
`Transform` comes from the same `AxisMap`
([ADR 0035](/adr/0035-cgns-plot3d-export)).

## AxisMap and the wiring invariant

Two faces on a seam need not use the same axes in the same direction. An
`i_lo` face may meet a partner's `j_hi` face whose `i` runs backwards. The
`AxisMap` records, for each along-axis of this face (in `FACE_ALONG_AXES`
order), the partner's axis and whether it is reversed:

```python
from shore.mesh.face import AxisMap

AxisMap.identity("i_hi")                       # along=(("j", False), ("k", False))
AxisMap.from_pairs(("i", True), ("k", False))  # j -> reversed i, k -> k
```

<Figure name="axis-map" ext="svg" caption="Two blocks sharing a face: this face's along-axes mapped onto the partner's, one of them reversed." />

**Wiring invariant.** `HexBlock._set_face_bc(name, bc, partner, axis_map)`
is the one place a face is wired. For `SHARED`, `DIRICHLET` and `PERIODIC` it
requires a partner and an `AxisMap`, and checks the map's partner axes
against the partner face's own along-axes; `WALL`, `FREE` and `SYMMETRY`
must have neither. A topology builder that wires a seam wrongly fails at
construction, not later as a wrong `cc.par`:

```python
import numpy as np

from shore.mesh.block import HexBlock
from shore.mesh.face import AxisMap, FaceBC

# one buffer of 9 x 5 x 5 nodes on [0, 2] x [0, 1] x [0, 1]
x, y, z = np.meshgrid(
    np.linspace(0, 2, 9), np.linspace(0, 1, 5), np.linspace(0, 1, 5), indexing="ij"
)
buffer = np.stack([x, y, z], axis=-1)

# two blocks cut from it: node plane i = 4 belongs to both (a view)
a = HexBlock("a", buffer[:5], periodic_j=False)
b = HexBlock("b", buffer[4:], periodic_j=False)
a._set_face_bc("i_hi", FaceBC.SHARED, b.face("i_lo"), AxisMap.identity("i_hi"))
b._set_face_bc("i_lo", FaceBC.SHARED, a.face("i_hi"), AxisMap.identity("i_lo"))

print(np.shares_memory(a.nodes, b.nodes))  # True

# the wrong map is refused at wiring time:
# AxisMap.from_pairs(("i", False), ("k", False)) on i_hi -> i_lo raises
# ValueError: partner axes ('i', 'k') do not match 'i_lo''s along-axes ('j', 'k')
```

`HexBlock` defaults to `periodic_j=True` (the legacy single-block O-grid);
pass `periodic_j=False` for an open block. Both sides of a seam are wired
explicitly, and `AxisMap.inverse` gives the partner side's map.

The face graph (BCs, partners, axis maps, associations) is what
`.adjacency.json` stores (version 3); `cc-par`, the CGNS writer and a
project loaded without replay all read it back.

## Seams: views and write-through

A seam is one set of points seen from two blocks. SHORE keeps it closed in
two ways.

- **Views.** Builders allocate one buffer and cut the blocks from it, so a
  `SHARED` seam is the same memory: writing one side writes the other.
  This works only where the seam's axis map is an identity or a flip in a
  buffer the builder laid out. Permuted seams (the cubed sphere's
  cap–equator edges) are `DIRICHLET` copies, and blocks read back from saved
  grids ([ADR 0033](/adr/0033-load-without-replay)) hold copies of every seam.
- **Point matching.** `shore.mesh.seams.seams(blocks, tol)` pairs each
  connected face once and locates the smaller face inside the larger one by
  its nodes, giving SIDS index ranges and a `Transform` on both sides. This
  is how a face shared with several partners is handled, such as the CH
  hybrid's `body.k_hi` under three channel blocks (`shared_faces` lists such
  faces). Nodes that do not coincide within `tol` raise `ParameterError`.
- **Write-through** ([ADR 0040](/adr/0040-seams-kept-closed)). `Aliases`
  groups the nodes that are one point across `SHARED` and `DIRICHLET` seams
  (a node where several blocks meet is one group). After an in-place edit,
  `sync` gives each group the value of its last member the edit changed, in
  block order, as views would ("last writer wins"); groups the edit left
  alone, or changed alike, are untouched. `Snap`, `FillFace`,
  `Redistribute`, `Orthogonalize` and `WallLayer` run inside this
  write-through, so a one-sided association or a loaded project does not
  open a seam. A seam not point-matched within $10^{-9}$ × the bounding
  diagonal makes the edit warn and run without closing seams.

`PERIODIC` seams do not coincide and are not aliased.

## Spacing per edge

Every `Edge` carries a `Spacing1D` (law `uniform`, `tanh`, `tanh2`,
`geometric` or `pinned`, with its parameters; see
[Spacing laws](/concepts/spacing-laws)). Defaults: `uniform` along `i` and
`j`, `geometric` with `ds = 1e-3`, `growth = 1.15` along `k`.

**Identity rule.** `HexBlock.__init__` gives all four edges of an axis the
*same* `Spacing1D` instance and remembers it. The wall-normal resolver
treats a k-edge whose `spacing is not` that instance as a user override:

- no k-edge replaced: the march uses the `(ds, growth)` it is given;
- all four replaced and equal by value: that law sets the k-steps;
- some replaced, or replaced with different values: `UserWarning`, default
  kept.

All blocks of one `Mesh` share a single k-step schedule, because k-coupled
seams must stay coincident layer by layer. Block 0 wins, with a
`UserWarning` if another block asks for a different one. Editing an existing
grid's spacing is a separate command, `redistribute`
([ADR 0041](/adr/0041-redistribute-edge)), which resamples the lines of an
axis and carries the change through the seams.

`Edge.first_step` / `last_step` request a C1 match at an end of the edge;
builders read them, and `Spacing1D(law="pinned")` is the law that honours
them.

## See also

- [ADR 0003](/adr/0003-association) (associations on vertices, edges, faces),
  [ADR 0035](/adr/0035-cgns-plot3d-export), [ADR 0036](/adr/0036-cc-par-as-xall-reads-it),
  [ADR 0040](/adr/0040-seams-kept-closed), [ADR 0041](/adr/0041-redistribute-edge)
- API: [`shore.mesh.block`](/reference/api/shore.mesh.block),
  [`shore.mesh.face`](/reference/api/shore.mesh.face),
  [`shore.mesh.edge`](/reference/api/shore.mesh.edge),
  [`shore.mesh.vertex`](/reference/api/shore.mesh.vertex),
  [`shore.mesh.seams`](/reference/api/shore.mesh.seams),
  [`shore.mesh.spacing`](/reference/api/shore.mesh.spacing)
- [Per-edge spacing](/cookbook/per-edge-spacing), [Redistribute an edge](/cookbook/redistribute)
- [Association](/concepts/association), [Grid quality](/concepts/quality)
