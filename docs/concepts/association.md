---
title: Association
---

# Association

An association links a block entity (a face, an edge or a vertex) to
geometry: the whole body, a set of CAD surfaces, a chain of curves, or a
point. The nodes of an associated entity are projected onto that geometry
whenever SHORE snaps the grid. It is the abstraction of ICEM CFD Hexa,
Salome HexaBlock and Pointwise, where an edge pinned to a leading-edge curve
stays on it however the grid is edited
([ADR 0003](/adr/0003-association), [ADR 0039](/adr/0039-association-editing)).

## The value

`shore.mesh.association.Association` is a frozen value:

| Field | Meaning |
|---|---|
| `geometry` | the geometry's name in the project (`"body"` by default) |
| `kind` | `BODY`, `SURFACE`, `CURVE` or `POINT` |
| `keys` | the `EntityKey` fingerprints of the entities (empty for `BODY`) |

A `SURFACE` association with several keys is a surface set; a `CURVE`
association with several keys is a chain. Associations name entities by
fingerprint, not by index, so they are re-found after the geometry is
re-read or edited ([Geometry: fingerprints](/concepts/geometry#fingerprints-and-matching)).

`Vertex`, `Edge` and `Face` each have an `association` slot (see
[Block model](/concepts/block-model)). The legacy flag `body_projected` is a
read-only property derived from it.

## Dimension rule

An entity may take geometry of its own dimension or higher:

| Block entity | May be associated with |
|---|---|
| face | the body, surfaces |
| edge | the body, surfaces, curves |
| vertex | the body, surfaces, curves, a point |

Associating a face with a curve would collapse it onto a line; the command
refuses it and the journal is unchanged.

## Snapping

`AssociationResolver(geometries, exact=False)` resolves each association to
the current entities (`match_keys`; a lost key raises `MeshInputError`
"no longer resolves; re-associate the entity") and projects points onto
them. For a set of entities, each point goes to the nearest one.

`snap_block(block, resolver)` projects in the order faces → edges →
vertices, so the lower-dimensional association wins at shared nodes: corners
pin edges, and edges pin face boundaries (ICEM's projection order).

`p.snap(exact=False)` runs this on every block of the current component:

- `exact=False` projects on the tier-1 tessellation of a CAD model;
  `exact=True` uses the exact geometry (tier 2). On STL both are exact.
- It resolves every project geometry by name. `"body"` means the mesh's own
  body; on a mesh loaded without replay, which has no `Body`, it means the
  single geometry its component was built on.
- It runs inside the seams' write-through
  ([ADR 0040](/adr/0040-seams-kept-closed)), so an association on one side of
  a seam moves the partner's nodes too.
- It is undoable: the node arrays are snapshot to the on-disk undo cache.

Snapping is closest-point projection, nothing more: it does not re-space the
nodes along the entity. Nodes of a long edge snapped to a short curve pile up
at the curve's ends. Use [`redistribute`](/cookbook/redistribute) to re-space
an edge after it is on its curve.

```python
import numpy as np

from shore.exceptions import ParameterError
from shore.project import Project
from shore.shapes import box

box((2.0, 1.0, 0.5)).export("part.stl")  # a 2 x 1 x 0.5 box centred at the origin

p = Project()
p.load_geometry("part.stl", name="part")
p.build_primitive("box", vmin=(-1.2, -0.7, -0.4), vmax=(1.2, 0.7, 0.4), ni=9, nj=7, nk=5)

# the STL's feature curve along z at x = -1, y = -0.5, and its lower corner
curves, points = p.geometries["part"].geometry.features()
c = next(i for i, xyz in enumerate(curves) if np.allclose(xyz[:, :2], [-1.0, -0.5]))
v = int(np.argmin(np.linalg.norm(points - [-1.0, -0.5, -0.25], axis=1)))

p.associate("box", "i_lo_j_lo", "part", "curve", (c,))  # the block's k-edge at i = j = 0
p.associate("box", "v0", "part", "point", (v,))  # its corner (0, 0, 0)
p.snap()
print(p.mesh.blocks[0].nodes[0, 0, :])  # the k-edge, now on the box edge x=-1, y=-0.5

try:
    p.associate("box", "k_lo", "part", "curve", (c,))  # a face onto a curve
except ParameterError as err:
    print(err)  # associate: a face cannot be associated with a curve (...)
```

The edge's end nodes (z = ±0.4) land on the curve's ends (z = ±0.25); the
interior nodes keep their z.

## Where associations come from

- **Topology builders.** A `Mesh` built on a body associates every block's
  `WALL` `k_lo` face, its four edges and four corners with the whole body
  (`associate_wall`, `Association.body()`). C-grid wake cuts and channel
  faces are `SHARED`, not walls, so they stay unassociated. Primitives and
  collars carry no associations.
- **The `associate` command.** `p.associate(block, entity, geometry="body",
  kind="body", indices=())` names the block entity as a face (`"k_lo"`), an
  edge (`"i_lo_k_lo"`) or a vertex (`"v0"` … `"v7"`), and the geometry
  entities by their current indices; it stores their fingerprints.
  `kind="none"` clears. Undo restores the previous association; nodes move
  only on `snap`. In the GUI the Association dock issues the same command.
- **Splits.** `split_blocks` gives each child entity the association of the
  lowest-dimensional source entity it lies on (vertex, then edge, then face);
  internal cut faces are unassociated.
- **Blockings.** `p.associate_blocking(block, entity, ...)` associates
  entities of the blocking skeleton ([Blocking](/concepts/blocking)), with the
  same validation. `premesh` projects them and copies them onto the
  `HexBlock`s it makes.

## Inheritance in a blocking

In a blocking ([ADR 0042](/adr/0042-blocking-skeleton-premesh)), an entity
without an association of its own takes one from the entities above it,
lowest dimension first, **when they agree**:

- an edge takes the association of its associated faces;
- a vertex takes that of its associated edges, else of its associated faces.

Associating the six faces of a hole around a body is then enough: its edges
and corners follow. An edge between two faces associated with *different*
geometry inherits nothing and stays straight until it gets an association
of its own, because the nearer of two surfaces is not their intersection
curve. Inheritance is evaluated at premesh time
(`effective_associations`) from the current state, not copied at split time,
so changing a face's association changes what its edges inherit.

## How associations are stored

- **In the journal** as the commands that set them:
  `p.associate('box', 'i_lo_j_lo', 'part', 'curve', (3,))`. A replay re-runs
  the command with those indices against the file the journal reads, and
  stores fresh fingerprints. The fingerprints matter when a geometry is
  replaced under its name (`load_geometry` again, with an edited file) while
  the blocks that name it stay: the next `snap` re-finds each entity by
  `match_keys`, or fails naming the lost ones.
- **In the face graph** `adjacency.json`, version 3: a deduplicated
  top-level `associations` table (`geometry`, `kind`, `keys`), and per block
  an `associated` map from face names, edge names and vertex indices to rows
  of that table. Version-2 sidecars still read. A project saves one face
  graph per component (`grids/<component>/adjacency.json`), so a project
  loaded without replay keeps its associations and can snap.
- **Not stored:** the blocking itself. It is rebuilt by replaying the
  journal; a project loaded without replay has its premeshed blocks (with
  their associations) but no blocking.

## Limits

- Fingerprint matching fails on large topological edits (a face split in
  CAD) and on large shape changes (+30 % span on a test wing lost some
  keys): re-associate the lost entities.
- Snap does not redistribute, and closest-point projection of a face onto a
  surface set picks the nearest surface per node, which is not an
  intersection curve: give such edges a curve association.
- The C-grid's trailing-edge columns are not yet associated with the TE
  curve (that needs a CAD-native TE detector).

## See also

- [ADR 0003](/adr/0003-association), [ADR 0039](/adr/0039-association-editing),
  [ADR 0040](/adr/0040-seams-kept-closed), [ADR 0042](/adr/0042-blocking-skeleton-premesh),
  [ADR 0049](/adr/0049-stl-feature-curves)
- API: [`Association`](/reference/api/shore.mesh.association#Association),
  [`AssociationResolver`](/reference/api/shore.mesh.association#AssociationResolver),
  [`snap_block`](/reference/api/shore.mesh.association#snap_block),
  [`effective_associations`](/reference/api/shore.mesh.blocking#effective_associations),
  [`Project`](/reference/api/shore.project.project#Project)
- [Associate with geometry](/cookbook/associate), [Blocking a wing–body](/tutorials/blocking),
  [Project scripts](/reference/run)
- [Geometry](/concepts/geometry), [Block model](/concepts/block-model)
