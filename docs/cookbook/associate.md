# Associate a block with geometry and snap

Use this when block nodes must lie on a CAD face, edge or corner (or on a whole
body): you record which geometry entity each block face, edge or vertex
belongs to, then `snap` projects the nodes onto it.

The script needs the `cad` extra and runs in about 3 s. Associations are a
Python-API and GUI feature; the CLI has no `associate` command.

## Find the geometry entities

<<< @/../examples/cookbook/associate.py#find_entities{python}

```text
surface 2, curve 5, point 1
```

- A CAD geometry exposes SURFACE, CURVE and POINT entities, indexed from 0 in
  `geometry.entities(kind)`. Indices follow the file's topology order, so look
  them up rather than hard-coding them: here by the entity closest to a probe
  point.
- An STL geometry has one BODY and, since [ADR 0049](/adr/0049-stl-feature-curves),
  feature CURVE and POINT entities; see [Feature curves of an STL](/cookbook/stl-curves).

## Associate a face, an edge and a vertex, then snap

<<< @/../examples/cookbook/associate.py#associate_snap{python}

```text
p.associate('box', 'v0', 'part', 'point', (1,))
v0         : [-1.  -0.5  0.5]
edge y, z  : [[-0.5  0.5]]
face z     : [0.5]
k_lo off the surface by at most 6.9e-15
```

- The entity is a face name (`k_lo`), an edge name (`j_lo_k_lo`) or a vertex
  `v0`..`v7` in `HexBlock.vertices()` order; see the
  [block model](/concepts/block-model).
- `kind` is `"body"`, `"surface"`, `"curve"`, `"point"`, or `"none"` to clear.
  An entity takes geometry of its own dimension or higher: a face a surface or
  the body, an edge also a curve, a vertex also a point.
- `associate` moves no node; `snap()` does. Where a node belongs to several
  associated entities, the lowest dimension wins: vertex, then edge, then
  face. Here `v0` goes to the corner (-1, -0.5, 0.5), the rest of the edge to
  the curve, the rest of the face to z = 0.5.
- The association stores fingerprints of the entities, not only their
  indices, so it survives a re-import of a modified CAD file
  ([ADR 0003](/adr/0003-association)).
- `snap(exact=True)` uses the exact projection on the CAD surfaces instead of
  their tessellation; on an STL both are the same.
- `snap` resolves every project geometry by name (`"part"` here). On a
  topology built from a geometry, every block's WALL face (`k_lo`, with its
  edges and corners) is already associated with the body, so `p.snap()` after
  an edit puts the wall back on it.
- Seam nodes shared with a neighbour are written through, so SHARED and
  DIRICHLET seams stay closed ([ADR 0040](/adr/0040-seams-kept-closed)).

## Clear or undo an association

<<< @/../examples/cookbook/associate.py#clear{python}

```text
None
True
```

- `kind="none"` clears; `undo()` restores the previous association.
- Associations are saved with the project and come back on load, with or
  without replay ([Save, load, replay](/cookbook/projects)).
- Refused with `ParameterError`, leaving project and journal unchanged: an
  unknown block, entity, geometry or kind; an index out of range; no index for
  a surface, curve or point; an index for `"body"`; a geometry entity of lower
  dimension than the block entity (a face to a curve).

## See also

- [Association](/concepts/association) and [Geometry](/concepts/geometry)
- [Rebuild a face on a CAD surface](/cookbook/cad-face)
- [`shore.mesh.association`](/reference/api/shore.mesh.association)
- [ADR 0003](/adr/0003-association), [ADR 0039](/adr/0039-association-editing)
