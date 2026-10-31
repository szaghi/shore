---
title: shore.mesh.blocking
outline: [2, 3]
---

# `shore.mesh.blocking`

Blocking: a block skeleton edited by commands, premeshed into blocks (ADR 0042).

A blocking is the topology of a multi-block grid without its nodes: skeleton
vertices (points), hexahedral blocks over them (8 vertex ids each, in
`HexBlock.vertices()` order: index `4 i + 2 j + k`), a node count on
every edge, and associations of vertices, edges and faces with geometry.
Blocks touching each other share vertex ids, so a face or an edge is the
same entity in every block that has it.

Editing keeps the blocking conforming:
- `split` cuts a block across one axis and carries the cut through every
  block sharing a cut edge (the whole sheet), at the same node;
- `delete` removes blocks (the hole they leave is a boundary).

`premesh` makes ordinary `HexBlock`s:
- vertices are projected onto their associations;
- each edge is a straight line of equally spaced nodes, projected onto its
  association and re-equidistributed along the projected polyline;
- each face is a TFI of its four edges, projected onto its association;
- each block is a TFI of its six faces.
Edges and faces are made once and copied into every block that has them,
so seams match bit for bit. Faces shared by two blocks become SHARED
seams; boundary faces with an association are WALL, the others FREE.

An entity without its own association takes the one of the entities above
it, lowest dimension first, when they agree: an edge the association of its
associated faces, a vertex that of its associated edges, else of its faces.

## `face_corners` <Badge type="info" text="function" /> {#face_corners}

```python
face_corners(face: str) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L60">source</a>

Local vertex indices (2, 2) of *face*, over its along-axes (FACE_ALONG_AXES).

## `Blocking` <Badge type="tip" text="dataclass" /> {#Blocking}

```python
Blocking(
    points: np.ndarray,
    blocks: list[tuple[str, tuple[int, ...]]],
    nodes: dict[tuple[int, int], int],
    associations: dict[EntityKey, Association] = field(default_factory=dict),
    next_label: int = 1,
    planes: tuple[Plane, ...] = (),
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L75">source</a>

A block skeleton (see the module docstring).

### `Blocking.copy` <Badge type="info" text="method" /> {#Blocking-copy}

```python
copy() -> Blocking
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L85">source</a>

### `Blocking.labels` <Badge type="info" text="method" /> {#Blocking-labels}

```python
labels() -> list[str]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L96">source</a>

### `Blocking.block` <Badge type="info" text="method" /> {#Blocking-block}

```python
block(label: str) -> tuple[int, ...]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L99">source</a>

### `Blocking.edge_count` <Badge type="info" text="method" /> {#Blocking-edge_count}

```python
edge_count(a: int, b: int) -> int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L105">source</a>

### `Blocking.entity` <Badge type="info" text="method" /> {#Blocking-entity}

```python
entity(label: str, name: str) -> tuple[EntityKey, str]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L108">source</a>

The key of face, edge or vertex *name* of block *label*, and which it is.

### `Blocking.shape` <Badge type="info" text="method" /> {#Blocking-shape}

```python
shape(label: str) -> tuple[int, int, int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L123">source</a>

### `Blocking.edges` <Badge type="info" text="method" /> {#Blocking-edges}

```python
edges() -> dict[tuple[int, int], list[tuple[int, str]]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L129">source</a>

Every edge (a &lt; b) -&gt; the (block index, edge name) that have it.

### `Blocking.polylines` <Badge type="info" text="method" /> {#Blocking-polylines}

```python
polylines() -> list[tuple[str, str, np.ndarray]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L138">source</a>

Every edge once, as `(block, edge name, nodes (n, 3))`: its nodes
equally spaced on the straight line between its skeleton vertices, in
the direction of the edge in the first block that has it.

### `Blocking.centres` <Badge type="info" text="method" /> {#Blocking-centres}

```python
centres(kind: str) -> list[tuple[str, str, np.ndarray]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L153">source</a>

Skeleton centres of every `"face"` (once, named by its first block),
`"block"` or `"vertex"` (its point, once, named `v0`..`v7` by
its first block): `(block, name, centre)`, the name `""` for a block.

### `Blocking.faces` <Badge type="info" text="method" /> {#Blocking-faces}

```python
faces() -> dict[tuple[int, ...], list[tuple[int, str]]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L171">source</a>

Every face (sorted vertex ids) -&gt; the (block index, face name) that have it.

## `box` <Badge type="info" text="function" /> {#box}

```python
box(vmin: Any, vmax: Any, ni: int, nj: int, nk: int, label: str = 'b0') -> Blocking
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L184">source</a>

One block between the corners *vmin* and *vmax*, with *ni* x *nj* x *nk* nodes.

## `split` <Badge type="info" text="function" /> {#split}

```python
split(blocking: Blocking, label: str, axis: str, fraction: float) -> Blocking
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L199">source</a>

Cut block *label* across *axis* at the node nearest *fraction* of it.

Every block sharing a cut edge is cut too, at the same node, so the
blocking stays conforming. The low part keeps the label, the high part
gets a new one (`b<n>`).

## `ogrid` <Badge type="info" text="function" /> {#ogrid}

```python
ogrid(
    blocking: Blocking,
    labels: list[str] | tuple[str, ...] = (),
    faces: list[tuple[str, str]] | tuple[tuple[str, str], ...] = (),
    exclude: list[tuple[str, str]] | tuple[tuple[str, str], ...] = (),
    factor: float = 0.5,
    nodes: int = 5,
) -> Blocking
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L304">source</a>

Insert an O-grid: a layer of blocks on faces of the block set *labels*.

The layered faces are *faces* (`(block, face)` pairs), or every face on
the boundary of the set, less *exclude*; they all lie on the set's
boundary. Boundary faces of the set that are not layered stay open: the
O-grid passes through them, and they must be faces of the blocking's
boundary. Every vertex of a layered face gets an inner copy, moved into
the set by `(1 - factor)` times half the sum of its edges that leave
the layered faces into the set (inside one block: the block scaled by
*factor* about its centre). The set's blocks take the inner copies; each
layered face gets a block between it and its copy, with *nodes* nodes
across. Labels: the set's blocks keep theirs, the layer blocks are new
`b<n>` in block and face order. With no *labels*, the set is every
block.

## `move` <Badge type="info" text="function" /> {#move}

```python
move(
    blocking: Blocking,
    label: str,
    vertex: str,
    xyz: Any,
    resolver: AssociationResolver | None = None,
) -> Blocking
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L452">source</a>

Move vertex *vertex* (`v0`..`v7`) of block *label* to *xyz*: onto
the projection of *xyz* on the vertex's association (own or inherited),
when it has one, so an associated vertex slides on its geometry.

## `delete` <Badge type="info" text="function" /> {#delete}

```python
delete(blocking: Blocking, labels: list[str] | tuple[str, ...]) -> Blocking
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L478">source</a>

Remove blocks *labels*; their faces left alone become boundary faces.

## `associate` <Badge type="info" text="function" /> {#associate}

```python
associate(
    blocking: Blocking,
    label: str,
    entity: str,
    association: Association | None,
) -> Blocking
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L490">source</a>

Set (or with None clear) the association of *entity* of block *label*.

## `symmetry` <Badge type="info" text="function" /> {#symmetry}

```python
symmetry(blocking: Blocking, point: Any, normal: Any) -> Blocking
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L506">source</a>

A copy of *blocking* with the symmetry plane through *point* (ADR 0053).

`premesh` makes the unassociated boundary faces lying on the plane
`SYMMETRY` faces, and `orthogonal_walls` keeps the lines leaving the
wall from the plane in it.

## `on_plane` <Badge type="info" text="function" /> {#on_plane}

```python
on_plane(points: np.ndarray, plane: Plane, tol: float) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L528">source</a>

Which of *points* (..., 3) lie within *tol* of *plane*.

## `plane_tolerance` <Badge type="info" text="function" /> {#plane_tolerance}

```python
plane_tolerance(points: np.ndarray) -> float
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L534">source</a>

1e-9 of the extent of *points*: their distance to a plane they lie on.

## `effective_associations` <Badge type="info" text="function" /> {#effective_associations}

```python
effective_associations(blocking: Blocking) -> dict[EntityKey, Association]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L540">source</a>

Every entity's association: its own, else inherited (module docstring).

## `tfi_face` <Badge type="info" text="function" /> {#tfi_face}

```python
tfi_face(
    i_lo: np.ndarray,
    i_hi: np.ndarray,
    j_lo: np.ndarray,
    j_hi: np.ndarray,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L607">source</a>

Coons patch (na, nb, 3) from its four boundary rows.

*i_lo* / *i_hi* are the rows at a = 0 / a = na - 1 (length nb), *j_lo* /
*j_hi* those at b = 0 / b = nb - 1 (length na). The boundary is copied.

## `tfi_volume` <Badge type="info" text="function" /> {#tfi_volume}

```python
tfi_volume(x: np.ndarray) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L625">source</a>

Fill the interior of *x* (ni, nj, nk, 3) by TFI of its six faces, in place.

## `premesh` <Badge type="info" text="function" /> {#premesh}

```python
premesh(blocking: Blocking, resolver: AssociationResolver | None) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/blocking.py#L652">source</a>

The blocking's blocks as `HexBlock`s (module docstring).
