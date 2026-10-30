---
title: shore.project.commands
outline: [2, 3]
---

# `shore.project.commands`

Commands: the only way to change a `Project` (ADR 0004).

A command is a frozen dataclass. `apply` changes the project and returns an
`Undo` (or `None` for output-only commands such as `WriteGeo`);
`to_source` renders the command as the one line of public Python API that
re-issues it, which is what the journal records.

Node-array changes (march, snap) are undone from snapshots that the
`UndoCache` keeps on disk, not in RAM: a 1e7-node grid is 240 MB.

A project holds named geometries and named components (ADR 0016); mesh
commands act on the current component. Commands that read files name those
fields in `inputs`: `Project.save` copies the files into the project, and
a replay resolves relative input paths against `input_root` (the saved
project's directory).

## `resolve_input` <Badge type="info" text="function" /> {#resolve_input}

```python
resolve_input(path: str) -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L72">source</a>

*path*, relative to `input_root` when that is set and *path* is relative.

## `Cancelled` <Badge type="tip" text="class" /> {#Cancelled}

```python
class Cancelled(ShoreError)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L97">source</a>

Raised inside a long command when its `CancelToken` was set.

## `CancelToken` <Badge type="tip" text="class" /> {#CancelToken}

```python
CancelToken()
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L101">source</a>

Cooperative cancellation flag, checked once per layer / sweep.

### `CancelToken.cancel` <Badge type="info" text="method" /> {#CancelToken-cancel}

```python
cancel() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L107">source</a>

### `CancelToken.check` <Badge type="info" text="method" /> {#CancelToken-check}

```python
check() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L110">source</a>

## `Context` <Badge type="tip" text="dataclass" /> {#Context}

```python
Context(progress: ProgressFn | None = None, cancel: CancelToken | None = None)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L116">source</a>

Per-execution hooks: progress reporting and cancellation.

### `Context.report` <Badge type="info" text="method" /> {#Context-report}

```python
report(fraction: float, message: str = '') -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L122">source</a>

### `Context.check` <Badge type="info" text="method" /> {#Context-check}

```python
check() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L126">source</a>

## `Undo` <Badge type="tip" text="dataclass" /> {#Undo}

```python
Undo(restore: Callable[[], None])
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L132">source</a>

How to revert one applied command.

## `UndoCache` <Badge type="tip" text="class" /> {#UndoCache}

```python
UndoCache(root: Path | None = None)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L138">source</a>

On-disk store of node-array snapshots, deleted with the project.

### `UndoCache.snapshot` <Badge type="info" text="method" /> {#UndoCache-snapshot}

```python
snapshot(mesh: Mesh) -> Callable[[], None]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L145">source</a>

Save every block's nodes; return a function restoring them in place.

### `UndoCache.close` <Badge type="info" text="method" /> {#UndoCache-close}

```python
close() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L159">source</a>

## `Command` <Badge type="tip" text="class" /> {#Command}

```python
class Command(ABC)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L188">source</a>

Base of all commands. Subclasses are frozen dataclasses.

### `Command.apply` <Badge type="info" text="method" /> {#Command-apply}

```python
apply(project: Project, ctx: Context) -> Undo | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L197">source</a>

### `Command.to_source` <Badge type="info" text="method" /> {#Command-to_source}

```python
to_source() -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L199">source</a>

`p.<method>(...)` with positional fields and non-default keywords.

## `LoadGeometry` <Badge type="tip" text="dataclass" /> {#LoadGeometry}

```python
LoadGeometry(
    path: str,
    anchor_strategy: str = 'auto',
    skip_hygiene: bool = False,
    name: str = 'body',
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L215">source</a>

Load a geometry (STL or STEP/IGES/BREP) under *name*.

Reloading a name replaces that geometry and discards the meshes of the
components built on it (ADR 0016).

### `LoadGeometry.apply` <Badge type="info" text="method" /> {#LoadGeometry-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L231">source</a>

## `UseComponent` <Badge type="tip" text="dataclass" /> {#UseComponent}

```python
UseComponent(name: str)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L250">source</a>

Make *name* the current component, creating it empty if new (ADR 0016).

### `UseComponent.apply` <Badge type="info" text="method" /> {#UseComponent-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L258">source</a>

## `BuildTopology` <Badge type="tip" text="dataclass" /> {#BuildTopology}

```python
BuildTopology(
    topology: str,
    nk: int,
    ds: float = 0.001,
    growth: float = 1.15,
    k_law: str = 'geometric',
    k_beta: float = 3.0,
    k_thickness: float | None = None,
    options: dict[str, Any] = field(default_factory=dict),
    params: dict[str, Any] = field(default_factory=dict),
    geometry: str = 'body',
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L272">source</a>

Build the topology's blocks on the body (`shore mesh` stage 2).

The k-spacing fields define the wall-normal law used both by the
topology (seam pinning) and by `March`. *options* are the topology
constructor's arguments (`ch`, `oh`); *params* go to its `build`.

### `BuildTopology.to_source` <Badge type="info" text="method" /> {#BuildTopology-to_source}

```python
to_source() -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L294">source</a>

### `BuildTopology.apply` <Badge type="info" text="method" /> {#BuildTopology-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L298">source</a>

## `BuildState` <Badge type="tip" text="dataclass" /> {#BuildState}

```python
BuildState(
    topology: str,
    nk: int,
    ds: float = 0.0,
    growth: float = 0.0,
    spacing: Any = None,
    writer: Callable[[Path], None] | None = None,
    marchable: bool = True,
    geometries: tuple[str, ...] = (),
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L324">source</a>

What later commands need from the command that made the mesh.

## `March` <Badge type="tip" text="dataclass" /> {#March}

```python
March(
    smooth_strength: float = 0.2,
    smooth_iters: int = 2,
    blend_normals_k: int = 0,
    smooth_mode: str = 'limited',
    march_scheme: str = 'algebraic',
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L345">source</a>

Extrude all blocks from the wall (`shore mesh` stage 3).

### `March.apply` <Badge type="info" text="method" /> {#March-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L356">source</a>

## `Snap` <Badge type="tip" text="dataclass" /> {#Snap}

```python
Snap(exact: bool = False)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L390">source</a>

Project associated nodes onto their geometry (ADR 0003).

### `Snap.apply` <Badge type="info" text="method" /> {#Snap-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L397">source</a>

## `Associate` <Badge type="tip" text="dataclass" /> {#Associate}

```python
Associate(
    block: str,
    entity: str,
    geometry: str = 'body',
    kind: str = 'body',
    indices: tuple[int, ...] = (),
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L470">source</a>

Associate a block face, edge or vertex with geometry (ADR 0003, 0039).

*entity* is a face (`"k_lo"`), an edge (`"i_lo_k_lo"`) or a vertex
(`"v0"` .. `"v7"`, `HexBlock.vertices()` order) of block *block* of
the current component. *kind* is `"body"`, `"surface"`, `"curve"`
or `"point"` with entity *indices* of geometry *geometry* (none for
`"body"`), or `"none"` to clear. An entity may be associated with
geometry of its own dimension or higher: a face with surfaces or the
body, an edge also with curves, a vertex also with a point. Nodes move
only on `snap`.

### `Associate.to_source` <Badge type="info" text="method" /> {#Associate-to_source}

```python
to_source() -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L492">source</a>

### `Associate.apply` <Badge type="info" text="method" /> {#Associate-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L496">source</a>

## `Redistribute` <Badge type="tip" text="dataclass" /> {#Redistribute}

```python
Redistribute(
    block: str,
    edge: str,
    law: str = 'uniform',
    beta: float = 3.0,
    cluster_at: float | None = None,
    first_step: float | None = None,
    last_step: float | None = None,
    interior: str = 'uniform',
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L571">source</a>

Redistribute the nodes along a block edge's axis (ADR 0041).

The spacing law (`uniform`, `tanh`, `tanh2` with *beta* and
*cluster_at*, `pinned` with *first_step* / *last_step* lengths and an
*interior* law) is taken over the arc length of *edge* of *block* (current
component) and applied, as arc-length fractions, to every grid line of the
block along that axis, and through SHARED / DIRICHLET seams to the blocks
sharing a face that contains the axis. Nodes stay on their lines' old
polylines; line ends do not move. Associated nodes are not re-snapped:
run `snap` after.

### `Redistribute.apply` <Badge type="info" text="method" /> {#Redistribute-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L596">source</a>

## `Orthogonalize` <Badge type="tip" text="dataclass" /> {#Orthogonalize}

```python
Orthogonalize(strength: float = 1.0)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L620">source</a>

Make the grid lines leave the walls along their normal (ADR 0048).

Every block of the current mesh with exactly one WALL face has its lines
across the wall reshaped (`shore.mesh.orthogonal`, tangent *strength*
times the chord); blocks whose faces the seams change are re-filled by
TFI. Refused, with nothing changed, if a block that was valid would be
inverted (lower *strength*).

### `Orthogonalize.apply` <Badge type="info" text="method" /> {#Orthogonalize-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L634">source</a>

## `WallLayer` <Badge type="tip" text="dataclass" /> {#WallLayer}

```python
WallLayer(first_step: float)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L678">source</a>

Cluster the nodes towards the walls (ADR 0050): on every line across
the wall of each block with one WALL face, geometric cells from the wall
with the first *first_step* long (an absolute length), on the line's own
polyline (after `orthogonalize`: along the wall normal). Refused if a
line cannot hold it or a valid block would invert.

### `WallLayer.apply` <Badge type="info" text="method" /> {#WallLayer-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L690">source</a>

## `FillFace` <Badge type="tip" text="dataclass" /> {#FillFace}

```python
FillFace(
    block: str,
    face: str,
    surface: int,
    smoothing: str = 'elliptic',
    geometry: str = 'body',
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L702">source</a>

Rebuild one block face from its boundary on a CAD surface (ADR 0007).

The face's four boundary rows stay as they are; the interior is filled by
`shore.surface.tfi.surface_grid` on SURFACE entity number *surface* of
the loaded geometry.

### `FillFace.apply` <Badge type="info" text="method" /> {#FillFace-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L719">source</a>

## `CreateBlocking` <Badge type="tip" text="dataclass" /> {#CreateBlocking}

```python
CreateBlocking(
    ni: int,
    nj: int,
    nk: int,
    geometry: str = 'body',
    margin: float = 0.0,
    vmin: tuple[float, float, float] | None = None,
    vmax: tuple[float, float, float] | None = None,
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L747">source</a>

Start the current component's blocking with one block (ADR 0042).

The block is the bounding box of geometry *geometry*, each side moved out
by *margin* times the box's extent along that axis, or the box between
*vmin* and *vmax*; it has *ni* x *nj* x *nk* nodes. Replaces the
component's blocking; its mesh changes only on `premesh`.

### `CreateBlocking.apply` <Badge type="info" text="method" /> {#CreateBlocking-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L767">source</a>

## `SplitBlocking` <Badge type="tip" text="dataclass" /> {#SplitBlocking}

```python
SplitBlocking(block: str, axis: str, fraction: float)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L795">source</a>

Cut block *block* of the blocking across *axis* at the node nearest
*fraction* of it, and every block sharing a cut edge with it (ADR 0042).
The low part keeps the label; the high part is `b<n>`, *n* counting up.

### `SplitBlocking.apply` <Badge type="info" text="method" /> {#SplitBlocking-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L807">source</a>

## `DeleteBlocks` <Badge type="tip" text="dataclass" /> {#DeleteBlocks}

```python
DeleteBlocks(blocks: tuple[str, ...])
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L815">source</a>

Remove blocks *blocks* from the blocking; the faces they leave become
boundary faces (ADR 0042).

### `DeleteBlocks.apply` <Badge type="info" text="method" /> {#DeleteBlocks-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L824">source</a>

## `Ogrid` <Badge type="tip" text="dataclass" /> {#Ogrid}

```python
Ogrid(
    blocks: tuple[str, ...] = (),
    faces: tuple[tuple[str, str], ...] = (),
    exclude: tuple[tuple[str, str], ...] = (),
    factor: float = 0.5,
    nodes: int = 5,
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L832">source</a>

Insert an O-grid in the blocking (ADR 0043): a layer of blocks on the
faces of the block set *blocks* (every block when empty), *faces* or all
the set's boundary faces less *exclude* (`(block, face)` pairs); see
`shore.mesh.blocking.ogrid`. *factor* sizes it, *nodes* across.

### `Ogrid.apply` <Badge type="info" text="method" /> {#Ogrid-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L847">source</a>

## `MoveVertex` <Badge type="tip" text="dataclass" /> {#MoveVertex}

```python
MoveVertex(block: str, vertex: str, xyz: tuple[float, float, float])
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L858">source</a>

Move a blocking vertex (ADR 0045): vertex *vertex* (`v0`..`v7`) of
block *block* to *xyz*, projected onto its association (own or inherited)
when it has one, so an associated vertex slides on its geometry.

### `MoveVertex.apply` <Badge type="info" text="method" /> {#MoveVertex-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L870">source</a>

## `AssociateBlocking` <Badge type="tip" text="dataclass" /> {#AssociateBlocking}

```python
AssociateBlocking(
    block: str,
    entity: str,
    geometry: str = 'body',
    kind: str = 'body',
    indices: tuple[int, ...] = (),
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L881">source</a>

Associate a face, edge or vertex of a blocking block with geometry
(ADR 0042): as `associate` (ADR 0039), on the blocking. Entities without
an association of their own inherit one on `premesh`.

### `AssociateBlocking.to_source` <Badge type="info" text="method" /> {#AssociateBlocking-to_source}

```python
to_source() -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L895">source</a>

### `AssociateBlocking.apply` <Badge type="info" text="method" /> {#AssociateBlocking-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L899">source</a>

## `Premesh` <Badge type="tip" text="dataclass" /> {#Premesh}

```python
Premesh(exact: bool = False)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L911">source</a>

Make the current component's mesh from its blocking (ADR 0042):
edges, faces and blocks by TFI, projected onto their associations
(*exact*: tier-2 CAD projection). Replaces the mesh.

### `Premesh.apply` <Badge type="info" text="method" /> {#Premesh-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L920">source</a>

## `WriteGeo` <Badge type="tip" text="dataclass" /> {#WriteGeo}

```python
WriteGeo(stem: str, per_block: bool = False)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L941">source</a>

Write `.geo` files (output only).

Naming: a primitive's own rule (`shore primitive`), unless
*per_block*; otherwise `shore mesh`'s rule: one `<stem>.geo` for a
single non-cubed-sphere block, else `<stem>_<label>.geo` per block
(always the latter with *per_block*, as `shore split` and
`shore cap-tip` write).

### `WriteGeo.apply` <Badge type="info" text="method" /> {#WriteGeo-apply}

```python
apply(project: Project, ctx: Context) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L958">source</a>

## `WriteAdjacency` <Badge type="tip" text="dataclass" /> {#WriteAdjacency}

```python
WriteAdjacency(path: str)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L975">source</a>

Write the `.adjacency.json` sidecar (output only).

### `WriteAdjacency.apply` <Badge type="info" text="method" /> {#WriteAdjacency-apply}

```python
apply(project: Project, ctx: Context) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L984">source</a>

## `WriteCgns` <Badge type="tip" text="dataclass" /> {#WriteCgns}

```python
WriteCgns(path: str, components: tuple[str, ...] = ())
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1008">source</a>

Write a CGNS file: a base per component, a zone per block (ADR 0035; output only).

### `WriteCgns.apply` <Badge type="info" text="method" /> {#WriteCgns-apply}

```python
apply(project: Project, ctx: Context) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1018">source</a>

## `WritePlot3d` <Badge type="tip" text="dataclass" /> {#WritePlot3d}

```python
WritePlot3d(path: str, components: tuple[str, ...] = ())
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1027">source</a>

Write a multi-block Plot3D grid file of the components' blocks (ADR 0035; output only).

### `WritePlot3d.apply` <Badge type="info" text="method" /> {#WritePlot3d-apply}

```python
apply(project: Project, ctx: Context) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1037">source</a>

## `BuildPrimitive` <Badge type="tip" text="dataclass" /> {#BuildPrimitive}

```python
BuildPrimitive(kind: str, params: dict[str, Any] = field(default_factory=dict))
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1047">source</a>

Build a parametric primitive (`shore primitive`); replaces the mesh.

### `BuildPrimitive.to_source` <Badge type="info" text="method" /> {#BuildPrimitive-to_source}

```python
to_source() -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1056">source</a>

### `BuildPrimitive.apply` <Badge type="info" text="method" /> {#BuildPrimitive-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1059">source</a>

## `BuildTipCap` <Badge type="tip" text="dataclass" /> {#BuildTipCap}

```python
BuildTipCap(
    body_geo: str,
    tip: str = 'jlo',
    le_frac: float = 0.05,
    te_frac: float = 0.95,
    span_axis: str = 'z',
    tip_inset: float | None = None,
    nk: int | None = None,
    n_cap_j: int = 12,
    cap_j_spacing: str = 'tanh',
    cap_j_beta: float = 4.0,
    auto_round_te: bool = True,
    te_thickness: float | None = None,
    te_thickness_frac: float | None = 0.002,
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1075">source</a>

Wing-tip cap from a C-grid body `.geo` (`shore cap-tip`); replaces the mesh.

### `BuildTipCap.apply` <Badge type="info" text="method" /> {#BuildTipCap-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1096">source</a>

## `BuildCollar` <Badge type="tip" text="dataclass" /> {#BuildCollar}

```python
BuildCollar(
    other: str,
    seed: tuple[float, float, float],
    n_loop: int,
    ds_surface: float,
    growth_surface: float,
    m_a: int,
    m_b: int,
    ds_wall: float,
    growth_wall: float,
    nk: int,
    label: str = 'collar',
    plane: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None,
    geometry: str = 'body',
    loop_law: str = 'uniform',
    loop_beta: float = 3.0,
    loop_cluster_at: float | None = None,
    loop_curvature: float = 0.0,
    split: bool = False,
    wake_length: float | None = None,
    volume: str | None = None,
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1111">source</a>

Collar block around the junction of the loaded body and *other* (ADR 0010).

Body A is the project's geometry, body B is loaded from the file *other*
(its path is what the journal records). Surface steps are
`ds_surface * growth_surface**j` (*m_a* on A, *m_b* on B), wall-normal
steps `ds_wall * growth_wall**k` (*nk* of them). With *plane*
(`(point, normal)`) the half-collar of a half model (ADR 0013). The
nodes along the junction curve follow ``Spacing1D(loop_law,
beta=loop_beta, cluster_at=loop_cluster_at)`` and *loop_curvature*
(ADR 0020). With *split*, `split_collar` (ADR 0024): the loop is cut at
sharp corners, with wakes of *wake_length* (default: the collar's width
on A), six blocks per piece (closed loops only). *volume* is
`CollarParams.volume` (ADR 0025): `"stations"` or `"smoothed"`; `None`,
the default (`"smoothed"`; `"stations"` when a project saved before
format 3 replays, ADR 0026).
Replaces the mesh with the collar's block(s), a separate Chimera
component.

### `BuildCollar.apply` <Badge type="info" text="method" /> {#BuildCollar-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1156">source</a>

## `LoadBlocks` <Badge type="tip" text="dataclass" /> {#LoadBlocks}

```python
LoadBlocks(geo_paths: tuple[str, ...], adjacency: tuple[str, ...] = ())
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1221">source</a>

Blocks from `.geo` files, wired from adjacency sidecars; replaces the mesh.

### `LoadBlocks.apply` <Badge type="info" text="method" /> {#LoadBlocks-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1231">source</a>

## `Split` <Badge type="tip" text="dataclass" /> {#Split}

```python
Split(splits: tuple[tuple[str, str, tuple[int, ...]], ...])
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1246">source</a>

Split blocks by `(label, axis, at)` directives (`shore split`).

### `Split.apply` <Badge type="info" text="method" /> {#Split-apply}

```python
apply(project: Project, ctx: Context) -> Undo
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1254">source</a>

## `WriteGrd` <Badge type="tip" text="dataclass" /> {#WriteGrd}

```python
WriteGrd(path: str, geo_paths: tuple[str, ...])
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1268">source</a>

Pack `.geo` files into an Xall `.grd` (`shore grd`, output only).

Packs from the `.geo` files, not from memory, so the `.grd` carries
the `.geo` print precision exactly as the file pipeline does.

### `WriteGrd.apply` <Badge type="info" text="method" /> {#WriteGrd-apply}

```python
apply(project: Project, ctx: Context) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1282">source</a>

## `WriteCcPar` <Badge type="tip" text="dataclass" /> {#WriteCcPar}

```python
WriteCcPar(
    adjacency: str,
    path: str,
    grd_filename: str,
    output_basename: str = 'cc.',
    free_face_bc: int = 40,
    multigrid_levels: tuple[int, int] = (4, 1),
    boundary_layer_thickness: float = -0.1,
    numerical_beach: float = -1.0,
    meta: dict[str, dict[str, Any]] = field(default_factory=dict),
    boxes: tuple[dict[str, Any], ...] = (),
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1289">source</a>

Xall `cc.par` from an adjacency sidecar (`shore cc-par`, output only).

*meta* maps block labels to `BlockMeta` overrides (`level`,
`group`, `priority`, `comment`, `free_face_bc`); *boxes* is a
list of marker-box JSON entries.

### `WriteCcPar.apply` <Badge type="info" text="method" /> {#WriteCcPar-apply}

```python
apply(project: Project, ctx: Context) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1312">source</a>

## `WriteXall` <Badge type="tip" text="dataclass" /> {#WriteXall}

```python
WriteXall(
    directory: str,
    grd: str = 'mesh.grd',
    np: int = 1,
    components: dict[str, dict[str, Any]] = field(default_factory=dict),
    output_basename: str = 'cc.',
    multigrid_levels: tuple[int, int] = (4, 1),
    boundary_layer_thickness: float = -0.1,
    numerical_beach: float = -1.0,
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1344">source</a>

The Xall files of the whole run, every component (ADR 0016, 0017; output only).

For each component with a mesh, in creation order, into
`directory/<component>/`: the `.geo` files (one per block), the
adjacency sidecar, `<component>.grd` and `cc.par`, written by the
same commands as `write_geo` / `write_adjacency` / `write_grd` /
`write_cc_par`. Then, in *directory*, *grd* and `cc.par` merged from
them (`merge_grd`, `merge_cc_par`) and one `proc.input` balanced
over the merged `.grd` on *np* ranks.

The `cc.par` header options apply to every component (`merge_cc_par`
requires them to agree). *components* maps a component to its options:
`free_face_bc` / `meta` / `boxes` for its `cc.par` (as
`write_cc_par`), `body` for all its blocks in `proc.input`, and
`proc_meta` (block label -&gt; `proc` / `group` / `body`) per block.
All of it is checked before anything is written.

### `WriteXall.apply` <Badge type="info" text="method" /> {#WriteXall-apply}

```python
apply(project: Project, ctx: Context) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1379">source</a>

## `WriteProcInput` <Badge type="tip" text="dataclass" /> {#WriteProcInput}

```python
WriteProcInput(
    grd: str,
    path: str,
    np: int = 1,
    meta: dict[str, dict[str, int]] = field(default_factory=dict),
    header_comment: str = 'generated by SHORE',
)
```
Bases: `Command`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1470">source</a>

Xall `proc.input` from a `.grd` (`shore proc-input`, output only).

*meta* maps 1-based block indices (as strings) to `proc` / `group`
/ `body` overrides.

### `WriteProcInput.apply` <Badge type="info" text="method" /> {#WriteProcInput-apply}

```python
apply(project: Project, ctx: Context) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/commands.py#L1487">source</a>
