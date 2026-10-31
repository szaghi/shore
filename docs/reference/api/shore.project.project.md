---
title: shore.project.project
outline: [2, 3]
---

# `shore.project.project`

Project: the stateful document the CLI, the Python API and the GUI drive.

Every change goes through a `Command` (`shore.project.commands`); the
public methods below build one and execute it. The project keeps the applied
commands as its history, so it can undo and redo them and render them as a
*journal*: a plain Python script of the same method calls that
`shore run` replays (ADR 0004).

A project holds named geometries and named components (ADR 0016): every
mesh command acts on the current component (`use_component`; the default
is `"main"`), and `write_xall` writes the Xall files of all of them.

On disk a project is a directory (format 3; 1 and 2 still load):

```
project.toml                 metadata: format, shore version, components
                             and their block labels, input hashes
journal.py                   the state-changing history, reading inputs/
inputs/                      copies of every file the journal reads
adjacency.json               the current component's face graph
grids/<component>/<n>.npz    node arrays, checked on load
grids/<component>/adjacency.json  each component's face graph (ADR 0033)

```
`Project.load` checks the input copies against their hashes, replays
`journal.py` with the inputs resolved against the project directory
(wherever it now is), and requires every grid to be bit-identical to the
saved one, so a load is also a determinism check. Format-1 projects (one
mesh, inputs at their journaled paths) still load. `load(replay=False)`
(ADR 0033) skips the rebuild: the journal is recorded (geometries loaded,
components selected), and each component's blocks are read from its grids and
rewired from its face graph.

## `Project` <Badge type="tip" text="class" /> {#Project}

```python
Project(undo_root: str | Path | None = None)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L112">source</a>

Meshing project: geometries, components and the command history that made them.

### `Project.mesh` <Badge type="info" text="property" /> {#Project-mesh}

```python
mesh: Mesh | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L130">source</a>

The current component's mesh.

### `Project.build` <Badge type="info" text="property" /> {#Project-build}

```python
build: BuildState | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L139">source</a>

What later commands need from the command that made the current mesh.

### `Project.blocking` <Badge type="info" text="property" /> {#Project-blocking}

```python
blocking: Blocking | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L148">source</a>

The current component's blocking (ADR 0042).

### `Project.blockings` <Badge type="info" text="property" /> {#Project-blockings}

```python
blockings: dict[str, Blocking | None]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L157">source</a>

Each component's blocking (`None` without one), in creation order.

### `Project.require_blocking` <Badge type="info" text="method" /> {#Project-require_blocking}

```python
require_blocking(command: str) -> Blocking
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L161">source</a>

### `Project.body` <Badge type="info" text="property" /> {#Project-body}

```python
body: Body | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L167">source</a>

The geometry named `"body"` (the default of every command).

### `Project.geometries` <Badge type="info" text="property" /> {#Project-geometries}

```python
geometries: Mapping[str, Body]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L172">source</a>

Loaded geometries by name (read-only).

### `Project.components` <Badge type="info" text="property" /> {#Project-components}

```python
components: dict[str, Mesh | None]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L177">source</a>

Each component's mesh (`None` while empty), in creation order.

### `Project.current_component` <Badge type="info" text="property" /> {#Project-current_component}

```python
current_component: str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L182">source</a>

### `Project.require_geometry` <Badge type="info" text="method" /> {#Project-require_geometry}

```python
require_geometry(name: str, command: str) -> Body
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L185">source</a>

### `Project.geometry_digest` <Badge type="info" text="method" /> {#Project-geometry_digest}

```python
geometry_digest(name: str) -> str | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L194">source</a>

SHA-256 of the file geometry *name* was loaded from (None if not loaded).

### `Project.set_geometry` <Badge type="info" text="method" /> {#Project-set_geometry}

```python
set_geometry(name: str, body: Body | None, digest: str | None = None) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L198">source</a>

Register (or, with `None`, remove) geometry *name*; for commands.

### `Project.discard_built_on` <Badge type="info" text="method" /> {#Project-discard_built_on}

```python
discard_built_on(name: str) -> dict[str, tuple[Mesh | None, BuildState | None]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L208">source</a>

Empty the components whose mesh was built on geometry *name*; return
their previous state for `restore_components`. For commands.

### `Project.restore_components` <Badge type="info" text="method" /> {#Project-restore_components}

```python
restore_components(old: dict[str, tuple[Mesh | None, BuildState | None]]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L222">source</a>

### `Project.select_component` <Badge type="info" text="method" /> {#Project-select_component}

```python
select_component(name: str) -> bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L227">source</a>

Make *name* current, creating it if new; True if it was created. For commands.

### `Project.drop_component` <Badge type="info" text="method" /> {#Project-drop_component}

```python
drop_component(name: str) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L235">source</a>

Remove component *name* (not the current one). For commands.

### `Project.load_geometry` <Badge type="info" text="method" /> {#Project-load_geometry}

```python
load_geometry(
    path: str | Path,
    anchor_strategy: str = 'auto',
    skip_hygiene: bool = False,
    *,
    name: str = 'body',
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L242">source</a>

Load geometry *name* from STL or STEP/IGES/BREP (ADR 0016).

### `Project.use_component` <Badge type="info" text="method" /> {#Project-use_component}

```python
use_component(name: str) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L253">source</a>

Make component *name* current, creating it empty if new (ADR 0016).

### `Project.build_topology` <Badge type="info" text="method" /> {#Project-build_topology}

```python
build_topology(
    topology: str,
    *,
    nk: int,
    ds: float = 0.001,
    growth: float = 1.15,
    k_law: str = 'geometric',
    k_beta: float = 3.0,
    k_thickness: float | None = None,
    options: dict[str, Any] | None = None,
    geometry: str = 'body',
    **params: Any,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L257">source</a>

Build the blocks of *topology* on *geometry* (see
`shore.pipeline.build_topology`).

### `Project.march` <Badge type="info" text="method" /> {#Project-march}

```python
march(
    smooth_strength: float = 0.2,
    smooth_iters: int = 2,
    blend_normals_k: int | None = None,
    smooth_mode: str = 'limited',
    march_scheme: str = 'algebraic',
    *,
    progress: ProgressFn | None = None,
    cancel: CancelToken | None = None,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L288">source</a>

Extrude the blocks from the wall; *progress* / *cancel* are hooks.

`blend_normals_k` is removed (ADR 0034): 0 warns, anything else is
refused when the command runs (a journal loaded without replay keeps it).

### `Project.snap` <Badge type="info" text="method" /> {#Project-snap}

```python
snap(exact: bool = False) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L312">source</a>

Project associated nodes onto their geometry (any project geometry, ADR 0039).

### `Project.redistribute` <Badge type="info" text="method" /> {#Project-redistribute}

```python
redistribute(
    block: str,
    edge: str,
    law: str = 'uniform',
    *,
    beta: float = 3.0,
    cluster_at: float | None = None,
    first_step: float | None = None,
    last_step: float | None = None,
    interior: str = 'uniform',
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L316">source</a>

Apply a spacing law along *edge* of *block* to the lines of its axis,
through the seams (`Redistribute`, ADR 0041); then `snap` if needed.

### `Project.associate` <Badge type="info" text="method" /> {#Project-associate}

```python
associate(
    block: str,
    entity: str,
    geometry: str = 'body',
    kind: str = 'body',
    indices: Sequence[int] = (),
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L334">source</a>

Associate face, edge or vertex *entity* of *block* (current component)
with entities *indices* of kind *kind* of geometry *geometry*, or clear
it (`kind="none"`); see `Associate` (ADR 0003, 0039).

### `Project.create_blocking` <Badge type="info" text="method" /> {#Project-create_blocking}

```python
create_blocking(
    ni: int,
    nj: int,
    nk: int,
    *,
    geometry: str = 'body',
    margin: float = 0.0,
    vmin: Sequence[float] | None = None,
    vmax: Sequence[float] | None = None,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L347">source</a>

Start the blocking with one block of *ni* x *nj* x *nk* nodes around
*geometry* (bounding box grown by *margin* x its extent per side) or
between *vmin* and *vmax* (`CreateBlocking`, ADR 0042).

### `Project.split_blocking` <Badge type="info" text="method" /> {#Project-split_blocking}

```python
split_blocking(block: str, axis: str, fraction: float) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L364">source</a>

Cut *block* across *axis* at *fraction*, through its sheet (`SplitBlocking`).

### `Project.delete_blocks` <Badge type="info" text="method" /> {#Project-delete_blocks}

```python
delete_blocks(blocks: Sequence[str]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L368">source</a>

Remove blocking blocks *blocks* (`DeleteBlocks`).

### `Project.ogrid` <Badge type="info" text="method" /> {#Project-ogrid}

```python
ogrid(
    blocks: Sequence[str] = (),
    *,
    faces: Sequence[Sequence[str]] = (),
    exclude: Sequence[Sequence[str]] = (),
    factor: float = 0.5,
    nodes: int = 5,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L372">source</a>

Insert an O-grid on the faces of blocking blocks *blocks* (all when
empty): *faces*, or the set's boundary less *exclude* (`(block, face)`
pairs), sized by *factor*, *nodes* across (`Ogrid`, ADR 0043).

### `Project.move_vertex` <Badge type="info" text="method" /> {#Project-move_vertex}

```python
move_vertex(block: str, vertex: str, xyz: Sequence[float]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L390">source</a>

Move blocking vertex *vertex* (`v0`..`v7`) of *block* to *xyz*,
onto its association if it has one (`MoveVertex`, ADR 0045).

### `Project.associate_blocking` <Badge type="info" text="method" /> {#Project-associate_blocking}

```python
associate_blocking(
    block: str,
    entity: str,
    geometry: str = 'body',
    kind: str = 'body',
    indices: Sequence[int] = (),
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L399">source</a>

Associate face, edge or vertex *entity* of blocking block *block* with
geometry, as `associate` does for mesh blocks (`AssociateBlocking`).

### `Project.symmetry_plane` <Badge type="info" text="method" /> {#Project-symmetry_plane}

```python
symmetry_plane(point: Sequence[float], normal: Sequence[float]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L413">source</a>

Declare a symmetry plane of the blocking through *point*, normal to
*normal* (`SymmetryPlane`, ADR 0053).

### `Project.premesh` <Badge type="info" text="method" /> {#Project-premesh}

```python
premesh(exact: bool = False) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L426">source</a>

Make the component's mesh from its blocking (`Premesh`, ADR 0042).

### `Project.orthogonalize` <Badge type="info" text="method" /> {#Project-orthogonalize}

```python
orthogonalize(strength: float = 1.0) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L430">source</a>

Make the grid lines leave the walls along their normal
(`Orthogonalize`, ADR 0048).

### `Project.wall_layer` <Badge type="info" text="method" /> {#Project-wall_layer}

```python
wall_layer(first_step: float) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L435">source</a>

Cluster the nodes towards the walls, the first cell *first_step* long
(`WallLayer`, ADR 0050).

### `Project.surface_grid` <Badge type="info" text="method" /> {#Project-surface_grid}

```python
surface_grid(
    block: str,
    face: str,
    *,
    surface: int,
    smoothing: str = 'elliptic',
    geometry: str = 'body',
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L440">source</a>

Rebuild *face* of *block* from its boundary on SURFACE entity *surface*
of *geometry* (TFI + elliptic smoothing in (u, v), ADR 0007).

### `Project.build_primitive` <Badge type="info" text="method" /> {#Project-build_primitive}

```python
build_primitive(kind: str, **params: Any) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L453">source</a>

Build a parametric primitive (`box`, `annulus`, `flat_caps`,
`h_grid_rect`, `h_grid_tube`); see `shore.pipeline.build_primitive`.

### `Project.build_collar` <Badge type="info" text="method" /> {#Project-build_collar}

```python
build_collar(
    other: str | Path,
    *,
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
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L458">source</a>

Collar block around the junction of *geometry* and *other*, a
geometry name or a file path (ADR 0010, 0016); with *plane* the
half-collar of a half model (ADR 0013); `loop_*` distribute the
nodes along the junction curve (ADR 0020); *split* cuts the collar
at sharp corners, with wakes of *wake_length* (ADR 0024); *volume*
the outer-boundary construction (ADR 0025). Replaces the current
component's mesh. See `BuildCollar`.

### `Project.build_tip_cap` <Badge type="info" text="method" /> {#Project-build_tip_cap}

```python
build_tip_cap(body_geo: str | Path, **options: Any) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L514">source</a>

Wing-tip cap from a C-grid body `.geo`; options as `BuildTipCap`.

### `Project.load_blocks` <Badge type="info" text="method" /> {#Project-load_blocks}

```python
load_blocks(
    geo_paths: list[str | Path],
    adjacency: list[str | Path] | tuple = (),
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L518">source</a>

Replace the mesh with blocks read from `.geo` files (+ sidecars).

### `Project.split` <Badge type="info" text="method" /> {#Project-split}

```python
split(splits: list[tuple[str, str, list[int] | tuple[int, ...]]]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L524">source</a>

Split blocks by `(label, axis, at)` directives.

### `Project.write_geo` <Badge type="info" text="method" /> {#Project-write_geo}

```python
write_geo(stem: str | Path, per_block: bool = False) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L530">source</a>

Write the `.geo` files (see `WriteGeo` for the naming rules).

### `Project.write_grd` <Badge type="info" text="method" /> {#Project-write_grd}

```python
write_grd(path: str | Path, geo_paths: list[str | Path]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L534">source</a>

Pack `.geo` files into an Xall `.grd`.

### `Project.write_cc_par` <Badge type="info" text="method" /> {#Project-write_cc_par}

```python
write_cc_par(adjacency: str | Path, path: str | Path, **options: Any) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L538">source</a>

Write an Xall `cc.par`; options as `WriteCcPar`.

### `Project.write_proc_input` <Badge type="info" text="method" /> {#Project-write_proc_input}

```python
write_proc_input(grd: str | Path, path: str | Path, **options: Any) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L546">source</a>

Write an Xall `proc.input`; options as `WriteProcInput`.

### `Project.write_adjacency` <Badge type="info" text="method" /> {#Project-write_adjacency}

```python
write_adjacency(path: str | Path) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L550">source</a>

Write the `.adjacency.json` sidecar.

### `Project.write_cgns` <Badge type="info" text="method" /> {#Project-write_cgns}

```python
write_cgns(path: str | Path, components: Sequence[str] = ()) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L554">source</a>

Write a CGNS file of *components* (all with a mesh when empty): a base
per component, a zone per block, 1-to-1 connections and BCs from the
face graph (ADR 0035). Needs `h5py` (`shore-mesh[cgns]`).

### `Project.write_plot3d` <Badge type="info" text="method" /> {#Project-write_plot3d}

```python
write_plot3d(path: str | Path, components: Sequence[str] = ()) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L560">source</a>

Write the blocks of *components* (all with a mesh when empty) as one
multi-block Plot3D grid file (ADR 0035).

### `Project.write_xall` <Badge type="info" text="method" /> {#Project-write_xall}

```python
write_xall(
    directory: str | Path,
    *,
    grd: str = 'mesh.grd',
    np: int = 1,
    components: dict[str, dict[str, Any]] | None = None,
    output_basename: str = 'cc.',
    multigrid_levels: tuple[int, int] = (4, 1),
    boundary_layer_thickness: float = -0.1,
    numerical_beach: float = -1.0,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L565">source</a>

Write the Xall files of every component, merged; per-component and
run-wide options as `WriteXall` (ADR 0016, 0017).

### `Project.execute` <Badge type="info" text="method" /> {#Project-execute}

```python
execute(
    command: Command,
    *,
    progress: ProgressFn | None = None,
    cancel: CancelToken | None = None,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L594">source</a>

Apply *command*, record it, and clear the redo stack.

While a project loads without replay (ADR 0033) only geometries and
component switches are applied; the other commands are recorded.

### `Project.undo` <Badge type="info" text="method" /> {#Project-undo}

```python
undo() -> Command
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L641">source</a>

Revert the last state-changing command (and any outputs after it).

### `Project.redo` <Badge type="info" text="method" /> {#Project-redo}

```python
redo() -> Command
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L657">source</a>

Re-apply the last undone command.

### `Project.history` <Badge type="info" text="property" /> {#Project-history}

```python
history: list[Command]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L667">source</a>

Applied commands, oldest first (undone ones excluded).

### `Project.require_mesh` <Badge type="info" text="method" /> {#Project-require_mesh}

```python
require_mesh() -> tuple[Mesh, BuildState]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L671">source</a>

### `Project.journal` <Badge type="info" text="method" /> {#Project-journal}

```python
journal() -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L678">source</a>

The history as a replayable Python script.

### `Project.save_journal` <Badge type="info" text="method" /> {#Project-save_journal}

```python
save_journal(path: str | Path) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L682">source</a>

### `Project.save` <Badge type="info" text="method" /> {#Project-save}

```python
save(directory: str | Path) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L688">source</a>

Write the project directory, inputs included (see module docstring).

### `Project.load` <Badge type="info" text="classmethod" /> {#Project-load}

```python
load(directory: str | Path, *, replay: bool = True) -> Project
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L773">source</a>

Replay a project directory's journal and verify its saved grids.

With `replay=False` (ADR 0033) the journal is only recorded:
geometries are loaded and components selected, but no mesh is rebuilt;
each component's blocks come from its saved grids and face graph. The
loaded history cannot be undone, and a loaded mesh cannot be marched.

### `Project.close` <Badge type="info" text="method" /> {#Project-close}

```python
close() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L864">source</a>

Delete the on-disk undo snapshots.

## `run_journal` <Badge type="info" text="function" /> {#run_journal}

```python
run_journal(path: str | Path, input_dir: str | Path | None = None) -> Project
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/project/project.py#L872">source</a>

Execute a journal script and return the `Project` it built (`p`).

With *input_dir*, relative input paths in the journal resolve against it
(a saved project's directory). A journal is Python: run only journals you
trust.
