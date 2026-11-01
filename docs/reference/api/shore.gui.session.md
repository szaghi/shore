---
title: shore.gui.session
outline: [2, 3]
---

# `shore.gui.session`

The GUI's owner of a project (ADR 0037): every action is a Project command.

`Session` holds no Qt object, so the GUI logic is tested headless. Widgets
call its methods only; each state change goes through `Project` (and so
into its journal, ADR 0004), and listeners are told what changed.

Long commands (march, builds) are run by the window on a worker thread, with
the ADR 0004 progress and cancel hooks; `Session.run` itself is
synchronous.

## `BlockKey` <Badge type="tip" text="dataclass" /> {#BlockKey}

```python
BlockKey(component: str, label: str)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L28">source</a>

A block of a component, as the viewport and the model tree name it.

## `SkeletonKey` <Badge type="tip" text="dataclass" /> {#SkeletonKey}

```python
SkeletonKey(component: str, block: str, edge: str)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L36">source</a>

A blocking entity of a component (ADR 0042, 0044): `edge` is an edge
or face name of blocking block `block`, or `""` for the block itself.

## `Session` <Badge type="tip" text="class" /> {#Session}

```python
Session(project: Project | None = None, directory: Path | None = None)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L45">source</a>

One project, its saved location and the GUI's change listeners.

### `Session.subscribe` <Badge type="info" text="method" /> {#Session-subscribe}

```python
subscribe(listener: Listener) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L57">source</a>

Call *listener* after every change of the project.

### `Session.run` <Badge type="info" text="method" /> {#Session-run}

```python
run(
    method: str,
    *args: Any,
    progress: ProgressFn | None = None,
    cancel: CancelToken | None = None,
    **kwargs: Any,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L67">source</a>

Call `Project.<method>`: one command, one journal line.

The project is unchanged if the command raises (the command layer
restores a failed or cancelled march, and a command that fails while
building leaves no state behind).

### `Session.run_cli` <Badge type="info" text="method" /> {#Session-run_cli}

```python
run_cli(command: Sequence[str]) -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L91">source</a>

Run a CLI command (`["mesh", "body.stl", "--nk", "6"]`) on this
project (ADR 0062) and return its output.

Its Project commands land in this project's history, each one
undoable. A command that fails (exit code not 0) is rolled back:
the project is as it was, and the error carries the command's output.

### `Session.run_journal` <Badge type="info" text="method" /> {#Session-run_journal}

```python
run_journal(path: str | Path) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L117">source</a>

Replace the project with the one a journal builds (`shore run`).

### `Session.view_file` <Badge type="info" text="method" /> {#Session-view_file}

```python
view_file(path: str | Path) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L126">source</a>

Show a file as `shore view` does: a `.geo` file loads as the
current component's blocks (`load_blocks`, one journal line); a
`.npz` surface (key `"k0"`) is shown only, outside the project.

### `Session.undo` <Badge type="info" text="method" /> {#Session-undo}

```python
undo() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L143">source</a>

### `Session.redo` <Badge type="info" text="method" /> {#Session-redo}

```python
redo() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L149">source</a>

### `Session.open` <Badge type="info" text="method" /> {#Session-open}

```python
open(directory: str | Path, *, replay: bool = True) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L155">source</a>

Replace the project with a saved one (ADR 0016; ADR 0033 without replay).

### `Session.new` <Badge type="info" text="method" /> {#Session-new}

```python
new() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L163">source</a>

### `Session.save` <Badge type="info" text="method" /> {#Session-save}

```python
save(directory: str | Path | None = None) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L169">source</a>

### `Session.journal_lines` <Badge type="info" text="method" /> {#Session-journal_lines}

```python
journal_lines() -> list[str]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L180">source</a>

The journal's command lines (without the header).

### `Session.blocks` <Badge type="info" text="method" /> {#Session-blocks}

```python
blocks() -> list[tuple[BlockKey, HexBlock]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L184">source</a>

Every block of every component with a mesh, components in project order.

### `Session.skeleton` <Badge type="info" text="method" /> {#Session-skeleton}

```python
skeleton() -> list[tuple[SkeletonKey, Any]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L193">source</a>

Every blocking edge of every component, with its nodes on the
straight line between its skeleton vertices (`Blocking.polylines`).

### `Session.skeleton_markers` <Badge type="info" text="method" /> {#Session-skeleton_markers}

```python
skeleton_markers(kind: str) -> list[tuple[SkeletonKey, Any]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L205">source</a>

A marker per blocking face (`kind="face"`: ``SkeletonKey(component,
block, face)`, the face named by its first block) or block (`"block"``:
`SkeletonKey(component, block, "")`), at its skeleton centre.

### `Session.tree` <Badge type="info" text="method" /> {#Session-tree}

```python
tree() -> list[tuple[str, bool, list[tuple[str, tuple[int, int, int], dict[str, str]]]]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L218">source</a>

`(component, is current, [(block label, shape, {face: BC name})])`.

### `Session.associations` <Badge type="info" text="method" /> {#Session-associations}

```python
associations(key: BlockKey) -> dict[str, str]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L232">source</a>

`{entity: label}` for the associated faces, edges and vertices
(`"v0"`..) of a block; labels like `"body"` or `"part: surface 2, 5"`.

### `Session.association_label` <Badge type="info" text="method" /> {#Session-association_label}

```python
association_label(association: Any) -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L253">source</a>

`"geometry"` for a body, else `"geometry: kind i, j"` (indices
re-found from the association's fingerprints in the current geometry).

### `Session.block` <Badge type="info" text="method" /> {#Session-block}

```python
block(key: BlockKey) -> HexBlock | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L270">source</a>

## `surface_display` <Badge type="info" text="function" /> {#surface_display}

```python
surface_display(nodes: Any) -> Any
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/session.py#L277">source</a>

A structured surface `(ni, nj, 3)` as a `geometry.base.Display`:
two triangles per quad, one body entity.
