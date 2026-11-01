---
title: shore.gui.window
outline: [2, 3]
---

# `shore.gui.window`

SHORE's main window (ADR 0006, 0037): model tree, viewport, command forms, console.

Every action goes through :class:`shore.gui.session.Session`, so it is a
Project command and a journal line (ADR 0004). Long commands run on a worker
thread with the ADR 0004 progress and cancel hooks; the window refreshes from
the session when they finish.

Layout (ADR 0006): model tree on the left, the viewport in the centre,
command forms on the right, the journal and the log at the bottom.

## `MainWindow` <Badge type="tip" text="class" /> {#MainWindow}

```python
MainWindow(
    session: Session | None = None,
    *,
    ask_path: AskPath = _ask_path,
    offscreen_viewport: bool = False,
)
```
Bases: `QMainWindow`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L92">source</a>

The SHORE GUI over one :class:`Session`.

### `MainWindow.select_skeleton` <Badge type="info" text="method" /> {#MainWindow-select_skeleton}

```python
select_skeleton(key: SkeletonKey, *, add: bool = False) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L266">source</a>

Select a blocking edge, face or block (`key.edge`: an edge or face
name, `""` for a block); *add* toggles it in the selection.

### `MainWindow.clear_blocking_selection` <Badge type="info" text="method" /> {#MainWindow-clear_blocking_selection}

```python
clear_blocking_selection() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L277">source</a>

### `MainWindow.ogrid_selection` <Badge type="info" text="method" /> {#MainWindow-ogrid_selection}

```python
ogrid_selection(factor: float = 0.5, nodes: int = 5) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L302">source</a>

Run `Project.ogrid` on the selection: the selected blocks are the
set (every block when none is selected), the selected faces the layered
faces (else the set's whole boundary), each named by its block in the
set.

### `MainWindow.delete_selection` <Badge type="info" text="method" /> {#MainWindow-delete_selection}

```python
delete_selection() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L341">source</a>

Run `Project.delete_blocks` on the selected blocks.

### `MainWindow.move_selection` <Badge type="info" text="method" /> {#MainWindow-move_selection}

```python
move_selection(xyz: Any) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L355">source</a>

Run `Project.move_vertex` for the one selected blocking vertex
(ADR 0045): to *xyz*, onto its association if it has one.

### `MainWindow.move_selection_to_pick` <Badge type="info" text="method" /> {#MainWindow-move_selection_to_pick}

```python
move_selection_to_pick() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L379">source</a>

Move the selected blocking vertex to the last picked geometry point.

### `MainWindow.associate_blocking_selection` <Badge type="info" text="method" /> {#MainWindow-associate_blocking_selection}

```python
associate_blocking_selection() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L386">source</a>

Associate each selected blocking face or edge with the geometry
selected in the Association dock (`Project.associate_blocking`).

### `MainWindow.split_selection` <Badge type="info" text="method" /> {#MainWindow-split_selection}

```python
split_selection() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L411">source</a>

Run `Project.split_blocking` through the selected skeleton node:
the picked edge's block, its axis and the node's fraction along it
(making the edge's component current first if needed).

### `MainWindow.edge_info` <Badge type="info" text="method" /> {#MainWindow-edge_info}

```python
edge_info() -> dict[str, Any] | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L435">source</a>

The selected edge: block, edge, axis, nodes, length, first and last cell.

### `MainWindow.redistribute_selection` <Badge type="info" text="method" /> {#MainWindow-redistribute_selection}

```python
redistribute_selection(law: str = 'uniform', **params: Any) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L469">source</a>

Run `Project.redistribute` on the selected edge (making its
component current first if needed).

### `MainWindow.select_form` <Badge type="info" text="method" /> {#MainWindow-select_form}

```python
select_form(method: str) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L609">source</a>

Show the form of *method* in the Command dock.

### `MainWindow.refresh` <Badge type="info" text="method" /> {#MainWindow-refresh}

```python
refresh() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L714">source</a>

Bring the tree, the journal and the viewport in line with the session.

### `MainWindow.set_visible` <Badge type="info" text="method" /> {#MainWindow-set_visible}

```python
set_visible(kind: str, visible: bool) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L781">source</a>

Show or hide all blocks, geometries, slices or the blocking (its
skeleton and markers).

### `MainWindow.add_slice` <Badge type="info" text="method" /> {#MainWindow-add_slice}

```python
add_slice(block: BlockKey, axis: str, index: int) -> SliceKey
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L791">source</a>

Show the constant-*axis* plane *index* of *block*.

### `MainWindow.remove_slice` <Badge type="info" text="method" /> {#MainWindow-remove_slice}

```python
remove_slice(key: SliceKey) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L807">source</a>

### `MainWindow.refresh_quality` <Badge type="info" text="method" /> {#MainWindow-refresh_quality}

```python
refresh_quality() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L859">source</a>

Colour boundary quads by the scaled Jacobian of their cells, or clear.

### `MainWindow.run_form` <Badge type="info" text="method" /> {#MainWindow-run_form}

```python
run_form(method: str, values: dict[str, Any] | None = None) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L906">source</a>

Run the form of *method* with *values* (field name -&gt; text or value).

### `MainWindow.cancel` <Badge type="info" text="method" /> {#MainWindow-cancel}

```python
cancel() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L982">source</a>

### `MainWindow.wait_idle` <Badge type="info" text="method" /> {#MainWindow-wait_idle}

```python
wait_idle(timeout: float = 600.0) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L986">source</a>

Process events until no command runs (tests and scripted use).

### `MainWindow.undo` <Badge type="info" text="method" /> {#MainWindow-undo}

```python
undo() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L997">source</a>

### `MainWindow.redo` <Badge type="info" text="method" /> {#MainWindow-redo}

```python
redo() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L1000">source</a>

### `MainWindow.new_project` <Badge type="info" text="method" /> {#MainWindow-new_project}

```python
new_project() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L1003">source</a>

### `MainWindow.open_project` <Badge type="info" text="method" /> {#MainWindow-open_project}

```python
open_project(directory: str | None = None, *, replay: bool = True) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L1006">source</a>

### `MainWindow.save_project` <Badge type="info" text="method" /> {#MainWindow-save_project}

```python
save_project(directory: str | None = None, *, ask: bool = False) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L1011">source</a>

### `MainWindow.pick_at` <Badge type="info" text="method" /> {#MainWindow-pick_at}

```python
pick_at(x: float, y: float, *, add: bool = False) -> Any
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L1021">source</a>

Pick at pixel (x, y) as a double click does (shift: *add*).

### `MainWindow.associate_selection` <Badge type="info" text="method" /> {#MainWindow-associate_selection}

```python
associate_selection() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L1108">source</a>

Associate the selected block entity with the selected geometry entities.

### `MainWindow.clear_association` <Badge type="info" text="method" /> {#MainWindow-clear_association}

```python
clear_association() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L1117">source</a>

## `main` <Badge type="info" text="function" /> {#main}

```python
main(project: str | None = None, *, replay: bool = True, smoke: bool = False) -> int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/window.py#L1133">source</a>

Run the GUI; *smoke* builds the window, shows it briefly and exits.
