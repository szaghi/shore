---
title: shore.gui.viewport
outline: [2, 3]
---

# `shore.gui.viewport`

The GUI's renderer interface (ADR 0006, 0037).

Only :mod:`shore.gui.viewport.pygfx_backend` imports pygfx, wgpu or
rendercanvas. This module holds the interface and the structured-grid
arithmetic the backends share: a block is drawn as its six boundary faces
(one quad per boundary cell face, in `FACES` order) and its twelve edges,
and a picked quad maps back to its (block, face, i, j, k) in O(1).

**Modules**

- [`shore.gui.viewport.pygfx_backend`](./shore.gui.viewport.pygfx_backend) — pygfx/wgpu implementation of the GUI viewport

## `SliceKey` <Badge type="tip" text="dataclass" /> {#SliceKey}

```python
SliceKey(block: Any, axis: str, index: int)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L22">source</a>

A constant-index plane of a block: `axis` in `"ijk"`, node `index`.

## `GeometryPick` <Badge type="tip" text="dataclass" /> {#GeometryPick}

```python
GeometryPick(geometry: str, entity: Any, xyz: np.ndarray)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L31">source</a>

A picked geometry entity and the float64 point on it.

## `SkeletonPick` <Badge type="tip" text="dataclass" /> {#SkeletonPick}

```python
SkeletonPick(key: Any, node: int, nodes: int, xyz: np.ndarray)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L40">source</a>

A picked blocking edge (ADR 0042): its key, the node nearest the pick
(of `nodes`, counted along the drawn polyline) and that node's point.

## `MarkerPick` <Badge type="tip" text="dataclass" /> {#MarkerPick}

```python
MarkerPick(key: Any, xyz: np.ndarray)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L51">source</a>

A picked blocking face or block marker (ADR 0044): its key and centre.

## `PickResult` <Badge type="tip" text="dataclass" /> {#PickResult}

```python
PickResult(key: Any, face: str, cell: tuple[int, int, int], xyz: np.ndarray)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L59">source</a>

A picked boundary cell face: its block, face and cell, and its centre.

## `Viewport` <Badge type="tip" text="class" /> {#Viewport}

```python
class Viewport(Protocol)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L68">source</a>

What the GUI needs from a renderer.

### `Viewport.set_blocks` <Badge type="info" text="method" /> {#Viewport-set_blocks}

```python
set_blocks(blocks: Sequence[tuple[Any, np.ndarray]]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L71">source</a>

Replace the scene with these `(key, nodes (ni, nj, nk, 3))` blocks.

### `Viewport.set_face_scalars` <Badge type="info" text="method" /> {#Viewport-set_face_scalars}

```python
set_face_scalars(
    values: dict[Any, np.ndarray] | None,
    clim: tuple[float, float],
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L74">source</a>

Colour each boundary quad by a value (`boundary_values` order), or clear.

### `Viewport.show_all` <Badge type="info" text="method" /> {#Viewport-show_all}

```python
show_all() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L79">source</a>

### `Viewport.look_at_point` <Badge type="info" text="method" /> {#Viewport-look_at_point}

```python
look_at_point(xyz: np.ndarray, distance: float) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L81">source</a>

### `Viewport.pick` <Badge type="info" text="method" /> {#Viewport-pick}

```python
pick(x: float, y: float) -> PickResult | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L83">source</a>

### `Viewport.project` <Badge type="info" text="method" /> {#Viewport-project}

```python
project(xyz: np.ndarray) -> tuple[float, float]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L85">source</a>

A world point (float64) to logical pixel coordinates.

### `Viewport.request_draw` <Badge type="info" text="method" /> {#Viewport-request_draw}

```python
request_draw() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L88">source</a>

## `face_grids` <Badge type="info" text="function" /> {#face_grids}

```python
face_grids(nodes: np.ndarray) -> list[np.ndarray]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L91">source</a>

The six boundary node grids `(na, nb, 3)` of a block, in `FACES` order.

## `quad_indices` <Badge type="info" text="function" /> {#quad_indices}

```python
quad_indices(na: int, nb: int, base: int = 0) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L103">source</a>

Quads `((na-1)*(nb-1), 4)` of an `na x nb` node grid, a-major.

## `quad_offsets` <Badge type="info" text="function" /> {#quad_offsets}

```python
quad_offsets(shape: tuple[int, int, int]) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L110">source</a>

Index of each face's first quad in a block's quad list, plus the total.

## `quad_cell` <Badge type="info" text="function" /> {#quad_cell}

```python
quad_cell(shape: tuple[int, int, int], quad: int) -> tuple[str, tuple[int, int, int]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L117">source</a>

The face and the (i, j, k) cell behind quad number *quad* of a block.

## `quad_centre` <Badge type="info" text="function" /> {#quad_centre}

```python
quad_centre(nodes: np.ndarray, face: str, cell: tuple[int, int, int]) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L136">source</a>

float64 centre of the boundary quad of *cell* on *face*.

## `boundary_values` <Badge type="info" text="function" /> {#boundary_values}

```python
boundary_values(cell_values: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L146">source</a>

Per boundary quad, the value of the cell behind it (`(ni-1, nj-1, nk-1)` input).

Quads follow `FACES` order and, within a face, `quad_indices` order.

## `edge_segments` <Badge type="info" text="function" /> {#edge_segments}

```python
edge_segments(nodes: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L156">source</a>

The 12 block edges as line-segment vertex pairs `(2m, 3)`.

## `colormap` <Badge type="info" text="function" /> {#colormap}

```python
colormap(values: np.ndarray, clim: tuple[float, float]) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L164">source</a>

RGBA float32 of *values* on a blue-white-red ramp over *clim* (NaN: grey).

## `slice_grid` <Badge type="info" text="function" /> {#slice_grid}

```python
slice_grid(nodes: np.ndarray, axis: str, index: int) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L176">source</a>

The node grid `(na, nb, 3)` of the constant-*axis* plane *index*.

## `slice_cell` <Badge type="info" text="function" /> {#slice_cell}

```python
slice_cell(
    shape: tuple[int, int, int],
    axis: str,
    index: int,
    quad: int,
) -> tuple[int, int, int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L182">source</a>

The cell of a slice quad: on the plane's + side, or its - side on the last plane.

## `slice_values` <Badge type="info" text="function" /> {#slice_values}

```python
slice_values(cell_values: np.ndarray, axis: str, index: int) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L195">source</a>

Per slice quad, the value of its cell (`slice_cell`), in quad order.

## `edge_nodes` <Badge type="info" text="function" /> {#edge_nodes}

```python
edge_nodes(nodes: np.ndarray, name: str) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L212">source</a>

The node polyline of block edge *name* (`HexBlock` edge naming).

## `face_edges` <Badge type="info" text="function" /> {#face_edges}

```python
face_edges(face: str) -> tuple[str, ...]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L220">source</a>

The four edges bounding *face*.

## `face_vertices` <Badge type="info" text="function" /> {#face_vertices}

```python
face_vertices(face: str) -> tuple[int, ...]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L225">source</a>

The four corners (vertex indices) of *face*.

## `nearest_edge` <Badge type="info" text="function" /> {#nearest_edge}

```python
nearest_edge(nodes: np.ndarray, face: str, xyz: np.ndarray) -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L241">source</a>

The edge of *face* closest to the point *xyz* (an edge pick).

## `nearest_vertex` <Badge type="info" text="function" /> {#nearest_vertex}

```python
nearest_vertex(nodes: np.ndarray, face: str, xyz: np.ndarray) -> int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/__init__.py#L246">source</a>

The corner of *face* closest to the point *xyz* (a vertex pick).
