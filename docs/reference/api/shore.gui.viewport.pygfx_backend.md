---
title: shore.gui.viewport.pygfx_backend
outline: [2, 3]
---

# `shore.gui.viewport.pygfx_backend`

pygfx/wgpu implementation of the GUI viewport (ADR 0006, 0037, 0038).

The only module of SHORE that imports pygfx, wgpu or rendercanvas.

Float precision: positions are kept in float64 on the CPU and uploaded as
`float32(p - origin)`; the origin follows the orbit target and moves
(re-uploading everything, "rebase") when the target drifts more than
`rebase_factor` view widths from it, or on `look_at_point`. Coordinates
reported by `pick` come from the float64 arrays (ADR 0006, spike item 3).

Drawn objects: per block its six boundary faces (quads) and twelve edges;
per slice plane its quads; per geometry its triangles, curves and points.

## `create_canvas` <Badge type="info" text="function" /> {#create_canvas}

```python
create_canvas(
    parent: Any = None,
    *,
    offscreen: bool = False,
    size: tuple[int, int] = (800, 600),
) -> Any
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L42">source</a>

A rendercanvas canvas: a Qt widget, or an offscreen canvas for tests.

## `PygfxViewport` <Badge type="tip" text="class" /> {#PygfxViewport}

```python
PygfxViewport(canvas: Any, *, rebase_factor: float = 100.0)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L98">source</a>

pygfx/wgpu implementation of `Viewport`.

### `PygfxViewport.set_blocks` <Badge type="info" text="method" /> {#PygfxViewport-set_blocks}

```python
set_blocks(blocks: Sequence[tuple[Any, np.ndarray]]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L147">source</a>

### `PygfxViewport.set_slices` <Badge type="info" text="method" /> {#PygfxViewport-set_slices}

```python
set_slices(slices: Sequence[tuple[SliceKey, np.ndarray]]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L161">source</a>

Replace the slice planes: `(SliceKey, the block's nodes)`.

### `PygfxViewport.set_geometries` <Badge type="info" text="method" /> {#PygfxViewport-set_geometries}

```python
set_geometries(geometries: Sequence[tuple[str, Any]]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L170">source</a>

Replace the geometries: `(name, shore.geometry.base.Display)`.

### `PygfxViewport.set_skeleton` <Badge type="info" text="method" /> {#PygfxViewport-set_skeleton}

```python
set_skeleton(edges: Sequence[tuple[Any, np.ndarray]]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L220">source</a>

Replace the blocking skeleton (ADR 0042): `(key, edge nodes (n, 3))`,
drawn over the blocks so it stays pickable.

### `PygfxViewport.set_markers` <Badge type="info" text="method" /> {#PygfxViewport-set_markers}

```python
set_markers(markers: Sequence[tuple[Any, np.ndarray]]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L242">source</a>

Replace the blocking face or block markers (ADR 0044): `(key, xyz)`,
points drawn over everything so each stays pickable.

### `PygfxViewport.marker_keys` <Badge type="info" text="property" /> {#PygfxViewport-marker_keys}

```python
marker_keys: list[Any]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L258">source</a>

### `PygfxViewport.skeleton_keys` <Badge type="info" text="property" /> {#PygfxViewport-skeleton_keys}

```python
skeleton_keys: list[Any]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L263">source</a>

### `PygfxViewport.show` <Badge type="info" text="method" /> {#PygfxViewport-show}

```python
show(kind: str, visible: bool) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L267">source</a>

Show or hide every `"block"`, `"slice"`, `"geometry"` or `"skeleton"` object.

### `PygfxViewport.visible` <Badge type="info" text="method" /> {#PygfxViewport-visible}

```python
visible(kind: str) -> bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L272">source</a>

### `PygfxViewport.set_face_scalars` <Badge type="info" text="method" /> {#PygfxViewport-set_face_scalars}

```python
set_face_scalars(
    values: dict[Any, np.ndarray] | None,
    clim: tuple[float, float] = (0.0, 1.0),
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L275">source</a>

Colour block and slice quads by value (keys: block keys, `SliceKey`), or clear.

### `PygfxViewport.face_colours` <Badge type="info" text="method" /> {#PygfxViewport-face_colours}

```python
face_colours(key: Any) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L291">source</a>

The uploaded per-quad RGBA of a block or slice (for tests).

### `PygfxViewport.keys` <Badge type="info" text="property" /> {#PygfxViewport-keys}

```python
keys: list[Any]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L297">source</a>

### `PygfxViewport.slice_keys` <Badge type="info" text="property" /> {#PygfxViewport-slice_keys}

```python
slice_keys: list[SliceKey]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L301">source</a>

### `PygfxViewport.geometry_names` <Badge type="info" text="property" /> {#PygfxViewport-geometry_names}

```python
geometry_names: list[str]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L305">source</a>

### `PygfxViewport.uploaded_vertices` <Badge type="info" text="method" /> {#PygfxViewport-uploaded_vertices}

```python
uploaded_vertices(key: Any) -> int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L308">source</a>

### `PygfxViewport.uploaded_positions` <Badge type="info" text="method" /> {#PygfxViewport-uploaded_positions}

```python
uploaded_positions(key: Any) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L312">source</a>

float64 world positions of what was uploaded for *key*'s first object.

### `PygfxViewport.view_width` <Badge type="info" text="method" /> {#PygfxViewport-view_width}

```python
view_width() -> float
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L333">source</a>

### `PygfxViewport.maybe_rebase` <Badge type="info" text="method" /> {#PygfxViewport-maybe_rebase}

```python
maybe_rebase() -> bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L338">source</a>

### `PygfxViewport.show_all` <Badge type="info" text="method" /> {#PygfxViewport-show_all}

```python
show_all() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L350">source</a>

### `PygfxViewport.look_at_point` <Badge type="info" text="method" /> {#PygfxViewport-look_at_point}

```python
look_at_point(
    xyz: np.ndarray,
    distance: float,
    view_dir: np.ndarray | None = None,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L363">source</a>

### `PygfxViewport.camera_position` <Badge type="info" text="method" /> {#PygfxViewport-camera_position}

```python
camera_position() -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L377">source</a>

The camera's float64 world position.

### `PygfxViewport.project` <Badge type="info" text="method" /> {#PygfxViewport-project}

```python
project(xyz: np.ndarray) -> tuple[float, float]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L381">source</a>

### `PygfxViewport.pixel_size` <Badge type="info" text="method" /> {#PygfxViewport-pixel_size}

```python
pixel_size(xyz: np.ndarray) -> float
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L388">source</a>

World size of one pixel at *xyz* (perspective).

### `PygfxViewport.pick` <Badge type="info" text="method" /> {#PygfxViewport-pick}

```python
pick(x: float, y: float) -> PickResult | GeometryPick | SkeletonPick | MarkerPick | None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L396">source</a>

### `PygfxViewport.highlight_polyline` <Badge type="info" text="method" /> {#PygfxViewport-highlight_polyline}

```python
highlight_polyline(points: np.ndarray) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L468">source</a>

Outline *points* (float64, (n, 3)): a selected edge or vertex.

### `PygfxViewport.on_pick` <Badge type="info" text="method" /> {#PygfxViewport-on_pick}

```python
on_pick(callback: Any) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L486">source</a>

Call `callback(pick result or None, shift)` on a double click in the view.

### `PygfxViewport.request_draw` <Badge type="info" text="method" /> {#PygfxViewport-request_draw}

```python
request_draw() -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L501">source</a>

### `PygfxViewport.draw_now` <Badge type="info" text="method" /> {#PygfxViewport-draw_now}

```python
draw_now() -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/viewport/pygfx_backend.py#L504">source</a>

Render synchronously and return the image (offscreen canvases).
