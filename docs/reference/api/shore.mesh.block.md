---
title: shore.mesh.block
outline: [2, 3]
---

# `shore.mesh.block`

Block ABC and HexBlock: structured grid block classes.

## `Block` <Badge type="tip" text="class" /> {#Block}

```python
class Block(ABC)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L15">source</a>

Abstract base for a structured multi-dimensional node block.

Subclasses define the concrete topology (number of vertices, edges,
faces) and provide named accessors.  The only shared state is `nodes`
(the coordinate array) and `label`.

### `Block.n_vertices` <Badge type="info" text="property" /> {#Block-n_vertices}

```python
n_vertices: int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L28">source</a>

### `Block.n_edges` <Badge type="info" text="property" /> {#Block-n_edges}

```python
n_edges: int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L32">source</a>

### `Block.n_faces` <Badge type="info" text="property" /> {#Block-n_faces}

```python
n_faces: int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L36">source</a>

### `Block.faces` <Badge type="info" text="method" /> {#Block-faces}

```python
faces() -> list[Face]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L39">source</a>

### `Block.edges` <Badge type="info" text="method" /> {#Block-edges}

```python
edges() -> list[Edge]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L42">source</a>

### `Block.vertices` <Badge type="info" text="method" /> {#Block-vertices}

```python
vertices() -> list[Vertex]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L45">source</a>

### `Block.wall_face` <Badge type="info" text="method" /> {#Block-wall_face}

```python
wall_face() -> Face
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L48">source</a>

The k=0 body-wall face.

### `Block.outer_face` <Badge type="info" text="method" /> {#Block-outer_face}

```python
outer_face() -> Face
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L53">source</a>

The k=nk-1 outer boundary face.

### `Block.face` <Badge type="info" text="method" /> {#Block-face}

```python
face(name: FaceName) -> Face
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L57">source</a>

Look up a face by canonical name.

### `Block.edge` <Badge type="info" text="method" /> {#Block-edge}

```python
edge(name: EdgeName) -> Edge
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L64">source</a>

Look up an edge by canonical name.

## `HexBlock` <Badge type="tip" text="class" /> {#HexBlock}

```python
HexBlock(
    label: str,
    nodes: np.ndarray,
    periodic_i: bool = False,
    periodic_j: bool = True,
    spacing_i: Spacing1D | None = None,
    spacing_j: Spacing1D | None = None,
    spacing_k: Spacing1D | None = None,
)
```
Bases: `Block`<br>
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L72">source</a>

Structured hexahedral block: (ni, nj, nk, 3) node array.

Topology
--------
- 8 vertices at the 8 corners.
- 12 edges: 4 along each axis direction.
- 6 faces: i_lo, i_hi, j_lo, j_hi, k_lo (wall), k_hi (outer).

The `nodes` array is owned externally (e.g. a slice of a topology's
ring allocation) so that adjacent blocks can share numpy views.

**Parameters**

- **`label`** — Human-readable block identifier (e.g. `'sub0'`, `'cap_north'`).
- **`nodes`** — Coordinate array, shape `(ni, nj, nk, 3)`.  May be a view.
- **`periodic_i, periodic_j`** — Periodicity flags forwarded to marching and smoothing.
- **`spacing_i, spacing_j, spacing_k`** — Default spacing laws for the three axis directions.  Individual edges may override by carrying a different `Spacing1D` instance.

### `HexBlock.n_vertices` <Badge type="info" text="property" /> {#HexBlock-n_vertices}

```python
n_vertices: int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L179">source</a>

### `HexBlock.n_edges` <Badge type="info" text="property" /> {#HexBlock-n_edges}

```python
n_edges: int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L183">source</a>

### `HexBlock.n_faces` <Badge type="info" text="property" /> {#HexBlock-n_faces}

```python
n_faces: int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L187">source</a>

### `HexBlock.faces` <Badge type="info" text="method" /> {#HexBlock-faces}

```python
faces() -> list[Face]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L190">source</a>

### `HexBlock.edges` <Badge type="info" text="method" /> {#HexBlock-edges}

```python
edges() -> list[Edge]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L193">source</a>

### `HexBlock.vertices` <Badge type="info" text="method" /> {#HexBlock-vertices}

```python
vertices() -> list[Vertex]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L196">source</a>

### `HexBlock.wall_face` <Badge type="info" text="method" /> {#HexBlock-wall_face}

```python
wall_face() -> Face
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L199">source</a>

### `HexBlock.outer_face` <Badge type="info" text="method" /> {#HexBlock-outer_face}

```python
outer_face() -> Face
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L202">source</a>

### `HexBlock.periodic` <Badge type="info" text="property" /> {#HexBlock-periodic}

```python
periodic: tuple[bool, bool]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/block.py#L275">source</a>
