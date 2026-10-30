---
title: shore.mesh.orthogonal
outline: [2, 3]
---

# `shore.mesh.orthogonal`

Grid lines leaving the walls along their normal (ADR 0048) and wall layers (ADR 0050).

After a premesh (ADR 0042) the lines across a block run from the wall to the
opposite face along the transfinite interpolation, so they leave a curved
wall obliquely. `orthogonal_walls` reshapes every line of each block with
exactly one WALL face into a cubic Hermite curve that leaves the wall node
along the wall's normal (with tangent `strength` times the line's chord)
and meets the opposite node with the line's old end direction; the nodes
keep their arc-length fractions, the wall and the opposite face stay put.

The normal of a wall node is the mean of the grid normals of every wall face
holding it (each oriented into its block), keyed by the node's exact
coordinates: blocks sharing a wall node use the same normal, and blocks
sharing a face build its lines identically, so seams stay closed.

## `grid_normals` <Badge type="info" text="function" /> {#grid_normals}

```python
grid_normals(surface: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/orthogonal.py#L34">source</a>

Unit normals (na, nb, 3) of a structured surface grid (na, nb, 3).

## `wall_normals` <Badge type="info" text="function" /> {#wall_normals}

```python
wall_normals(blocks: list[HexBlock]) -> dict[bytes, np.ndarray]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/orthogonal.py#L40">source</a>

Per wall node, by its coordinates' bytes: the mean grid normal of the
wall faces holding it, each oriented into its block.

## `hermite_lines` <Badge type="info" text="function" /> {#hermite_lines}

```python
hermite_lines(
    lines: np.ndarray,
    normals: np.ndarray,
    strength: float,
    samples: int = 64,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/orthogonal.py#L63">source</a>

Lines (m, n, 3), wall node first, reshaped to leave along *normals* (m, 3).

Each line becomes the cubic Hermite curve from its wall node (tangent
`strength * chord * normal`) to its last node (tangent `chord` times
the old last segment's direction); its nodes keep their arc-length
fractions, measured on the old line and placed on the curve.

## `orthogonal_walls` <Badge type="info" text="function" /> {#orthogonal_walls}

```python
orthogonal_walls(blocks: list[HexBlock], strength: float = 1.0) -> list[str]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/orthogonal.py#L103">source</a>

Reshape, in place, the lines of every block with exactly one WALL face
to leave the wall along its normal; the labels of the blocks changed.

## `layer_fractions` <Badge type="info" text="function" /> {#layer_fractions}

```python
layer_fractions(length: float, n: int, first_step: float) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/orthogonal.py#L141">source</a>

Arc-length fractions (n,) of geometric cells from the wall, the first
*first_step* long, over a line of *length* (one line at a time, so a
line's answer does not depend on the others).

## `wall_layer` <Badge type="info" text="function" /> {#wall_layer}

```python
wall_layer(blocks: list[HexBlock], first_step: float) -> list[str]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/orthogonal.py#L157">source</a>

Redistribute, in place, the nodes of every line across the wall of each
block with one WALL face: geometric cells from the wall, the first
*first_step* long, on the line's own polyline; the labels changed.
