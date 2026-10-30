---
title: shore.mesh.redistribute
outline: [2, 3]
---

# `shore.mesh.redistribute`

Redistributing the nodes of a block edge (ADR 0041).

The spacing law of one block edge is applied, as a distribution of arc-length
fractions, to every grid line of the block along that axis, and through the
seams to the blocks that share a face containing that axis (reversed across
a flipped seam), transitively. Every node stays on its line's old polyline
and every line keeps its end nodes, so faces across the axis do not move.

The law is taken over the picked edge's arc length: tanh / tanh2 cluster by
`beta` (`cluster_at` as a fraction), `pinned` holds `first_step` and
`last_step` (lengths along the picked edge) around an `interior` law.

## `edge_axis` <Badge type="info" text="function" /> {#edge_axis}

```python
edge_axis(edge: str) -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/redistribute.py#L31">source</a>

The axis an edge runs along: the letter its name does not fix.

## `fractions` <Badge type="info" text="function" /> {#fractions}

```python
fractions(spacing: Spacing1D, nodes: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/redistribute.py#L38">source</a>

Arc-length fractions (n,) of *spacing* over the polyline *nodes* (n, 3).

## `resample` <Badge type="info" text="function" /> {#resample}

```python
resample(lines: np.ndarray, t: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/redistribute.py#L52">source</a>

Lines (..., n, 3) resampled at arc-length fractions *t* along themselves.

Collapsed lines (zero length: polar axes) are kept; every line keeps its
two end nodes exactly.

## `plan` <Badge type="info" text="function" /> {#plan}

```python
plan(
    blocks: list[HexBlock],
    block: str,
    edge: str,
    t: np.ndarray,
) -> list[tuple[HexBlock, str, np.ndarray]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/redistribute.py#L83">source</a>

Every (block, axis, fractions) the redistribution reaches through seams.

## `redistribute` <Badge type="info" text="function" /> {#redistribute}

```python
redistribute(
    blocks: list[HexBlock],
    block: str,
    edge: str,
    spacing: Spacing1D,
) -> list[tuple[str, str]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/redistribute.py#L132">source</a>

Apply *spacing* along *edge* of *block*, in place; the (block, axis) changed.

New positions are computed from the old ones for every block before any
is written, so blocks sharing seam nodes as views are not resampled twice.
