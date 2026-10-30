---
title: shore.quality
outline: [2, 3]
---

# `shore.quality`

Structured-grid quality metrics (roadmap phase 4).

All functions take one block's node array `(ni, nj, nk, 3)` and follow the
marcher's orientation convention: a valid cell has a positive Jacobian
`(e_i x e_j) . e_k` with `k` pointing away from the wall.

- `scaled_jacobian`: per cell, the minimum over its 8 corners of
  `det[e_i, e_j, e_k] / (|e_i| |e_j| |e_k|)` with the edges leaving that
  corner (the Verdict / CUBIT definition). 1 for a rectangular cell, the
  cosine of the shear angle for a sheared one, &lt;= 0 for a folded corner.
- `k_orthogonality_deg`: per node of layer `k`, the angle between the
  k-segment to layer `k + 1` and the layer's normal `t_i x t_j`
  (central-difference tangents). 0 is orthogonal; &gt; 90 means the segment
  points back through the layer.
- `k_angle_deg`: the same angle against given directions (e.g. the body's
  surface normals), unbiased at block-edge nodes.
- `k_spacing` / `k_growth`: k-segment lengths and their ratios.
- `aspect_ratio`: per cell, longest over shortest of its 12 edges.
- `summary`: the extremes of all of these for a list of blocks.

## `scaled_jacobian` <Badge type="info" text="function" /> {#scaled_jacobian}

```python
scaled_jacobian(nodes: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/quality.py#L36">source</a>

Minimum corner scaled Jacobian per cell, shape (ni-1, nj-1, nk-1).

## `k_orthogonality_deg` <Badge type="info" text="function" /> {#k_orthogonality_deg}

```python
k_orthogonality_deg(nodes: np.ndarray, k: int = 0) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/quality.py#L68">source</a>

Angle (deg) between the k-segment `k -> k+1` and layer `k`'s normal, (ni, nj).

## `k_angle_deg` <Badge type="info" text="function" /> {#k_angle_deg}

```python
k_angle_deg(nodes: np.ndarray, directions: np.ndarray, k: int = 0) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/quality.py#L77">source</a>

Angle (deg) between the k-segment `k -> k+1` and given *directions* (ni, nj, 3).

With the body's outward surface normals at layer 0 this is the wall
orthogonality against the true surface: unlike `k_orthogonality_deg`
it is not biased at block-edge nodes, where a one-sided layer normal is
not the surface normal across the seam.

## `k_spacing` <Badge type="info" text="function" /> {#k_spacing}

```python
k_spacing(nodes: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/quality.py#L97">source</a>

Length of every k-segment, shape (ni, nj, nk-1).

## `k_growth` <Badge type="info" text="function" /> {#k_growth}

```python
k_growth(nodes: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/quality.py#L102">source</a>

Ratio of consecutive k-segment lengths, shape (ni, nj, nk-2).

## `aspect_ratio` <Badge type="info" text="function" /> {#aspect_ratio}

```python
aspect_ratio(nodes: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/quality.py#L108">source</a>

Longest over shortest of each cell's 12 edges, shape (ni-1, nj-1, nk-1).

## `summary` <Badge type="info" text="function" /> {#summary}

```python
summary(blocks: Iterable[Any]) -> dict[str, float]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/quality.py#L122">source</a>

Extremes over *blocks* (objects with `.nodes`, or node arrays).
