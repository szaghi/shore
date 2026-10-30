---
title: shore.volume.jacobian
outline: [2, 3]
---

# `shore.volume.jacobian`

Jacobian field computation for structured volume grids.

All public functions operate on grids of shape (ni, nj, nk, 3).
The j-direction is treated as periodic (O-grid topology).

## `jacobian_field` <Badge type="info" text="function" /> {#jacobian_field}

```python
jacobian_field(grid: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/jacobian.py#L12">source</a>

Compute the per-cell Jacobian scalar field for a structured volume.

Each cell Jacobian is the scalar triple product of the three edge
vectors from the cell's (i, j, k) corner:

    J = (e_i x e_j) . e_k

where e_i, e_j, e_k are the edge vectors in the i-, j-, and k-directions
respectively. J &gt; 0 for a valid cell; J ≤ 0 for an inverted cell.

The j-direction is periodic: the cell at j = nj-1 uses j = 0 as its
j+1 neighbour. The i-direction is not periodic: cells only exist for
i ∈ [0, ni-2].

**Parameters**

- **`grid`** — Structured volume, shape (ni, nj, nk, 3).

**Returns**

- `numpy.ndarray, shape (ni-1, nj, nk-1)` — Per-cell Jacobian values. There are nj cells in the j-direction (one per column, with the last column wrapping to j=0).
