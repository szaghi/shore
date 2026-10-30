---
title: shore.surface.equator
outline: [2, 3]
---

# `shore.surface.equator`

Equator band construction from cap seam boundaries.

Given the two cap k=0 layers (north and south), the equator's k=0 layer is
built so its first/last rows match the cap seams bit-exactly, and its
interior follows great-circle meridians from each south seam node to the
corresponding north seam node.

The cap is built first (gnomonic flat-square seed); the equator inherits the
cap-edge node distribution along its seam (this is how `j_spacing` flows
from cap to equator).  The equator's `i_spacing` controls only the
distribution along each meridian, so the cap and equator are decoupled in
the i-direction and coupled in the j-direction.

Public API
----------
build_equator(cap_south, cap_north, ni, nj, projector, i_spacing, i_beta)
    Return `(ni, nj, 3)` equator k=0 layer.

## `build_equator` <Badge type="info" text="function" /> {#build_equator}

```python
build_equator(
    cap_south: np.ndarray,
    cap_north: np.ndarray,
    ni: int,
    theta_cap: float,
    projector: object,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/equator.py#L118">source</a>

Return the `(ni, nj, 3)` equator k=0 layer.

Row `i=0` is the south seam (cap_south boundary, traversed in sub-block
order); row `i=ni-1` is the north seam.  Interior rows are built by
slerp along great-circle meridians from each south seam node to the
corresponding north seam node.  Each meridian is then ray-cast through
the cached Projector so its interior nodes lie on the body surface.

**Parameters**

- **`cap_south, cap_north`** — Cap k=0 layers, shape `(nc, nc, 3)` each.
- **`ni`** — Number of equator latitude rows (including both seam rows).
- **`theta_cap`** — Polar exclusion angle in radians used to build the cap blocks.  Drives the meridional slerp interval `[theta_cap, pi - theta_cap]` for the i-spacing law.
- **`projector`** — `Projector` with cached BVH and anchor.
- **`i_spacing`** — Latitude spacing law (`'uniform'`, `'tanh'`, `'tanh2'`).
- **`i_beta`** — Clustering strength for tanh/tanh2.

**Returns**

- `numpy.ndarray, shape ``(ni, nj, 3)`` with ``nj = 4*(nc-1)```
