---
title: shore.mesh.spacing
outline: [2, 3]
---

# `shore.mesh.spacing`

Spacing1D: a single edge's spacing law and parameters.

## `Spacing1D` <Badge type="tip" text="dataclass" /> {#Spacing1D}

```python
Spacing1D(
    law: Literal['uniform', 'tanh', 'tanh2', 'geometric', 'pinned'],
    beta: float = 3.0,
    ds: float | None = None,
    growth: float = 1.15,
    total_thickness: float | None = None,
    cluster_at: float | None = None,
    first_step: float | None = None,
    last_step: float | None = None,
    interior: Literal['uniform', 'tanh', 'tanh2'] = 'uniform',
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/spacing.py#L21">source</a>

One-dimensional spacing law for a block edge.

**Parameters**

- **`law`** — Spacing algorithm: `'uniform'`, `'tanh'`, `'tanh2'`, `'geometric'`, or `'pinned'`.
- **`beta`** — Clustering strength for tanh/tanh2 laws (also used as the *interior* law strength when `law='pinned'` and `interior` is tanh-based).
- **`ds`** — First-layer thickness for geometric law.
- **`growth`** — Geometric growth ratio (geometric law only).
- **`total_thickness`** — Total extent required for tanh/tanh2 laws.
- **`cluster_at`** — Asymmetric clustering target for tanh2 (longitude axis use case).
- **`first_step, last_step`** — Required first and last interval lengths for `law='pinned'`.  Both must be &gt; 0.  When set on a non-pinned law they're ignored by `Spacing1D.build` directly but read by topology/primitive builders as a request for C1 endpoint matching at multi-block seams.
- **`interior`** — Interior law for `law='pinned'` (`'uniform'`, `'tanh'`, or `'tanh2'`).  Defaults to `'uniform'`.  Ignored for other laws.

### `Spacing1D.build` <Badge type="info" text="method" /> {#Spacing1D-build}

```python
build(n: int, start: float = 0.0, stop: float = 1.0) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/spacing.py#L60">source</a>

Return n coordinates on [start, stop] according to the spacing law.

For `'geometric'` this returns step sizes (n-1,), not coordinates.
For all other laws it returns coordinates (n,).

**Parameters**

- **`n`** — Number of points (or layers for geometric).
- **`start, stop`** — Interval endpoints.

**Returns**

- `numpy.ndarray` — Shape `(n,)` for uniform/tanh/tanh2/pinned; `(n-1,)` for geometric.

### `Spacing1D.build_steps` <Badge type="info" text="method" /> {#Spacing1D-build_steps}

```python
build_steps(nk: int) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/spacing.py#L105">source</a>

Return (nk-1,) wall-normal step sizes.

Convenience wrapper for the k-direction; delegates to
`volume.spacing.build_k_steps` for all laws.

**Parameters**

- **`nk`** — Total number of layers including k=0.
