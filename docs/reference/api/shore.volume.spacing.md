---
title: shore.volume.spacing
outline: [2, 3]
---

# `shore.volume.spacing`

1D grid spacing laws for surface angle redistribution and wall-normal extrusion.

All functions are pure numpy — no scipy dependency.

Public API
----------
build_i_coords         : latitude parameter redistribution
build_j_coords         : longitude parameter redistribution (periodic axis)
build_k_steps          : wall-normal step sizes
pinned_distribution    : 1D coords with first/last interval pinned to caller-given
                         step sizes (used for C1 matching at multi-block seams)

Low-level primitives (reusable):
uniform, tanh_one_sided, tanh_two_sided, geometric

## `uniform` <Badge type="info" text="function" /> {#uniform}

```python
uniform(n: int, start: float, stop: float) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/spacing.py#L28">source</a>

Uniform distribution: n points in [start, stop].

**Parameters**

- **`n`** — Number of points.
- **`start, stop`** — Interval endpoints (inclusive).

**Returns**

- `numpy.ndarray, shape (n,)`

## `tanh_one_sided` <Badge type="info" text="function" /> {#tanh_one_sided}

```python
tanh_one_sided(n: int, start: float, stop: float, beta: float) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/spacing.py#L45">source</a>

Roberts one-sided tanh clustering toward *start*.

Dense near `start`, coarse near `stop`.

Uses the reversed Roberts formula so that the fine spacing is at
`tau=0` (the wall / cap boundary):

    t(tau) = 1 - tanh(beta \* (1 - tau)) / tanh(beta)

**Parameters**

- **`n`** — Number of points.
- **`start, stop`** — Interval endpoints.
- **`beta`** — Clustering strength (&gt; 0).  Larger values → tighter clustering. beta=1 is mild; beta=5 is aggressive.

**Returns**

- `numpy.ndarray, shape (n,)`

## `tanh_two_sided` <Badge type="info" text="function" /> {#tanh_two_sided}

```python
tanh_two_sided(
    n: int,
    start: float,
    stop: float,
    beta: float,
    cluster_at: float | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/spacing.py#L74">source</a>

Roberts two-sided tanh clustering.

**Parameters**

- **`n`** — Number of points.
- **`start, stop`** — Interval endpoints.
- **`beta`** — Clustering strength (&gt; 0).
- **`cluster_at`** — If None: symmetric clustering at both `start` and `stop`. If a value in `[start, stop]`: clustering at that *interior* point — cells are densest at `cluster_at` and coarsen monotonically away from it on both sides (a single dense band, not the two-ended default). Used for j-direction longitude clustering (e.g. leading-edge longitude on an airfoil).  When `cluster_at` lands on an endpoint the result reduces to one-sided clustering at that end.

**Returns**

- `numpy.ndarray, shape (n,)`

## `build_i_coords` <Badge type="info" text="function" /> {#build_i_coords}

```python
build_i_coords(
    ni: int,
    theta_start: float,
    theta_stop: float,
    spacing: str = 'uniform',
    beta: float = 3.0,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/spacing.py#L146">source</a>

Return *ni* latitude coordinates in `[theta_start, theta_stop]`.

**Parameters**

- **`ni`** — Number of latitude points.
- **`theta_start, theta_stop`** — Colatitude range in radians (e.g. `theta_cap` to `pi - theta_cap`).
- **`spacing`** — `'uniform'` — equal angular spacing (default, backward-compatible). `'tanh'`    — one-sided Roberts: cluster near `theta_start` (the cap-boundary row; highest-quality zone for caps). `'tanh2'`   — two-sided Roberts: cluster at both `theta_start` and `theta_stop`.
- **`beta`** — Clustering strength for tanh/tanh2 (ignored for uniform).

**Returns**

- `numpy.ndarray, shape (ni,)`

## `build_j_coords` <Badge type="info" text="function" /> {#build_j_coords}

```python
build_j_coords(
    nj: int,
    spacing: str = 'uniform',
    beta: float = 3.0,
    cluster_lon: float | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/spacing.py#L183">source</a>

Return *nj* longitude coordinates in `[0, 2*pi)`, periodic.

**Parameters**

- **`nj`** — Number of longitude points (periodic: endpoint excluded).
- **`spacing`** — `'uniform'` — equal angular spacing (default, backward-compatible). `'tanh2'`   — two-sided Roberts: cluster toward `cluster_lon` (defaults to `pi`, the midpoint, if not supplied). Note: `'tanh'` (one-sided) is not meaningful for a periodic axis.
- **`beta`** — Clustering strength for tanh2 (ignored for uniform).
- **`cluster_lon`** — Target longitude in radians for tanh2 clustering.  Only used when `spacing='tanh2'`.  Default `pi`.

**Returns**

- `numpy.ndarray, shape (nj,)` — Values in `[0, 2*pi)` with endpoint excluded.

**Notes**

For tanh2, the algorithm generates nj+1 symmetric tanh-spaced points on
`[0, 2*pi]`, then rolls the array so the dense region is centred at
`cluster_lon`.  The last (repeated) point is dropped to give a periodic
grid.

## `build_k_steps` <Badge type="info" text="function" /> {#build_k_steps}

```python
build_k_steps(
    nk: int,
    ds: float,
    growth: float = 1.15,
    spacing: str = 'geometric',
    beta: float = 3.0,
    total_thickness: float | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/spacing.py#L239">source</a>

Return `(nk-1,)` array of wall-normal step sizes.

**Parameters**

- **`nk`** — Total number of layers including k=0.
- **`ds`** — First-layer thickness (geometric law) or ignored for tanh/tanh2.
- **`growth`** — Geometric growth ratio (geometric law only).
- **`spacing`** — `'geometric'` — `ds * growth**k` (default, backward-compatible). `'tanh'`      — one-sided Roberts from 0 to `total_thickness`: cluster near the wall. `'tanh2'`     — two-sided Roberts from 0 to `total_thickness`: cluster at both wall and far-field.
- **`beta`** — Clustering strength for tanh/tanh2.
- **`total_thickness`** — Total wall-normal extent.  Required for tanh/tanh2; ignored for geometric.

**Returns**

- `numpy.ndarray, shape (nk-1,)` — Physical step sizes, all positive.

## `pinned_distribution` <Badge type="info" text="function" /> {#pinned_distribution}

```python
pinned_distribution(
    n: int,
    start: float,
    stop: float,
    first_step: float,
    last_step: float,
    spacing: str = 'uniform',
    beta: float = 3.0,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/spacing.py#L298">source</a>

Return `(n,)` coords with pinned first and last intervals.

`out[1] - out[0] == first_step` and `out[-1] - out[-2] == last_step`;
the interior `n-2` coords lie on `[start + first_step, stop - last_step]`
distributed by the chosen law.

Used for C1 spacing matching at multi-block seams: each edge whose two
endpoints share a step size with neighbouring blocks needs the first and
last cells pinned to those neighbour steps, while the interior is free
to follow a clustering law.

**Parameters**

- **`n`** — Number of points (&gt;= 2).
- **`start, stop`** — Interval endpoints (inclusive).
- **`first_step, last_step`** — Required first and last interval lengths (both &gt; 0).
- **`spacing`** — Interior law: `'uniform'`, `'tanh'`, or `'tanh2'`.
- **`beta`** — Clustering strength for tanh/tanh2 (ignored for uniform).

**Returns**

- `numpy.ndarray, shape ``(n,)```

**Raises**

- `ParameterError` — If `n < 2`, `first_step` or `last_step` is non-positive, or the two pinned intervals don't fit inside `[start, stop]`.

**Notes**

Edge cases:
    - `n == 2`: law is unused; `out = [start, stop]` and the caller
      must ensure `first_step == stop - start` (no constraint on last).
    - `n == 3`: interior collapses to a single coord at
      `start + first_step`, requiring
      `first_step + last_step == stop - start`.
