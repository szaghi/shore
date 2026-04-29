# shore.volume.spacing

Pure-numpy 1D distribution kernels: low-level laws (`uniform`,
`tanh_one_sided`, `tanh_two_sided`) plus three domain builders
(`build_i_coords`, `build_j_coords`, `build_k_steps`) and the
endpoint-pinned variant `pinned_distribution`.

For the higher-level [`Spacing1D`](./api-spacing) wrapper used by
edges and topologies, see `shore.mesh.spacing`. For the algorithmic
overview, see [Spacing laws](/algorithm/spacing-laws).

::: tip Module split
The functions on this page are pure numerical kernels — no
SHORE types, no scipy. The user-facing dataclass that selects
between these laws is [`Spacing1D`](./api-spacing) in
`shore.mesh.spacing`.
:::

## Low-level laws

### `uniform`

```python
def uniform(n: int, start: float, stop: float) -> np.ndarray
```

`n` equally spaced coordinates in `[start, stop]`, both endpoints
included. Wraps `np.linspace`.

### `tanh_one_sided`

```python
def tanh_one_sided(n: int, start: float, stop: float, beta: float) -> np.ndarray
```

`n` Roberts-clustered coordinates with cells dense at `start`,
sparse at `stop`. `beta > 0` is the clustering strength
(typical: 1 = mild, 3 = moderate (default), 5 = aggressive).

### `tanh_two_sided`

```python
def tanh_two_sided(
    n: int,
    start: float,
    stop: float,
    beta: float,
    cluster_at: float | None = None,
) -> np.ndarray
```

`n` Roberts-clustered coordinates symmetric about `cluster_at`
(default: midpoint). Dense at both ends when `cluster_at = 0.5`
(in normalised parameter); shifts asymmetrically when `cluster_at`
is offset.

## Domain builders

### `build_i_coords`

```python
def build_i_coords(
    ni: int,
    theta_start: float,
    theta_stop: float,
    spacing: str = "uniform",
    beta: float = 3.0,
) -> np.ndarray
```

`ni` latitude coordinates in `[theta_start, theta_stop]` (radians).
`spacing` ∈ `"uniform"` / `"tanh"` / `"tanh2"`.

### `build_j_coords`

```python
def build_j_coords(
    nj: int,
    spacing: str = "uniform",
    beta: float = 3.0,
    cluster_lon: float | None = None,
) -> np.ndarray
```

`nj` longitude coordinates in `[0, 2*pi)` with the endpoint
excluded (periodic). `spacing` ∈ `"uniform"` / `"tanh2"`. `"tanh"`
(one-sided) is not meaningful on a periodic axis. `cluster_lon`
specifies the longitude where the tanh2 dense zone is centred
(default: `pi`).

### `build_k_steps`

```python
def build_k_steps(
    nk: int,
    ds: float,
    growth: float = 1.15,
    spacing: str = "geometric",
    beta: float = 3.0,
    total_thickness: float | None = None,
) -> np.ndarray
```

Wall-normal *step sizes* — array of shape `(nk-1,)`, all positive.

| `spacing` | Required parameters | Behaviour |
|---|---|---|
| `"geometric"` (default) | `ds`, `growth` | `step[k] = ds * growth**k`; total thickness floats. |
| `"tanh"` | `total_thickness` | `nk` coordinates by one-sided Roberts on `[0, total_thickness]`, then `np.diff`. Cluster near the wall. |
| `"tanh2"` | `total_thickness` | Two-sided Roberts; cluster at wall and far-field. |

`build_k_steps` returns step sizes (the form `Mesh.march()`
consumes), not coordinates. Coordinates would be
`np.cumsum(np.r_[0.0, steps])`.

## `pinned_distribution`

```python
def pinned_distribution(
    n: int,
    start: float,
    stop: float,
    first_step: float,
    last_step: float,
    spacing: str = "uniform",
    beta: float = 3.0,
) -> np.ndarray
```

`n` coordinates on `[start, stop]` whose **first and last
intervals exactly equal** `first_step` and `last_step`. The
remaining `n-2` interior points are distributed by the chosen
interior law (`"uniform"` / `"tanh"` / `"tanh2"`) on the reduced
interval `[start + first_step, stop - last_step]`.

This is the kernel behind [`Spacing1D(law="pinned")`](./api-spacing).
The cubed-sphere equator builder uses it directly to enforce C1
seam matching at every cap-equator seam vertex.

| Parameter | Description |
|---|---|
| `n` | Total number of points (≥ 3, since the first and last are pinned). |
| `start`, `stop` | Interval endpoints. |
| `first_step`, `last_step` | Required first / last interval lengths. Both must be > 0 and their sum must be < `stop - start`. |
| `spacing` | Interior law: `"uniform"` (default), `"tanh"`, or `"tanh2"`. |
| `beta` | Clustering strength for tanh / tanh2. |

## Examples

```python
from shore.volume.spacing import (
    build_k_steps, pinned_distribution, tanh_one_sided
)

# Geometric wall-normal: ds=1e-3, growth=1.15, 40 layers → 39 steps
steps = build_k_steps(40, ds=1e-3, growth=1.15)
print(steps[0], steps[-1])

# Tanh wall-normal: fixed total thickness 2.0
steps = build_k_steps(40, ds=0, growth=0, spacing="tanh",
                     total_thickness=2.0, beta=4.0)
print(steps.sum())  # 2.0

# Pinned 1D distribution
coords = pinned_distribution(20, 0.0, 1.0,
                             first_step=0.05, last_step=0.10,
                             spacing="tanh", beta=4.0)
print(coords[1] - coords[0],   coords[-1] - coords[-2])
# 0.05 0.10
```

## See also

- [`shore.mesh.spacing.Spacing1D`](./api-spacing) — the per-edge
  wrapper.
- [Spacing laws](/algorithm/spacing-laws) — algorithmic overview.
- [Per-edge spacing](/guide/per-edge-spacing) — task-oriented
  walkthrough.
