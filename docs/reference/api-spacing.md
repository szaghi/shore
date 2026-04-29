# shore.mesh.spacing

`Spacing1D` — the per-edge spacing law data class — and the
companion [`Edge`](#edge) descriptor that wires a `Spacing1D` onto
each of a block's 12 axis-parallel edges. Together they are the
public Python knobs for surface and wall-normal node distribution.

For the underlying numerics (`uniform`, `tanh`, `tanh2`,
`build_k_steps`, etc.) see [`shore.volume.spacing`](./api-volume-spacing).
For an algorithmic overview, see [Spacing laws](/algorithm/spacing-laws);
for the per-edge usage pattern,
[Per-edge spacing](/guide/per-edge-spacing).

## `Spacing1D`

```python
@dataclass
class Spacing1D:
    law:              Literal["uniform", "tanh", "tanh2", "geometric", "pinned"]
    beta:             float = 3.0
    ds:               float | None = None
    growth:           float = 1.15
    total_thickness:  float | None = None
    cluster_at:       float | None = None
    first_step:       float | None = None
    last_step:        float | None = None
    interior:         Literal["uniform", "tanh", "tanh2"] = "uniform"
```

A single edge's spacing law and parameters. Five laws are
supported:

| Law | Returns | Required parameters |
|---|---|---|
| `"uniform"` | `n` equally spaced coordinates on `[start, stop]`. | — |
| `"tanh"` | `n` one-sided Roberts-clustered coordinates. | `beta` (clustering strength). |
| `"tanh2"` | `n` two-sided Roberts-clustered coordinates. | `beta`; `cluster_at` (optional asymmetric target). |
| `"geometric"` | `n-1` step sizes growing as `ds * growth^k`. **Returns step sizes, not coordinates.** | `ds`; `growth`. |
| `"pinned"` | `n` coordinates with caller-chosen first and last interval lengths and an interior law. Used for C1 endpoint matching at multi-block seams. | `first_step`, `last_step`; `interior` (default `"uniform"`); `beta` (when `interior` is tanh-based). |

Tanh / tanh2 *also* accept `total_thickness` when used as a
wall-normal step generator via `build_steps()`.

### Methods

#### `Spacing1D.build`

```python
def build(self, n: int, start: float = 0.0, stop: float = 1.0) -> np.ndarray
```

Generic 1D distribution on `[start, stop]`. Returns shape `(n,)`
for `uniform` / `tanh` / `tanh2` / `pinned`; shape `(n-1,)` for
`geometric` (step sizes, the natural form for marching).

#### `Spacing1D.build_steps`

```python
def build_steps(self, nk: int) -> np.ndarray
```

Convenience wrapper for the wall-normal direction. Returns
`(nk-1,)` step sizes regardless of `law`. For `tanh` / `tanh2` it
internally requires `total_thickness`.

This is the form `Mesh.march()` consumes when a block carries a
user-overridden k-spacing on its edges; see
[Per-edge spacing](/guide/per-edge-spacing).

### Examples

**Geometric wall-normal spacing** — the default for `nk` layers:

```python
from shore.mesh.spacing import Spacing1D

sp = Spacing1D(law="geometric", ds=1e-3, growth=1.15)
steps = sp.build_steps(40)   # (39,) step sizes
```

**Tanh wall-normal spacing** — fixed total thickness, cluster at
the wall:

```python
sp = Spacing1D(law="tanh", beta=4.0, total_thickness=2.0)
steps = sp.build_steps(40)
# sum(steps) == 2.0; first step << last step
```

**Pinned 1D distribution** — first / last interval lengths fixed,
interior tanh-clustered:

```python
sp = Spacing1D(
    law="pinned",
    first_step=0.05,
    last_step=0.10,
    interior="tanh",
    beta=4.0,
)
coords = sp.build(20, start=0.0, stop=1.0)
# coords[1] - coords[0]  == 0.05
# coords[-1] - coords[-2] == 0.10
```

This is the law the cubed-sphere equator builder uses internally
to match the cap edge's first cell at every seam vertex (C1 seam
pin). The `pinned` law makes that machinery available as public
API for any topology or primitive that needs C1 endpoint matching.

## `Edge`

```python
@dataclass
class Edge:
    index:           int
    name:            EdgeName
    axis:            Literal["i", "j", "k"] | None
    vertex_indices:  tuple[int, int]
    n_nodes:         int
    spacing:         Spacing1D
    body_projected:  bool = False
    first_step:      float | None = None
    last_step:       float | None = None
```

One edge of a structured block. A `HexBlock` has **12 edges**: 4
along each of the three axes. The edge's `name` encodes the *two
constant axes* (e.g. `"j_lo_k_lo"` runs along i with j=0, k=0); the
`axis` field gives the free direction.

| Edge name | Axis | Block slice | Endpoints |
|---|---|---|---|
| `j_lo_k_lo` | i | `nodes[:, 0, 0, :]` | sub-block south-pole, j=0, k=0 |
| `j_hi_k_lo` | i | `nodes[:, -1, 0, :]` | south-pole, j=nj-1, k=0 |
| `j_lo_k_hi` | i | `nodes[:, 0, -1, :]` | south-pole, j=0, k=nk-1 |
| `j_hi_k_hi` | i | `nodes[:, -1, -1, :]` | south-pole, j=nj-1, k=nk-1 |
| `i_lo_k_lo` | j | `nodes[0, :, 0, :]` | i=0, k=0 |
| `i_hi_k_lo` | j | `nodes[-1, :, 0, :]` | i=ni-1, k=0 |
| `i_lo_k_hi` | j | `nodes[0, :, -1, :]` | i=0, k=nk-1 |
| `i_hi_k_hi` | j | `nodes[-1, :, -1, :]` | i=ni-1, k=nk-1 |
| `i_lo_j_lo` | k | `nodes[0, 0, :, :]` | the 4 wall-normal edges |
| `i_hi_j_lo` | k | `nodes[-1, 0, :, :]` |  |
| `i_lo_j_hi` | k | `nodes[0, -1, :, :]` |  |
| `i_hi_j_hi` | k | `nodes[-1, -1, :, :]` |  |

`edge.nodes(block)` returns a `(n_nodes, 3)` view into the block's
coordinate array.

### Per-edge spacing override

`HexBlock.__init__` assigns the *same* `Spacing1D` instance to all
edges of a given axis (`spacing_i`, `spacing_j`, `spacing_k`). To
override the spacing on a subset of edges, replace the
`edge.spacing` attribute post-construction:

```python
from shore.mesh.spacing import Spacing1D

new_sp_k = Spacing1D(law="tanh", beta=4.0, total_thickness=2.0)
for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
    block.edge(name).spacing = new_sp_k
```

`Mesh.march()` checks block 0's k-edges by **identity** against the
block's stored default. If all 4 k-edges were replaced with the
same instance, the new schedule is used; partial replacement emits
a `UserWarning` and the block falls back to the default. See
[`Mesh.march`'s k-step precedence](./api-mesh#mesh-march) and the
[per-edge spacing guide](/guide/per-edge-spacing) for the full
recipe.

### `first_step` / `last_step`

`Edge.first_step` and `Edge.last_step` request **C1 endpoint
matching** on a per-vertex basis: when set, the first / last
interval along the edge has the chosen physical length, with the
interior shaped by `edge.spacing`. Topology and primitive builders
read these fields and apply a precedence rule at corners where
multiple incident edges request different lengths (cap / body
geometry wins; the overridden edges are listed in a
`UserWarning`).

Setting these does *not* require `spacing.law == "pinned"`: the
consumer is free to read the steps and switch laws as needed.

## See also

- [`shore.mesh.Mesh`](./api-mesh) — `Mesh.march`'s k-step precedence.
- [`shore.mesh.topology`](./api-topology) — `HexBlock`, `Edge`'s home.
- [`shore.volume.spacing`](./api-volume-spacing) — the lower-level
  numerics: `uniform`, `tanh_one_sided`, `tanh_two_sided`,
  `pinned_distribution`, `build_k_steps`, `build_i_coords`,
  `build_j_coords`.
- [Spacing laws](/algorithm/spacing-laws) — algorithmic overview.
- [Per-edge spacing](/guide/per-edge-spacing) — task-oriented
  walkthrough.
