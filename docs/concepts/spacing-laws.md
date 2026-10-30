# Spacing laws

A spacing law places $n$ nodes on an interval. SHORE uses the same small
set of laws for every direction: the surface i and j directions of the
body-fitted topologies, the wall-normal step schedule, the primitives,
the channel blocks of the hybrids, and edge redistribution. The 1D kernels
are in [`shore.volume.spacing`](/reference/api/shore.volume.spacing)
(pure numpy); [`Spacing1D`](/reference/api/shore.mesh.spacing#Spacing1D)
is the object every block edge carries.

<Figure name="spacing-laws" ext="svg" caption="Node distributions of uniform, tanh, tanh2 and geometric on one axis, with the cell size along the axis below." />

## The laws

Write $\tau = i/(n-1) \in [0, 1]$ and map $t(\tau) \in [0,1]$ onto
$[x_0, x_1]$ as $x = x_0 + (x_1 - x_0)\,t$.

### `uniform`

$t = \tau$. Equal cells; the default for the surface directions and the
primitives.

### `tanh` (one-sided Roberts)

Cells cluster toward the **start** of the interval:

$$
t(\tau) = 1 - \frac{\tanh\bigl(\beta(1 - \tau)\bigr)}{\tanh\beta} .
$$

$\beta > 0$ sets the strength (default 3.0): about 1 is mild, 5 is
aggressive.

### `tanh2` (two-sided Roberts)

Cells cluster toward **both ends**:

$$
t(\tau) = \frac12\left(1 + \frac{\tanh\bigl(\beta(\tau - \tfrac12)\bigr)}{\tanh(\beta/2)}\right).
$$

With `cluster_at` set to an interior point $c$, the law instead builds
one dense band at $c$: two one-sided `tanh` segments, $[x_0, c]$ dense at
its right end and $[c, x_1]$ dense at its left end, with the nodes split
in proportion to the two lengths (at least 2 each). A `cluster_at` on an
end point reduces to one-sided clustering at that end.

### `geometric`

Wall-normal only: step sizes $\Delta s_k = \Delta s_0\, r^k$,
$k = 0 \dots n-2$. The total thickness follows from $\Delta s_0$, $r$
and $n$; it is not prescribed.

### `pinned`

The first and last intervals are fixed to given lengths and the interior
$n - 2$ nodes follow an interior law (`uniform`, `tanh` or `tanh2`) on
$[x_0 + \Delta s_{first},\ x_1 - \Delta s_{last}]$:

$$
x_1 - x_0 = \Delta s_{first}, \qquad x_{n-1} - x_{n-2} = \Delta s_{last}.
$$

This is the construction behind C1 spacing across block seams: the
cubed-sphere equator pins each meridian's end intervals to the cap edge
([k=0 construction](/concepts/k0-construction#per-j-c1-seam-pinning)).
`ParameterError` is raised for $n < 2$, non-positive steps, or two pins
that do not fit in the interval.

```python
from shore.mesh.spacing import Spacing1D

sp = Spacing1D(law="pinned", first_step=0.05, last_step=0.10, interior="tanh", beta=4.0)
x = sp.build(8, start=0.0, stop=1.0)  # x[1] - x[0] == 0.05, x[-1] - x[-2] == 0.10
```

## Which law where

| Direction | `shore mesh` flag | Laws | Default |
|---|---|---|---|
| Surface i (meridian, cubed sphere and O-grid) | `--surface-i-spacing`, `--surface-i-beta` | `uniform`, `tanh`, `tanh2` | `uniform`, 3.0 |
| Surface j (cap edge / longitude) | `--surface-j-spacing`, `--surface-j-beta` | `uniform`, `tanh2` | `uniform`, 3.0 |
| Wall-normal k | `--volume-k-spacing`, `--volume-k-beta`, `--volume-k-thickness` | `geometric`, `tanh`, `tanh2` | `geometric`, 3.0 |
| Channel blocks (OH, CH) | `--k-law-channel`, `--k-beta-channel` | `uniform`, `tanh`, `tanh2` | `uniform`, 3.0 |

### Surface i

The meridian runs from the south cap seam (`i = 0`) to the north cap seam
(`i = ni − 1`) in colatitude. `tanh` clusters toward the south seam only;
`tanh2` toward both seams. On the cubed sphere the end intervals are
pinned to the cap edge when the steps are compatible, and the law shapes
the interior.

### Surface j

On the cubed sphere the j law sets the distribution along the cap edges
on $[-\tan\theta_{cap},\ +\tan\theta_{cap}]$, and the equator inherits it
through the seam. `tanh2` clusters toward the cap corners (the 3-block
junctions). `tanh` is not offered: the law must be symmetric across the
four quadrants. (`--j-cluster-lon` is deprecated and ignored.)

### Wall-normal k

`k = 0` is the wall, `k = nk − 1` the outer layer.

| Law | Fine cells | Fixes | Needs |
|---|---|---|---|
| `geometric` | at the wall | the first cell | `--ds` (0.001), `--growth` (1.15) |
| `tanh` | at the wall | the total thickness | `--volume-k-thickness` |
| `tanh2` | at the wall and the outer layer | the total thickness | `--volume-k-thickness` |

With `tanh`/`tanh2` the first cell follows from `nk`, the thickness and
$\beta$; `--ds` and `--growth` are ignored. Choose `geometric` when the
first cell (the $y^+$ target) matters most and `tanh` when the outer
distance must be fixed. All blocks of a mesh share one k schedule
([Marching](/concepts/hyperbolic-marching#step-schedule)). Recipes:
[Set the wall spacing](/cookbook/wall-spacing).

## Per-edge spacing

By default a block applies one law per axis to its four parallel edges.
Each of the twelve edges of a [`HexBlock`](/concepts/block-model) can
instead carry its own `Spacing1D`. The block stamps the same instance on
the four edges of an axis, so an override is detected by **identity**
(`edge.spacing is not block._default_spacing_*`).

When the four parallel edges of axis i disagree, the primitives blend
them by transfinite interpolation over the transverse face:

$$
t(i, j, k) = (1-\eta)(1-\zeta)\,t_{00}(i) + \eta(1-\zeta)\,t_{10}(i)
+ (1-\eta)\zeta\,t_{01}(i) + \eta\zeta\,t_{11}(i),
$$

$\eta = j/(n_j-1)$, $\zeta = k/(n_k-1)$. With four equal edges this is the
plain 1D distribution.

| Where | i / j | k |
|---|---|---|
| `BoxMesh` (`edge_pins=`) | TFI | TFI |
| `AnnularCylinderMesh` (`edge_pins=`) | TFI on r; on θ only for a partial sector | TFI |
| `FlatCapsCylinder` | no (the radial and arc directions are fixed by the construction) | no |
| O-grid | no (j is periodic) | k-edge override read by `Mesh.march` |
| Cubed sphere | no (the seam pin is geometric, not user-driven) | k-edge override read by `Mesh.march`, one schedule for all 6 blocks |

`edge_pins` is a Python-API feature; the CLI does not expose it. For
existing grids, [`Project.redistribute`](/cookbook/redistribute) resamples
the lines along one edge with any law and carries the change through
connected faces ([ADR 0041](/adr/0041-redistribute-edge)). Worked
examples: [Per-edge spacing](/cookbook/per-edge-spacing).

## Choosing $\beta$

| $\beta$ | Effect |
|---|---|
| 1–2 | mild clustering |
| 3 (default) | moderate |
| 4–5 | strong: the cells at the coarse end grow large |
| > 5 | usually too strong: coarse cells degenerate |

Check the result with [`shore check`](/reference/check) and the
[quality metrics](/concepts/quality): a sharp drop in the scaled Jacobian
where fine and coarse cells meet means $\beta$ is too large.

## See also

- [`shore mesh`](/reference/mesh): every spacing flag
- [k=0 construction](/concepts/k0-construction): how the i and j laws meet at the cap seam
- [Marching](/concepts/hyperbolic-marching): the step schedule in use
- [`shore.volume.spacing`](/reference/api/shore.volume.spacing), [`Spacing1D`](/reference/api/shore.mesh.spacing#Spacing1D)
