# Intersection curves

A collar grid starts from the curve where two components meet: a wing
root on a fuselage, a strut on a plate, an appendage on a hull. SHORE
traces that curve from a seed point using nothing but closest-point
projection onto the two bodies, so the same code works on STL meshes,
CAD models and analytic surfaces.

There are two kernels in `shore.surface.intersect`:

| Function | Curve | Use |
|---|---|---|
| `intersection_loop(body_a, body_b, seed, n)` | closed loop | whole bodies; the j direction of a closed collar |
| `intersection_curve(body_a, body_b, seed, n, *, plane=)` | open curve, both ends on a symmetry plane | half models; the j direction of a half-collar |

Both return an `(n, 3)` array of nodes on both bodies, ordered along
$\mathbf n_b \times \mathbf n_a$ (the bodies' outward normals). That order
is what the collar needs: `march_surface(A, curve)` marches on A away
from B, and `march_surface(B, curve[::-1])` on B away from A.

<Figure name="intersection-curve" caption="Two bodies and their intersection loop, traced from a seed point." />

## Tracing

1. **Onto the curve.** Alternating projections onto A and B take the
   seed onto the intersection. If they do not converge to within the
   on-curve tolerance ($10^{-9}$ of the model diagonal), there is no
   intersection near the seed and the call is refused. When the bodies
   meet in several loops, the seed selects one, and the loop's node 0 is
   where the seed lands.
2. **Predictor–corrector.** The tracer steps along the tangent
   $\mathbf t = \mathbf n_b \times \mathbf n_a$ by $4 \cdot 10^{-3}$ of the
   model diagonal, then corrects each step back onto both bodies by
   alternating projections (traced points converged to $10^{-10}$).
   A loop is traced until it returns to its start.
3. **Resampling.** The traced polyline is resampled to $n$ nodes (see
   *Node distribution*) and every node is polished back onto both bodies
   to $10^{-13}$ of the diagonal, all nodes at once.

Where the normals are parallel there is no tracing direction: grazing
(tangential) intersections are out of scope and raise `ParameterError`.
A sphere resting on a plane is refused within seconds, because
alternating projections converge too slowly at a tangent contact to reach
the on-curve tolerance.

On the four analytic junctions of the acceptance spec (plane, tilted
plane, cylinder and sphere against a cylinder, $n = 64$) the nodes lie on
both bodies within $10^{-12}$, chords are uniform to 1.0015, and the
length matches the exact curve to $-0.040\,\%$, which is the inscribed
64-gon's own deficit $(\pi/64)^2/6$.

### Open curves on a symmetry plane

On a half model the bodies stop at the symmetry plane, and the junction
runs from the plane back to it. `intersection_curve` traces from the seed
in both directions until a step reaches the plane **within the on-curve
tolerance**: a projection onto a half body never goes past the plane, so
a strict "has crossed" test would never fire and the trace would stall.
Each end is then put exactly on the plane by Newton steps along the
curve tangent, polished onto both bodies. The plane is
`(point, normal)`; the side of the seed is kept, so the normal's sign
does not matter, and the bodies may be whole or cut at the plane.

The curve is not required to meet the plane at a right angle (on truly
mirror-symmetric bodies it does). The half-collar warns when it does not;
see [Collar blocks](/concepts/collars).

## Node distribution

By default the nodes are equally spaced in arc length. Two controls
change that, and they compose:

- **`spacing`**, a `Spacing1D` law applied to a computational coordinate:
  `uniform` or `tanh2` (with an optional `cluster_at` in $[0, 1]$) on a
  loop; `uniform`, `tanh` or `tanh2` on an open curve. `tanh` is not
  allowed on a loop: it would jump in spacing at node 0.
- **`curvature`**, a weight $\alpha \ge 0$ that maps the computational
  coordinate to arc length by equidistributing the monitor

  $$
  m(s) = 1 + \alpha\,\kappa(s)\,\frac{L}{2\pi},
  $$

  where $\kappa$ is the curvature and $L$ the curve length. At
  $\alpha = 0$ the nodes are equally spaced; on a circle $m$ is constant,
  so a circular loop stays uniform for any $\alpha$. The weight is
  dimensionless. Only with $\alpha = 0$ does the law hold exactly by arc
  length.

To make the curvature meaningful the traced polyline is first refined:
segments next to a vertex that turns by more than 2° are bisected and the
new points polished onto both bodies, down to $10^{-5}$ of the diagonal.
Curvature is then averaged over half a uniform cell, never taken at a
single vertex: on a faceted STL every facet edge is a spike, and vertex
curvature folded collars at $\alpha = 3$. On an open curve the curvature
is measured on the curve together with its mirror image, so the ends on
the plane see the bend the plane cuts through.

Adjacent cells are limited to a ratio of 1.2. Measured on a plate pierced
by a $1 \times 0.2$ elliptic post (curvature ratio 125 between ends and
flanks), 128 nodes:

| | uniform | $\alpha = 1$ | $\alpha = 3$ | $\alpha = 10$ |
|---|---|---|---|---|
| largest turn per cell | 42.6° | 5.5° | 4.4° | 4.0° |
| collar minimum scaled Jacobian | 0.008 | 0.23 | 0.25 | 0.27 |

$\alpha = 1$ is the suggested starting point for rounded sections. With too
few nodes the 1.2 ratio caps how much clustering is possible.

The defaults (`spacing=None`, `curvature=0.0`) keep the plain arc-length
resampling bit for bit, so saved projects replay unchanged;
`Spacing1D("uniform")` is the same.

## Corners

A junction crossing an edge of either body (a trailing edge through a
plate, a prismatic strut, a body pierced across an edge) has a corner,
where the tangent jumps. A plain step overshoots the edge and projection
pulls it back, so a naive tracer stalls. SHORE detects corners from the
trace itself, without CAD edges or STL crease detection:

- **Failed step.** A step fails when it advances less than half a step, or
  its direction or the tangent turns by 30° or more (`corner_deg`).
- **Adaptive walk.** A failed step is redone as a walk of 1/16 steps,
  recursively. A failure that persists at $10^{-6}$ of the diagonal
  (`corner_scale`) is a corner. A walk that finds none was a tight smooth
  bend, and the original step is kept, so smooth traces are bit-identical
  to the tracer without corner handling.
- **Leaving the corner.** The curve's points at distance $r$ are found by
  alternating projections onto A, B and the sphere of radius $r$, from 26
  directions; the way out is the one not on the incoming branch.
- **Locating it.** The tangent lines of the two branches are crossed at
  points ever closer to the corner ($r$ divided by 8 per round). The
  corner is exact to round-off: $2 \cdot 10^{-15}$ on analytic bodies, the
  STL vertex exactly, $2 \cdot 10^{-11}$ on a CAD spline wing.

Every corner becomes a fixed node, like node 0 (or the two plane ends).
With `curvature = 0`, each segment between fixed nodes gets cells in
proportion to its length (at least one, by largest remainders) and the law
applies within it. Too few nodes for the segments raises `ParameterError`
naming the corners.

### Grading across corners

With `curvature > 0` the spacing is graded over the whole curve, corners
included. The target spacing $g(s) = C / m(s)$ is capped by cones of slope
$\ln 1.2$ in both directions (wrapping round for a loop), and $C$ is found
by bisection so that the curve holds exactly the required number of cells.
Nodes equidistribute $1/g$, integrated on a sampling at a quarter of the
smallest spacing; segments get cells by their share of that integral. On
the STL NACA 0012 wings through a plate (sharp and blunt TE), adjacent
cells differ by at most 1.198 and the corners stay exact. Without grading
the blunt base (0.0025 wide) was one tiny cell between cells 8.7× larger.

Curves without corners keep the earlier ratio limiter, which limits the
cells and then rescales to the length, so their output is unchanged.

## Parameters

`TraceParams` holds the tracer's settings; lengths are relative to the
model diagonal. They are not per-case tuning knobs.

| Field | Default | Meaning |
|---|---|---|
| `step` | 4e-3 | predictor step |
| `tolerance` | 1e-13 | convergence of the output nodes |
| `trace_tolerance` | 1e-10 | convergence of the traced points |
| `on_curve` | 1e-9 | a point is on the curve this close to both bodies |
| `polish_iterations` | 200 | alternating projections per polish |
| `polish_period` | 12 | longest exact cycle of polish iterates detected |
| `max_steps` | 200 000 | trace step limit |
| `refine_turn_deg` | 2 | refine where a vertex turns more |
| `refine_min` | 1e-5 | shortest refined segment |
| `curvature_window` | 0.5 | curvature averaged over this fraction of a uniform cell |
| `max_cell_ratio` | 1.2 | adjacent cells of the curvature distribution |
| `corner_deg` | 30 | a step turning more is checked for a corner |
| `corner_scale` | 1e-6 | a step failing at this size is a corner |
| `ring_iterations` | 200 | projections onto the curve and a sphere |

**Polish cycles.** On sliver STL triangles, projection noise can trap a
polished point in an exact cycle (period 2, 5, ...) that never meets the
tolerance. The polish keeps its last 12 iterates; when a new one repeats
one of them bit for bit, it returns what the remaining iterations would
have returned. Results are unchanged; a 1000-point-section wing loop went
from 48 s to 8 s. `polish_period = 0` turns the detection off.

## Example

```python
import numpy as np

from shore.geometry.stl import StlGeometry
from shore.mesh.spacing import Spacing1D
from shore.shapes import cylinder, sphere
from shore.surface.intersect import intersection_curve, intersection_loop

ball = StlGeometry(sphere(radius=1.0, subdivisions=4))
post = StlGeometry(cylinder(radius=0.4, height=3.0, sections=64))
seed = np.array([0.4, 0.0, 0.9])

loop = intersection_loop(ball, post, seed, 96, curvature=1.0)
plane = (np.zeros(3), np.array([0.0, 1.0, 0.0]))  # y = 0; the seed's side is kept
half = intersection_curve(
    ball, post, np.array([0.0, 0.4, 0.9]), 49, plane=plane, spacing=Spacing1D("tanh2", beta=2.0)
)
print(loop.shape, half.shape, half[[0, -1], 1])  # ends on y = 0
```

From the CLI and projects the same controls are `--loop-law`,
`--loop-beta`, `--loop-cluster-at` and `--loop-curvature` on
`shore collar`, and `loop_*` on `Project.build_collar`.

## Limits

- **Grazing intersections** (parallel normals) raise; they are not traced.
- **One loop per seed.** There is no automatic search for every loop.
- **One plane** for open curves: curves ending on two different planes
  (a pylon, a tail) or on the open edge of a body are out of scope. Curves
  tangent to the plane at their ends raise; curves meeting the plane more
  than twice on the seed's side are not handled.
- **Corners** under 30°, cusps, a seed on a corner and a curve meeting the
  plane at a corner are out of scope.
- **Cost on slivers.** Tracing on sliver STL triangles spends seconds in
  trimesh's broad phase (5–8 s for a 1000-point-section wing loop).
- Corners are traced exactly, but a single collar block still folds at
  sharp ones; see [Collar blocks](/concepts/collars).

## See also

- ADRs: [0009 — Intersection loops](/adr/0009-intersection-curves),
  [0011 — Open curves on a plane](/adr/0011-open-intersection-curves),
  [0020 — Junction-curve spacing](/adr/0020-junction-curve-spacing),
  [0021 — Corners](/adr/0021-junction-curves-with-corners),
  [0022 — Grading across corners](/adr/0022-collars-at-blunt-corners),
  [0031 — Polishing cycles](/adr/0031-polishing-exact-cycles).
- [`shore.surface.intersect`](/reference/api/shore.surface.intersect) —
  [`intersection_loop`](/reference/api/shore.surface.intersect#intersection_loop),
  [`intersection_curve`](/reference/api/shore.surface.intersect#intersection_curve),
  [`TraceParams`](/reference/api/shore.surface.intersect#TraceParams).
- [Surface marching](/concepts/surface-marching), [Collar blocks](/concepts/collars),
  [Spacing laws](/concepts/spacing-laws).
- [`shore collar`](/reference/collar), [Collar tutorial](/tutorials/collar),
  [Half models](/cookbook/half-models).
