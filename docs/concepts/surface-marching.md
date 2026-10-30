# Surface marching

Surface marching grows a structured grid over a body from a curve on it,
one row at a time, the way hyperbolic volume marching grows layers off a
wall. It is the kernel of SHORE's collar grids: from the junction curve,
a surface grid is marched away from the other body on each component.
The scheme is the SURGRD type of Chan & Steger (1992) and Chan & Buning
(1995): a hyperbolic system for orthogonality and cell area, linearised
and solved once per layer, with every new layer projected onto the body.

`shore.surface.march.march_surface(geometry, curve, steps, *,
periodic=False, entity=BODY, plane=None)` returns an `(n, m + 1, 3)`
array: row 0 is the curve, row $j + 1$ lies on the
$\mathbf n \times \mathbf t$ side of row $j$, every node on the surface.
Only `Geometry.project` (with normals) is used, so it works on STL and
CAD bodies alike.

<Figure name="surface-march" caption="A surface grid marched over a body from a curve on it." />

## The linearised system

Let $\mathbf r$ be the current row, $\mathbf t = \mathbf r_\xi$ its
tangent (a central difference along the row), $\mathbf n$ the surface
normal and $\Delta\mathbf r$ the increment to the next row. About the
local solution $\Delta\mathbf r_0 = d\,(\mathbf n \times \mathbf t) /
|\mathbf n \times \mathbf t|$, the three conditions are

$$
\begin{aligned}
\mathbf t \cdot \Delta\mathbf r + \Delta\mathbf r_0 \cdot D(\Delta\mathbf r) &= 0
&&\text{(lines orthogonal to the row)}\\
(\mathbf n \times \mathbf t) \cdot \Delta\mathbf r
+ (\Delta\mathbf r_0 \times \mathbf n) \cdot D(\Delta\mathbf r) &= \Delta A
&&\text{(prescribed cell area)}\\
\mathbf n \cdot \Delta\mathbf r &= 0
&&\text{(stay in the tangent plane)}
\end{aligned}
$$

with $D$ the central difference along the row. Implicit dissipation
$-\varepsilon_i \nabla^2(\Delta\mathbf r)$ and explicit tangential
dissipation $\varepsilon_e\,T(\nabla^2\mathbf r)$ are added, and the
block-tridiagonal system (cyclic for a periodic curve) is solved with one
sparse LU factorisation per row. The new row is projected onto the
surface; that projection also returns the next row's normals, so the
geometry is queried once per node and row.

## Stabilisation

Every term is ramped by $\rho = \min(1, \Delta s / h)$, the step over the
local node spacing, so thin rows march along the local normal and only
thick ones feel the smoothing:

- **Area smoothing.** The prescribed area $\Delta A = \Delta s\,h$ is
  blended towards its neighbours (`area_smoothing` 0.5 × $\rho$,
  2 passes). Fronts that collide at a concave corner share their area
  instead of crossing. The marching distance is kept within
  $[\Delta s / 2, 2\Delta s]$.
- **Contraction sensor.** Where the area-smoothed pure-normal trial row
  contracts, the dissipation is raised:
  $\varepsilon_i = \varepsilon_{i,0}\,\rho\,(1 + 16\,c)$ and
  $\varepsilon_e = \varepsilon_{e,0}\,\rho^2\,(1 + 16\,c)$, with $c \in [0, 1]$
  the contraction.

The same stabilisation lengthens the steps of converging fronts: marching
north from a sphere's equator overshoots the schedule by 3.5 % (accepted
bound 5 %). Pure normal marching, the alternative, folds at concave
corners (2 cells on the spec's concave L, where this scheme gives a
minimum angle of 23.9°).

Measured on exact surfaces (25 rows, steps $0.01 \cdot 1.08^j$): straight
lines on a plane and rings on a cylinder are exact to $10^{-15}$; a
generator marched round a cylinder is within 0.056 % of the geodesic
distance; a helix on a cylinder keeps orthogonality within 0.04°; a convex
L keeps 45°.

## Ends of open curves

- **Free** (`plane=None`): the end nodes march along their own
  $\mathbf n \times \mathbf t$.
- **On a symmetry plane** (`plane=(point, normal)`): each end is solved as
  an interior node whose missing neighbour is a mirror ghost, the image of
  its inner neighbour across the plane, both in the differences (tangent,
  Laplacian, contraction sensor) and in the linear system (the ghost's
  increment is the reflected increment). After the projection, alternating
  projections put the end on both the body and the plane (to $10^{-14}$ of
  the diagonal). On mirror-symmetric bodies the result equals the mirrored
  periodic march of the whole curve to $3 \cdot 10^{-14}$ of the first
  step, and the end tangent is within 0.12° of the plane normal. Snapping
  freely marched ends onto the plane afterwards gives up to 5° instead.
  A curve whose ends are off the plane, or `plane` with `periodic=True`,
  is refused.

STL vertex normals are interpolated barycentrically, so on a coarsely
triangulated flat face they tilt (up to 14° on a coarse half box), and a
half mesh's normals are not the mirror of the whole mesh's. Marches on
such STL half models keep their ends on the plane and stay valid, but do
not match the mirrored whole-body march.

## Creases and corners

Junction loops with corners (blunt and sharp trailing edges, struts) put
nodes of the front on edges of the surface. Three mechanisms, each
triggered only where it is needed, so marches that do not trigger them
are bit-identical:

| Situation | Trigger | Treatment |
|---|---|---|
| **Crease**: a front corner on a surface edge; the node's normal leans along $\mathbf t$ | $\lvert\mathbf n \cdot \hat{\mathbf t}\rvert > 0.1$ (`crease_tolerance`) | the marching frame uses $\mathbf n - (\mathbf n \cdot \hat{\mathbf t})\hat{\mathbf t}$, normalised |
| **Sharp crease**: the front folds back on itself (a sharp TE) | at a crease, $\lvert\mathbf t\rvert$ below 0.5 × the node's arc spacing (`sharp_crease`) | $\mathbf t$ lengthened to the arc spacing |
| **Clustered corner**: the step dwarfs the front spacing | step > 24 × the smallest spacing of the nodes not closer than on the starting curve (`max_step_ratio`) | the row is marched in $\lceil \text{ratio} / 2 \rceil$ sub-steps, each projected (`substep_ratio` 2) |
| **Fold**: converging fronts overtake each other | a cell reversed against the majority, or a front segment pointing backwards | nodes within 2 of the fold smoothed along the front (1-2-1) and re-projected, up to 200 times (`untangle_*`) |

Why each is needed:

- **Creases.** The scheme assumes $\mathbf n \perp \mathbf t$. At a
  crease the area row divides by $|\mathbf n \times \mathbf t| <
  |\mathbf t|$, the marching distances vary from node to node, and the
  orthogonality row turns that into tangential drift: up to 0.13 chord
  marching up a prismatic blunt wing, with 10–73 flipped cells. After the
  fix: 0 flips, drift 0.011–0.060.
- **Sharp creases.** Across a 16° fold the central difference is about
  $\sin 8°$ of the spacing, and the prescribed area gives an outsized step
  that drags the neighbours. The periodic march of a sharp-TE wing with
  curvature clustering flipped 66 cells; after the fix none.
- **Sub-steps.** With curvature clustering the corner spacing is
  $10^{-4}$–$3 \cdot 10^{-5}$ while the first step is $5 \cdot 10^{-3}$.
  Counting only nodes not closer than on the starting curve keeps
  converging fronts (whose spacing shrinks legitimately) from triggering
  it.
- **Untangling.** At a concave corner, sub-stepping does not converge;
  the corner cells shrink until nodes overtake each other along the
  front.

After these, collar walls (single and split) no longer fold on any of the
tested junctions: L-shape plates, blunt bases and sharp creases at loop
sizes 96–384 with curvature 0–3. Checked with a four-corner cell sine
(each corner of each surface cell): at least 0.06 everywhere.

## Parameters

`MarchParams` holds the defaults; they are not per-case tuning knobs.

| Field | Default | Meaning |
|---|---|---|
| `eps_implicit` | 1.0 | implicit dissipation at $\rho = 1$ |
| `eps_explicit` | 0.5 | explicit tangential dissipation at $\rho = 1$ |
| `area_smoothing` | 0.5 | area blend towards the neighbours at $\rho = 1$ |
| `area_iters` | 2 | area smoothing passes |
| `convergence_gain` | 16 | extra dissipation per unit of front contraction |
| `distance_clip` | 2 | distance kept in $[\Delta s / 2, 2\Delta s]$ |
| `crease_tolerance` | 0.1 | $\lvert\mathbf n \cdot \hat{\mathbf t}\rvert$ above it: crease |
| `sharp_crease` | 0.5 | $\lvert\mathbf t\rvert$ below this × arc spacing: lengthened |
| `max_step_ratio` | 24 | sub-step trigger |
| `substep_ratio` | 2 | sub-steps of at most this × the spacing |
| `untangle_iterations` | 200 | smoothing passes at a fold |
| `untangle_reach` | 2 | nodes smoothed each side of a folded cell |
| `end_tolerance` | 1e-14 | ends on surface and plane (× diagonal) |
| `end_iterations` | 100 | alternating projections at an end |

## Example

```python
import numpy as np

from shore.geometry.stl import StlGeometry
from shore.shapes import cylinder, sphere
from shore.surface.intersect import intersection_loop
from shore.surface.march import march_surface

ball = StlGeometry(sphere(radius=1.0, subdivisions=4))
post = StlGeometry(cylinder(radius=0.4, height=3.0, sections=64))
loop = intersection_loop(ball, post, np.array([0.4, 0.0, 0.9]), 64)

steps = 0.01 * 1.12 ** np.arange(15)
on_ball = march_surface(ball, loop, steps, periodic=True)  # away from the post
on_post = march_surface(post, loop[::-1], steps, periodic=True)  # away from the ball
print(on_ball.shape, on_post.shape)  # (64, 16, 3) twice
```

## Limits

- **Converging fronts overshoot** the step schedule by a few per cent
  (3.5 % on the sphere): the price of the concave-corner robustness.
- **Concave notches.** Marching deeper into a concave notch than the
  notch is wide concentrates the front's loss of length on the cells at
  the bisector; they become long and thin (aspect about 2700 on the
  L-shape plate). The surface grid stays valid, but a collar volume
  built on it folds there. Spreading bunched fronts was tried and not
  adopted (ADR 0029, 0030).
- **Run time.** Sub-stepping makes strongly clustered loops (384 nodes
  with curvature) take minutes: up to 150 s for a blunt-wing collar,
  uncapped by decision.
- **Not fixed:** the sharp STL wing at 384 loop nodes with curvature
  (spacing $3 \cdot 10^{-5}$ on sliver triangles) still flips a few cells.
- **Two planes.** `march_surface` holds both ends on one plane; the split
  collar uses a private variant with one plane per end.

## See also

- ADRs: [0008 — Surface marching](/adr/0008-surface-marching),
  [0012 — Ends on a plane](/adr/0012-march-ends-on-a-plane),
  [0023 — Creases](/adr/0023-surface-marching-at-creases),
  [0024 — Sharp creases](/adr/0024-split-collar-at-sharp-corners),
  [0027 — Clustered corners](/adr/0027-surface-marching-at-clustered-corners),
  [0029 — Converging fronts](/adr/0029-spreading-converging-fronts),
  [0030 — Concave notches](/adr/0030-concave-notches).
- [`shore.surface.march`](/reference/api/shore.surface.march) —
  [`march_surface`](/reference/api/shore.surface.march#march_surface),
  [`MarchParams`](/reference/api/shore.surface.march#MarchParams).
- [Intersection curves](/concepts/intersection-curves),
  [Collar blocks](/concepts/collars),
  [Marching schemes](/concepts/hyperbolic-marching) (the volume counterpart).
