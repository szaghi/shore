---
title: Grid quality
---

# Grid quality

SHORE measures grids with a small, fixed set of metrics in
`shore.quality`. They are **frozen**: the acceptance specs under
`tests/spec/` state their thresholds in these metrics, so the definitions
cannot change without a deliberate re-freeze
([ADR 0005](/adr/0005-marching-schemes)). Every function takes one block's
node array `(ni, nj, nk, 3)` and assumes the marcher's orientation: a valid
cell is right-handed, $(\mathbf{e}_i \times \mathbf{e}_j) \cdot \mathbf{e}_k
> 0$, with $k$ pointing away from the wall.

## The metrics

### Scaled Jacobian

`scaled_jacobian(nodes)`, shape `(ni-1, nj-1, nk-1)`: per cell, the minimum
over its 8 corners of

$$
J_s = \frac{\det[\mathbf{e}_i, \mathbf{e}_j, \mathbf{e}_k]}
           {|\mathbf{e}_i|\,|\mathbf{e}_j|\,|\mathbf{e}_k|}
$$

with $\mathbf{e}_i, \mathbf{e}_j, \mathbf{e}_k$ the three cell edges leaving
that corner, each oriented along $+i$, $+j$, $+k$. This is the Verdict /
CUBIT definition. $J_s = 1$ for a rectangular cell, the cosine of the shear
angle for a sheared one, and $\le 0$ when a corner is folded. A corner with
a zero-length edge scores $-1$. Taking all 8 corners matters: a cell can
have a positive Jacobian at one corner and be folded at another.

<Figure name="quality-jacobian" caption="A grid coloured by the scaled Jacobian of each cell." />

### Wall orthogonality

`k_orthogonality_deg(nodes, k=0)`, shape `(ni, nj)`: the angle between the
k-segment from layer $k$ to $k+1$ and the layer's normal
$\mathbf{t}_i \times \mathbf{t}_j$, with tangents by central differences
(one-sided on the layer's edges). 0° is orthogonal; above 90° the segment
points back through the layer.

The one-sided tangents bias it on block-edge nodes, where the true surface
normal spans the seam: unsmoothed normal extrusion, exactly orthogonal in
the interior, scored 9.4° mean / 13.9° max on a coarse cubed sphere's edge
nodes. `k_angle_deg(nodes, directions, k=0)` measures the same angle
against given directions instead; with the body's surface normals at the
wall it is unbiased, and it is what the marching specs use (taken
orientation-free, $\min(a, 180° - a)$).

### Wall-normal spacing and growth

- `k_spacing(nodes)`, shape `(ni, nj, nk-1)`:
  $|\mathbf{x}_{k+1} - \mathbf{x}_k|$ along every k-line. Its first layer is
  the first cell height, to compare with the requested `ds`.
- `k_growth(nodes)`, shape `(ni, nj, nk-2)`: the ratio of consecutive
  k-segments. A geometric law with `growth = 1.15` reads 1.15.

### Aspect ratio

`aspect_ratio(nodes)`, shape `(ni-1, nj-1, nk-1)`: longest over shortest of
each cell's 12 edges; infinite for a cell with a zero-length edge. Boundary
layers are stretched by design (wall cells of $10^{-3}$ under surface cells
of $10^{-1}$ read 100), so it is a diagnostic, not a pass criterion.

### Summary

`summary(blocks)` reduces a block list (objects with `.nodes`, or arrays)
to extremes:

| Key | Value |
|---|---|
| `scaled_jacobian_min` | minimum $J_s$ over all cells |
| `folded_cells` | cells with $J_s \le 0$ |
| `wall_orthogonality_max_deg`, `wall_orthogonality_mean_deg` | `k_orthogonality_deg` at $k = 0$ |
| `aspect_ratio_max` | maximum aspect ratio |
| `k_growth_max` | maximum k-growth (NaN when no block has `nk > 2`) |

<Figure name="quality-histogram" ext="svg" caption="Distribution of the scaled Jacobian and of the wall orthogonality for the sphere tutorial grid." />

## An example

The same cubed-sphere grid, marched with the default algebraic scheme and
with the hyperbolic scheme:

```python
import warnings

from shore import quality
from shore.exceptions import MarchingError
from shore.project import Project
from shore.shapes import sphere

warnings.simplefilter("ignore", UserWarning)  # the cap/equator C1 fallback notice
sphere(radius=1.0, subdivisions=4).export("ball.stl")

for scheme in ("algebraic", "hyperbolic"):
    p = Project()
    p.load_geometry("ball.stl")
    p.build_topology("cubed_sphere", nk=21, ds=0.005, growth=1.15, ni=33, nj=64)
    try:
        p.march(march_scheme=scheme)
    except MarchingError as err:
        print(f"{scheme:10s} refused: {str(err).split(' Cubed')[0]}")
        continue
    s = quality.summary(p.mesh.blocks)
    print(f"{scheme:10s} min SJ {s['scaled_jacobian_min']:.3f}  folded {s['folded_cells']}"
          f"  wall angle mean/max {s['wall_orthogonality_mean_deg']:.2f}"
          f"/{s['wall_orthogonality_max_deg']:.2f}")
```

```text
algebraic  refused: Folded cell (scaled Jacobian -0.089) in cap_north at layer k=19, cell (15, 15), 4 folded in the layer.
hyperbolic min SJ 0.488  folded 0  wall angle mean/max 0.40/4.22
```

With the default 9° cap the algebraic march folds the corner cells of the
polar caps (the 3-block junctions) in its outer layers. Since
[ADR 0052](/adr/0052-no-silent-cap-folds) the cubed-sphere march judges
every new cell by its scaled Jacobian, corner cells included, and refuses
the grid; before, its one-corner guard left those cells out and returned
16 folded cells here. A larger cap (`theta_cap` 30° or more) or the
hyperbolic scheme marches this body cleanly. The other topologies still
use the one-corner guard: run `summary` on every grid you hand to a
solver.

## What the acceptance specs require

Thresholds are set per case, before implementation, from measured
references. They are not universal targets, but they show what the
algorithms are held to (measured values in brackets where the ADR gives
them):

| Grid (spec) | Required |
|---|---|
| Hyperbolic march, cubed sphere 60 × 60 × 10 (`test_hyperbolic_marching.py`, H2) | min $J_s \ge 0.70$; wall angle to the body normal mean ≤ 2.0°, max ≤ 4.6° (algebraic reference 2.05° / 4.6°, $J_s$ 0.754) |
| Hyperbolic march, coarse sphere (H1) | wall-angle mean ≤ 1.15 × pure normal extrusion (3.49°), max ≤ 10°; first cell within 5 % of `ds` |
| Hyperbolic march, concave prism (H3) | 15 layers, 0 folded cells |
| Collar block (`test_collar.py`) | 0 folded cells, min $J_s \ge 0.05$ (prototype 0.081) |
| Collar volume at corners (`test_collar_volume.py`) | min $J_s$ from ≥ 0.01 (blunt TE, extreme clustering; 0.023) to ≥ 0.2 (split lens collar; 0.437) |
| Blocking premesh (`test_blocking.py`) | a split box: $> 0.999$ (Cartesian); 26 blocks around a sphere: $> 0.5$ (0.71); the same journal on a stretched sphere: $> 0.3$ (0.46) |
| O-grid insertion (`test_ogrid.py`) | $> 0.5$, or $> 0.3$ on the harder cases; seams closed (gap 0) |
| Orthogonalize and wall layer, sphere and wing blockings (`test_orthogonalize.py`, `test_wall_layer.py`) | min $J_s$ no lower than before the edit; $> 0.2$ to $> 0.3$ by case (a 40° wall angle brought under 25°) |
| Phase-7 wing–body, 90 blocks (`test_wing_body.py`) | min $J_s > 0.1$ (0.17, edited model 0.135); wall angle median < 12°, p95 < 35° (9.0° / 25.3°); first cells exactly 2e-3 |

The wing–body values are honest for an interactive blocking around a
wing–body junction: sharp edges (trailing edge, tip) keep wall angles near
45°, and the minimum $J_s$ sits at the concave junction
([ADR 0051](/adr/0051-wing-body-acceptance)).

## `shore check`

`shore check FILE.geo [--warn-threshold 1e-6]` validates one `.geo` block
layer by layer, and exits with code 1 if any layer has a non-positive
Jacobian. It uses the marcher's guard, not the frozen metrics: the
dimensional triple product $(\mathbf{e}_i \times \mathbf{e}_j) \cdot
\mathbf{e}_k$ at each cell's $(i, j, k)$ corner only, so
`--warn-threshold` is an absolute cell volume, and a fold at another corner
can pass. `Mesh.check()` does the same per block in memory. For the
corner-wise metrics, call `shore.quality` on the blocks
([Check grid quality](/cookbook/quality)).

## See also

- [ADR 0005](/adr/0005-marching-schemes) (the frozen metrics and the re-freeze of
  wall orthogonality), [ADR 0051](/adr/0051-wing-body-acceptance)
- API: [`shore.quality`](/reference/api/shore.quality),
  [`summary`](/reference/api/shore.quality#summary),
  [`shore.volume.jacobian`](/reference/api/shore.volume.jacobian)
- [Check grid quality](/cookbook/quality), [`shore check`](/reference/check),
  [A sphere, end to end](/tutorials/sphere)
- [Marching schemes](/concepts/hyperbolic-marching), [Block model](/concepts/block-model)
