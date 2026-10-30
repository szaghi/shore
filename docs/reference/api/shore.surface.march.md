---
title: shore.surface.march
outline: [2, 3]
---

# `shore.surface.march`

Hyperbolic surface marching: the collar-grid kernel (ADR 0008).

`march_surface` grows a structured grid over a surface from a curve, one
layer per step (Chan & Steger 1992, SURGRD type). Layer `j + 1 = r + dr`
is marched from layer `r` to the `n x t` side (`n` the surface normal,
`t = r_xi`), with `dr` from the linearised system

    t . dr          + dr0 . D(dr)        = 0          (k-lines orthogonal to the layer)
    (n x t) . dr    + (dr0 x n) . D(dr)  = dA         (prescribed cell area)
    n . dr                               = 0          (stay in the tangent plane)

about the local solution `dr0 = dist * (n x t) / |n x t|` (`D` the
central difference along the layer), plus implicit dissipation
`-eps_i Lap(dr)` and explicit tangential dissipation `eps_e T(Lap r)`.
One sparse solve per layer (SuperLU, deterministic); the new layer is then
projected onto the surface, which also gives the next layer's normals, so
the geometry is queried once per node and layer.

Stabilisation, ramped by `rho = min(1, step / h)` (`h` the node spacing)
so thin layers march along the local normal:

- area smoothing: `dA` is blended toward its neighbours (fronts that
  collide at a concave corner share their area instead of crossing);
- a convergence sensor raises the dissipation where the pure-normal trial
  layer contracts.

The same stabilisation lengthens the steps of converging fronts: marching
north from a sphere's equator overshoots the schedule by ~3.5 % (ADR 0008).
Open curves march their end nodes along their own `n x t`, or, with
`plane`, hold them on a symmetry plane (ADR 0012): each end is solved as
an interior node whose missing neighbour is the mirror image of its inner
neighbour across the plane, then put on the surface and the plane by
alternating projections. On mirror-symmetric bodies this is the symmetric
half of the periodic march of the mirrored closed curve.

At a crease the node's normal is made orthogonal to the front tangent
(ADR 0023). At a sharp crease, where the front folds back on itself (the
sharp trailing edge of a wing), the central difference is much shorter than
the node spacing; its tangent is then lengthened to the mean arc spacing of
the node, so the prescribed area does not turn into an outsized step
(ADR 0024). The split collar (ADR 0024) also holds the two ends of an open
curve on two different planes (`_Ends`).

At corners where the loop is clustered (ADR 0027) a layer whose step exceeds
`max_step_ratio` x the front spacing (counting only nodes not closer than
on the starting curve, so converging fronts do not trigger it) is marched in
sub-steps; and a layer that folds against the previous one (a cell
reversed, or a front segment pointing backwards: converging fronts at a
concave corner) is smoothed along the front round the fold and re-projected
until it does not. Marches where neither happens are unchanged.

## `MarchParams` <Badge type="tip" text="dataclass" /> {#MarchParams}

```python
MarchParams(
    eps_implicit: float = 1.0,
    eps_explicit: float = 0.5,
    area_smoothing: float = 0.5,
    area_iters: int = 2,
    convergence_gain: float = 16.0,
    distance_clip: float = 2.0,
    crease_tolerance: float = 0.1,
    sharp_crease: float = 0.5,
    max_step_ratio: float = 24.0,
    substep_ratio: float = 2.0,
    untangle_iterations: int = 200,
    untangle_reach: int = 2,
    end_tolerance: float = 1e-14,
    end_iterations: int = 100,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/march.py#L66">source</a>

Defaults of the surface marcher (not per-case tuning knobs).

## `march_surface` <Badge type="info" text="function" /> {#march_surface}

```python
march_surface(
    geometry: Geometry,
    curve: np.ndarray,
    steps: np.ndarray,
    *,
    periodic: bool = False,
    entity: EntityRef = BODY,
    params: MarchParams = DEFAULT,
    plane: tuple[np.ndarray, np.ndarray] | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/march.py#L88">source</a>

March a surface grid from *curve* over *entity* of *geometry*.

**Parameters**

- **`geometry`** — Any backend; only `project` (with normals) is used.
- **`curve`** (np.ndarray, shape (n, 3)) — The starting curve, n &gt;= 3, on the surface. With *periodic* a closed loop whose last node is not repeated.
- **`steps`** (np.ndarray, shape (m,)) — Positive marching distance of each layer.
- **`periodic`** — Whether *curve* is closed.
- **`entity`** — The entity marched over (default: the whole body).
- **`params`** — Marcher settings.
- **`plane`** — `(point, normal)`: hold both ends of an open *curve* on this plane with a symmetry condition (ADR 0012). The curve's ends must lie on it. `None`: the ends march freely.

**Returns**

- `np.ndarray, shape (n, m + 1, 3)` — `[:, 0]` is *curve*; layer j + 1 lies on the `n x t` side of layer j, every node on the surface.

**Raises**

- `ParameterError` — *curve* not (n, 3) with n &gt;= 3; *steps* not a non-empty 1-D array of positive finite values; *curve* off the surface by more than 1e-6 of the model diagonal; *plane* not a point and a non-zero normal, given with *periodic*, or with the curve's ends off it by more than that tolerance.
