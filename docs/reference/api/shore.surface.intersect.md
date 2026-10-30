---
title: shore.surface.intersect
outline: [2, 3]
---

# `shore.surface.intersect`

Intersection loop of two bodies, the start curve of a collar grid (ADR 0009).

`intersection_loop` finds the closed curve where two bodies meet, near a
seed point, with `Geometry.project` only (so STL, CAD and analytic bodies
work alike):

1. alternating projections onto A and B take the seed onto the curve;
2. the curve is traced with predictor steps along `t = n_b x n_a` (the
   bodies' outward normals), each corrected back onto both bodies by
   alternating projections, until the trace returns to its start;
3. the traced polyline is resampled to *n* nodes equally spaced in arc
   length, and every node is polished back onto both bodies.

The `n_b x n_a` order makes `march_surface(A, loop)` march on A away
from B, and `march_surface(B, loop[::-1])` on B away from A: the two halves
of a collar (ADR 0008). Alternating projections converge fast where the
bodies cross at a clear angle; grazing intersections are out of scope.

`intersection_curve` is the open variant for half models (ADR 0011): the
curve on the seed's side of a plane, traced both ways from the seed until
it reaches the plane (within the on-curve tolerance: half-model bodies stop
at the plane, so the trace lands on it without crossing), each end put
exactly on the plane by Newton steps along the curve tangent.

Both take a node distribution along the curve (ADR 0020): a `Spacing1D`
law by arc length and/or `curvature` clustering, which equidistributes
`1 + alpha * kappa * L / (2 pi)` on the traced polyline (refined first
where it turns, so that the curvature and the arc lengths are resolved)
before the law is applied, with adjacent cells limited to
`TraceParams.max_cell_ratio`. The defaults keep the plain arc-length
resampling above, bit for bit.

## `TraceParams` <Badge type="tip" text="dataclass" /> {#TraceParams}

```python
TraceParams(
    step: float = 0.004,
    tolerance: float = 1e-13,
    trace_tolerance: float = 1e-10,
    on_curve: float = 1e-09,
    polish_iterations: int = 200,
    polish_period: int = 12,
    max_steps: int = 200000,
    refine_turn_deg: float = 2.0,
    refine_min: float = 1e-05,
    curvature_window: float = 0.5,
    max_cell_ratio: float = 1.2,
    corner_deg: float = 30.0,
    corner_scale: float = 1e-06,
    ring_iterations: int = 200,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/intersect.py#L50">source</a>

Tracer settings (not per-case tuning knobs); lengths relative to the model diagonal.

## `intersection_loop` <Badge type="info" text="function" /> {#intersection_loop}

```python
intersection_loop(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    n: int,
    *,
    spacing: Spacing1D | None = None,
    curvature: float = 0.0,
    entity_a: EntityRef = BODY,
    entity_b: EntityRef = BODY,
    params: TraceParams = DEFAULT,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/intersect.py#L74">source</a>

The closed intersection loop of *body_a* and *body_b* nearest *seed*.

**Parameters**

- **`body_a, body_b`** — The two bodies; only `project` (with normals) is used.
- **`seed`** (np.ndarray, shape (3,)) — A point near the junction; selects the loop when there are several.
- **`n`** — Number of nodes, &gt;= 3.
- **`spacing`** — Node distribution by arc length from node 0 round to node 0 again: `uniform` or `tanh2` (default: uniform).
- **`curvature`** — Clustering weight *alpha* &gt;= 0 where the loop bends (ADR 0020).
- **`entity_a, entity_b`** — The entities intersected (default: the whole bodies).
- **`params`** — Tracer settings.

**Returns**

- `np.ndarray, shape (n, 3)` — Nodes distributed by *spacing* and *curvature* (by default equally spaced in arc length), on both bodies, ordered along `n_b x n_a`; `[0]` is where the seed lands; the last node is not repeated.

**Raises**

- `ParameterError` — n &lt; 3; *seed* not a 3-vector; a law other than `uniform` or `tanh2`, or a bad *curvature*; no intersection near *seed*; a trace that does not close.

## `intersection_curve` <Badge type="info" text="function" /> {#intersection_curve}

```python
intersection_curve(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    n: int,
    *,
    plane: tuple[np.ndarray, np.ndarray],
    spacing: Spacing1D | None = None,
    curvature: float = 0.0,
    entity_a: EntityRef = BODY,
    entity_b: EntityRef = BODY,
    params: TraceParams = DEFAULT,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/intersect.py#L173">source</a>

The open intersection curve of *body_a* and *body_b* ending on *plane*.

**Parameters**

- **`body_a, body_b`** — The two bodies; only `project` (with normals) is used. They may end at the plane (a half model).
- **`seed`** (np.ndarray, shape (3,)) — A point near the junction, off the plane; selects the curve and the side of the plane that is kept.
- **`n`** — Number of nodes, &gt;= 3.
- **`plane`** — `(point, normal)`, two 3-vectors; the normal's sign does not matter.
- **`spacing`** — Node distribution by arc length from the first end to the last: `uniform`, `tanh` or `tanh2` (default: uniform).
- **`curvature`** — Clustering weight *alpha* &gt;= 0 where the curve, or its mirror image through the plane, bends (ADR 0020).
- **`entity_a, entity_b`** — The entities intersected (default: the whole bodies).
- **`params`** — Tracer settings.

**Returns**

- `np.ndarray, shape (n, 3)` — Nodes distributed by *spacing* and *curvature* (by default equally spaced in arc length), on both bodies and on the seed's side, ordered along `n_b x n_a`; `[0]` and `[-1]` on the plane.

**Raises**

- `ParameterError` — n &lt; 3; *seed* not a 3-vector; *plane* not a point and a non-zero normal; a law other than `uniform`, `tanh` or `tanh2`, or a bad *curvature*; *seed* on the plane; no intersection near *seed*; a curve that closes, or meets the plane tangentially, or does not reach it.
