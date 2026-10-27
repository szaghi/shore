# 0021 — Junction curves with corners

- **Status:** Accepted
- **Date:** 2026-10-27

## Context

The tracer of ADR 0009 follows the intersection curve with steps along
`t = n_b x n_a`, each step polished back onto both bodies. Where the curve
crosses an edge of either body (a crease), the tangent jumps:
- the step overshoots the edge;
- projection pulls the point back onto it;
- the trace stalls there until `max_steps` runs out.

Common junctions fail this way:
- a sharp or blunt trailing edge through a plate;
- any prismatic strut;
- a body pierced across the edge of another.

ADR 0020 found this. Prototyped, all of these failed: a NACA 0012 wing
through a plate (sharp and blunt TE), and analytic lens, rectangle,
L-shape and box-edge junctions.

## Investigation (2026-10-27, prototypes in the session scratchpad; T1)

**A corner from `Geometry.project` only.** Everything below uses
projection only, so the same code works for STL, CAD and analytic bodies.

- **Failed step.** A plain step fails when it advances less than half a
  step, or leaves the tangent by 30° or more, or the tangent turns by 30°
  or more (`TraceParams.corner_deg`, the STL crease angle of ADR 0014).
- **Adaptive walk.** A failed step is redone as a walk with a sixteenth of
  the step, recursively. A step that still fails at `1e-6` × the diagonal
  (`corner_scale`) is a corner. A walk that finds none means a tight smooth
  bend, such as a leading edge that turns 40° per step, and the original
  step is kept. This keeps every smooth trace bit for bit (below).
- **The way out of a corner.** The curve's points at distance *r* from a
  point are found by alternating projections onto A, B and the sphere of
  radius *r* (`_ring`), from 26 directions. The way out is the point that
  is not on the incoming branch.
- **Locating the corner.** The tangent lines of the two branches are
  crossed at points ever closer to the corner (the radius divided by 8
  each round). The corner is exact to round-off: 2e-15 on the analytic
  bodies, and the STL TE vertex exactly. If the branches turn less than
  30° at the smallest scale, there is no corner.
- **Two corners closer than a step.** The blunt TE's base corners are
  1/10 of a step apart. The first version located corners directly
  between a failed step's ends, saw the straight base between them,
  declared a smooth bend and stalled. The adaptive walk resolves both.
- **A step that lands on the corner.** Projection snapped one step
  exactly onto the L-shape's vertex. The corner then replaces that point
  instead of duplicating it; a duplicate made the curvature NaN.

**Fixed nodes and segments.** Node 0 (or the plane ends) and every corner
are fixed nodes. Each segment between them gets cells in proportion to its
length, or to its ADR 0020 monitor integral with curvature, at least one
each by largest remainders. ADR 0020's law then applies within the
segment.

| Junction (n = 64) | Corners | Corner error | Trace |
|---|---|---|---|
| plane × biconvex lens (tips 22.8°) | 2 convex | 2e-15 | 0.04 s |
| plane × rectangle 1 × 0.1 (blunt base) | 4 × 90° | 0 | 0.05 s |
| plane × L-shape | 5 convex, 1 concave | 0 | 0.08 s |
| box × sphere on its edge (corners on A) | 2 | 0 | 0.06 s |
| STL plate × NACA 0012, sharp / blunt TE (n = 96) | 1 / 2 | 0 / 0 | 0.5 / 0.8 s |
| CAD plate × spline NACA 0012, sharp / blunt TE | 1 / 2 | 2e-11 / 0 | 0.6 / 0.8 s |

**More notes:**
- **Smooth curves keep their output.** They come out bit for bit as before,
  including tight bends: the elliptic junction of ADR 0020, its STL post,
  and a faceted sphere × cylinder.
- **Laws per segment.** On the rectangle's straight segments, the default
  law and `tanh2` hold to 1e-15. With curvature 1, the lens's cells next to
  its tips shrink to 0.13-0.35 × the uniform ones.
- **Collar at corners (measured, out of scope).** `build_collar` on these
  loops folds at every corner:
  - lens 1784 cells, rectangle 13, L-shape 109, box edge 19;
  - STL sharp wing 1818, blunt 1354;
  - curvature clustering makes it worse.

  A single block cannot wrap a trailing edge. That needs its own ADR.

## Decision

User decisions of 2026-10-27:
- **Scope:** the tracer and corner nodes, with the collar at corners only
  measured.
- **Corners handled:** sharp and blunt convex, concave, and corners on
  either body.
- **Detection:** from the trace.
- **Acceptance:** corners exact; smooth curves unchanged; spacing per
  segment; STL and CAD wings.

- `intersection_loop` and `intersection_curve` trace through corners,
  with unchanged signatures. Every corner is a node, and segments get cells
  and laws as above. `ParameterError` names the corners when *n* leaves a
  segment without a cell.
- `TraceParams` gains `corner_deg = 30`, `corner_scale = 1e-6` and
  `ring_iterations = 200`.
- `_Tracer.trace` and `trace_to_plane` return the polyline and its corner
  indices. `refine` keeps fixed vertices.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Corners from CAD edges, STL crease detection | Two code paths, and no analytic bodies. The trace finds corners on any `Geometry`. |
| Locate the corner directly between a failed step's ends | Fails when two corners are closer than a step (the blunt TE). |
| Always adaptive steps (shrink on tight bends) | Changes tight smooth bends, so saved 0.2.x projects with such collars would not replay. |
| Spread the turning angle into the curvature only (no fixed nodes) | No node on the corner: the collar's grid line misses the edge. |

## Consequences

- Junction curves of prismatic struts and of wings with sharp or blunt
  trailing edges are now available, with exact corner nodes. Spacing per
  segment clusters at corners with `curvature`.
- Collars on them still fold at the corners. The next step is a collar
  that wraps corners: a block split at the corner nodes, or a C/H-type
  layout.

## Verification

### Acceptance spec (frozen 2026-10-27, before implementation)

**Specs.**
- `tests/spec/test_loop_corners.py` runs in the main CI job.
- `tests/spec/test_loop_corners_cad.py` runs in the CAD job.
- `tests/spec/data/loop_bends_default.npz` holds the tight-bend curves.
  It was written by the code before this ADR (HEAD, loaded from
  `git show`) and is hashed in `SHA256SUMS`.

**Verdict.**
- Main: 22 passed, 0 skipped, 0 xfailed. At freeze: 10 passed (the golden
  K3 cases) and 12 xfailed.
- CAD: 2 passed. At freeze: 2 xfailed.

**Stall budget.** The corner cases pass `TraceParams(max_steps=2000)`.
The traces need fewer than 300 steps, and the old tracer would otherwise
spin through 200 000 steps (minutes) before failing.

**Harness validation.**
- With the prototype behind the interface: 22 + 2 passed.
- At freeze, K4's "too few nodes" case first XPASSed, because the old
  tracer also raises `ParameterError` (the trace never closes). It now
  requires the message to name the corners.

| Case | Required |
|---|---|
| K1 (× 4) | loops through corners (lens, rectangle, L-shape, box edge): corners on distinct nodes ≤ 1e-12 × diagonal; nodes on both bodies |
| K2 (× 2) | open curves (rectangle at x = 0, lens at x = 0.25): corners on nodes, ends on the plane |
| K3 (× 10) | ADR 0020's golden curves and the tight-bend golden curves, within 1e-12 |
| K4 (× 4) | per segment: uniform and `tanh2` within 1e-9, cell counts within 1 of the share; curvature: corner cells ≤ 0.5 × uniform, ratio ≤ 1.21; too few nodes raises naming the corners |
| K5 (× 2) | STL NACA 0012 wing through a plate, sharp and blunt: TE corners on nodes ≤ 1e-9 × diagonal, ≤ 5 s |
| K6 (× 2, CAD job) | CAD spline wing, sharp and blunt: TE points on nodes ≤ 1e-7 × diagonal, ≤ 10 s |

**Out of scope:**
- the collar at corners;
- a seed on a corner;
- a curve meeting the plane at a corner;
- corners under 30°;
- cusps.

### Review (2026-10-27, user decisions)

1. Corner node indices are not exposed yet: the collar-at-corners ADR adds
   the API when its consumer exists.
2. The threshold stays a fixed 30° in `TraceParams.corner_deg` (CAD and
   analytic bodies have no crease angle).
3. Next step: a collar that wraps corners.

## Acceptance (2026-10-27): junction curves with corners spec met

**Kernel.** `_Tracer.trace` and `trace_to_plane` check each plain step
(`_good`) and redo a failed one with `_fine`, a recursive walk of 1/16
steps.
- **Corner:** a failure at `corner_scale` is a corner. `_turn` leaves it
  through `_ring`, and `_locate` places it.
- **Tight bend:** a walk without a corner keeps the plain step, so smooth
  traces are bit-identical.
- **Dedup:** `_add_corner` replaces a traced point that sits on the corner.

`refine` renumbers fixed vertices. `_resample` and `_allocate` give the
segments their cells and laws.

**Half models.** The first full-suite run failed ADR 0018's L2, a
half-collar on STL half bodies: one step short of the plane the trace
stalls on the cut edge of the open half body. It was taken for a corner
with no way out. The fix:
- a fine walk stops when the plane is within its step;
- a stalled main step with the plane within one step ends on the plane,
  as the old tracer eventually did.

A step that turns more than 30° but still advances is not a stall: the
half ellipse's tip is next to the plane, and its golden curve stays
bit-identical.

**Verdicts (T2; `scripts/check_spec.py` reports all files ok):**
- `tests/spec/test_loop_corners.py`: **22 passed, 0 skipped, 0 xfailed**.
- `tests/spec/test_loop_corners_cad.py`: **2 passed** (run locally with
  `shore[cad]`; the CAD job runs it).
- All 7 `@PENDING` lines are removed.

**Unit tests** (`tests/test_intersection.py`):
- line crossing, including the parallel fallback;
- allocation with a minimum of one cell;
- corner dedup;
- `refine` keeps fixed vertices;
- `_ring` finds both branches of a circle.

**Docs.** `docs/reference/collar.md` states that corners are traced and
that collars still fold there.
