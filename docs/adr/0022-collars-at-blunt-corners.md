# 0022 — Graded spacing across corners; collars at blunt corners

- **Status:** Accepted
- **Date:** 2026-10-27

## Context

ADR 0021 traces junction curves through corners and makes every corner a
node, but the collar block (ADR 0010) built on such loops folds at the
corners. This step asks what the single collar block can do at corners.
The user decided to cover corners of up to about 90° first, in the single
block. Sharp trailing edges need a different topology. The corner fan is
ruled out, because Xnavis does not accept collapsed wall cells. That leaves
a split, multi-block collar.

## Investigation (2026-10-27, prototypes in the session scratchpad; T1)

**Where collars fold at corners.** Measured on ADR 0021's junctions, with
`n_loop` 96, surface steps `0.005 * 1.12**j` and the collar thickness in
units of the junction spacing *h* (0.005):

| Junction | 1.9 × h | 4.7 × h | 11.8 × h |
|---|---|---|---|
| rectangle prism, convex 90° corners | valid (min sJ 0.19) | valid (0.05) | 13 folded, outer layer only |
| box edge × sphere r 0.5 (corners on A) | valid (0.21) | valid (0.15) | 6 folded |
| L-shape, one concave corner | 4 folded | 8 folded | 109 folded |
| biconvex lens, sharp tips | — | — | 1784 folded, walls included |

- **Convex 90° corners** work in the single block up to a thickness limit,
  like ADR 0010's thick-collar limit.
- **Concave corners** fold even in thin collars.
- **More hyperbolic layers** make it much worse, from 400 to 2500 folded
  cells.
- **Smoothing across stations** near the corners, in j, moves whole i-lines
  against fixed walls and folds thousands of cells. That approach is
  rejected.
- An earlier "box edge" case used a sphere of radius 1. It touched the
  box's bottom edge tangentially and produced one zero-Jacobian cell; that
  came from the test geometry, not from the corner.

**Grading across corners.** With `curvature > 0`, ADR 0021 distributes
each segment between corners on its own. The ratio limit then never sees
the corners. On the blunt NACA 0012 STL wing, whose base is 0.0025 wide,
adjacent loop cells differed by 8.7× without curvature and 6.6× with it.

Two repairs failed:
- **Snapping corners to the nearest node of a whole-loop distribution:**
  ratio 28.
- **Iterating the limiter with each segment rescaled to its length:** it
  cannot converge when a segment has too few cells to grade.

Both failures led to a flaw in **ADR 0020's limiter**. ADR 0020 limits the
*cells*, then rescales them to the curve length. When the limiter binds
hard, the rescale moves the clusters. On the blunt wing, the cluster meant
for the TE (s = 1.537) ended up at s = 1.898, and the TE got the largest
cells.

**The fix: a graded spacing function.**
- The spacing g(s) = C / m(s) is capped by cones of slope ln 1.2, so
  adjacent cells differ by at most 1.2.
- C is found by bisection so that the curve holds exactly the number of
  cells it needs.
- Nodes equidistribute 1/g, integrated on a sampling finer than the cells.
  Integrating only at the polyline's vertices made cells jump wherever the
  vertices were a trace step apart (ratio 2.3).
- Segments between corners get cells in proportion to the integral of 1/g
  and are placed along it.

Results:
- **Adjacent cells:** at most 1.198 on both STL wings, sharp and blunt, at
  curvature 1 and 3, and at most 1.207 on the analytic corner loops.
- **Corners:** still exact.
- **Thin collars:** the rectangle and box-edge collars with curvature 1 are
  valid at 1.9 × h (min sJ 0.28 and 0.22).

**The STL wing collar** still folds at every thickness: 524 cells at
1.9 × h uniform, about 6000 with curvature. The cause is not the corners.
ADR 0008's hyperbolic **surface** marcher flips cells at the clustered
leading edge once the surface step reaches 5-15 × the front spacing. This
happens even when marching straight up a prismatic wing, and the uniform
loop flips there too, being too coarse.

## Decision

User decision of 2026-10-27, the bounded step:

- **Corner curves with `curvature > 0`:** `intersection_loop` and
  `intersection_curve` place their nodes by the graded spacing function
  over the whole curve. Adjacent cells stay within `max_cell_ratio` (1.2)
  everywhere, including across corners and node 0.
- **Unchanged:**
  - curves with corners and `curvature = 0` keep ADR 0021's uniform cells
    per segment;
  - curves without corners keep ADR 0020's distribution, pinned bit for
    bit by ADR 0021's golden file.
- **The single collar block at convex corners up to 90°** is accepted
  within the measured thickness: 4.7 × h with a uniform loop and 1.9 × h
  with curvature 1. These are frozen as regression guards.
- **Recorded for later decisions:**
  - ADR 0020's limiter flaw on smooth curves (a re-freeze changes saved
    projects with curvature collars);
  - the surface marcher at a clustered leading edge;
  - concave and sharp corners;
  - thick collars at corners.

## Alternatives considered

| Alternative | Why not |
|---|---|
| j-smoothing near corners | Moves whole i-lines against fixed walls; thousands of folds. |
| More hyperbolic layers | 400-2500 folds. |
| Corner fan (single block) | Collapsed wall cells, which Xnavis does not accept. |
| Snap corners to a global distribution | Ratio 28: rounding is a whole cell where a segment is one cell long. |
| Limiter iterated with per-segment rescale | Infeasible when a segment has too few cells to grade. |
| Replace ADR 0020's limiter on smooth curves too | A re-freeze of ADR 0021's golden; saved projects with curvature collars would not replay. Left to the user. |

## Consequences

- Loops with corners and curvature are graded smoothly. A blunt trailing
  edge's base is no longer one tiny cell between large ones.
- Collars at convex 90° corners are valid, but thin: under 5 × h, or
  under 2 × h with curvature. Thicker ones, concave corners and sharp
  trailing edges need the split collar.
- The wing-plate collar is blocked by surface marching at a clustered
  leading edge. That is the next problem to solve.

## Verification

### Acceptance spec (frozen 2026-10-27, before implementation)

**Spec.** `tests/spec/test_corner_grading.py` runs in the main CI job.

**Verdict.** It must report 13 passed, 0 skipped, 0 xfailed. At freeze it
reports 5 passed and 8 xfailed. The 5 are regression guards: the box-edge
loop, which is graded already because it has no segment much shorter than
its cells, and the four collar cases. Validated with the prototype behind
the interface: 13 passed.

| Case | Required |
|---|---|
| G1 (× 4) | STL NACA 0012 wing through a plate, blunt and sharp, curvature 1 and 3: adjacent chords ≤ 1.22, TE corners on nodes |
| G1 (× 4) | lens, rectangle, L-shape (pending) and box edge (guard), curvature 1: ≤ 1.22, corners exact |
| G1 | open rectangle at x = 0, curvature 1: ≤ 1.22, corners exact |
| G2 (× 4) | rectangle and box edge × sphere r 0.5 collars: uniform at 4.7 × h, min sJ ≥ 0.025; curvature 1 at 1.9 × h, min sJ ≥ 0.1 |

**Out of scope:**
- thicker collars at corners;
- concave and sharp corners;
- the STL wing collar;
- ADR 0020 on smooth curves.

### Review (2026-10-27, user decisions)

1. ADR 0020's limiter on smooth curves is replaced later, in its own
   re-freeze ADR. That ADR first measures how far clusters move on real
   cases.
2. Next step: surface marching at a clustered leading edge.

## Acceptance (2026-10-27): corner-grading spec met

**Kernel.** `_resample` with `curvature > 0` builds the monitor over the
whole curve. `_graded_spacing` caps `C / m` by cones of slope `ln 1.2` in
both directions (tripled for a loop), with C found by bisection on the
cell count. Segments get cells by the integral of 1/g on a sampling at
0.25 of the smallest spacing and are placed along it. With
`curvature = 0`, ADR 0021's uniform cells per segment are unchanged.

**Verdict for `tests/spec/test_corner_grading.py`: 13 passed, 0 skipped,
0 xfailed** (T2; `scripts/check_spec.py` reports all files ok). The 3
`@PENDING` lines are removed. ADR 0020's and 0021's specs pass unchanged.

**Unit test** (`tests/test_intersection.py`): on open and closed curves,
the graded spacing's slope is at most `ln 1.2`, the cell count is exact,
and the cluster stays at the corner.

**Docs.** `docs/reference/collar.md` gives the thickness limits at blunt
corners.
