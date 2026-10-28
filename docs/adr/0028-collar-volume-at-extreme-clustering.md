# 0028 — The collar volume at extreme clustering (investigation)

- **Status:** Proposed (investigation; no change adopted)
- **Date:** 2026-10-28

## Context

Since ADR 0027 the surface grids no longer flip at clustered loop corners,
and collar walls are valid. With curvature clustering at n_loop 192-384,
though, the collar volumes still fold, while the collar is 100-300 × as
thick as the corner spacing. ADR 0027 recorded this as the volume's limit,
and the user chose it as the next step (2026-10-28).

The user was away for this step and asked for it to go ahead on my
recommendations:
- the cases are ADR 0027's targets;
- diagnose first;
- bit-identical where nothing triggers;
- no folds where reachable.

## Investigation (2026-10-28, scripts in the session scratchpad `xv/`; T1)

**Where the folds are** (thickness 1.9 × h, smoothed volume):

| Collar | Folded cells | Where |
|---|---|---|
| L-shape, n_loop 192 curvature 1 | 67 | plate side, from layer 1 up, far edge near the concave bisector |
| L-shape, n_loop 384 uniform | 100 | same |
| blunt wing, n_loop 192 curvature 1 / 3 | 4 / 5 | last layer, B side |
| blunt wing, n_loop 384 curvature 1 | 61 | layers 9-13, junction and B side |
| split lens, n_loop 192 curvature 1 | 81 | B blocks, at the stations next to the cut |
| split sharp wing, n_loop 192 curvature 1 | 22 | B blocks, at the stations next to the cut |
| split 60° wedge, n_loop 192 curvature 1 | 742 | B blocks, at the stations next to the cut |

### L-shape: the plate grid itself is degenerate near the concave bisector

The first fold is in a near-wall hyperbolic layer (k0 = 9, so these are
not the outer part). Its bottom quad at the plate's far edge:
- is 0.038 long and 9e-5 wide;
- has its front direction 18° off the marching direction;
- has one corner close to 180°: the nodes are almost in a line, bunched
  where the fronts from the two faces met.

Near-wall layers move such nodes sideways by about 1e-5, enough to make the
quad non-convex. A repair that replaces folded near-wall nodes by a pure
normal offset fires but cannot help: offsetting a nearly flat quad still
tips it over.

ADR 0027's untangling keeps the front in order but does not spread the
bunched nodes. Its spec metric (the cell sine at one corner per cell) misses
these flat corners.

### Split collars: the volume at the cut

- **Where:** the folds run through the whole B height at the stations next
  to the sharp corner, from layer 4 or 9 upward, where cells are 40 × longer
  than wide.
- **Why the smoothing never acts there:** along a split piece the cut turns
  the wall by only 8-30°, below ADR 0026's 60° trigger, so the stations
  next to the cut are never smoothed.

Tried:

| Change | wedge / lens / sharp wing (n_loop 192, curvature 1) |
|---|---|
| current | 742 / 81 / 22 |
| near-wall layers by pure normal offset | 565 / 761 / 354 |
| offset smoothing ×2 or ×4 | unchanged (it never triggers) |
| `backstep_ratio` 0.8 | 620 / 57 / 16 |
| smoothing forced everywhere | 402 / 49 / 9 |
| trigger lowered to 5° | 402 / 49 / 9 |

None reaches zero.

### Blunt wing

At n_loop 192 only 4-5 cells fold, on the last layer. At n_loop 384 there
are 61, in layers 9-13 at the junction. Not investigated further.

## Decision

No change. Two root causes are identified, and both remedies need the
user's decision:
1. **Spread converging fronts in the surface march**, not only keep them
   in order. Every converging march would change, including ADR 0008's
   concave-L golden that ADR 0027 just re-froze.
2. **A volume construction that holds when the collar is 20-300 × as thick
   as the local spacing.** ADR 0025's levers are exhausted at the cut of a
   split collar. Candidates:
   - a hyperbolic or elliptic 3-D volume with a free outer boundary;
   - collars that thin out where the spacing is extreme;
   - coarser clustering of the loop at sharp corners for split collars.

Also recorded: ADR 0027's spec metric should check all four corners of a
surface cell (the flat corners above pass its one-corner sine).

## Alternatives considered

The rows of the tables above, all short of zero folds.

## Consequences

- **What works:** collars at extreme clustering have valid walls (ADR
  0027), and their volumes fold far less than before ADR 0025-0027. They
  are not valid yet.
- **Practical guidance** (documented in the collar reference): keep the
  collar thinner than about 10 × the smallest loop spacing at corners, or
  use less curvature clustering there.

## Verification

None: an investigation record, with no code and no spec. The measurements
are T1, reproducible with the scratchpad scripts named above.

## Correction (2026-10-28, ADR 0029)

The four-corner sine shows that the L-shape plate quads are convex (minimum
0.06-0.20), not flat: they are bunched, with adjacent front cells up to
9.5 × apart. Spreading them (ADR 0029) was not adopted.
