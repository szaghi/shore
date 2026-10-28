# 0029 — Spreading converging fronts (not adopted); a four-corner surface metric

- **Status:** Rejected (the spreading); Accepted (the four-corner metric)
- **Date:** 2026-10-28

## Context

ADR 0028 traced the L-shape collar's remaining folds to the plate's surface
grid at the concave corner, where converging fronts bunch the nodes.

User decisions of 2026-10-28:
- spread converging fronts where they bunch, keeping every golden
  bit-identical;
- freeze acceptance on four-corner cell quality, valid L-shape collars, no
  regression and a spacing-ratio bound;
- re-freeze ADR 0027's M1 to a four-corner metric.

## Investigation (2026-10-28, scripts in the session scratchpad `sp/`; T1)

**The four-corner sine** is the sine at each corner of a surface cell, of
the two edges leaving it, against the surface normal. It corrects ADR 0028:
the L-shape plate quads are convex (minimum 0.06-0.20), not flat. What they
are is bunched. Adjacent front cells differ by up to 9.5 ×, against ≤ 1.9 ×
on every other junction and ≤ 2.32 on the golden marches.

**Spreading.** Where adjacent front cells differ by more than 3 ×, the nodes
round them are smoothed along the front (1-2-1) and re-projected, after
ADR 0027's untangling.
- **Absolute trigger** (ratio > 3):
  - bunching is bounded at 3 ×;
  - the corner sine at n_loop 192 rises from 0.20 to 0.66;
  - L-shape collar folds 67 → 42, 100 → 38, 1169 → 896.

  But it changes ADR 0024's golden blunt-wing collars: that loop places a
  0.0025 base cell next to much larger ones on purpose.
- **Relative trigger** (ratio > 3 × the ratio on the starting curve):
  - every golden is bit-identical (full suite: 1677 passed);
  - collar folds 67 → 68, 100 → 38, and **1169 → 1551** (L-shape at
    n_loop 384 with curvature).

  That last case is a regression.

The residual L-shape folds sit at the concave corner station, in the
plate's far-edge row, from layer 6: still the near-wall hyperbolic march,
whose boundary rows take the scheme's local solution. Spreading the front
does not reach them.

## Decision

User decision of 2026-10-28: the spreading is **not adopted**.
- The version that keeps the goldens regresses one collar.
- The version that helps the collars changes frozen goldens.

The **four-corner metric** is adopted as a guard. ADR 0027's M1 is
re-frozen to it with bound 0.05; every case holds it today (minimum 0.0626,
L-shape at n_loop 384 uniform).

## Alternatives considered

| Alternative | Why not |
|---|---|
| Spreading with an absolute ratio trigger | Changes ADR 0024's golden collars (intentional base-cell clustering). |
| Spreading relative to the starting curve | Regresses the L-shape collar at n_loop 384 with curvature (1169 → 1551). |
| ADR 0027's one-corner sine | Misses nearly flat corners; replaced by the four-corner sine. |

## Consequences

- **No code change.** Concave-corner collars at n_loop 192-384 still fold
  in their volume.
- **Next lever:** the far-edge row of the near-wall hyperbolic march at the
  concave corner station (its boundary-row treatment), for a later ADR.
- **The surface-quality guard is stricter:** four corners per cell.

## Verification

**Re-frozen spec.** `tests/spec/test_march_corners.py` (M1, four-corner
sine ≥ 0.05): 15 passed, 0 skipped, 0 xfailed, on the current code.
