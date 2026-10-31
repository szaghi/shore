# 0056 — Refuse folded cells in the O-grid and C-grid marches

- **Status:** Accepted
- **Date:** 2026-10-31

## Context

ADR 0052 gave the cubed-sphere march a guard that judges every new cell by
its minimum corner scaled Jacobian (`shore.quality`). The other marched
topologies kept only the old check, `_check_jacobian`: one corner of each
cell, compared with an absolute tolerance of −1e-6. It lets folds through:
- a fold at another corner of the cell;
- an inverted cell whose volume is below 1e-6.

**The owner asked for the fix (2026-10-31).**

## Investigation (2026-10-31, `.night/2026-10-31-bcs-guards-handedness/proto/`; T1)

With the new check disabled, a search over O-grid (capsule, box) and
sharp-TE C-grid (NACA 0012 wing) marches found folds the old check accepts:

| Case | Old check | Folded cells | Min scaled Jacobian |
|---|---|---|---|
| O-grid, box 2 × 1 × 0.6, ni 25, nj 33, nk 20, ds 0.01 | accepted | 4 | −0.004 |
| C-grid, NACA 0012 wing, nk 20, ds 2e-3, growth 1.2 | accepted | 4 | −0.98 |

The specs' standard cases are clean in every marched block: O-grid on a
sphere, C-grid with growth 1.1, CH's wrap.

## Decision

1. **The same guard everywhere.** `sphere_cap._check_folds` is called after
   the old check, which stays, in three marchers:
   - `_march_ogrid`, with `periodic_j=True`, so the cells across the j wrap
     count;
   - `_march_cgrid_three`, on its upper, body and lower blocks; CH's wrap
     uses it;
   - `_march_cgrid`, the 1-block legacy grid.

   A cell at or below 0 raises `MarchingError`. The message names the block,
   layer and cell, and the remedy: the hyperbolic scheme, a smaller growth,
   more smoothing.
2. **Marches without a fold are unchanged, bit for bit.** The guard only
   raises.
3. CH's channel and outflow blocks are filled after the march, not marched,
   so the guard does not apply to them (see ADR 0057).

## Consequences

- No marched topology hands out a folded grid without an error any more:
  "SHORE never writes a folded grid" now holds for every marcher.
- C-grids marched deep, about a quarter chord with growth 1.2, are refused
  where they used to fold; use the hyperbolic scheme.

## Verification

`tests/spec/test_marcher_folds.py` (main job), frozen before the
implementation:
- G1–G3: the O-grid on the box and the C-grid wing are refused, from the
  API and the CLI.
- G4 (×3): clean O-grid, C-grid and CH cases march with every marched cell
  positive.
- G5: the clean O-grid and C-grid equal their grids from the base commit
  (`data/marcher_folds_clean.npz`).

At freeze: 4 passed, 3 xfailed.
