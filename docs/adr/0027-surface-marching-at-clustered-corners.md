# 0027 — Surface marching at clustered corners

- **Status:** Accepted
- **Date:** 2026-10-28

## Context

ADR 0025 found collars whose wall itself folds: the plate's surface march
flips cells with curvature clustering at n_loop 192, at the L-shape's
concave corner (10 cells) and at the blunt base (4). No volume construction
can repair a folded wall.

User decisions of 2026-10-28:
- cover concave corners and the blunt base, after a broader sweep;
- diagnose first;
- marches bit-identical where the fix does not trigger;
- acceptance: no flipped cells, collar valid, drift and spacing bounds,
  goldens unchanged.

## Investigation (2026-10-28, prototypes in the session scratchpad; T1)

**Sweep.** Periodic marches on both bodies:
- junctions: rectangle, L-shape, blunt and sharp STL wings, lens, box
  edge, 90° and 60° wedges;
- n_loop 96 / 192 / 384, curvature 0 / 1 / 3.

| Regime | Flipped cells |
|---|---|
| R1 plate at a concave corner (L-shape, n_loop 192-384) | 4-34, growing with refinement |
| R2 plate at the blunt base (curvature, n_loop 192-384; cells 1.6e-4 to 3e-5) | 2-24 |
| R3 plate round sharp tips (single block) | the split collar's job (ADR 0024) |
| R4 along sharp creases on B (curvature, n_loop 192-384) | 94-1033; ADR 0024's fix holds at n_loop 96 only |
| R5 sharp STL wing, n_loop 384 with curvature | singular solver |
| rectangle, box edge, 90° wedge | 0, up to n_loop 384 |

The split collar shows the same wall folds with curvature at n_loop
192-384: the lens folds 259-586 cells, and the sharp wing hits the singular
solver.

**Two mechanisms.**
- **Convex fans and sharp creases (R2, R4, R5).** Curvature clustering
  makes the corner spacing 1e-4 to 3e-5, while the first surface step is
  5e-3: the step is 20-150 × the spacing. Uniform sub-stepping settles it:

  | March | m = 1 | m = 4 | m = 16 | m = 64 |
  |---|---|---|---|---|
  | blunt base | 4 | 0 | 0 | 0 |
  | sharp-wing crease | 179 | 48 | 0 | 0 |
  | lens crease | 94 | 0 | 0 | 0 |

- **Converging fronts (R1).** The shortening concentrates on the corner
  cells. At n_loop 384, uniform, the corner cell goes from 1.04e-2 to
  3.8e-4 in 8 layers while the step grows to 1.2e-2, until nodes overtake
  each other along the front: a segment reverses, and the next layer's
  cells invert. Sub-stepping does not converge here (3, 177, 6 flips for
  m = 4, 16, 64).

**Triggers that leave valid marches alone.**
- **Sub-stepping.** A step-to-spacing trigger over all nodes would fire on
  ADR 0008's golden concave-L march (ratio 372, from its converged corner).
  Counting only nodes not closer together than on the starting curve, the
  goldens peak at 10.8 (`ellipse_b`) and ADR 0023's crease cases at 21.
  - A trigger of 12 changes a frozen collar (ADR 0025's blunt case, min
    0.023 → -0.008).
  - A trigger of 16 still dips below that case's bound.
  - A trigger of 24 fixes every failing march and leaves that collar
    unchanged (0.0231).
- **Untangling.** It fires only on an actual fold: a new layer's cell
  inverted, or a front segment reversed against the previous layer.
  Checking only inverted cells missed the first reversal. With both checks
  it fixes every concave case.

**ADR 0008's golden concave-L march folds.** Its last three layers have 16
reversed segments, which ADR 0008's spec never checked and ADR 0023 froze
bit for bit. The repair changes it by at most 1.6e-4 and removes them. It
is re-frozen (user decision).

**Results of the prototype.**

| March | flipped (before → after) | other |
|---|---|---|
| L-shape plate, n_loop 192 curvature 1 / 384 uniform / 384 curvature 1 | 10 / 4 / 34 → 0 | min cell sine 0.67 / 0.18 / 0.24 |
| blunt plate, n_loop 192 curvature 1 / 384 curvature 3 | 4 / 24 → 0 | sine 0.72 |
| sharp crease, n_loop 192 curvature 1 | 179 → 0 | drift 0.0995 → 0.0107 |
| lens crease, n_loop 192 curvature 1 / 384 curvature 3 | 94 / 1033 → 0 | drift 0.058 / 1.0 → 0.009 / 0.006 |
| wedge crease, n_loop 384 curvature 1 | 182 → 0 | drift 0.225 → 0.005 |

- **Collar walls no longer fold anywhere:**
  - before: 4-586 cells per case, plus singular-solver errors;
  - after: every wall layer is valid, and the split sharp wing runs.
- **Collar volumes still fold** at these clusterings: the collar is 100-300
  × thicker than the corner spacing, which is the volume's limit (user
  decision: record, do not chase here). Examples:
  - L-shape n_loop 192 curvature 1: 320 → 67 folded cells;
  - blunt wing n_loop 192: 135 → 4 (curvature 1), 49 → 5 (curvature 3);
  - split lens: 5313 → 81.
- **Not fixed:** the sharp STL wing at n_loop 384 with curvature (spacing
  3e-5 on sliver triangles) still flips 2-119 cells.
- **Cost:** sub-stepping takes up to 150 s (blunt collar, n_loop 384,
  curvature 3) and 280 s (split sharp wing, n_loop 384). Accepted
  uncapped (user decision).
- **Full suite with the prototype:** only ADR 0023's C4 for the concave-L
  golden fails, as expected.

## Decision

- **Sub-stepping.** `MarchParams.max_step_ratio` is 24 and
  `substep_ratio` is 2. A layer whose step exceeds 24 × the smallest
  spacing of the nodes not closer than on the starting curve is marched in
  `ceil(ratio / 2)` sub-steps, each projected (and its ends put on their
  planes).
- **Untangling.** `MarchParams.untangle_iterations` is 200 and
  `untangle_reach` is 2. Where a new (sub-)layer folds against the
  previous one (`_inverted`: a cell reversed against the majority, or a
  front segment pointing backwards), the new nodes within 2 of the folded
  cells are smoothed along the front (1-2-1) and re-projected, until none
  folds.
- **Re-freeze.** ADR 0023's C4 for `m_concave_l_plane` becomes C4b: the
  repaired march (`data/march_default.npz` rewritten for that entry only)
  and no reversed segment.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Uniform sub-stepping, step / min spacing over all nodes | Fires on ADR 0008's converging golden (ratio 372). |
| Trigger 12 or 16 | Changes ADR 0025's frozen blunt collar (min 0.023 → -0.008 / 0.008). |
| Sub-stepping for concave corners | Does not converge (3 / 177 / 6 flips). |
| Stronger convergence dissipation everywhere | Changes every converging march, golden included. |
| Untangling on inverted cells only | Misses the first reversal; L-shape 34 → 21. |
| Cap the sub-steps | May leave the extreme cases flipped (user decision: no cap). |

## Consequences

- **Surface grids no longer flip** at clustered convex corners, sharp
  creases and concave corners. Collar walls, single and split, are valid
  there.
- **Runs that trigger sub-stepping are slower,** up to minutes at n_loop
  384 with curvature.
- **Still open:**
  - the collar volume at extreme clustering (thickness / corner spacing
    100-300);
  - the sharp STL wing at n_loop 384 with curvature.

## Verification

### Acceptance spec (frozen 2026-10-28, before implementation)

**Spec.** `tests/spec/test_march_corners.py` runs in the main CI job,
together with the re-frozen `test_march_creases.py` (C4b) and
`data/march_default.npz`.

**Verdicts.**
- `test_march_corners.py`: 15 passed, 0 skipped, 0 xfailed. At freeze: 2
  passed, 13 xfailed.
- `test_march_creases.py`: 22 passed. At the re-freeze: 21 passed, 1
  xfailed.

Validated with the prototype behind the interface: 37 passed.

| Case | Required |
|---|---|
| M1 (× 5) | plate marches at concave corners and the blunt base: no flipped cell, cell sine ≥ 0.1 |
| M2 (× 4) | sharp-crease marches: no flipped cell, drift ≤ 0.02 |
| M3 (× 3) | collar walls (L-shape, blunt wing, split lens; n_loop 192, curvature 1) without folds |
| M4 (× 3) | untriggered collars valid (guards); the parameters |
| C4b | ADR 0008's concave-L march equals its re-frozen golden; no reversed segment |

### Review (2026-10-28, user decisions)

1. The trigger stays at 24 and the target at 2.
2. Next step: the collar volume at extreme clustering.

## Acceptance (2026-10-28): clustered-corner marching spec met

**Kernel.** As decided. `_march` runs `_substeps` sub-steps per layer, each
followed by `_untangle` (which uses `_inverted`).

**Verdicts** (T2; `scripts/check_spec.py` reports all files ok):
- `tests/spec/test_march_corners.py`: **15 passed, 0 skipped, 0 xfailed**;
- `tests/spec/test_march_creases.py`: **22 passed, 0 skipped, 0 xfailed**.

All 14 `@PENDING` lines are removed.

**Unit tests** (`tests/test_surface_march.py`):
- a run of clustered cells gets `ceil(ratio / 2)` sub-steps;
- the same cells after a front converged there get none;
- `_inverted` finds two nodes that overtake each other.
