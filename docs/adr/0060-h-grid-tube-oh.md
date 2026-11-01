# 0060 — An O-H cross-section for the H-grid tube

- **Status:** Accepted
- **Date:** 2026-11-01

## Context

`HGridTubeMesh` (`shore primitive h-grid-tube`) filled a straight circular
tube with one block on a polar grid: i along the stream, j round the wall
(periodic), k from the centreline to the wall. Measured on 2026-11-01 (T1):

1. **The centreline was a wall.** The k = 0 row collapses every ring onto
   one point, and its face was marked `WALL`. So `cc.par` wrote a viscous
   wall (code 1) and CGNS a `BCWall` on a line inside the fluid. ADR 0022
   records that Xnavis does not accept collapsed wall cells. ADR 0055 left
   the face as it was.
2. **The first ring of cells was degenerate.** Each cell next to the axis
   has a face collapsed to a line, and `shore.quality` scores every one of
   them −1.
3. **The block was left-handed for two of the three axes.** The cells off
   the axis scored −0.966 for `--streamwise-axis x` (the default) and `z`,
   and +0.966 for `y`. The cross-section axes for `y`, `("x", "z")`, were
   not cyclic like the other two.

The module's justification did not hold either:
- the O-grid has no polar axis: it stops at its cap angle;
- nothing verified that Xall handles the singularity.

**The owner chose the O-H cross-section (option B, 2026-11-01).**

## Decision

1. **Five blocks: `ring_0` .. `ring_3` and `core`.**
   - The core is a square whose corners lie at `core_frac` (0.5) of the
     radius, with a Cartesian grid inside.
   - Ring `q` joins the core's side to a quarter of the wall, from −45 + 90
     q to 45 + 90 q degrees. Its lines run straight from the side to the
     wall, so at the core's corners they follow the 45° diagonals: the
     three blocks there meet at 90 / 135 / 135°, and the worst cell is
     about sin 135° = 0.707.
   - Every block is right-handed, for every streamwise axis: the cross axes
     are cyclic and each block's j is reversed when needed.
2. **Faces:**
   - the rings' `k_hi` (the wall) is `WALL`;
   - every block's `i_lo` / `i_hi` (inflow and outflow) is `FREE`;
   - the 16 faces between blocks are `SHARED`, with axis maps from their
     corner nodes (`shore.mesh.seams.wire_by_corners`, generalised from the
     tip caps' wiring).

   There is no wall on the axis, and there is no axis.
3. **Parameters:**
   - `n_circ` is the cell count round the wall and must now be a multiple
     of 4;
   - `n_radial` is the rings' node count from the core to the wall;
   - `k_law="tanh"` clusters at the wall, as the docstring always promised;
     the polar build clustered at the centreline;
   - `label` is ignored, with a `DeprecationWarning`;
   - `core_frac` is new.

   `blocks[0]` is a ring, so ADR 0055's B4 (the tube keeps its wall at
   `blocks[0].k_hi`) holds unchanged.
4. **The polar layout is removed, not kept as an option:** it is invalid
   on all three counts above.

## Alternatives considered

- **Option A:** make the axis `FREE` (or an axis code) and fix the
  handedness. The degenerate ring stays and still fails every quality
  check.
- **A bulged core** (sides curved towards the wall): better cells near the
  corners. The straight square already gives 0.707; this is left for later.

## Consequences

- Breaking for the tube's users:
  - one block becomes five, with new labels;
  - `n_circ` must be a multiple of 4;
  - `k_law="tanh"` now clusters at the wall.
- Projects saved with a tube before ADR 0060 replay to the new blocks, so
  their saved grids do not match. Load them without replay.
- `shore check` and `quality.summary` no longer report the tube as inverted.

## Verification

`tests/spec/test_h_grid_tube_oh.py` (main job), frozen before the
implementation:
- T1 (×3 streamwise axes): five blocks, every cell above 0.5.
- T2 (×3): walls only on the rings' `k_hi`, on the radius; inflow and
  outflow FREE; 16 SHARED faces matched through their axis maps.
- T3: the cross-section's area equals the `n_circ`-gon's.
- T4: the wall nodes are evenly spaced (a guard).
- T5: the stream and its law (a guard).
- T6: `k_law="tanh"` clusters at the wall.
- T7: `n_circ` must be a multiple of 4.
- T8: the CLI writes five valid files, and the face graph's walls are the
  four ring faces.

At freeze: 2 passed, 10 xfailed.
