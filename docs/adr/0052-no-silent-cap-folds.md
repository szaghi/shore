# 0052 — Refuse folded cells in the cubed-sphere march, corner cells included

- **Status:** Accepted
- **Date:** 2026-10-30

## Context

Writing the documentation (night of 2026-10-30) showed that the cubed-sphere
march can return a folded grid without a word. The guard after each layer,
`sphere_cap._check_jacobian`, has two blind spots:

1. **It skips the four corner cells of every block** (`exclude_corners=True`
   for the four equatorial blocks and both caps). Its docstring gives the
   reason: noise from the "degenerate winding" at the block junctions.
2. **It tests one corner of each cell**, with the triple product compared to
   an absolute tolerance of `-1e-6`. That is not scale-free: a thin wall cell
   of volume 1e-7 could invert and still pass.

The folds sit exactly in the skipped cells, at the 3-block junctions of the
cap corners. Neither the march nor `shore check` reports them;
`quality.summary` does. **The owner asked for the fix with an ADR and a
frozen spec (2026-10-30).**

## Investigation (2026-10-30, `.night/2026-10-30-folds-symmetry/proto/`; T1)

A sweep over two bodies (icosphere with 4 subdivisions; capsule with radius
0.5 and height 2), three sizes ((ni, nj, nk, ds) = (33, 64, 25, 0.01),
(41, 64, 31, 1e-3), (25, 48, 21, 0.02)), cap angles 9–50° and both marching
schemes, growth 1.15. Folded means minimum corner scaled Jacobian <= 0:

| theta_cap | algebraic | hyperbolic | C1 fallback warning |
|---|---|---|---|
| 9° (the default, 0.05π) | **folds in all 6 cases**: 4–104 cells, min −0.08 to −0.99 | no fold | yes |
| 20° | no fold (min ≥ 0.17) | no fold | yes |
| 30° | no fold (min ≥ 0.57) | no fold (min ≥ 0.57) | no |
| 35–40° | no fold | no fold | no |
| 45–50° | capsule: no fold at 45°, `MarchingError` at 50° (33/64/25) | **capsule folds**: 8–24 cells, min −0.006 to −0.32 | no |

Every fold is at a cap corner. The existing unit test
`test_cap_gnomonic.py::TestMarchingNoNegativeJacobian` (default cap, nk 20)
had a folded grid (min −0.42) and passed, because it asserted only that the
march did not raise.

## Decision

1. **Every cell counts.** After each layer, `sphere_cap._check_folds` takes
   the minimum corner scaled Jacobian (`shore.quality.scaled_jacobian`, the
   frozen metric) of every cell between the two layers of every block, the
   corner cells included. A cell at or below 0 raises `MarchingError`.
   - The message names the block, the layer, the worst cell and the number
     of folds, and the remedies: a larger cap (30° or more), the hyperbolic
     scheme, a smaller growth.
   - The scaled Jacobian is scale-free, so thin wall cells are judged like
     thick ones.
2. **The old guard stays as it is.** It can only raise where it raised
   before, so the new guard adds refusals and moves no node. A march with no
   fold is bit-identical to before.
3. **Scope:** the cubed-sphere march (`Mesh._march_caps` and the OH wrap,
   which uses it). The other marchers keep their guards, see Consequences.
4. **The default cap angle is not changed here.** 30° is clean in every
   swept case, but the topology default backs a frozen golden: ADR 0034's R8,
   `data/blend_removal_default.npz`. Changing it is a re-freeze and is left to
   the owner, with the evidence above. With the 9° default, many default
   marches are now refused with a message that says what to change, instead
   of returning folded grids.

## Alternatives considered

- **Mend the folds** by smoothing or untangling the cap corners: a new
  algorithm with its own spec; refusing is the honest minimum.
- **Include the corner cells in the old single-corner check:** the folds
  are visible at other corners of those cells, and the absolute tolerance
  stays scale-dependent.
- **Change the default to 30° in the CLI only**: CLI, config and Python
  defaults would disagree.

## Consequences

- `shore mesh body.stl` with the default cap angle now fails on many
  bodies, at the cap corners, with the remedy in the message. Tutorials and
  recipes already use `--theta-cap-deg 35`.
- One unit test was changed: `TestMarchingNoNegativeJacobian` now marches
  with a 30° cap and asserts every cell positive. The test's premise was
  wrong, not the guard.
- **Left open:** the other marchers (O-grid, C-grid, CH) still use one
  corner and an absolute tolerance. Proposed: the same scaled guard, in its
  own ADR once measured on their test grids.
- **Left to the owner:** a 30° default for the cubed sphere (re-freeze of
  ADR 0034 R8's cubed-sphere goldens).

## Verification

`tests/spec/test_cap_folds.py`, frozen before the implementation, in the main job:
- F1–F4: folding cases refused, from the API and the CLI; failed marches
  are not recorded and leave the mesh as it was.
- F5: clean cases march, every cell positive.
- F6: a clean case is bit-identical to the grid marched at the base commit
  (`data/cap_folds_clean.npz`).
- F7: the OH hybrid still marches.

At freeze: 6 passed, 4 xfailed.
