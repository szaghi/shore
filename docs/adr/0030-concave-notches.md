# 0030 — Collars in concave notches (investigation)

- **Status:** Rejected (no change adopted)
- **Date:** 2026-10-28

## Context

ADR 0029 left the L-shape collar's folds at its concave corner station, in
the plate's far-edge row. User decisions of 2026-10-28:
- L-shape collars;
- diagnose first;
- bit-identical where nothing triggers;
- acceptance: valid collars, measured reduction, no regression.

Midway the user allowed two more attempts: a concave cut, and a
redistribution that would re-freeze ADR 0008's concave-L golden.

## Investigation (2026-10-28, scripts in the session scratchpad `fe/`, `sp/`; T1)

**The mechanism is geometric, not the boundary rows.**
- The fold cell at the plate's far edge on the bisector is 0.037 long and
  1.4e-5 wide (aspect about 2700).
- Interior rows creep along the bisector by about 3e-6 per near-wall layer,
  while the boundary row stays put. Heights are uniform; that shift alone
  folds such a cell.
- Marching into a notch, the front loses about 2D of length at distance D
  from a 90° concave corner. Nodes march along their normals, so the whole
  loss lands on the cells at the bisector.

**Three remedies, none adopted.**

| Remedy | L-shape collar folds (n_loop 192 curvature 1 / 384 uniform / 384 curvature 1); before 67 / 100 / 1169 |
|---|---|
| Spread bunched cells (ADR 0029) | 68 / 38 / 1551 |
| Concave cut: one piece, both ends held on the bisector plane | 414 / 112 / 1507 (+23 wall folds); the two cut faces differ by up to 0.06 |
| Share the contraction along the front, full move per layer | 2801 / 51 / 6064, with surface and wall folds |
| Same, partial move β = 0.05 | 19 / 40 / 986, with surface folds |
| Same, partial move β = 0.2 | 557 / 3 / 1739, with surface folds |

- Holding the ends on the bisector plane lets the end node slide, but its
  neighbours still converge toward it.
- Moving nodes along the front far from the node below them makes the cells
  between layers shear over.
- Every partial move fixes one case and breaks another, re-introducing
  surface folds that ADR 0027 had eliminated.

## Decision

No change. Collars whose plate march runs deeper into a concave notch than
the notch is wide (here a 0.186 march into a 0.5 notch, at n_loop 192-384)
fold at the bisector. The documented limit: keep the collar's width on A
small at concave corners, or use n_loop 96 there, where they are valid up
to 1.9 × (ADR 0025).

A fix needs nodes that leave the front where it converges, which a single
structured block cannot do: for example a grid whose j-lines end at the
bisector (an H-type block in the notch). That is a topology for a later ADR.

## Consequences

- **No code change.** ADR 0008's concave-L golden is not re-frozen.
- **The collar reference** states the notch limit.

## Verification

None (an investigation record). The measurements are T1.
