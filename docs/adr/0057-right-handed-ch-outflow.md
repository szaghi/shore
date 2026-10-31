# 0057 — Right-handed CH outflow blocks (and the tip caps, recorded)

- **Status:** Accepted for the CH outflow blocks; the wing-tip caps are
  recorded, not changed
- **Date:** 2026-10-31

## Context

Every cell of three CH outflow blocks had a negative scaled Jacobian:
`chan_outflow_top_upper`, `chan_outflow_wake_upper`, and the blunt TE's
`chan_outflow_wake_base`. Their channel axes run i along the chord, j along
the span and k upwards. With the chord along x, the span along z and up
along y, that frame is left-handed, since x̂ × ẑ = −ŷ while k points +ŷ. The
lower-side blocks, whose k points down, are right-handed.

- `quality.summary` counted those cells as folded.
- `shore check` reported them as "Inverted".
- Whether Xall accepts left-handed blocks was never checked; solvers
  commonly assume right-handed ones.

**The owner asked for the fix (2026-10-31).**

## Decision

1. **The left-handed outflow blocks run j against their partners.**
   `fill_channel_blocks` plants their inner face as the partner's face with
   j reversed (`inner_face[:, ::-1]`). Their nodes are the old ones with j
   reversed; every other block is unchanged, bit for bit.
2. **The axis maps say so.** Each seam between a flipped and an unflipped
   block maps j with a flip:
   - upper outflow blocks ↔ the wrap's `upper.i_lo` and `chan_upper.i_lo`;
   - the sharp wake cut, `wake_upper.k_lo` ↔ `wake_lower.k_lo`;
   - the blunt `wake_base` ↔ `wake_lower` and ↔ the wrap's `wake_base`.

   Seams between two flipped blocks stay identity: `wake_upper` ↔
   `top_upper` and `wake_upper` ↔ `wake_base`. The `cc.par` orientation
   codes and the CGNS transforms follow from the maps.

## The wing-tip caps: recorded, not changed

The C-grid's conformal tip caps (`cap_tips`, and `shore cap-tip`, which
uses the same builder) have more than a handedness problem, measured on a
NACA 0012 wing with both tips:

| Blocks | Cells | Problem |
|---|---|---|
| `cap_jlo_lower`, `cap_jlo_wake_base`, `cap_jhi_upper` | all between −1 and −0.99 | mirrored |
| `cap_jlo_le` | all ≤ 0 | degenerate |
| `cap_jhi_le` | 121 of 605 negative, the rest positive | twisted |
| the four `cap_*_te_*` blocks | none: one node along i | empty |

- The mirroring comes from one shared construction for both sides and both
  tips: the upper and lower caps take j from the chord cut to their own
  wall, and k from the tip wall inwards.
- A flip would fix the three mirrored blocks. The LE and TE blocks need a
  new construction, which an axis flip cannot give.
- **Proposed:** a separate ADR to redesign the cap (an H-cap with a proper
  LE block and no empty TE blocks), with a frozen spec on validity and
  seams.

## Consequences

- CH grids are right-handed in every block; `shore check` and
  `quality.summary` stop reporting the outflow blocks.
- Code that indexed the upper outflow blocks against the wrap's j now
  reverses it. One unit test (`test_ch_topology`'s wake-seam test) was
  changed accordingly.

## Verification

`tests/spec/test_ch_handedness.py` (main job), frozen before the
implementation, for the sharp and the blunt TE:
- H1: every cell positive.
- H2: every SHARED face matches its partner through its axis map (a guard;
  a mutation check with one map back to the identity fails it).
- H3: the grid equals the base commit's (`data/ch_handedness_before.npz`),
  with the flipped blocks' j reversed.

At freeze: 2 passed, 4 xfailed.
