# 0025 — The collar volume at loop corners

- **Status:** Proposed
- **Date:** 2026-10-28

## Context

ADR 0024 split the collar at sharp corners and found that every other fold
at a loop corner is in the volume, not the wall:
- thick collars at convex corners;
- concave corners;
- the blunt trailing edge.

ADR 0010 builds the volume above the near-wall layers per loop station, as
a 2-D (i, k) sheet:
- the station's i-line is offset in its best-fit plane;
- the swallowtail at the concave crease is trimmed;
- the nodes are redistributed by arc length;
- TFI and Winslow fill the sheet.

The user asked for this as its own ADR (decision of 2026-10-27). It was
done overnight on the recommendations, as the user asked.

## Investigation (2026-10-27/28, prototypes in the session scratchpad; T1)

**Dissection of folded cells.** For each folded cell, at its worst corner,
the edge vectors were examined. Four mechanisms showed up:
1. **The outer line's out-of-plane coordinate came from the wrong node.**
   `_outer_line` interpolates the out-of-plane part at an index that is
   linear in arc fraction. With geometric surface steps (1.12^j), a node
   takes another node's offset. Where a station's i-line curves out of its
   best-fit plane (near loop corners) that is a sideways shift. On the
   L-shape: 0.009, ten times the 0.0009 j-spacing.
2. **The best-fit plane twists the stations next to corners.** The
   plate's i-lines curve in plan there (plane residual 0.17 against 0.05
   away from corners). Flattening them moves the outer points sideways
   along the loop.
3. **One station alone cannot fan round a convex corner, nor stop crossing
   at a concave one.** Neighbouring stations' outer boundaries interfere
   once the collar is thicker than the loop spacing.
4. **Near-wall layers at two close corners (the blunt base).** The central
   difference across unequal cells tilts the corner normal.

**What was tried, on 9 benchmark cases (323 folded cells before):**

| Change | Folded cells |
|---|---|
| layer smoothing in (i, j), plain or tangential, ramped in k | 500-15000: it shears the thin k cells |
| station offsets rotated to a smoothed loop tangent | 290-610 |
| full 3-D Winslow with the outer boundary fixed | rectangle 20 → 12, 300 s |
| index fix alone | 399 |
| 3-D offset polyline (no flattening) alone | 341 |
| removing the loop-direction part of the offsets | 282 |
| corner-aware normals | 309 |
| **station offsets smoothed along the loop (σ = 0.5 thickness, one kernel for all i)** | **130** |
| + 3-D offset polyline + local smoothing where the boundary steps back | 49-53 |
| + each station's line trimmed again | 20-24 |

**Validation on 49 cases.** The cases were not used for tuning:
- rectangle, L-shape and blunt wing at n_loop 96 and 192, 1.9 / 4.7 /
  11.8 × h, uniform and curvature 1;
- box edge;
- ellipse;
- split sharp wing and lens.

| | `"stations"` (ADR 0010) | `"smoothed"` |
|---|---|---|
| all 49 cases | 5332 folded cells | 1843 |
| without the 6 whose wall folds | 2367 | 72 |
| rectangle, all | 284 | 0 |
| L-shape, uniform | 1144 | 12 |
| split collars, curvature, 4.7 × h | 351 | 0 |

- **Six cases fold in the wall itself:** the L-shape and the blunt wing
  with curvature at n_loop 192. There the plate's surface march flips 10
  and 4 cells, so no volume construction can help. That is a marcher
  problem.
- **What remains in the volume:**
  - concave corners at 4.7-11.8 × h: 2-19 cells;
  - the blunt base at n_loop 96: 2-15 cells;
  - one regression: the blunt wing at 1.9 × h with curvature, 6 → 10.
- **No harm where `"stations"` is valid:**
  - box edge 0.161 (0.168 before);
  - ellipse 0.270 (0.301);
  - ADR 0022's rectangle at 4.7 × h 0.048 (0.040);
  - half ellipse 0.309 (0.315).
- **Cost:** the same (7.4 s against 7.6 s for the rectangle at n_loop 192).

## Decision

**`CollarParams.volume`.**
- `"stations"` is the default: ADR 0010, bit for bit.
- `"smoothed"` is the new construction. It applies to `build_collar`
  (closed and half), `split_collar` and `extrude_collar`.
- It reaches the project and CLI as `BuildCollar.volume`,
  `Project.build_collar(volume=)` and `shore collar --collar-volume`.

**`"smoothed"` (`collar._smoothed_outer`).**
1. **Offset line.** Per station, the outer line lies on the station's 3-D
   offset polyline (`_offset_line`). The swallowtail is found in the
   best-fit plane, as before. The crossing becomes the midpoint of the two
   3-D crossing points (`_trim_origins` keeps both segments' fractional
   origins). The nodes are placed by arc length in 3-D, with the crease on
   the junction row. Nothing is flattened, which removes mechanisms 1 and 2.
2. **Smoothing along the loop.** The offsets `outer - inner` are smoothed
   along the loop by a Gaussian of width `offset_smoothing` (0.5) × the
   thickness, in junction-row arc length, with one kernel for every i.
   Rows stay aligned in i; stations agree with their neighbours.
3. **Local repair.** Where some row of the outer boundary advances along
   the loop by less than `backstep_ratio` (0.5) of the inner row, the two
   stations' offsets are smoothed again (1-2-1), until none does (at most
   `backstep_iterations`, 200).
4. **Second trim.** Each station's line is trimmed and redistributed again.
   Averaging offsets of differently turned stations can reorder the nodes
   round the crease image: that caused the rectangle's 4.7 × h regression
   (0 → 4 cells).

**Not adopted:** corner-aware normals. They are worth ≤ 10 cells on the
blunt base and cost the rectangle 4.

**Default.** `"stations"` stays the default until the user decides. Making
`"smoothed"` the default changes every collar, so ADR 0024's golden
collars would have to be re-frozen. The evidence above is for that
decision.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Split blocks at convex corners | A corner block's wall collapses onto the edge (ADR 0022, 0024). |
| 3-D Winslow over the volume | Fixed-boundary folds stay; 40 × slower. |
| (i, j) smoothing of the layers | Shears the thin k cells: thousands of folds (also ADR 0022). |
| Per-row redistribution of the outer surface in j | Breaks the alignment in i (worse). |
| Outer boundary with the loop-direction part removed | 282 folded cells; the stations still disagree. |
| Make it the default now | Re-freezes ADR 0024's golden collars without the user. |

## Consequences

- **Thick collars at convex corners are valid with `"smoothed"`:**
  prismatic struts, rectangular bases and the split collars with curvature
  clustering. Concave corners are much better at 1.9-4.7 × h.
- **Still open:** concave corners and the blunt base at larger thickness,
  and plate surface marches that flip with curvature at n_loop 192 (a
  marcher ADR).
- **Default:** the user decides whether to switch to `"smoothed"`; that
  needs a re-freeze.

## Verification

### Acceptance spec (frozen 2026-10-28, before implementation)

**Spec.** `tests/spec/test_collar_volume.py` runs in the main CI job and
reuses ADR 0024's golden `data/collar_default.npz`.

**Verdict.** It must report 25 passed, 0 skipped, 0 xfailed. At freeze it
reports 3 passed and 22 xfailed. The 3 are regression guards: the default
collars. Validated with the prototype behind the interface (a copy without
`@PENDING`): 25 passed.

| Case | Required |
|---|---|
| V1 (× 11) | `"smoothed"` valid where `"stations"` folds, with bounds at about half the measured minimum (rectangle, L-shape, blunt wing, split sharp wing and lens) |
| V2 (× 5) | no harm where `"stations"` is valid (box edge, ellipse, ADR 0022's rectangle at its 0.025 bound, split wing, half ellipse on its plane) |
| V3 (× 7) | default `"stations"`; default and explicit `"stations"` equal to ADR 0024's golden collars |
| V4 | `ParameterError` for an unknown volume |
| V5 | `Project.build_collar(volume="smoothed")` equals the kernel and is journalled |

**Out of scope:**
- concave corners at 4.7 × h and thicker;
- the blunt base at n_loop 96;
- loops whose surface grid folds;
- the default switch.
