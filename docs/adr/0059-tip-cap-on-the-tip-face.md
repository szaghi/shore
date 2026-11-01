# 0059 — Wing-tip caps on the tip face: an O-H overset cap

- **Status:** Accepted
- **Date:** 2026-10-31

## Context

The C-grid's tip caps (`cap_tips`, `shore cap-tip`) were wrong in more
than their orientation (ADR 0057). Since ADR 0058 they are refused.
- **The cap was a slab of the solid wing.** The C-grid's tip stations sit
  2 % of the span inside the tips. The cap filled the section between the
  station and a plane `tip_inset` (1e-3 chord) further out, with walls on
  its outer plane and on the section's sides. All of that lies inside the
  wing. The tip face itself had no grid, and neither did the strip of side
  wall between the station and the tip.
- **The blocks did not hold together.**
  - The upper and lower blocks shared one construction, so half of them
    were mirrored.
  - The LE block wrapped the LE between two cuts onto a line, which made it
    degenerate.
  - The TE blocks were one station long, so they had no cells.
  - A wake-base cap block meshed fluid, not wall.

**The owner asked for the redesign (2026-10-31).**

## Investigation (2026-10-31, `.night/2026-10-31-tip-cap/proto/`; T1)

An O-H footprint, extruded along the tip plane's normal, was prototyped on
an extruded NACA 0012 wing:
- sharp (rounded 0.2 %) and blunt TE;
- `ni_body` 40, 41, 80, 120 and 201;
- `le_frac` 0.02–0.1;
- flare 0, 30 and 45°.

The worst cell was 0.30 (flare 45°). At the defaults the worst is 0.41 on
NACA 0012 and 0.41 on a swept (20°), tapered (0.6) NACA 4412.

A TE ring that also took surface stations (a TE cut at 95 % chord) fell to
0.15 at the cut. So the TE block spans the TE base alone.

## Decision

1. **The cap is an overset component on the tip face.** Each tip gets five
   blocks, whatever the TE:
   - `core`, inside the section, inset by `core_frac` (0.5) of the local
     half-thickness from both sides;
   - a ring joining the core to the section's wall: `upper` (the C-grid's
     i = 0 side), `le` (round the LE, between the LE cuts at `le_frac`),
     `lower`, and `te` (across the TE base).

   Ring lines run straight from the core to the wall nodes. The core's
   corners sit on 45° diagonals from the LE cuts and the TE corners. So
   each corner where three blocks meet splits about 90 / 135 / 135°, and
   no corner is 180°.
2. **The wall is the section, on the tip face.** `k_lo` of every cap block
   is WALL. Together the five faces cover the section exactly, and its
   edge nodes are the C-grid's own tip-station wall nodes, bit for bit.
   - With `cap_tips` the C-grid's first and last stations are put on the
     tips (`TIP_STATION_MARGIN` = 1e-5 of the span, instead of 2 %). So the
     wrap covers the whole side wall and the cap the tip face.
   - `shore cap-tip` reads a body made without caps. Its default
     `tip_inset` moves the wall out to the tip face: `0.02 / 0.96` of the
     distance between the first and last stations.
3. **k leaves the wall along the plane's outward normal.** The layers use
   the wrap's wall-normal schedule: `spacing_k` on the C-grid, the body's
   mean k-line lengths for `shore cap-tip`. Ring lines lean outwards by up
   to `flare_deg` (default 30°), growing linearly from the core to the
   edge, so the cap overhangs the side wall it must overlap.
4. **Faces:** every block is right-handed.
   - The 8 faces between blocks are SHARED (16 SHARED faces per tip).
     `wire_cap_blocks` finds their axis maps from the corner nodes.
   - Every other face is FREE: a fringe for the overset assembly.
   - No `wake_base` cap block: the wake base is fluid.
5. **Interface changes:**
   - New parameters: `cap_core_frac` / `cap_flare_deg` on the C-grid,
     `--cap-core-frac` / `--cap-flare-deg` on `shore mesh`, `core_frac` /
     `flare_deg` on `CGridTipCap` and `shore cap-tip`.
   - The TE cut (`cap_te_frac`, `--te-frac`, `te_frac`) is ignored, with a
     `DeprecationWarning`.
   - `tip_inset` is the wall's offset from the station, not a slab
     thickness.
   - On `CGridTipCap`, `is_blunt`, `wake_base_cross_section` and the
     TE-cut indices are ignored, with a warning.
6. **ADR 0058's guard stays** after construction.

## Alternatives considered

- **A conformal cap SHARED with the wrap's end faces:** the tip wall and
  the wrap's end face are coplanar along the tip edge, so the cells there
  would have 180° corners. A real conformal closure needs the wrap to turn
  round the tip (an O-O topology), which is a new topology.
- **Keeping the TE cut at 95 % chord:** worst cell 0.15 at the cut.
- **Flipping the old blocks:** fixes only the mirrored ones, and the cap
  stays inside the solid.

## Consequences

- `--cap-tips` and `shore cap-tip` work again.
  - Capped C-grids keep their 4-block wrap (a sharp TE is still rounded).
  - Their stations move to the tips.
  - The cap blocks are new: 5 per tip, renamed.
- Projects saved with caps before ADR 0059 do not replay to their saved
  grids. Load them without replay (`--no-replay`).
- **Left open:** rounded or non-planar tips, and a conformal (O-O) tip
  closure. The flare does not guarantee the overlap that an Xall assembly
  needs at every tip; check the hole cutting.

## Verification

`tests/spec/test_tip_cap.py` (main job), frozen before the implementation.
Unless marked, the cases run on a sharp and a blunt NACA 0012 wing.
- C1: five blocks per tip, every cell's scaled Jacobian above 0.2, left
  alone by the march.
- C2: 16 SHARED faces per tip, matched through their axis maps.
- C3: `k_lo` is the only WALL face.
- C4: the wall is planar, on the STL's tip, and holds the tip-station
  nodes.
- C5: the wall's area equals the section's.
- C6: the core's layers follow the wrap's schedule along the normal.
- C7: the flare is 30° at the edge, or 0° when set to 0, leaning out of
  the section.
- C8: a swept, tapered, cambered wing.
- C9: `Project.build_tip_cap` puts the wall on the tip face.
- C10: `shore mesh --cap-tips` and `shore cap-tip` write valid files.
- C11: save and load, bit for bit.
- C12: `cap_te_frac` is deprecated and changes nothing.

At freeze: 0 passed, 19 xfailed.
