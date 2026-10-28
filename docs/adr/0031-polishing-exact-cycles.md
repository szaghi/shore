# 0031 — Polishing: stop at exact cycles

- **Status:** Accepted
- **Date:** 2026-10-28

## Context

An open follow-up since ADR 0021: tracing a junction loop on sliver STL
triangles can be slow. The case is a NACA 0012 wing extruded from
1000-point sections with 2 stations, through a plate: 48 s. The fix must
keep every other result bit-identical. User decision of 2026-10-28: this
follow-up first, while the collar work pauses.

## Investigation (2026-10-28; T1)

**Where the time goes.**
- Profile: 98 % of the 50 s is in trimesh's `closest_point`, nearly all of
  it in one call of the tracer's `polish`, within `_resample`.
- That call runs its full 200 iterations. 90 of its 93 points converge by
  iteration 4. One is trapped in an exact 2-cycle by projection noise
  (3.4e-11 per step, against a 5.7e-13 tolerance). Another repeats with an
  exact period 5 in its last bits.
- Every iteration still projects all 93 points. The extrusion's
  triangles are slivers, so the broad phase returns nearly every face
  (68 million candidate ids).

**Detecting the cycle.**
- A whole-array 2-cycle test never fires: the period-5 point keeps the
  array changing.
- An exact repeat of the whole array at any period up to 12 (here 10)
  does fire. The rest of the loop is then known: the first step within
  the cycle that meets the tolerance, otherwise the iterate the remaining
  iterations end on. Returning it gives the same result bit for bit.

| Loop (NACA 0012, 2 stations) | before | after | result |
|---|---|---|---|
| 1000-point sections, blunt TE | 48.4 s | 8.2 s | bit-identical |
| 1000-point sections, sharp TE | 5.5 s | 5.7 s | bit-identical |
| 200-point sections, blunt TE | 1.4 s | 1.4 s | bit-identical |

The 5-8 s that remain are trimesh's broad phase on slivers, a cost of
every projection.

## Decision

- `_Tracer.polish` keeps its last `TraceParams.polish_period` (12)
  iterates.
- When a new iterate equals one of them bit for bit (`_period`),
  `_cycle_end` returns what the remaining iterations would return.
- `polish_period = 0` turns the detection off.

## Consequences

- **Faster:** junction loops on sliver STLs no longer spend minutes in
  cycling polish.
- **Bit-identical:** every result is unchanged (full suite and all golden
  files).
- **Still open:** the broad phase on slivers. A better spatial index for
  `StlGeometry.project` would need its own measurements.

## Verification

- **Unit tests** (`tests/test_intersection.py`):
  - period detection;
  - the cycle end, including a step inside the cycle that meets the
    tolerance;
  - a polish over bodies whose projections cycle with period 3, never
    converging: same result bit for bit with and without detection.
- **Full suite:** passes unchanged (the golden loops and marches, within
  their tolerances).
