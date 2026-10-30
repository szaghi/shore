# 0047 — Wall entities projected from across the block (rejected)

- **Status:** Rejected
- **Date:** 2026-10-29

## Context

ADR 0046 option 2: initialise each wall entity by projecting its partner
across the block (the same entity on the opposite face) onto the geometry,
instead of the TFI of the projected edges. The aims were:
- edges between vertices on a body that wrap round it instead of cutting
  through it (ADR 0046's finding 3);
- grid lines leaving the wall along its normal, since a closest-point
  projection is along the normal.

**The user asked for all four ADR 0046 options to be done** (night of
2026-10-29). This one was prototyped as an opt-in premesh mode, so that ADR
0042's frozen premesh results would not change.

## Investigation (2026-10-29, scratchpad `night47/option1-rejected.patch`, `wing46/`; T1)

**Prototype:** `premesh(walls="normal")`.
- Each wall vertex, edge and face is projected from its partner across the
  block; entities on several walls with different partners are left alone.
- Everything else is rebuilt from the new wall entities.
- ADR 0042–0045's specs were unchanged with the default mode: 56 passed.

**The wing on a plane** of ADR 0046, route B, without the root cap:

| Variant | Min scaled Jacobian | Wall angle median / p95 | Wall spacing ratio |
|---|---|---|---|
| `tfi` (ADR 0042) | 0.44 | 35.9° / 53.3° | 1.2 |
| `normal`, plain projection | −1 | 21.5° / 74° | 10¹¹ (coincident nodes) |
| `normal`, nodes spread evenly along the projected lines | 0.44 | 35.8° / 53.4° | 1.7 |

- **Plain projection from outside crowds nodes onto the sharp TE and the
  LE.** Many outer nodes have the same closest point, so cells collapse.
- **Spreading the nodes again** (resampling the projected lines evenly)
  undoes the normality: the angles return to the TFI result.
- **Wall nodes evenly spaced and grid lines along the normals cannot both
  come from where the wall nodes are placed.** One must come from moving
  the interior: smoothing (ADR 0046 option 1) or marching (option 3).

**Finding 3 was wrong** (see ADR 0046's correction): vertices moved onto
the wing do not fold edges under `tfi`.
- The `normal` mode also overrides moved vertices with their partners'
  projections, which defeats ADR 0045's vertex moves.

## Decision

**Rejected.** No code change; the prototype patch is kept in the session
scratchpad. ADR 0046's option 2 is closed, and option 1 (smoothing after
premesh) is the next step.

## Consequences

- Premesh stays as ADR 0042 defines it.
- Wall quality on wings is a smoothing (or marching) problem: ADR 0048 on.
