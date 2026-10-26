# 0011 — Open intersection curves ending on a symmetry plane

- **Status:** Accepted
- **Date:** 2026-10-26

## Context

Wing–body collars are usually built on half models: the bodies end at the
symmetry plane, and the junction is an open curve from the plane back to
the plane. `intersection_loop` (ADR 0009) traces closed loops only. The
ADR 0009 review moved open curves here.

## Decision

`shore.surface.intersect.intersection_curve(body_a, body_b, seed, n, *,
plane, entity_a=BODY, entity_b=BODY) -> (n, 3)`. User decisions of
2026-10-26:

- **Symmetry plane only.** The curve ends on a plane given as
  `(point, normal)`. The seed's side of the plane is the side that is
  kept, so the normal's sign does not matter.
- **Curve kernel only.** The surface marcher's end constraint (layer ends
  held on the plane) and the half-collar block come in later specs.
- **A new function.** `intersection_loop` and its frozen spec stay as they
  are.

Method (prototype in the session scratchpad):
1. Alternating projections take the seed onto the curve, as in ADR 0009.
2. The ADR 0009 tracer runs from there in both directions along
   `±(n_b × n_a)`, until a step reaches the plane **within the on-curve
   tolerance**.
3. Newton steps along the curve tangent, each polished back onto both
   bodies, put the end exactly on the plane.
4. The two halves are joined, ordered along `n_b × n_a`, and resampled by
   arc length. The interior nodes are polished, and the two ends are the
   Newton points.

The "within tolerance" test in step 2 is necessary. On a half model the
bodies stop at the plane, so a projection never goes past it: the trace
lands *on* the plane and stalls there. A strict "has crossed" test never
fires. On the STL half model, that variant ran into its step limit (all
cases: 278 s), while the analytic half bodies did not show the problem.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_open_curve.py` runs in the main CI job.

**Verdict.** It must report 20 passed, 0 skipped, 0 xfailed. At freeze it
reports 0 passed, 20 xfailed, every case failing only with the
`ImportError` of the missing `intersection_curve`.

**Cases.** n = 33 throughout; the analytic bodies of ADR 0009's spec, plus
`Half`, the y >= 0 part of a body symmetric about y = 0, with exact
projection onto it.

| Case | Required | Prototype |
|---|---|---|
| O1 half models (× 4: the ADR 0009 pairs cut at y = 0) | on both bodies <= 1e-9 × diagonal; ends on the plane; every node at y >= 0; chords max/min <= 1.01; length within 0.5 % of the exact half curve; the seed's branch | 7e-13; 0; 1.0003; -0.044 % |
| O2 whole bodies | the seed's side is kept (seed at y < 0); a plane that is not a symmetry plane (x = 0.5, arc 2π/3) | -0.040 %, -0.018 % |
| O3 STL half model (half box + half post, 96 sections) | on both meshes <= 1e-9 × diagonal; ends at (±0.5, 0, 0); within 5 s (re-freeze) | 2e-16; exact; 0.13 s |
| O4 orientation (× 4 half models) | one `march_surface` layer on A outside B; on B, from the reversed curve, outside A | all |
| O5 interface (9) | deterministic; `ParameterError` for n < 3, a bad seed (× 2), a bad plane (× 2), a seed on the plane, no intersection, a closed curve that never meets the plane (within 5 s) | |

**Harness validation (T1).** The prototype was run behind the frozen
interface. It passes 15 of the 20 cases and fails exactly the 5 input
checks it does not make: n < 3, both bad seeds, the 2-D plane, and the
seed on the plane. The zero-normal case passes only by accident (the
prototype raises for the wrong reason); the implementation must check it.
With the strict-crossing test, the prototype also fails O3, so that case
discriminates.

**Measured for the next spec.** On the half models, ends marched freely by
`march_surface` (2 steps of 0.02) leave the plane by 0.9-1.2e-3, about 5 %
of a step. The half-collar needs the marcher to hold them on the plane.

**Out of scope:**
- ends on anything other than a plane (a third body, the open edges of a
  body);
- curves tangent to the plane at their ends;
- curves that meet the plane more than twice on the seed's side;
- the marcher's end constraint;
- spacing laws other than uniform;
- run time.

**For review (your call before implementation starts):**
1. One plane only. A wing–body half model with a pylon or a tail can need
   a curve that ends on two different planes. Leave that out?
2. The curve is not required to meet the plane at right angles. On a true
   symmetry plane it does, by symmetry. Should the kernel check it and
   warn or raise when the bodies are not mirror-symmetric there, where a
   collar's i-face will not be planar?
3. The time bound applies only to the closed-curve error case (5 s). The
   strict-crossing trap shows that a stall costs minutes. Add a run-time
   bound to O3?

## Review decisions and re-freeze (2026-10-26, user)

1. **One plane only.** Curves ending on two different planes (a pylon or
   tail junction) are deferred until a real case needs them.
2. **The right-angle check belongs to the half-collar spec.** That is where
   a curve meeting the plane obliquely matters: the bodies are not
   mirror-symmetric there, and the collar's i-face will not be planar.
   The curve kernel stays as frozen.
3. **Re-freeze:** O3 also requires the curve within 5 s (prototype 0.13 s).
   A trace that stalls at the plane now fails fast and explicitly, not
   after the step limit (minutes).

Verdict after re-freeze: still 0 passed, 20 xfailed; the prototype behind
the interface still passes 15 of 20, failing the same 5 input checks.

## Acceptance (2026-10-26): open-intersection-curve spec met

`intersection_curve` in `src/shore/surface/intersect.py` implements the
decision, with two new `_Tracer` methods: `trace_to_plane` (reaches or
crosses within the on-curve tolerance, detects a curve that closes first)
and `end_on_plane` (Newton along the tangent). Beyond the prototype it
checks the arguments (n, seed, plane, a seed on the plane, a junction point
off the seed's side) and raises `ParameterError` where the curve meets the
plane tangentially.

One defect found while measuring, and fixed: when the trace lands on the
plane within tolerance (always on half bodies), the end point had been
polished only to the trace tolerance (1e-10). It sat 7.3e-10 off the half
sphere, inside the spec's 6.9e-9 bound but 1,000 times the prototype's
error. Now it is polished to the output tolerance first. A unit test pins
the end points below 1e-11; it fails without the fix.

**Verdict for `tests/spec/test_open_curve.py`: 20 passed, 0 skipped,
0 xfailed** (T2; `scripts/check_spec.py` reports all files ok). All 12
`@PENDING` lines are removed. Full suite: 1347 passed, 0 xfailed.

Measured (T1, n = 33, half models cut at y = 0):

| Pair | On A | On B | Ends off the plane | Chords max/min | Length vs exact |
|---|---|---|---|---|---|
| plane × cylinder | 0 | 0 | 0 | 1.00001 | -0.040 % |
| tilted plane × cylinder | 1.4e-13 | 1.2e-16 | 0 | 1.00026 | -0.041 % |
| cylinder × cylinder | 1.8e-13 | 7.9e-17 | 0 | 1.00016 | -0.044 % |
| sphere × offset cylinder | 6.7e-13 | 2.8e-17 | 0 | 1.00024 | -0.040 % |
| STL half box + half post | 2.3e-16 | 2.2e-16 | 6e-17 | 1.0005 | |

Each curve takes 0.01-0.02 s (analytic) and 0.13 s (STL).

**Unit tests** (`tests/test_intersection.py`): a faceted icosphere and
cylinder cut by a plane; a plane tangent to the curve raises; the end
points' accuracy.

Next for half-model collars: the surface marcher's end constraint, then
the half-collar block (with the right-angle check, review decision 2).
