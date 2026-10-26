# 0009 — Intersection loops of two bodies, for collar grids

- **Status:** Proposed
- **Date:** 2026-10-26

## Context

A collar grid starts from the curve where two components meet, for
example a wing root on a fuselage or an appendage on a hull. ADR 0008
marches a surface grid from a given curve. This ADR provides that curve.

## Decision

`shore.surface.intersect.intersection_loop(body_a, body_b, seed, n, *,
entity_a=BODY, entity_b=BODY) -> (n, 3)`, all user decisions of
2026-10-26:

- **Protocol-only tracing.** Alternating projections from the seed onto A
  and B converge onto the curve. The tracer then takes predictor steps
  along `t = n_b × n_a`, each followed by an alternating-projection
  corrector, until the loop closes. Only `Geometry.project` (with
  normals) is used, so STL, CAD and analytic bodies work alike.
- **User seed.** The caller gives a point near the junction, as a GUI pick
  would. When the bodies meet in several loops, the seed selects one.
- **Resampled loop.** The loop has `n` nodes, equally spaced in arc length
  and polished back onto both bodies. It is ordered along `n_b × n_a`, so
  that `march_surface(A, loop)` and `march_surface(B, loop[::-1])` march
  the two halves of the collar away from the other body.

## Alternatives considered

- **Backend-native intersection** (OCCT `BRepAlgoAPI_Section` for CAD,
  trimesh mesh–mesh intersection for STL). Rejected: two code paths, and
  exact checks only in the CAD job.
- **Automatic search for every loop.** Rejected for now: it needs a global
  search whose completeness is hard to guarantee.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_intersection.py` runs in the main CI job.

**Verdict.** It must report 13 passed, 0 skipped, 0 xfailed. At freeze it
reports 0 passed, 13 xfailed, every case failing only with
`ModuleNotFoundError`.

**Rules.** As in ADR 0005: per-case `@PENDING` markers, with the hashes in
`tests/spec/SHA256SUMS`. The bodies are exact analytic surfaces defined in
the spec, and every loop case uses n = 64.

| Case | Required | Prototype |
|---|---|---|
| X1-X4: plane × cylinder (circle), tilted plane × cylinder (ellipse), cylinder × cylinder and sphere × offset cylinder (two loops each; the seed's is required) | on both bodies <= 1e-9 × diagonal; chord lengths max/min <= 1.01; length within 0.5 % of the exact curve; only the seed's loop; `loop[0]` at the seed | 1.2e-15; 1.0003 |
| X5: each pair with `march_surface` | one layer on A from the loop lies outside B; on B, from the reversed loop, outside A | |
| X6: interface | deterministic; `ParameterError` for n < 3, a bad seed, and no intersection near the seed | |

**Harness validation (T1).** The scratch prototype was run behind the
frozen interface, with the accepted `march_surface` running the X5
compositions.
- It passes all 9 geometric cases: the four loops, the four collar-side
  compositions, and determinism.
- It fails exactly the 4 input-error cases, because it has no input
  validation.

So the bounds are attainable, the orientation composes with the marcher,
and the error cases discriminate.

**For review (your call before implementation starts):**
1. Only uniform arc-length spacing is required. Collars usually want
   clustering, for example near a trailing edge. Add a spacing law now, or
   in the collar-block spec?
2. Grazing (tangential) intersections, where the bodies' normals are
   parallel, are out of scope. Alternating projections slow down there, and
   the tracing direction `n_b × n_a` vanishes. Should the kernel detect
   this case and raise?
3. Open intersection curves, which end on a boundary such as a symmetry
   plane, are out of scope. Wing–body collars on half models need them.

## Acceptance (2026-10-26): intersection-loop spec met

`src/shore/surface/intersect.py` implements the decision. Its settings
(`TraceParams`, lengths relative to the model diagonal):
- predictor step 4e-3;
- traced points converged to 1e-10;
- output nodes polished to 1e-13, all nodes at once.

It also checks what the prototype did not:
- the arguments (integer `n >= 3`, a finite 3-vector seed);
- that the seed lands on both bodies ("no intersection near the seed");
- a trace that does not close;
- parallel normals, where there is no tracing direction (the grazing case
  is out of scope, but it fails cleanly).

**Verdict for `tests/spec/test_intersection.py`: 13 passed, 0 skipped,
0 xfailed** (T2; `scripts/check_spec.py` reports all files ok). All 6
`@PENDING` lines are removed; the parametrised cases share them. Full
suite: 1285 passed, 0 xfailed.

Measured (T1):

| Pair | On both bodies | Chords max/min | Length vs exact | Branch, start |
|---|---|---|---|---|
| plane × cylinder | 0 | 1.00004 | -0.040 % | ok, node 0 |
| tilted plane × cylinder | 1.5e-13 | 1.00026 | -0.041 % | ok, node 0 |
| cylinder × cylinder | 5.9e-14 | 1.00153 | -0.044 % | ok, node 0 |
| sphere × offset cylinder | 6.7e-13 | 1.00151 | -0.040 % | ok, node 0 |
| required | <= 6.9e-9 | <= 1.01 | within 0.5 % | |

The length shortfall is the inscribed 64-gon's own deficit, (π/64)²/6 ≈
0.040 %, not a tracing error.

**Performance.** Tracing with the output tolerance cost 6.3 s per loop on
two STL meshes, from about 12,600 single-point trimesh queries. Converging
the traced points only to 1e-10, with a 2.7× longer step, brings that to
1.8 s. The cost is chord uniformity, 1.0002 → 1.0015.

**Unit tests** (`tests/test_intersection.py`):
- a faceted icosphere pierced by a faceted cylinder: the loop lies on both
  meshes (< 1e-8), only the seed's (upper) loop is traced, and the chords
  are even to 5 %;
- the collar half marched on the sphere leaves the cylinder.

The review questions above are still open. Next for collars: the collar
block (two-sided surface grids extruded into a volume block).

## Review decisions and re-freeze (2026-10-26, user)

1. **Spacing laws along the loop move to the collar-block spec.** Where to
   cluster depends on the block topology; the loop can be redistributed
   along its arc length there.
2. **Re-freeze:** new case **X7** requires a sphere tangent to a plane to
   raise `ParameterError` within 5 s. The case discriminates: the
   unguarded prototype silently returns a meaningless "loop". In the
   kernel, the "no intersection near the seed" guard fires first, because
   alternating projections converge too slowly at a tangent contact to
   reach the on-curve tolerance. The parallel-normals guard covers
   tangency met during tracing.
3. **Open intersection curves move to the collar-block spec**, with the
   marcher's end constraint (ADR 0008, decision 2).

Verdict after re-freeze: 14 passed, 0 skipped, 0 xfailed.
