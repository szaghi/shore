# 0019 — Deprecate the marcher's normal blending (`blend_normals_k`)

- **Status:** Proposed
- **Date:** 2026-10-26

## Context

`blend_normals_k` (default 0) blends the body's surface normal into the
legacy marchers' direction over the first *K* layers, at full weight at the
wall. ADR 0015 restored the old vertex normals for it, because the
crease-aware normals (ADR 0014) jump at sharp edges and fold the first
layer. ADR 0015 left open a direction that is both correct on faces and
smooth across edges.

The goal was to find such a normal source. The investigation instead found
that blending itself has no demonstrated benefit. No commit records why it
was introduced; the tests only guard that it does not damage the
cubed-sphere seam.

## Investigation (2026-10-26, prototypes in the session scratchpad; T1)

Five normal sources were tried:
- the old vertex normals;
- the crease-aware normals;
- an *edge band*: crease-aware normals turning to the edge bisector within
  0.05 or 0.15 of a sharp edge;
- *grid smoothing*: crease-aware normals Laplacian-smoothed over the
  layer's nodes, 5 or 20 sweeps.

Each was blended over 5 layers in cubed-sphere, O-grid and C-grid marches.
The direction error is the angle between the first layer's direction and
the exact surface normal, at nodes more than 0.2 from any sharp edge.

| Case | No blend | Old normals | Crease-aware | Edge band 0.15 | Grid smoothing 5 |
|---|---|---|---|---|---|
| cylinder, O-grid / cubed | 10° / 12° | 41° / 44° | 10° / fold | 10° / 12° | 15° / 17° |
| box, O-grid / cubed | 15° / 15° | 55° / 54° | fold / fold | fold / 15° | 34° / 33° |
| cylinder, coarse 9 × 16 | 23° | 46° | fold | 18° (min sJ 0.37 vs 0.51) | 56° |
| box, coarse 9 × 16 | 42° | 53° | fold | 18° (min sJ **0.03**) | 52° |
| concave "flower", smooth | 15-20° | the same or 2° better | the same | the same | 25-28° |
| NACA 0012 C-grid | marches | **fold** | **fold** | **fold** | **fold** |

- **With the marcher's Laplacian smoothing off, the unblended directions
  are already exact off the edges** (0.0°) and never fold. The 10-15°
  above comes from that smoothing, which no normal source can touch.
- **On a C-grid, blending folds with every source,** the old one included.
  The blend gives the wake-cut nodes, which are off the body, the normal of
  the nearest trailing-edge point. That is a defect of the blend itself,
  not of the normals.
- **Blending improved the direction in one setting only:** coarse grids
  with the edge band, and there it brings the box to the edge of folding.
  It is never both better and safe. Across edges it smooths neighbouring
  directions (the old normals: 8-28° between neighbours, against 49-76°
  unblended) only by being wrong elsewhere.

## Decision

User decision of 2026-10-26: **deprecate `blend_normals_k`.**
- `blend_normals_k > 0` emits a `DeprecationWarning` naming the option and
  this ADR, and keeps today's behaviour (ADR 0015's normals) for one
  release; it is removed in a later one.
- `shore mesh` prints a deprecation notice on the console for
  `--blend-normals-k` and the `[march] blend_normals_k` config key, because
  Python hides `DeprecationWarning` from library code by default.
- On C-grid topologies (`cgrid`, and `ch`, which marches its wrap the same
  way), `blend_normals_k > 0` raises `ParameterError` before marching: it
  cannot produce a grid there.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_blend_deprecation.py` runs in the main CI job.

**Verdict.** It must report 8 passed, 0 skipped, 0 xfailed. At freeze it
reports 1 passed and 7 xfailed:
- B5 passes: it is a regression guard (no warning without blending);
- the 7 pending cases fail for the missing behaviour: no warning, no CLI
  notice, and on the C-grid the fold itself (`MarchingError` at k = 1)
  instead of the required `ParameterError`.

| Case | Required |
|---|---|
| B1 (× 2) | `Mesh.march` with `blend_normals_k > 0` warns (cubed sphere, O-grid) |
| B2 | `Project.march` and a journal replay warn; the replay is bit-identical |
| B3, B4 | `--blend-normals-k` and the config key print a notice; exit 0 |
| B5 | no warning at 0 (regression guard) |
| B6 | C-grid with blending raises `ParameterError`; the mesh is unchanged |
| B7 | the C-grid CLI run exits 1 with the reason |

**Out of scope:**
- removing the option (a later release);
- CH beyond sharing the C-grid code path (a unit test).
