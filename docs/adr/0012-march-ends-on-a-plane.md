# 0012 — Surface marching with ends on a symmetry plane

- **Status:** Proposed
- **Date:** 2026-10-26

## Context

`intersection_curve` (ADR 0011) gives the open junction curve of a half
model, from the symmetry plane back to it. `march_surface` (ADR 0008)
marches the ends of an open curve freely. On the half models, the free
ends leave the plane by 0.017-0.056 over the 15 collar layers: the
half-collar would have no face on the plane.

## Decision

`march_surface(..., plane=(point, normal))`, user decisions of 2026-10-26:
- **One plane for both ends**, as `intersection_curve` produces. Curves
  ending on two planes stay deferred (ADR 0011, review decision 1).
- **Mirror ghost node.** This is the classical symmetry-plane condition.
  Each end node is solved as an interior node; its missing neighbour is
  the mirror image of its inner neighbour across the plane, both in the
  differences (tangent, Laplacian, contraction sensor) and in the linear
  system (the ghost's increment is the reflected increment of the inner
  node). After the projection onto the body, alternating projections put
  each end on both the body and the plane.
- **Orthogonality and an unchanged interior** are required, on exactly
  mirror-symmetric bodies.

`plane=None` keeps the accepted free-end march.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_march_ends.py` runs in the main CI job.

**Verdict.** It must report 17 passed, 0 skipped, 0 xfailed. At freeze it
reports 0 passed, 17 xfailed, every case failing only with the `TypeError`
of the missing `plane` argument.

**Cases.** Curves from `intersection_curve` with n = 33, and the collar
schedule of ADR 0010: 15 layers `0.01 × 1.12^j`.

| Case | Required | Mirror ghost (prototype) | March, then snap the ends |
|---|---|---|---|
| E1 half models (4 pairs × 2 sides) | nodes on the body, ends on the plane (<= 1e-9 × diagonal); y >= 0 | 7e-13; 0 | |
| | end tangent (second-order, one-sided) within 0.5° of the plane normal | 0.013-0.119° (the same as the whole-body grid) | 0.07-5.0° |
| | equal to the mirrored whole-body periodic march (<= 1e-9 × first step) | 2-3e-14 | 0.01-1.6 |
| E2 oblique plane x = 0.5 (× 2 sides) | ends on the plane, nodes on the body; valid cells | 2e-14; valid | |
| E3 STL half box + half post (× 2) | ends on the plane, nodes on the meshes; y >= 0; valid cells | 5e-16; valid | |
| E4 interface (5) | deterministic; `ParameterError` for a bad plane (× 2), `periodic=True` with a plane, curve ends off the plane | | |

**Harness validation (T1).** The prototype was run behind the frozen
interface. It passes 13 of 17 cases and fails exactly the 4 input checks it
does not make. The "march, then snap" variant fails all 8 E1 cases: the
right angle and the mirror match discriminate.

**Finding: STL normals are not mirror-symmetric at a cut.** `StlGeometry`
interpolates vertex normals barycentrically (ADR 0002, bit-identical to
`Projector`). On the coarse half box, the normal on the flat top face is
tilted by up to 14° (-0.24, 0, 0.97), and the half mesh's normals are not
the mirror image of the whole mesh's. Marched on the half box, the end
tangent drifts to 5.1° and the grid differs from the mirrored whole-mesh
march by 1.5 first steps (post: 0.9°, 0.085). That is the geometry, not
the end condition, so E3 requires only the ends on the plane and valid
cells. The same interpolation affects collars on any coarsely triangulated
flat STL face; it is outside this spec.

**Out of scope:**
- ends on two different planes, or one free end and one on a plane;
- curves that meet the plane obliquely on mirror-symmetric bodies;
- the half-collar block;
- run time.

**For review (your call before implementation starts):**
1. The finding above: should the marcher (or `StlGeometry`) use face
   normals, or snap the end normals into the plane? It is a separate
   decision that affects ADR 0002's bit-identity with `Projector`.
2. The mirror-match bound (1e-9 × the first step) is five orders above the
   prototype's 3e-14. It requires the symmetric solution, not a particular
   floating-point order. Keep it?
