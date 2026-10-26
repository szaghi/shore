# 0013 — The half-collar block on a symmetry plane

- **Status:** Proposed
- **Date:** 2026-10-26

## Context

ADR 0011 traces the open junction curve of a half model and ADR 0012
marches its surface halves with their ends on the symmetry plane. This ADR
extrudes them into the half-collar volume block, the last piece for
half-model collars.

## Decision

User decisions of 2026-10-26:

- **`build_collar(..., plane=None)`.** With a plane, it traces
  `intersection_curve`, marches both halves with `march_surface(plane=)`,
  and returns an open block (`periodic_j` false) whose two j-faces lie on
  the plane. `Project.build_collar` gains the same option.
- **A new `FaceBC.SYMMETRY`.** The j-faces carry it into
  `.adjacency.json`, and `write_cc_par` emits it as the natural-BC code
  **2**. In cc.par a natural BC is written as its positive index *p* in
  `[1, 19]`; the overset preprocessor maps it to its internal code
  `rcc(p)`, -2 for symmetry (`BC_NATURAL_SIMMETRY` in
  `overset_exploded_block_object.F90`). This is the same convention as the
  existing `9` for extrapolation. The table in `docs/reference/cc-par.md`
  lists the internal negative codes and does not say so; that needs a note.
- **Mirror match and faces on the plane** are required on mirror-symmetric
  half models.

Method: the half-collar is the half of the ADR 0010 volume extruded from
the **mirrored closed wall**. The wall is reflected across the plane into a
j-periodic wall with `2 (n_loop - 1)` stations, extruded by the accepted
volume stage, and the first `n_loop` stations are kept, with the j-faces
put on the plane exactly. This reuses the accepted hyperbolic kernel
(ADR 0005) unchanged, where a symmetry condition inside it would have
changed an accepted kernel. The cost is extruding twice the stations.

The volume stage becomes public as `shore.collar.extrude_collar(wall,
steps_k, junction)`, so the spec can build the mirror reference
independently of `build_collar`.

When the curve meets the plane obliquely (end tangent more than 1° off the
plane normal), `build_collar` warns (ADR 0011, review decision 2). The
mirrored wall is then not a valid surface, and the block's quality is not
guaranteed.

## Why not compare with a whole-body collar

A whole-body loop seeded on the plane is resampled from a different traced
polyline, and its nodes differ from the half curve's by 2e-6 to 2.7e-5. An
exact comparison would be dominated by that resampling noise. The mirror
reference is therefore built from the half-collar's own wall, through the
public `extrude_collar`.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_half_collar.py` runs in the main CI job.

**Verdict.** It must report 23 passed, 0 skipped, 0 xfailed. At freeze it
reports 0 passed, 23 xfailed, every case failing on the missing interface
(`plane` argument, `extrude_collar`, the command option).

**Cases.** `n_loop` = 33 on the ADR 0011 half models; 15 surface layers
`0.01 × 1.12^j` on each side, 30 wall-normal layers `1e-4 × 1.2^k`.

| Case | Required | Prototype |
|---|---|---|
| H1 wall (× 4) | the composition of the accepted kernels (<= 1e-12) | exact |
| H2 mirror (× 4) | equal to `extrude_collar` of the mirrored wall, first `n_loop` stations (<= 1e-9 × thickness) | by construction |
| H3 validity (× 4) | j-faces on the plane; y >= 0; j-lines within 0.5° of the plane normal; 0 folds, min scaled Jacobian >= 0.05, first cell 0.95-1.05, outer layer within 1 % | 0; 0.014-0.119°; 0.081-0.660; 0.984-1.001; 0.997-1.002 |
| H4 topology (× 4) | shape, `periodic_j` false, j-faces SYMMETRY, k_lo WALL, others FREE | |
| H5 oblique plane | warns; j-faces on the plane | warns; 2 folded cells (not required) |
| H6 interface (2) | `extrude_collar` reproduces the accepted whole collar bit for bit; a bad plane raises | |
| H7 output (2) | SYMMETRY in the sidecar; cc.par code 2, family 0 | |
| C1, C2 command (STL half box + half post) | one valid block with SYMMETRY j-faces; journal, replay | 0 folds, min sJ 0.512, faces exact, 0.9 s |

**Harness validation (T1).** The prototype was run behind the frozen
interface, standing in for `build_collar(plane=)` and `extrude_collar`. It
passes 15 of 23 cases and fails exactly the 8 that need what it does not
have: `FaceBC.SYMMETRY` (H4 × 4, H7 × 2) and the command option (C1, C2).

**Out of scope:**
- two planes;
- the quality of half models that are not mirror-symmetric at the plane;
- the STL vertex-normal symmetry (ADR 0012, a separate decision);
- run time.

**For review (your call before implementation starts):**
1. H2 compares the block with `extrude_collar` of its own mirrored wall.
   For this implementation that is true by construction; it binds any
   future one (for example a symmetry condition inside the hyperbolic
   kernel) to the same result. H3's right angle is the independent check.
   Enough, or also compare with a whole-body collar at a loose bound (the
   loops differ by up to 2.7e-5)?
2. SYMMETRY faces are not connections. `split` keeps them as they are, and
   `balance` may cut through them like FREE faces. Is that right for Xall?
3. Fix the misleading cc.par code table in `docs/reference/cc-par.md` in
   the implementation commit?
