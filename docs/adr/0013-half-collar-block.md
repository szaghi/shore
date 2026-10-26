# 0013 — The half-collar block on a symmetry plane

- **Status:** Accepted
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

## Review decisions (2026-10-26, user)

1. **H2 plus H3 is enough.** The self-mirror reference binds any future
   implementation to the symmetric solution, and H3's right angle is the
   independent check.
2. **SYMMETRY faces split like FREE faces.** A sub-face of a symmetry face
   is still a symmetry face with patch code 2, and there is nothing to
   pair. `split` copies face BCs generically, so it needs no change.
3. **The cc.par code table is fixed in this change:**
   `docs/reference/cc-par.md` now says that patches carry the positive
   natural-BC index, and that the negative codes are internal.

## Acceptance (2026-10-26): half-collar spec met

- `shore.collar.build_collar(..., plane=)` builds the half-collar as
  decided.
- `extrude_collar(wall, steps_k, junction)` is the public volume stage, and
  it checks the wall, the junction row and the steps.
- `FaceBC.SYMMETRY` is new. `write_cc_par` emits it as code 2 with family 0,
  and the adjacency sidecar round-trips it by name.
- `BuildCollar` / `Project.build_collar` gain `plane`, journaled as plain
  float tuples.
- The oblique-plane warning threshold is `CollarParams.oblique_deg` = 1°.

The closed collar is bit-identical to the committed version, checked on the
tilted-plane and sphere junctions; the ADR 0010 spec still reports 26
passed.

**Verdict for `tests/spec/test_half_collar.py`: 23 passed, 0 skipped,
0 xfailed** (T2; `scripts/check_spec.py` reports all files ok). All 11
`@PENDING` lines are removed. Full suite: 1392 passed, 0 xfailed.

Measured (T1, the spec's settings; the same numbers as the prototype, which
uses the same construction):

| Half model | Folded | Min scaled Jacobian | First cell | Outer layer | j-line vs plane normal |
|---|---|---|---|---|---|
| plane × cylinder | 0 | 0.491 | 0.992-1.000 | 1.000 | 0.014° |
| tilted plane × cylinder | 0 | 0.081 | 0.984-1.000 | 0.998-1.000 | 0.065° |
| cylinder × cylinder | 0 | 0.532 | 0.985-1.000 | 0.998-1.002 | 0.119° |
| sphere × offset cylinder | 0 | 0.660 | 0.990-1.001 | 0.997-1.001 | 0.048° |
| STL half box + half post (command) | 0 | 0.512 | 0.988-1.000 | | |

**Unit tests** (`tests/test_collar.py`): `extrude_collar` argument errors;
the end-angle measure on a half circle; a malformed plane in the command.

Half-model collars are complete: the open curve (ADR 0011), the marcher's
ends (ADR 0012) and the block (this ADR). The STL vertex-normal symmetry
(ADR 0012, review decision 1) is still open.
