# 0007 — Surface grids: TFI and elliptic smoothing in (u, v) per face

- **Status:** Proposed
- **Date:** 2026-10-26

## Context

SHORE builds every k = 0 surface by casting rays from an interior anchor
(`surface/remesh.py`, `surface/cap.py`, `surface/equator.py`), so it only
handles star-shaped bodies. CAD input (ADR 0002) and interactive blocking
(roadmap phases 6 and 7) need the opposite approach, as in ICEM and
Pointwise: the user or a topology fixes a block face's four edge
distributions on CAD curves, and the mesher fills the face on its CAD
surface.

The standard method is algebraic transfinite interpolation (TFI) followed
by elliptic smoothing. Plain TFI is exact on tensor-product data but folds
on strongly curved boundaries. Elliptic (Winslow) smoothing repairs folds,
but on its own it also washes out the clustering that CFD walls need.

ADR 0002 measured tier-1 projection at 2.5e4–1.3e5 points/s. Smoothing on
the surface with re-projection every sweep (1e5 nodes × 100 sweeps) would
take about 100 s per face. It decided to smooth in (u, v) per face and to
evaluate on the surface only at the end.

## Decision

1. **Kernel** `shore.surface.tfi.surface_grid(geometry, surface, boundary,
   smoothing="elliptic") -> SurfaceGrid(points, params)`:
   - the boundary rows of `boundary (ni, nj, 3)` are projected once onto the
     SURFACE entity to get their (u, v);
   - TFI of those parameters;
   - with `"elliptic"`, Winslow smoothing in (u, v) using the surface metric
     (the Christoffel terms of the parametrisation), with control functions
     of Thomas–Middlecoff type, taken from the boundary distributions, so
     that clustering is kept;
   - one evaluation on the surface at the end.

   The metric comes from a bounded number of geometry queries per node (for
   example, a parameter lattice sampled once), never per sweep.
2. **Command** `Project.surface_grid(block, face, *, surface,
   smoothing="elliptic")` rebuilds one block face from its boundary on a
   SURFACE entity of the loaded geometry. It is journaled and undoable like
   `Snap`.
3. Out of this first version:
   - faces whose (u, v) range crosses a periodic seam or a singular point;
   - block faces spanning several CAD faces;
   - boundary-orthogonality control (Sorenson / Hilgenstock).

## Alternatives considered

- **TFI on physical points, then projection.** Rejected: chords cut through
  curved surfaces. On the spec's cylinder, physical TFI is off by 0.046
  where (u, v) TFI is exact to 1e-15.
- **Smoothing on the surface with re-projection every sweep.** Rejected:
  ADR 0002 throughput.
- **Plain Winslow, without control functions.** Rejected: on the clustered
  rectangle the first cell grows 1.7–13× and the grid moves up to 33 local
  spacings.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

There are two modules:

| Module | Job | Verdict for acceptance | At freeze |
|---|---|---|---|
| `tests/spec/test_surface_tfi.py` | main (no CAD extra) | 23 passed, 0 skipped, 0 xfailed | 0 passed, 23 xfailed |
| `tests/spec/test_surface_tfi_cad.py` | CAD | 4 passed, 0 skipped, 0 xfailed | 0 passed, 4 xfailed |

The rules are those of ADR 0005: per-case `@PENDING`, and hashes in
`tests/spec/SHA256SUMS`. The kernel cases use exact analytic surfaces
written in numpy in the spec itself: the identity plane, a plane and a
cylinder with a distorted parametrisation, and a sphere patch. Every
reference is therefore exact and the work count is observable.

| Case | Required | Evidence |
|---|---|---|
| S1 TFI, tensor data | exact (1e-12 × diag); done in (u, v) | (u, v) TFI 1.1e-15; physical TFI 0.046 |
| S2 invariants, 4 surfaces × 2 modes | boundary = input, nodes on the surface, `evaluate(params) = points` (1e-12 × diag) | surfaces invert to ≤ 7.5e-15 |
| S3 metric, uniform physical boundaries | elliptic ≤ 0.02 h from the uniform grid; TFI ≥ 0.3 h | prototype (u, v) Winslow 0.0007 h (plane), 0.0008 h (cylinder); TFI 0.65 h, 1.03 h |
| S4 clustering, tanh tensor rectangle (ds 8.0e-4) | first cell 0.95–1.05 ds; ≤ 0.05 h | Thomas–Middlecoff 1.000, 0.0000 h; plain Winslow 1.66–13.3, 33 h |
| S5 fold repair, S-shaped overhang 81×61 | TFI folds; elliptic 0 folded, min corner angle ≥ 2° | TFI 29 folded; prototype (plain and Thomas–Middlecoff) 0 folded, 5.4° |
| S6 sphere patch | 0 folded | validity only |
| S7 work | ≤ 10 geometry points per node | |
| S8 interface | default `"elliptic"`, deterministic, `ParameterError` cases | |
| CAD: command on a STEP cylinder, annulus-sector inner face, clustered | reproduces the exact tensor face (≤ 0.05 h), on the CAD face (≤ 1e-9 × diag), rest of block bit-identical; journal, replay, undo; rejects unknown block or face and STL | reference face on the CAD face to 2.2e-16 |

**Harness validation (T1).**
- The analytic surfaces invert to round-off, and an offset along the
  normal reads back exactly.
- The fold metric reads +1 on valid grids and −1 on a mirrored grid.
- Each case's own "discriminates" check holds on the harness's (u, v) TFI.
- Every pending case fails only on the missing kernel or command:
  `ModuleNotFoundError` for all 23 kernel cases, and the missing
  `Project.surface_grid` for all 4 CAD cases.

**The S5 domain was changed twice while validating.**
1. A first overhang with a near-horizontal zig-zag gave an almost
   self-touching boundary, where no grid can be valid.
2. A second overhang had nodes that were not uniform in arc length. Plain
   Winslow passed it, but Thomas–Middlecoff folded 3 cells, because the
   control functions copy that spacing. No clip on P and Q separates the two
   cases: the overhang's |P| (0.21) lies inside the range the clustering
   case needs (0.15 / 0.25). With arc-length-uniform nodes P = 0, and the
   case is fair.
3. At 41×31 every overhang that folds TFI also folds discrete Winslow in
   the boundary cell row at the concave turn, where it is under-resolved. At
   81×61 the window is clean: a = 0.16 and c = 0.3 fold TFI in 29 cells,
   and elliptic smoothing in none.

**For review (your call before implementation starts):**
1. Seam-crossing and singular-point faces are out of scope, yet many CAD
   faces are full cylinders or spheres. Acceptable for v1, with the command
   rejecting them?
2. Clustered boundaries on a strongly concave boundary (S4 combined with
   S5) are out of scope. The prototype shows Thomas–Middlecoff folding there
   when the spacing is irregular. Should the spec require it, which would
   make limiting the control functions part of the method?
3. S6 (sphere) checks validity only. Should a quality bound be added (for
   example, a minimum angle no worse than TFI's)?
4. Is the work budget of 10 points per node right? The measured need is
   about 1 for the evaluation, plus the metric lattice.
