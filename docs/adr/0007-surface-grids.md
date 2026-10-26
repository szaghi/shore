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

## Kernel accepted (2026-10-26)

`src/shore/surface/tfi.py` implements Decision 1. **Verdict for
`tests/spec/test_surface_tfi.py`: 23 passed, 0 skipped, 0 xfailed** (T2,
run independently; `scripts/check_spec.py` reports all files ok). All 13
`@PENDING` lines in the module are removed; the parametrised cases share
them. The CAD half is pending: it needs `Project.surface_grid`.

Method:
- the boundary is projected once (`exact=True`) to get its (u, v);
- TFI in (u, v);
- Winslow in (u, v) with Christoffel terms and linearly blended
  Thomas-Middlecoff control functions (`|P|, |Q| <= 1.5`), solved by Picard
  iteration with one SuperLU factorisation per step, shared by u and v;
- one evaluation at the end, with the boundary rows copied from the input.

The metric comes from a parameter lattice twice the grid's resolution
(`lattice_refine = 2`), evaluated once. Its finite-difference derivatives
are interpolated bilinearly.

Measured margins (T1, `EllipticParams` defaults; the same values for every
case):

| Case | Measured | Required | Prototype |
|---|---|---|---|
| S3 distorted plane / cylinder | 0.00068 h / 0.00114 h (13 / 14 Picard steps) | <= 0.02 h | 0.0007 h / 0.0008 h |
| S4 clustering | first cell 1.00000; 7e-13 h (1 step: exact) | 0.95-1.05; <= 0.05 h | exact |
| S5 overhang | 0 folded, 5.42 deg (47 steps) | 0 folded, >= 2 deg | 5.4 deg |
| S6 sphere patch | 0 folded; min corner angle 68.3 deg | 0 folded | (lat-lon TFI: 88.1 deg) |
| S7 work | 5.19 points per node | <= 10 | |

Unit checks in `tests/test_surface_tfi.py`:
- the lattice derivatives (all five) match the exact sphere within 2e-3;
- TFI is exact on bilinear data;
- the control function equals the discrete `-x_ss / x_s` on a tanh row,
  and the clamp holds;
- geometry work is the same for 1 and 400 iterations.

S6 note (review question 3): on the sphere, Winslow lowers the minimum
corner angle from the lat-lon TFI's 88 deg to 68 deg. Winslow equalises
the grid; it does not preserve angles. The spec asks only for validity
there.

## Acceptance (2026-10-26): surface-grid spec met

`Project.surface_grid(block, face, *, surface, smoothing="elliptic")` is
the command `FillFace` (`shore/project/commands.py`). It works like this:
- It validates the block, the face name and the SURFACE index (an STL body
  has none), and runs the kernel on the face's node view.
- Only then does it snapshot the mesh, so a rejected call leaves no undo
  entry and the nodes untouched.
- It writes the result through the view, so a SHARED face aliased by
  another block stays consistent.
- It journals as `p.surface_grid('annulus', 'i_lo', surface=0)`.

**Verdicts** (T2, with `scripts/check_spec.py` reporting all files ok):

| Module | Job | Result |
|---|---|---|
| `test_surface_tfi.py` | main | 23 passed, 0 skipped, 0 xfailed |
| `test_surface_tfi_cad.py` | CAD | 4 passed, 0 skipped, 0 xfailed |

All `@PENDING` lines are removed; removing them is the only edit made to
either module since the freeze. Full local suite: 1244 passed.

Measured on the CAD case (T1):

| Check | Measured | Required |
|---|---|---|
| Distance from the primitive's exact clustered tensor face | 7.7e-4 local spacing | <= 0.05 |
| Distance from the CAD cylinder (exact projection) | 2.2e-16 | <= 3.8e-9 |
| Rest of the block | bit-identical | bit-identical |
| Run time | 0.01 s | |

The review questions above are still open; changing any answer is a
re-freeze.
