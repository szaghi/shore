# 0036 — cc.par connections as Xall reads them

- **Status:** Accepted
- **Date:** 2026-10-28

## Context

ADR 0035's CGNS writer had to work around a face shared with several
partners: CH's `body.k_hi` is shared with three channel blocks, but the
face graph names one partner per face.

Checking how `cc.par` handles that face turned up two defects, both
predating 2026-10-28. The user was away, under the night protocol, so the
choices below are mine and are listed for review.

1. **CH meshes cannot be written as an Xall run.** `write_xall` and
   `shore cc-par` fail with "asymmetric connection": the channel blocks'
   `k_lo` faces point at `body.k_hi`, whose single patch points back at only
   one of them.
2. **Connection codes on lo faces are wrong for Xall.** This affects every
   multi-block topology.

## Investigation (2026-10-28, scripts in the session scratchpad `cg/`, `cc36/`; T1)

**Who reads `cc.par`.** Only the overset pre-processor
(`~/fortran/Xnavis-projects/overset/src`); the Xall and Xnavis solvers
read its output. I read the code:
- `readpar.f90:140-230` reads the patches:
  - extents are clamped to the block;
  - a connection's `family` must point back at it (reciprocity);
  - extents must be divisible by `2**(levels-1)`;
  - no check compares the sizes of two partner patches.
- `metrix.f90:685-860` (`adiaface`, `tangind`, `normind`) maps a connection
  patch's ghost nodes onto its partner's nodes:
  - **Along the face:** digit `l` names the partner axis `(l + 1) // 2`.
    An odd digit runs up from the partner patch's lower extent, an even one
    down from its upper extent.
  - **Across the face:** with `m = 1` on a lo face and `0` on a hi face,
    an odd digit puts the partner face at `m × upper` and an even one at
    `(1 - m) × upper`. So **odd pairs opposite sides (lo ↔ hi), even pairs
    the same sides (lo ↔ lo, hi ↔ hi)**, whichever side this face is on.
- Several patches on one face are allowed (each patch is processed on its
  own). That covers CH.

**The interpreter is checked on a real case.** I wrote an interpreter of
those rules and ran it on a solver-run INSEAN case
(`~/fortran/insean/xnavis/xnavis_c/test/CC`: `assetto.grd`, 223 blocks, and
its `cc.par`). It finds **0 errors in 834 connection patches**.
- In that file, a `j_lo` ↔ `j_hi` seam carries code 135 on both sides (89
  pairs); `k_lo` ↔ `k_hi` likewise (98 pairs).

**SHORE's codes.** On the same seams SHORE writes 135 on the hi side and
145 / 136 on the lo side. `_orientation_code_from_axis_map` picked the
normal digit from "is the partner face at its maximum" alone. That is
Xall's rule on hi faces only: on lo faces the parity is reversed.

The interpreter on SHORE's `cc.par`:
- **Wrong:** cubed sphere, sharp and blunt C-grids, OH, split cubed
  sphere, flat caps and the split collar all have lo-face patches whose
  ghost nodes Xall takes from the wrong side of the partner block (nodes
  off by up to the block size).
- **Unaffected:** the O-grid, box, annulus, H-grids and collars, which
  have no connections.
- **With the parity fixed:** every one of them maps onto coincident nodes.

Only the lo-face codes change. In `tests/test_cc_par.py`, 10 unit tests pin
the old lo-face codes (a "historical table"); nothing shows that table was
ever checked against Xall.

**CH.** With one patch per partner on `body.k_hi`, each over the partner's
sub-range, the cc.par of CH (sharp and blunt) is reciprocal and maps onto
coincident nodes.
- The sub-ranges need node coordinates, which `cc_par` (it reads only the
  adjacency sidecar) does not have.
- So the sidecar gains an optional `patches` list on such faces only: per
  partner, a 0-based node `range`, the `partner` and the `axis_map`.
- It is written by `write_adjacency_json` from the blocks' nodes. The seam
  location is ADR 0035's, moved from `io/cgns.py` to `mesh/seams.py` and
  shared.
- Faces with one partner are not touched, so every other sidecar stays
  byte-identical, and blocks reloaded from rounded `.geo` files are never
  matched.

**Also found:** `write_xall` of a one-component project failed (`merge_grd:
need at least 2 input files`). It is fixed separately, before this ADR:
`write_xall` merges through internal helpers that accept one input, and the
public `merge_*` functions keep their two-input contract.

## Decision

- **Normal digit:** odd when the two faces are on opposite sides of their
  blocks, even when on the same side (Xall's `normind`). Along-face digits
  are unchanged.
- **Faces with several partners:**
  - `write_adjacency_json` adds `patches` to them (the face's own
    `partner` / `axis_map` stay, for older readers);
  - `write_cc_par` writes a patch per entry;
  - a single-patch face pointing at such a face gets the sub-patch that
    points back.
- **`mesh/seams.py`:** the seam location is shared by `cgns` and the
  sidecar (`seams`, `shared_faces`, `transform_vector`).
- **Unit tests:** `tests/test_cc_par.py`'s lo-face codes are updated to
  Xall's convention.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Split CH's `body` block into three | Changes CH's block count and every CH output, `.geo` files included. |
| Infer CH's sub-ranges in `cc_par` from patch sizes | The order of the three channel blocks along `body.k_hi` is not in the sidecar. |
| Version 4 sidecar | The key is additive and optional; readers that ignore it see today's sidecar. |
| Keep the codes, document them | Xall builds ghost cells from the wrong nodes on every lo-face seam. |

## Consequences

- **Every multi-block `cc.par` changes** on its lo-face connection patches
  (the normal digit only). Runs prepared with earlier SHORE versions should
  be regenerated.
- **CH meshes can be run in Xall.**
- **Not checked:** running Xall itself, and Xall's cell-centre mapping
  (`calcorn.f90`). Extents must still be divisible by
  `2**(multigrid levels - 1)`, or Xall stops; that is a mesh-size choice.
  SHORE does not check it.

## Verification

### Acceptance spec (frozen 2026-10-28, before implementation)

`tests/spec/test_cc_par_xall.py` runs in the main CI job, with
`data/xall_before.zip`: the `cc.par` and sidecar of nine block sets,
written before this ADR.

**Verdicts.**
- Required: 45 passed, 0 skipped, 0 xfailed.
- At freeze: 31 passed (guards: Q0, Q1 for the sets without connections, Q2
  for the symmetric sets, Q3, Q5), 14 xfailed.

| Case | Required |
|---|---|
| Q0 (× 2) | the interpreter accepts the real file's codes (135 / 135) on hand-built j and k seams and rejects SHORE's old lo-side codes |
| Q1 (× 11) | every connection, as Xall reads it, lands on coincident nodes (1e-12 × diagonal); reciprocal, equal lengths, inside the block |
| Q2 (× 11) | each block's patches tile its six faces once |
| Q3 (× 9) | against the writer before: only lo-face connection codes differ, only in the normal digit |
| Q4 (× 2) | CH sharp and blunt: `write_xall` works; `body.k_hi` gets three patches and is the only face with `patches` |
| Q5 (× 9) | the other sidecars are byte-identical to before |
| Q6 | a CH project saves, loads without replay and saves the same files |

### Acceptance (2026-10-28): spec met

**Implementation.** As decided:
- `mesh/seams.py` (`seams`, `shared_faces`, `transform_vector`), shared by
  `io/cgns.py`;
- `write_adjacency_json` writes `patches` on shared faces;
- `write_cc_par`: `_sub_patch` and the family resolution for shared faces;
  `_orientation_code_from_axis_map`'s normal digit by side parity;
- `tests/test_cc_par.py`: 10 lo-face codes updated (three of them, 415,
  135 and 325, appear with the same face pairs in the INSEAN file);
- `docs/reference/cc-par.md` rewritten on the codes.

**Verdicts** (T2 locally; `scripts/check_spec.py` reports all files ok):
- `tests/spec/test_cc_par_xall.py`: **45 passed, 0 skipped, 0 xfailed**.
  All `@PENDING` lines are removed.
- `tests/spec/test_exports.py` still 61 passed (the seam code moved).

### Review questions (for the user)

1. Is the INSEAN case (`assetto`) the right reference for the code
   convention, or does a newer pre-processor read `cc.par` differently?
2. Have SHORE-made multi-block runs been used in Xall? If so, their ghost
   cells on lo-face seams were wrong, and those runs should be redone.
