# 0051 — Phase 7 acceptance: a wing–body blocking replayed after a CAD edit

- **Status:** Accepted
- **Date:** 2026-10-30

## Context

Roadmap phase 7 (GUI v2, interactive blocking) ends on an acceptance test:
*a wing–body blocking built interactively, replayed from its journal after a
CAD edit.* ADR 0042–0050 delivered the pieces:
- blocking, split and premesh;
- O-grids;
- selection in the view;
- vertex moves;
- wall-normal lines;
- STL curves;
- wall layers.

This ADR runs the acceptance itself.

**Decided unattended (night of 2026-10-30), for review:**
- **Geometry:** a CAD model (build123d, STEP). No boolean engine is
  installed for trimesh, so a fused STL was not possible, and a CAD edit is
  the acceptance's own wording.
- **Model:** a half model with a symmetry plane at z = 0, like ADR 0046's
  wing on a plane, plus a fuselage:
  - the fuselage is a capsule along x (radius 0.5, cylinder length 5);
  - the wing is a blunt-TE NACA 0012 (chord 1) along +z to z = 3, fused
    with the fuselage;
  - 8 faces, 17 edges.
- **Edit:** span 2.7 and fuselage radius 0.55. The edit must stay inside
  the blocking's hole: a longer span would pierce the wing's block. The
  radius change moves the junction, which tests re-finding its edges by
  fingerprint.
- **No implementation is pending:** the specs pass at freeze. They record
  that the existing commands meet the acceptance, and they guard it.

## The blocking (T1, scratchpad `wb51/`; T2 by the specs)

Every step is a journaled command:
1. **Box and cuts:** a box [-6, 7] × [-4, 4] × [0, 6] (53 × 33 × 25
   nodes), cut at x = −3, −0.3, 1.3, 3, y = ±0.8, ±0.3, z = 0.8, 3.3.
2. **Hole:** the 9 blocks round the fuselage and the block round the wing.
   - Their 25 faces that bound the hole (not on the symmetry plane) are
     associated with the body.
   - The 4 edges at the base of the wing's block are associated with the
     3 junction edges (the CAD edges on both the fuselage and the wing;
     ADR 0039's curve association).
   - The 10 blocks are deleted.
3. **O-grid** on the hole's 25 faces (factor 0.7, 17 nodes across; ADR
   0043).
4. **`premesh(exact=True)`:** tier-2 projection. Tier 1 left wall nodes on
   the junction up to 4e-3 off the body: the CAD edge's tessellation and
   the faces' do not coincide.
5. **`orthogonalize(0.5)`** (ADR 0048): strength 1 is refused, because the
   wing's side layers at the concave junction would invert; there the
   normal is the bisector.
6. **`wall_layer(2e-3)`** (ADR 0050).

**Results:**

| | Original model | Edited model (replayed) |
|---|---|---|
| blocks | 90 | 90 |
| min scaled Jacobian | 0.17 | 0.135 |
| seams | closed (gap 0) | closed |
| wall nodes off the body (exact) | 9.5e-10 | 2.5e-9 |
| first cells | 2e-3 | 2e-3 |
| wall angle median / p95 (off plane, TE, tip, junction) | 9.0° / 25.3° | 9.5° / 26.4° |
| junction ring on the junction curves | < 1e-8 | < 1e-8 |

- **Before `orthogonalize`** the wall angle was 24.7° median and 68.6° p95.
- **Without the junction associations,** no wall node lands on the
  junction (all more than 1e-3 off).
- **Mutation:** a zero wall tangent in `orthogonalize` fails 2 of the 6
  headless cases.

## Interactive

The GUI spec builds the same blocking through the window:
- forms for the geometry, the box, the splits, premesh, orthogonalize and
  the wall layer;
- the Blocking dock's face, edge and block selections (ADR 0044), with the
  Association dock's geometry for the body and the junction curves, Delete
  and O-grid.

The session's journal, replayed headless, gives the session's grid bit for
bit, and on the edited model a valid grid on the new body.

## Decision

Phase 7's acceptance is met. The specs are frozen as guards:
- `tests/spec/test_wing_body.py` (CAD job);
- `tests_gui/spec/test_gui_wing_body.py` (GUI job, `shore[gui,cad]`).

## Consequences

- **Phase 7 is closed on this acceptance**, with recorded limits:
  - a half model;
  - an edit inside the hole;
  - concave junctions need `orthogonalize` below strength 1;
  - sharp edges (TE, tip) keep wall angles near 45°;
  - wall layers are algebraic (hyperbolic marching of a blocking's walls
    deferred, ADR 0050).
- **Open:** full-span and multi-component models, blockings that adapt to
  edits leaving the hole, 3-D elliptic smoothing.

## Verification

**Verdicts** (T2 locally; `scripts/check_spec.py` reports all files ok):
- `tests/spec/test_wing_body.py`: **6 passed, 0 skipped, 0 xfailed**;
- `tests_gui/spec/test_gui_wing_body.py`: **2 passed, 0 skipped, 0
  xfailed**.

| Case | Required |
|---|---|
| A1 | the model: 8 faces, 3 junction edges |
| A2 | 90 blocks; Jacobian > 0.1; seams closed; walls within 1e-8 of the body; first cells 2e-3; wall angle median < 12°, p95 < 35°; junction ring within 1e-8 of the junction |
| A3 | journal replay bit for bit |
| A4 | the edited model: as A2 on the new body and junction |
| A5 | without the junction associations no wall node reaches the junction |
| A6 | the journal's command sequence (as the GUI builds it) |
| I1 | built in the GUI: valid, closed, on the body, first cells; its journal replays bit for bit |
| I2 | the GUI's journal on the edited model: valid, on the new body |
