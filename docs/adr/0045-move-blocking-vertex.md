# 0045 — Moving blocking vertices (phase 7, slice 4)

- **Status:** Proposed (spec frozen 2026-10-29)
- **Date:** 2026-10-29

## Context

Roadmap phase 7 lists "vertex drag on associated entities". The blocks of
ADR 0042–0044 are boxes cut by splits and O-grids, and their vertices sit
where those operations put them. Shaping a blocking to a body (pulling a
corner to a wing tip, a vertex onto a curve) needs vertices moved one by
one.

**Decided unattended (night of 2026-10-29), for review:**
- **No mouse drag.** The command moves one vertex to a point. The GUI sends
  typed coordinates or the last geometry pick. Mouse drag is out of scope:
  it cannot be tested headless, and the command is what the journal records.
- **An associated vertex slides on its geometry.** The stored point is the
  projection of the target onto the vertex's association (own or
  inherited, ADR 0042).

## Investigation (2026-10-29, prototype in the session scratchpad `night43/`; T1)

- `move(blocking, label, vertex, xyz, resolver)` replaces one skeleton
  point. Blocks share vertex ids, so every block having the vertex follows.
- **Projection at move time, not only at premesh.** The skeleton then shows
  the vertex where premesh will put it, on the geometry. After a geometry
  edit, replay re-projects the stored point (premesh projects vertices
  anyway): a vertex moved onto the old sphere ends on the edited one.
- **Results (prototype):**
  - a vertex shared by two blocks of a row moves in both;
  - premesh corners equal the target, and the shared face matches bit for
    bit;
  - minimum scaled Jacobian 0.84;
  - an inherited association (from an associated face) projects a target
    of (0.3, 0.9, 0.9) onto the sphere's STL.
- **GUI:** the Blocking dock's vertex mode adds a marker per vertex (ADR
  0044). **Move vertex** takes "x, y, z"; **Move vertex to geometry pick**
  takes the last geometry pick.
- **Mutation the spec catches:** no projection (1 case fails).

## Decision

- **Model:** `shore.mesh.blocking.move`; `Blocking.centres("vertex")`.
- **Command:** `MoveVertex` / `p.move_vertex(block, vertex, xyz)`. It is
  journaled; undo restores the point.
- **GUI:**
  - vertex markers (`Session.skeleton_markers("vertex")`);
  - `move_selection(xyz)` (panel field `move_xyz`) and
    `move_selection_to_pick()` (`last_geometry_point`);
  - a *Move blocking vertex* form.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Mouse drag in the view | Not testable headless; the journal records the end point anyway. Later, on top of this command. |
| Store the raw target, project only at premesh | The skeleton would show the vertex off its geometry. |
| Refuse moves that invert a block | Moves are often made through invalid intermediate states; premesh quality is reported by the existing quality view. |

## Consequences

- **A blocking can be shaped by hand,** with associated vertices staying on
  their geometry.
- **Not done:**
  - mouse drag;
  - moving several vertices, or edges and faces, at once;
  - checking validity after a move.

## Verification

### Acceptance specs (frozen 2026-10-29, before implementation)

**Verdicts.**
- `tests/spec/test_move_vertex.py` (main job): required 9 passed, 0 skipped,
  0 xfailed. At freeze: 9 xfailed. With the prototype: 9 passed.
- `tests_gui/spec/test_gui_move_vertex.py` (GUI job): required 5 passed, 0
  skipped, 0 xfailed. At freeze: 5 xfailed. With the prototype: 5 passed.

| Case | Required |
|---|---|
| M1 | a shared, unassociated vertex: point, journal line, premesh corners, shared face |
| M2 | an inherited association: onto the STL; premesh keeps it |
| M3 | the edited body replayed: re-projected onto it |
| M4 | undo; replay bit for bit |
| M5 (× 5) | refusals: block, vertex, entity kind, two coordinates, NaN |
| X1 | vertex markers, each vertex once |
| X2 | a picked vertex moved from the panel: journal, skeleton redrawn |
| X3 | a vertex moved to a picked geometry point |
| X4 (× 2) | refusals: two entities selected, a bad coordinate text |
