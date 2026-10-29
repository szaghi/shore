# 0044 — Selecting blocking faces and blocks in the GUI (phase 7, slice 3)

- **Status:** Proposed (spec frozen 2026-10-29)
- **Date:** 2026-10-29

## Context

ADR 0043's O-grid names its block set and faces as labels and `(block,
face)` pairs, typed into a form. Its Consequences named that as the pain
point: the faces of a hole around a body belong to the blocks around it,
and their labels come from the split history. ADR 0042's GUI picks only
skeleton edges, for splitting.

**Decided unattended (night of 2026-10-29), for review:** this slice comes
before vertex moves. Vertex moves are the next roadmap item, but every later
blocking operation (O-grid, delete, associate, and the wing–body acceptance)
needs faces and blocks picked in the view. The design below was chosen
without the usual scope questions.

## Investigation (2026-10-29, prototype in the session scratchpad `night43/`; T1)

**Markers, not translucent faces.**
- A picked translucent face is the one nearest the eye, so the faces of a
  hole behind a lattice's outer faces could not be picked.
- Instead, face mode and block mode draw a marker at each face's or block's
  skeleton centre (pickable points, drawn over everything, like the CAD
  vertices of ADR 0039). Any face can then be picked from a view in which
  its marker is not covered by another.
- Each face is drawn once, named by its first block.

**One selection.**
- Edges (skeleton picks), faces and blocks go into one list of
  `SkeletonKey`s.
- A pick replaces the selection; shift adds the entity, or removes it if
  already selected (ICEM-like).
- The Blocking dock's selection mode switches the markers; the selection
  survives mode changes, so blocks and faces can be combined.
- After a command, keys that no longer exist are dropped.

**Actions on the selection.**
- **O-grid:** the selected blocks are the set (every block when none is
  selected). The selected faces are the layered faces, each renamed to its
  block inside the set, since a face is named by its first block, which may
  lie outside the set.
- **Delete blocks:** the selected blocks.
- **Associate:** each selected face or edge with the geometry selected in
  the Association dock.
- One command per entity, so one journal line each. The selection's
  component is made current first.

**Results (prototype):**
- a split row shows 11 face markers and 2 block markers;
- picking, shift-adding and toggling select as described;
- a picked block's O-grid journals `p.ogrid(('b1',), factor=0.6,
  nodes=4)`;
- the six hole faces of a 3 × 3 × 3 lattice, selected, give one
  `p.ogrid((), faces=(...))` and 32 blocks.

**Mutation the spec catches:** shift without the toggle (1 case fails).

## Decision

- **Model:** `Blocking.centres(kind)`.
- **Session:** `Session.skeleton_markers(kind)`. `SkeletonKey.edge` holds
  an edge or face name, or `""` for a block.
- **Viewport:** `PygfxViewport.set_markers`, `marker_keys`, and
  `MarkerPick(key, xyz)`.
- **Window, Blocking dock:**
  - `blocking_mode` (edge / face / block);
  - `blocking_selection`, `select_skeleton(key, add=)`,
    `clear_blocking_selection()`;
  - `ogrid_selection(factor, nodes)` (panel fields `ogrid_factor`,
    `ogrid_nodes`), `delete_selection()`,
    `associate_blocking_selection()`;
  - the Display dock's blocking toggle also hides the markers.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Translucent pickable faces | The nearest face wins: a hole's faces hide behind the outer ones. |
| Select a face by two of its edges (ICEM) | Two picks per face; markers are one. |
| Rubber-band (box) selection | Mouse-drag handling not testable headless; later. |
| Vertex moves first (roadmap order) | Every blocking operation needs faces and blocks picked; done first. |

## Consequences

- **Every blocking command can now be driven from the view:** create (form)
  → split (edge pick) → O-grid / delete / associate (selection) → premesh
  (form).
- **Not done:** rubber-band selection; markers hidden by depth (all
  markers are always drawn); vertex moves.

## Verification

### Acceptance spec (frozen 2026-10-29, before implementation)

`tests_gui/spec/test_gui_blocking_select.py` (GUI job): required 9 passed,
0 skipped, 0 xfailed. At freeze: 9 xfailed. With the prototype behind it:
9 passed.

| Case | Required |
|---|---|
| K1 | markers per mode, at the centres; hidden with the blocking |
| K2 (× 2) | picking markers: select, add, toggle, replace; an edge pick selects too |
| K3 | an O-grid of a picked block from the panel: the journal line |
| K4 | a lattice's hole: six faces selected, one O-grid, 32 blocks |
| K5 | associate selected faces with the selected geometry; delete a picked block |
| K6 (× 3) | refusals: nothing selected, faces given to delete, no geometry |
