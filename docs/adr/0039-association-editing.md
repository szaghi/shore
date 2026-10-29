# 0039 — Association editing (GUI slice 3)

- **Status:** Proposed
- **Date:** 2026-10-29

## Context

ADR 0003 introduced associations: a block face, edge or vertex linked to
geometry (the body, CAD surfaces, curves or points), whose nodes `snap`
projects onto it with vertex → edge → face priority. The topology builders
associate every wall with its body.

Until now nothing else set an association:
- no command, no journal line;
- `Mesh.snap` resolved only `"body"` (the mesh's own body).

GUI slice 2 (ADR 0038) made CAD faces, edges and vertices pickable. Slice 3
lets the user associate a block entity with them.

The user was away; under the night protocol the choices below are mine and
are listed for review.

## Investigation (2026-10-29, prototype in the session scratchpad `as39/`; T1)

**The command.** An association must be a project command so it is
journaled, undoable and saved. `Associate(block, entity, geometry, kind,
indices)` names:
- the block entity: face `k_lo`, edge `i_lo_k_lo`, vertex `v0`..`v7` in
  `HexBlock.vertices()` order;
- the geometry by its project name, and entity indices.

It stores `Association.of(...)`, ADR 0003's fingerprints, so the
association is re-found after a re-import. Saved projects already keep
associations in their face graphs (sidecar v3), with and without replay.

**Dimensions.** An entity may take geometry of its own dimension or higher:
- a face: surfaces or the body;
- an edge: also curves;
- a vertex: also points.

Projecting a face onto a curve or a point would collapse it. The first
prototype had the comparison reversed; the spec's refusal cases caught
it.

**Snap across geometries.** Snap now resolves every project geometry by
name.
- `"body"` still means the mesh's own body.
- A mesh loaded without replay (ADR 0033) has no `Body`. There `"body"` is
  the single geometry its component was built on. **Before this ADR such a
  project could not snap at all** ("association refers to geometry 'body',
  which this mesh lacks"), a defect from ADR 0033.
- Wall-only snaps give the same nodes, bit for bit, as resolving `"body"`
  alone (A7).

**Seam aliasing after a load without replay.** Built topologies share seam
nodes between blocks as numpy views. Blocks rebuilt from saved grids do
not. So a snap after such a load projects a seam node once per block, and
a replayed project projects it twice. The two differ at round-off (7e-18
on a cubed sphere).
- This ADR does not change it; A6 accepts 1e-12.
- After a load without replay, seams are coincident copies, not shared
  memory. Later in-place edits could separate them.
- This is listed for review.

**The GUI.**
- A selection mode (face, edge, vertex) turns a block pick into the picked
  face, its nearest edge or its nearest corner (`nearest_edge`,
  `nearest_vertex`; checked against `HexBlock`'s edge names, edge slices
  and vertex order).
- A geometry pick selects an entity; shift adds one of the same geometry
  and kind.
- **Associate** runs the command; for a block of another component it runs
  `use_component` first (two journal lines).
- The model tree gains an association column, and the associated edges and
  vertices follow the six faces as rows (ADR 0037's spec pins only the
  first six rows and two columns).

**Prototype run** (headless):
- the three selection modes picked `k_hi`, `j_lo_k_hi` and `v5`;
- a top and a side surface were selected with shift;
- **Associate** wrote `p.associate('box', 'k_hi', 'part', 'surface', (2,
  3))` and the tree showed it;
- undo, redo and snap through the form all worked (the face landed on the
  surfaces);
- a face → curve association was refused with the journal unchanged.

The core and GUI specs pass with the prototype: 17 and 9 cases.

## Decision

- **Command `Associate`:** `Project.associate(block, entity,
  geometry="body", kind="body", indices=())`.
  - Kinds `body | surface | curve | point | none`; the dimension rule
    above.
  - Undo restores the previous association. Nodes move only on `snap`.
  - Journal `p.associate('box', 'k_lo', 'part', 'surface', (2,))`, with
    empty indices left out.
- **Snap** (`Snap.apply`) resolves the project's geometries; `"body"` as
  above.
- **GUI:**
  - the selection-mode combo, `pick_at(x, y, add=)`;
  - an Association dock (the selections, **Associate**, **Clear
    association**);
  - the tree's association column with labels `"body"` / `"part: surface
    2, 3"`, with indices re-found from the fingerprints;
  - edge and vertex selections highlighted.
- **CI:** the core spec needs CAD and runs in the `test-cad` job; the GUI
  spec in `test-gui`.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Associations as view state | They change what `snap` does: project data. |
| Indices in the association | Not persistent across re-imports; ADR 0003 chose fingerprints. |
| A face associated with a curve | Collapses the face onto a line. |
| Refuse associating a block of another component | `use_component` is one journaled line; the user picked that block. |

## Consequences

- **Association editing is journaled.** A GUI association session replays
  headless.
- **Projects loaded without replay can snap.**
- **Not done:** association of slice cells, multi-block selection,
  dragging associated vertices (ADR 0006's vertex drag), tier-2 (exact)
  snap from the GUI (`snap(exact=True)` exists in the form).

## Verification

### Acceptance specs (frozen 2026-10-29, before implementation)

**Verdicts.**
- `tests/spec/test_associate.py` (CAD job): required 17 passed, 0
  skipped, 0 xfailed; at freeze 1 passed (A7, a guard), 16 xfailed.
- `tests_gui/spec/test_gui_associate.py` (GUI job): required 9 passed, 0
  skipped, 0 xfailed; at freeze 9 xfailed.

| Case | Required |
|---|---|
| A1 (× 3) | face → surface, edge → curve, vertex → point: journal line, fingerprinted association, nodes on the entity after snap (1e-9) |
| A2 | vertex → edge → face priority at shared nodes |
| A3 | undo, redo, clear |
| A4 (× 8) | bad block, entity, geometry, kind, indices; dimension rule; nothing changes |
| A5 (× 2) | saved and loaded with / without replay: associations and snap equal |
| A6 | a cubed sphere loaded without replay snaps (within 1e-12 of the replayed snap) |
| A7 | wall-only snap unchanged, bit for bit (guard) |
| U1 (× 3) | face, nearest edge, nearest corner selection |
| U2 | geometry selection, shift adds |
| U3 | associate through the window: journal, tree column, undo |
| U4 | snap through the form after associating |
| U5 | refusals logged, journal unchanged |
| U6 | wall associations shown as `body` |
| U7 | a block of another component: `use_component` first |

### Review questions (for the user)

1. **The dimension rule.** Should a face associated with a curve also be
   refused in ICEM's way? (It is refused.)
2. **Seam aliasing after a load without replay.** Rebuild the shared views
   from the face graph when loading (so seams stay one memory), or accept
   coincident copies?
