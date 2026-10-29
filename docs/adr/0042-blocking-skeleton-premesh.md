# 0042 — Blocking: a block skeleton, split and premeshed (phase 7, slice 1)

- **Status:** Accepted
- **Date:** 2026-10-29

## Context

Roadmap phase 7 (GUI v2) is interactive blocking in the style of ICEM CFD:
build a block topology around the geometry by editing it (split, O-grid
insertion, collapse, merge, vertex moves), associate it with the geometry,
and mesh it. Its acceptance is a wing–body blocking built interactively and
replayed from its journal after a CAD edit.

Until now every SHORE mesh came from a builder: a topology, a primitive, a
collar, blocks read from files. The user can edit the result (associate,
snap, redistribute: ADR 0039–0041) but cannot make a block topology of
their own.

**User decisions (2026-10-29):**
- model: **skeleton + premesh**. A blocking holds no nodes. Commands edit
  it, and a premesh command makes ordinary `HexBlock`s, so the rest of the
  pipeline (snap, redistribute, Xall, CGNS, the GUI) works unchanged;
- first slice: **create, split, premesh**.

## Investigation (2026-10-29, prototype in the session scratchpad `block42/`; T1)

**The blocking.**
- Skeleton vertices (points), and blocks of 8 vertex ids in
  `HexBlock.vertices()` order (index `4 i + 2 j + k`). Blocks that touch
  share vertex ids, so a face or an edge is one entity in every block that
  has it.
- A node count per edge. A block's four parallel edges always share a
  count: `create` gives them one, and a split cuts them all at the same
  node.
- Associations keyed by entity: a vertex id, an edge's vertex pair, a
  face's four vertex ids.

**Split.**
- `split(block, axis, fraction)` cuts at node `round(fraction (n - 1))`, so
  counts stay integers. The new vertices lie on the straight skeleton edges.
- The cut is carried through every block that shares a cut edge (the
  sheet). Carrying it through shared edges, not only shared faces, keeps the
  blocking conforming where blocks meet only along an edge.
- An edge cut twice, or a block reached along two axes or nodes, is refused.
- The low part keeps the label; each cut block's high part is `b<n>`, with
  `n` counting up over the blocking in block order. Labels are therefore
  deterministic on replay.
- Associations of cut edges and faces pass to both halves; new vertices and
  edges get none of their own (they inherit, below).

**Delete.** Removing the centre of a 3 × 3 × 3 lattice leaves a hole whose
six faces are boundary faces. ICEM builds external blockings this way
before O-grids. Delete was not in the chosen slice, but a blocking around a
body needs it; it is a few lines and is included.

**Premesh.**
- Vertices are projected onto their associations.
- Each edge is computed once: equally spaced nodes on the straight line,
  then, if associated, projected and re-equidistributed along the projected
  polyline (three passes, `redistribute.resample`).
- Each face is computed once: a Coons patch (TFI) of its four edges, its
  interior then projected onto its association.
- Each block is a TFI of its six faces.
- Edges and faces are copied into every block that has them, so seams match
  bit for bit. They are oriented by the eight corner permutations.
- Shared faces become SHARED seams, with axis maps from the corner ids.
  Boundary faces are WALL when associated and FREE otherwise.
- The blocks carry the associations, so `snap` after premesh moves nodes by
  at most round-off (1e-16).

**Inheritance.** An entity without its own association takes one from the
entities above it, lowest dimension first, when they agree:
- an edge, from its associated faces;
- a vertex, from its associated edges, else from its associated faces.

Associating the six faces of the hole with the body is then enough: the
hole's edges and corners follow. An edge between two faces associated with
different geometry stays straight until it gets an association of its own,
because nearest-of-two is not the intersection curve.

**Results (prototype):**
- a split box premeshes to the Cartesian grid within 1e-12;
- 26 blocks around the unit sphere (a 3 × 3 × 3 lattice of [-3, 3]³, cut at
  ±1, centre deleted, its faces associated with the body):
  - minimum scaled Jacobian 0.71;
  - wall nodes on the STL;
  - every seam closed (gap 0);
  - 50 ms;
- the journal replays bit for bit, and save/load and the Xall and CGNS
  writers work;
- the same journal replayed with an ellipsoid in place of the sphere
  re-projects the walls and stays valid.
- A rotated, right-handed neighbour block gets axis maps that match its
  nodes. Slice 1's commands never make one; O-grid insertion will.

**GUI.**
- The skeleton is drawn as orange lines over the blocks, so it stays
  pickable, like CAD curves (ADR 0039). Each edge is drawn with its node
  count.
- A pick returns the nearest node (pygfx's `vertex_index`). **Split here**
  cuts the picked edge's block at that node.
- The blocking commands also have forms.

**Spec pitfalls found while validating** (not prototype defects):
- STL files store float32, so distances must be measured to the file as
  read, not to the in-memory mesh it was written from (3e-8 apart).
- An ellipsoid poking out of the hole cube (1.3 × along x) pushes through
  the first layer of blocks (Jacobian 0.27). The edited body is now
  1.15 × 0.85 × 0.95 (0.46).
- A box corner belongs to three edges, so a corner pick may name any of
  them.

**Mutations the spec catches:** inheritance off (4 cases fail); a cut not
carried through the sheet (6 fail).

## Decision

- **`shore.mesh.blocking`:** `Blocking` (`points`, `blocks`, `nodes`,
  `associations`; `labels`, `block`, `shape`, `entity`, `edges`, `faces`,
  `polylines`), `box`, `split`, `delete`, `associate`,
  `effective_associations`, `premesh`, `tfi_face`, `tfi_volume`.
- **Commands** (journaled, undo restores the previous blocking or mesh),
  with `Project.blocking` (the current component's) and `Project.blockings`:
  - `CreateBlocking` / `p.create_blocking(ni, nj, nk, *, geometry, margin,
    vmin, vmax)`;
  - `SplitBlocking` / `p.split_blocking(block, axis, fraction)`;
  - `DeleteBlocks` / `p.delete_blocks(labels)`;
  - `AssociateBlocking` / `p.associate_blocking(block, entity, geometry,
    kind, indices)`, validated like `associate` (shared `_association`);
  - `Premesh` / `p.premesh(exact=False)`: replaces the component's mesh;
    `BuildState("blocking", ..., marchable=False, geometries=<associated>)`.
- **Persistence:** the blocking is rebuilt by replaying the journal. A
  project loaded without replay (ADR 0033) has its mesh but no blocking.
- **GUI:**
  - `Session.skeleton()` and `SkeletonKey`;
  - `PygfxViewport.set_skeleton` and `skeleton_keys`;
  - `SkeletonPick`;
  - the Display dock's Blocking toggle;
  - a Blocking dock with **Split here** (`MainWindow.split_selection`);
  - forms for the five commands.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Blocking = live `HexBlock`s edited in place | Every edit would remesh and re-wire views; the user chose skeleton + premesh. |
| Carry cuts through shared faces only | Leaves hanging vertices on blocks that share only an edge with the sheet. |
| Split at a free fraction (not a node) | Needs non-integer counts or a later re-count; node-aligned splits keep counts exact. |
| Inherit associations at split time (copy to new entities) | Stale copies after the parent association changes; inheritance at premesh follows the current state. |
| Edges between differently associated faces projected on both | Nearest-of-two is not the intersection curve; such edges need a curve association. |
| Leave delete for a later slice | No external blocking (around a body) is possible without it. |

## Consequences

- **A user-made topology** around any body, without a builder. The
  premeshed blocks go through the existing pipeline unchanged.
- **Quality is TFI's:** the hole edges of a Cartesian lattice around a body
  give 60° cells (minimum scaled Jacobian 0.7 around a sphere). O-grid
  insertion (a later slice) is the remedy.
- **Not done:**
  - merge, collapse, O-grid insertion, vertex moves;
  - spacing laws on blocking edges (use `redistribute` after premesh);
  - elliptic smoothing;
  - the blocking saved apart from the journal.
- **Next slices:** O-grid insertion, then vertex moves along associated
  entities. Then the roadmap's wing–body acceptance.

## Verification

### Acceptance specs (frozen 2026-10-29, before implementation)

**Verdicts.**
- `tests/spec/test_blocking.py` (main job): required 17 passed, 0 skipped,
  0 xfailed. At freeze: 17 xfailed. With the prototype behind it: 17
  passed.
- `tests_gui/spec/test_gui_blocking.py` (GUI job): required 6 passed, 0
  skipped, 0 xfailed. At freeze: 6 xfailed. With the prototype: 6 passed.

| Case | Required |
|---|---|
| B1 | create around a geometry: corners from its bounds and the margin |
| B2 | create between corners; the journal line |
| B3 (× 2) | split: counts, cut position, labels; rounding to a node |
| B4 | a cut carried through a row of three blocks; premesh valid, seams closed |
| B5 | premesh without associations: the Cartesian grid within 1e-12; SHARED / FREE |
| B6 | 26 blocks around a sphere: 6 WALL, 54 FREE, 96 SHARED; wall on the STL within 1e-9; scaled Jacobian > 0.5; seam gap 0 |
| B7 | snap after premesh moves no node more than 1e-12 |
| B8 | journal replay and save/load bit for bit; Xall with 26 blocks |
| B9 | the journal replayed on an edited body: wall on it, scaled Jacobian > 0.3, seams closed |
| B10 | inheritance: one associated face, disagreeing faces (straight edge), an own association |
| B11 (× 5) | refusals, project and journal unchanged |
| B12 | undo of split, delete and premesh |
| G1 | the skeleton drawn: 12 edges, their nodes |
| G2 | a picked node, Split here: the journal line, 20 edges |
| G3 (× 2) | refusals: nothing picked; an end node |
| G4 | premesh drawn; the skeleton picked over the blocks; hidden |
| G5 | the delete and associate forms |

## Acceptance (2026-10-29): both specs met

**Implementation.** As decided, plus:
- unit tests (`tests/test_blocking_unit.py`):
  - edge and face tables;
  - entity keys;
  - split counts and associations;
  - refusals;
  - TFI of a face and of a sheared volume;
  - the axis map of a rotated neighbour, checked against its nodes;
  - inheritance;
  - polylines;
- `docs/reference/run.md` (Blocking) and `docs/reference/gui.md` (the
  Blocking dock).

**Verdicts** (T2 locally; `scripts/check_spec.py` reports all files ok):
- `tests/spec/test_blocking.py`: **17 passed, 0 skipped, 0 xfailed**;
- `tests_gui/spec/test_gui_blocking.py`: **6 passed, 0 skipped, 0
  xfailed**.

All `@PENDING` lines are removed. Full suites: main 1896 passed, GUI 72
passed.
