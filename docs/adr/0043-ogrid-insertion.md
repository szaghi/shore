# 0043 — O-grid insertion in a blocking (phase 7, slice 2)

- **Status:** Proposed (spec frozen 2026-10-29)
- **Date:** 2026-10-29

## Context

ADR 0042's blocking around a body is a Cartesian lattice with a hole. At the
hole's edges three blocks share the body's smooth surface, so the grid lines
leave the wall up to 41° off its normal (26 blocks around a sphere). ICEM's
remedy is the O-grid: a layer of blocks wrapped around the body, so that the
lines leaving the wall are nearly normal to it.

**User decisions (2026-10-29):**
- one command for both cases: an O-grid inside a block set (a core and a
  ring of blocks, optionally open on excluded faces), and a layer on the
  faces of a hole;
- the size given by a scale factor;
- the GUI through a form only (picking blocks and faces later).

## Investigation (2026-10-29, prototype in the session scratchpad `ogrid43/`; T1)

**The operation.** A block set S (default: every block) and its layered
faces W. W is either given, or S's boundary faces (faces of S blocks not
shared with another S block) less the excluded ones.
- Every vertex `v` of a W face gets an inner copy `v'`.
- The S blocks take the copies, shrinking away from W.
- Each W face gets a layer block between itself and its copy, oriented like
  its S block (right-handed), with the face's node counts along it and
  `nodes` across.
- S's boundary faces outside W stay open: the layer passes through them,
  and the copies on them stay in their plane. This is the half-O at a
  symmetry plane.
- A face that touches the layer must be one of:
  - inside S;
  - layered;
  - open on the blocking's boundary.

  Otherwise the blocking would stop conforming, and the operation is
  refused.

**Where the copies go.** `v' = v + (1 − factor) · ½ Σ e`, summed over the
distinct across edges at `v`. An across edge is the edge of an S block from
`v` that leaves a W face into its block.
- **Inside one block,** the across edges are the block's three edges at
  `v`, so `v'` is the block scaled by `factor` about its centre (ICEM's
  factor).
- **Two blocks in a row:** the middle vertices move only across, so the
  shared face stays plane.
- **An open face** gives no across edge, so copies on it stay in its plane.
- **The faces of a hole:** the across edges point away from the hole, so the
  layer grows into the blocks around it.
- The obvious alternative was rejected. Moving `v` toward the mean centroid
  of its S blocks cancels at a hole corner: 7 octant blocks give a 0.14
  shift where the across edges give 1.

**Associations.**
- A layered face keeps its own, now on the layer's outer face.
- A face or edge that the copies replace passes its association on.
- An open face also passes its association to the layer faces lying on it.

**Results (prototype; wall angle = largest angle between a wall node's first
grid segment and the sphere's normal):**

| Blocking around the unit sphere | blocks | min scaled Jacobian | wall angle |
|---|---|---|---|
| lattice [-3, 3] cut at ±1 (ADR 0042) | 26 | 0.70 | 40.8° |
| same + hole layer (factor 0.5) | 32 | 0.44 | 5.7° |
| lattice [-4, 4] cut at ±2 (ADR 0042) | 26 | 0.40 | 24.6° |
| same, centre O-gridded, core associated and deleted | 32 | 0.58 | 4.8° |

- The global minimum Jacobian is the O-grid's own diagonal corners away
  from the wall: 1/√3 for a cube's ring.
- One block, a half-O and a two-block row all premesh valid (≥ 0.577), with
  seam gap 0.
- **The centre-block route fails** when the centre box is not larger than
  the body (±1 around the unit sphere: Jacobian −0.13). The core's faces
  project out to the sphere, which touches the old box. The hole-face mode
  has no such limit.

**Spec pitfalls found while validating:**
- The half-O's top corners move to z = 0.5, not 0: every corner has a
  layered face above it.
- A half-O also has 32 skeleton edges, because every corner lies on a
  layered side face.
- The hole faces must be named before the centre block is deleted.

**Mutations the spec catches:**
- centroid placement instead of across edges: 2 cases fail;
- no association transfer: 1 fails;
- no conformity check: 1 fails.

## Decision

- **`shore.mesh.blocking.ogrid(blocking, labels=(), faces=(), exclude=(),
  factor=0.5, nodes=5)`.**
- **Command:** `Ogrid` / `p.ogrid(blocks=(), *, faces=(), exclude=(),
  factor=0.5, nodes=5)`. Faces are `(block, face)` pairs. It is journaled,
  and undo restores the blocking.
- **Labels:** S blocks keep theirs; layer blocks are `b<n>` in block and
  face order.
- **GUI:** an O-grid form, with `blocks`, `faces` and `exclude` as literals.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Copies toward the mean centroid of their blocks | Cancels at hole corners: thin layers. |
| Two commands (inside / around) | The same operation; the user chose one. |
| Absolute layer thickness | The user chose a scale factor (ICEM's). |
| O-grid through faces shared with blocks outside the set (propagate into them) | Non-conforming unless those blocks are split too; refused, add them to the set. |

## Consequences

- **Body-fitted layers around a hole or inside a block set,** with grid
  lines within 6° of the wall's normal around a sphere, where ADR 0042's
  lattice gave 41°.
- **Not done:**
  - picking blocks and faces in the GUI;
  - O-grids through faces shared with blocks outside the set;
  - smoothing of the layer;
  - collapse, merge and vertex moves.
- **Naming the hole faces** takes `(block, face)` pairs. The GUI picking
  slice will make this easier.

## Verification

### Acceptance specs (frozen 2026-10-29, before implementation)

**Verdicts.**
- `tests/spec/test_ogrid.py` (main job): required 13 passed, 0 skipped, 0
  xfailed. At freeze: 13 xfailed. With the prototype behind it: 13 passed.
- `tests_gui/spec/test_gui_ogrid.py` (GUI job): required 3 passed, 0
  skipped, 0 xfailed. At freeze: 3 xfailed. With the prototype: 3 passed.

| Case | Required |
|---|---|
| O1 | one block: core at factor scale, labels, counts, 7 blocks valid, seams closed |
| O2 | an open face: core in its plane, its association on the 5 faces on it, those faces WALL |
| O3 | two blocks, all of them: the across-edge placement, the shared face kept |
| O4 | the hole of ADR 0042's sphere blocking: 32 blocks, wall on the STL, Jacobian > 0.3, wall angle < 10° (> 30° without) |
| O5 | the centre block O-gridded, its core associated and deleted: wall angle < 10° |
| O6 | a layered face keeps its association on the layer only |
| O7 | journal line, replay bit for bit, undo |
| O8 (× 6) | refusals, project and journal unchanged |
| W1–W3 | the GUI form: an O-grid, a half-O, a refusal |
