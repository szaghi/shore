# 0041 — Redistributing the nodes along a block edge (GUI spacing editor)

- **Status:** Proposed
- **Date:** 2026-10-29

## Context

Roadmap phase 6 lists spacing control. A block edge carries a `Spacing1D`,
but changing the spacing of an existing mesh meant rebuilding it with other
parameters. That is possible only for what the builder exposes, and only
for the builder's own axes.

**User decisions (2026-10-29):**
- a generic command that redistributes a picked edge, propagated to the
  block and through the seams;
- seams kept closed by ADR 0040's write-through, done first.

## Investigation (2026-10-29, prototype in the session scratchpad `redist41/`; T1)

**The law and the lines.**
- The law (`uniform`, `tanh`, `tanh2`, `pinned`; `geometric` is the
  march's) is built over the picked edge's arc length, which gives
  arc-length fractions `t`.
- Every grid line of the block along that axis is resampled at `t` along
  its own old polyline: one batched linear interpolation, vectorized over
  all lines.
- Nodes stay on their lines' polylines and line ends stay put. The two
  faces across the axis do not move, and the block stays structured.

**Through seams.**
- A face that contains the axis and is SHARED or DIRICHLET passes the
  distribution to its partner block, along the partner axis the AxisMap
  names, reversed when that axis is flipped. This repeats until no new
  (block, axis) is reached.
- A (block, axis) reached twice with different fractions is refused, as are
  PERIODIC faces and faces shared with several blocks (CH's `body.k_hi`).
- New positions are computed for every block before any is written, so
  blocks sharing seam memory are not resampled twice. ADR 0040's sync then
  closes the seams.

**Results (prototype):**
- a box edge takes the law exactly;
- a cubed sphere redistributed along i changes the four equatorial blocks,
  leaves the caps bit for bit, and keeps every seam closed;
- the C-grid's wake cut passes the reversed fractions to the lower block;
- CH's shared face is refused.

**Spec pitfalls found while validating** (not prototype defects):
- Resampled lines are measured along their **old** polylines. The new
  nodes' own chords cut corners, so fractions measured on them differ by
  ~1e-3.
- A first D5 checked the wake cut itself. ADR 0040's sync makes the seam
  coincide whatever the distribution, so D5 passed even with the flip
  handling removed. It now checks a line off the seam and fails under that
  mutation.

## Decision

- **`shore.mesh.redistribute`:** `redistribute(blocks, block, edge,
  Spacing1D)`, `plan`, `resample`, `fractions`, `edge_axis`.
- **Command:** `Redistribute` / `Project.redistribute(block, edge,
  law="uniform", *, beta, cluster_at, first_step, last_step, interior)`.
  - Journaled; undo restores the nodes.
  - Associated nodes are not re-snapped: `snap` after if needed.
- **GUI:** a Spacing dock acts on an edge selected in edge mode
  (ADR 0039).
  - It shows the edge's axis, nodes, length, and first and last cell.
  - It sets the law and its parameters; **Redistribute** runs the command,
    making the edge's component current first if needed.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Re-run the builder with other spacing | Only what each builder exposes, and only on its own axes. |
| Redistribute the picked edge only, TFI the faces | Breaks the structured lines through the block and its seam partners. |
| Same fractions on the partner without the flip | Mis-orders nodes across flipped seams (C-grid wake cut). |
| Snap associated nodes automatically | Associations move nodes only on `snap` (ADR 0039). |

## Consequences

- **Any mesh's spacing along a block axis can be changed afterwards,**
  including topologies whose builders do not expose it.
- **Not done:**
  - smoothing after redistribution;
  - redistributing through multi-partner or periodic faces;
  - geometric (march) laws;
  - a per-block-only option.

## Verification

### Acceptance specs (frozen 2026-10-29, before implementation)

**Verdicts.**
- `tests/spec/test_redistribute.py` (main job): required 15 passed, 0
  skipped, 0 xfailed; at freeze 15 xfailed.
- `tests_gui/spec/test_gui_spacing.py` (GUI job): required 5 passed, 0
  skipped, 0 xfailed; at freeze 5 xfailed.

| Case | Required |
|---|---|
| D1 (× 4) | each law on a straight edge; pinned holds its end cells |
| D2 | a curved O-grid: new nodes on the old polylines at fractions `t` |
| D3 | faces across the axis unchanged |
| D4 | cubed sphere: four blocks change, caps do not, seams closed |
| D5 | C-grid wake cut: reversed fractions on a line off the seam |
| D6 (× 5) | refusals, project and journal unchanged |
| D7 | undo, journal line |
| D8 | a project loaded without replay: seams closed |
| P1–P4 (× 5) | the GUI's edge info, redistribution, refusals and panel |
