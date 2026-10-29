# 0040 — Seams kept closed by in-place edits

- **Status:** Proposed
- **Date:** 2026-10-29

## Context

SHARED and DIRICHLET seams are one set of points on both sides.
- Built topologies share SHARED seam nodes as numpy views: one buffer cut
  into blocks.
- DIRICHLET seams (the cubed sphere's cap–equator edges, whose axes are
  permuted) are copies even after a build.
- Blocks read back from saved grids (ADR 0033) hold copies of every seam.

ADR 0039 found the consequence and listed it for review. An in-place edit
that changes one side of a seam does not change the other unless they share
memory. A snap with an association on one side of a DIRICHLET seam, or of
any seam after a load without replay, opens the seam.

**The user chose (2026-10-29):**
- rebuild the shared views;
- then, once rebuilding views turned out impossible in general (below),
  write-through seams;
- this first, the redistribute command after.

## Investigation (2026-10-29, prototype in the session scratchpad `seam40/`; T1)

**Shared views cannot be rebuilt in general.** A view shares memory only
where the seam's axis map is an identity or a flip, in a buffer laid out by
the topology builder. A face graph read from disk may hold permuted axis
maps and cycles, for which no buffer layout exists. That is why builders
make the permuted cap–equator seams DIRICHLET copies.

**Write-through gives views' behaviour without the layout.**
- **Alias groups:** `Aliases` joins the nodes that are one point across
  SHARED and DIRICHLET seams into groups, with the point matching of ADR
  0035/0036 (`seams(...)`) and connected components. A node where several
  blocks meet is one group.
- **Sync:** after an in-place edit, each group takes the value of its last
  member the edit changed, in block order, as views would ("last writer
  wins").

**Why after the edit, and only changed members.** A first version copied
each block's changed nodes to its partners as the snap went. On a built
cubed sphere that changed wall snaps by round-off (≤ 4e-17): each DIRICHLET
side used to be projected on its own, from identical inputs to identical
outputs. Syncing once at the end, from changed members only, leaves alone
any group whose members all changed alike. Wall-only snaps of built meshes
are therefore bit for bit as before (S2: cubed sphere, O-grid, C-grid, CH,
OH).

**Prototype results:**
- one-sided associations across a SHARED seam (`sub0.j_hi`) and a
  DIRICHLET seam (`cap_north.i_lo`) leave every seam closed (gap 0), on a
  replayed project and on one loaded without replay;
- the two agree within 1e-17.

Before this ADR the DIRICHLET seam opened even on the replayed project.

## Decision

- **`shore.mesh.seams`:**
  - `Aliases(blocks, tol)` and `Aliases.sync(changed)`;
  - `seams(..., bcs=)` selects the connection BCs (SHARED and DIRICHLET
    here; PERIODIC seams do not coincide).
- **Commands:** `Snap` and `FillFace` run their edits inside
  `_seams_written_through(blocks)`.
- **Unmatched seams:** if a seam is not point-matched within 1e-9 × the
  bounding diagonal, the edit warns (`UserWarning`) and runs without
  closing seams.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Rebuild shared views on load | No buffer layout exists for permuted axis maps or cycles (above). |
| Write through block by block, during the edit | Changes built meshes' wall snaps by round-off. |
| Average the members | Not what views do; moves nodes that one side alone should set. |
| Accept the open seams | One-sided associations and loaded projects would leave gaps in the grid. |

## Consequences

- **Seams stay closed** under snap and face rebuilds, after a build or a
  load, with associations on one side.
- **ADR 0039's review question 2 is answered:** coincident copies plus
  write-through.
- **Cost:** each snap matches the seams once (the CGNS writer's point
  matching).
- **The redistribute command (next) relies on it:** moving the nodes of an
  edge must move its seam partners.

## Verification

### Acceptance spec (frozen 2026-10-29, before implementation)

`tests/spec/test_seam_sync.py` runs in the main job.

**Verdicts.**
- Required: 9 passed, 0 skipped, 0 xfailed.
- At freeze: 6 passed (S2 × 5, S4: guards), 3 xfailed.
- With the prototype behind it: 9 passed.

| Case | Required |
|---|---|
| S1 (× 2) | one-sided associations across SHARED and DIRICHLET seams: every seam closed after snap, replayed and loaded without replay; the two within 1e-12 |
| S2 (× 5) | wall-only snaps of built meshes bit-identical to projecting each block on its own |
| S3 | an unmatched seam: snap warns and still snaps |
| S4 | undo restores both sides |
