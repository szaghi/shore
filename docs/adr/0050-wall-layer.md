# 0050 — Wall layers with an absolute first cell (algebraic, after premesh)

- **Status:** Proposed (spec frozen 2026-10-29)
- **Date:** 2026-10-29

## Context

ADR 0046 option 3 was to march the wall layer: the blocking for the far
field, SHORE's hyperbolic marching for the layer next to the wall, which is
what viscous walls need. **The user asked for all of ADR 0046's options
(night of 2026-10-29).**

A hyperbolic march of a blocking's wall needs a marcher that works on any
set of structured wall patches. Its normals must be seam-aware across
patches joined by arbitrary axis maps, with ghost rows from neighbours. Its
layers must be stitched to the blocks they replace.
- SHORE has such marchers only per topology: the cubed sphere's
  `_march_caps` and the C-grid's `_march_cgrid_three`.
- Generalising them is a large solver, and a risky one to build unattended.
  It is deferred, with this ADR as its algebraic substitute.

What a viscous wall needs from the layer:
1. lines leaving the wall along its normal (ADR 0048);
2. a prescribed first cell height;
3. smooth growth away from the wall.

(1) is done. This ADR adds (2) and (3) on top of it.

## Investigation (2026-10-29, scratchpad `layer50/`; T1)

**Geometric cells per line.**
- On every line across the wall of each block with one WALL face, the cells
  grow geometrically from the wall. The first is exactly `first_step` (an
  absolute length); the ratio comes from the line's own length (Newton,
  `topology.oh_channel._solve_geometric_ratio`).
- The nodes are placed on the line's polyline (ADR 0041's `resample`), so
  after `orthogonalize` they stay on the wall-normal curves. The wall node
  and the far node do not move.
- A first step not below a line's uniform spacing `L / (n − 1)` is refused:
  the cells would have to shrink away from the wall.

**Seams bit for bit, and why that took two fixes.** Lines on a face shared
by two wall blocks are computed in each block, so their results must not
depend on the batch:
- **Newton's iteration count:** a vectorised solve iterates until every
  entry converges, so the count depends on the batch, and the same line
  gave different last bits in two blocks. Fixed by solving one line at a
  time.
- **Summation order:** `sum(axis=1)` over a strided view sums in a
  layout-dependent order. Fixed with a sequential `cumsum` on a contiguous
  copy.

Seam gaps went from 2–6e-16 to 0. The command also runs inside ADR 0040's
write-through, as a backstop.

**Results:**

| Case | First step | Growth ratio | Min scaled Jacobian | Wall angle |
|---|---|---|---|---|
| sphere hole layer, 17 nodes across, orthogonalized | 1e-2 | 1.14–1.27 | 0.44 (kept) | max 1.9° |
| same | 1e-3 | 1.39–1.52 | 0.44 | max 1.9° |
| wing on a plane (route B), 25 nodes across, orthogonalized | 1e-4 | 1.41–1.45 | 0.44 | median 3.3°, p95 8.9° |
| ADR 0042's lattice (no layer) | 1e-3 | — | refused: 8 corner blocks would invert | — |

- **On the lattice,** the face blocks' clustered side faces are shared
  with edge and corner blocks. Those have no wall and are re-filled by TFI
  between faces of different spacing, which inverts them.
- **Wall layers need a layer topology** (an O-grid around the body, ADR
  0043), and the refusal enforces it.

**Mutation the spec catches:** clustering at the far end of lines whose
wall is on the high side of their block (3 cases fail).

## Decision

- **`shore.mesh.orthogonal`:** `layer_fractions`, `wall_layer`.
- **Command:** `WallLayer` / `p.wall_layer(first_step)`. It shares
  `Orthogonalize`'s frame (`_reshape_wall_blocks`): write-through, re-fill,
  and refusal on inversion.
- **GUI:** a *Wall layer* form.

## Alternatives considered

| Alternative | Why not (now) |
|---|---|
| Hyperbolic marching of the blocking's wall patches | Needs a general seam-aware multi-patch marcher; deferred (see Context). |
| `redistribute` (ADR 0041) with a `pinned` first step | Fractions of the picked edge's length: the first cell differs from line to line. |
| A growth-ratio limit (e.g. 1.2) with the layer thickness free | Needs the layer's node count chosen from it; the blocking fixes the count first. Later. |

## Consequences

- **A viscous-ready wall layer from a blocking:**
  `premesh → orthogonalize → wall_layer(ds)`. Lines leave the wall along
  its normal, with the first cell prescribed (wing: median 3.3°).
- **Not done:**
  - hyperbolic marching;
  - growth-ratio limits;
  - blocks with several walls;
  - wall layers without a layer topology (refused).
- **The roadmap's wing–body acceptance** now has the blocking and the wall
  quality it needs. What remains is a wing–body case itself: one STL or CAD
  model with the junction as a curve (ADR 0049).

## Verification

### Acceptance specs (frozen 2026-10-29, before implementation)

**Verdicts.**
- `tests/spec/test_wall_layer.py` (main job): required 6 passed, 0 skipped,
  0 xfailed. At freeze: 6 xfailed. With the prototype: 6 passed.
- `tests_gui/spec/test_gui_wall_layer.py` (GUI job): required 2 passed, 0
  skipped, 0 xfailed. At freeze: 2 xfailed. With the prototype: 2 passed.

| Case | Required |
|---|---|
| L1 | the sphere's hole layer: first cells 1e-3 (rtol 1e-9), growing, ends fixed, seams closed, Jacobian kept, wall angle < 3° |
| L2 | the wing: first cells 1e-4, Jacobian > 0.3, seams closed, median wall angle < 5°, p95 < 12° |
| L3 (× 3) | refusals: 0, above a line's uniform spacing, a lattice without a layer |
| L4 | journal, replay, undo |
| W1–W2 | the GUI form: first cells; a refusal |
