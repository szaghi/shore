# 0048 — Grid lines leaving the walls along their normal

- **Status:** Proposed (spec frozen 2026-10-29)
- **Date:** 2026-10-29

## Context

ADR 0046 measured the wing on a plane: premesh's TFI grid lines leave the
wall 30° off its normal (median, smooth region) and up to 45°. ADR 0047
showed that where the wall nodes go cannot fix it; the interior must move.

ADR 0046 listed two interior options: elliptic smoothing (option 1) and
wall-normal marching (option 3). **The user asked for all of ADR 0046's
options (night of 2026-10-29).** This ADR is the first, targeted step of
option 1: the measured defect is the lines' direction at the wall, not the
interior's smoothness.

## Investigation (2026-10-29, scratchpad `ortho48/`; T1)

**Wall-normal Hermite lines.**
- In every block with exactly one WALL face, each grid line across the wall
  becomes a cubic Hermite curve.
  - At the wall node its tangent is `strength · chord · normal`.
  - At the opposite node it keeps the line's old end direction.
- The nodes keep their arc-length fractions (64 samples per curve).
- The wall and the opposite face do not move.
- **Normals:**
  - from the wall grid itself (the cross product of its tangents), not the
    STL facets, which are piecewise constant;
  - averaged per wall node over every wall face holding it, keyed by the
    node's exact coordinates.

  Blocks sharing a node therefore share its normal, and blocks sharing a
  face build its lines identically. Seams stay closed bit for bit (gap 0).
- **Neighbours:** blocks without a wall whose faces the seams change (ADR
  0040's write-through) are re-filled by TFI.
- **Refusal:** a strength that would invert a block that was valid is
  refused, with nothing changed.

**Results:**

| Blocking | Wall angle before → after (strength 1) | Min scaled Jacobian |
|---|---|---|
| sphere, ADR 0043 hole layer | max 5.7° → 1.6° | 0.44 → 0.44 |
| sphere, ADR 0042 lattice | max 40.8° → 19.0° | 0.70 → 0.26 |
| sphere lattice, strength 2 | refused (8 corner blocks would invert) | — |
| wing on a plane, ADR 0046 route B (9 nodes across the layer); smooth region | median 29.8° → 7.8°, p95 43.7° → 14.4° | 0.44 → 0.44 |

- The wing's TE and tip edges stay at about 45°. A sharp edge has no single
  normal, and one block per wall face there cannot follow two; that needs a
  cut at the edge.
- Combined with `redistribute` (ADR 0041), every node stays on its
  orthogonal line. That gives wall clustering along wall-normal lines, the
  algebraic counterpart of marching.

**Spec pitfalls found while validating:**
- The wall angle is measured on the first chord, so it depends on the
  number of nodes across the layer: 5 nodes give a median of 13.6°, 9 give
  7.8°. The spec fixes 9.
- An angle threshold after `redistribute` did not show that the lines were
  kept: 4.8°, against 5.7° before any reshaping. The spec checks instead
  that every node stays on its line.

**Mutation the spec catches:** a zero wall tangent (no reshaping), 4 cases
fail.

## Decision

- **`shore.mesh.orthogonal`:** `grid_normals`, `wall_normals`,
  `hermite_lines`, `orthogonal_walls`.
- **Command:** `Orthogonalize` / `p.orthogonalize(strength=1.0)`, after
  premesh. It is journaled; undo restores the nodes.
- **GUI:** an *Orthogonal wall lines* form.

## Alternatives considered

| Alternative | Why not (now) |
|---|---|
| Wall nodes projected from across the block | ADR 0047: crowds nodes on sharp edges, or gives back the TFI angles. |
| STL facet normals | Piecewise constant: lines would jitter from facet to facet. |
| 3-D Winslow with Thomas–Middlecoff orthogonality control | A larger multi-block solver; it addresses interior smoothness, not yet measured as a defect. Left for a later ADR. |
| Hyperbolic marching of the wall layer | ADR 0046 option 3; its own ADR. |

## Consequences

- **Premeshed walls get wall-normal lines** in one command, valid and with
  seams closed.
- **Not done:**
  - sharp wall edges;
  - blocks with two or more walls;
  - interior elliptic smoothing.

## Verification

### Acceptance specs (frozen 2026-10-29, before implementation)

**Verdicts.**
- `tests/spec/test_orthogonalize.py` (main job): required 7 passed, 0
  skipped, 0 xfailed. At freeze: 7 xfailed. With the prototype: 7 passed.
- `tests_gui/spec/test_gui_orthogonalize.py` (GUI job): required 2 passed,
  0 skipped, 0 xfailed. At freeze: 2 xfailed. With the prototype: 2 passed.

| Case | Required |
|---|---|
| R1 | the sphere's hole layer: max wall angle > 5° before, < 3° after; Jacobian kept; walls and opposite faces fixed; seams closed |
| R2 | ADR 0042's lattice: > 40° → < 25°; valid (> 0.2); seams closed |
| R3 (× 2) | refusals: an inverting strength, 0; mesh and journal unchanged |
| R4 | the wing (route B): median > 25° before; median < 10° and p95 < 20° after; Jacobian > 0.3; seams closed |
| R5 | `redistribute` after: every node on its reshaped line |
| R6 | undo; replay bit for bit |
| G1–G2 | the GUI form; a refusal |
