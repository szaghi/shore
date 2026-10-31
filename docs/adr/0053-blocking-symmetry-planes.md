# 0053 — Symmetry planes of a blocking

- **Status:** Accepted
- **Date:** 2026-10-30

## Context

ADR 0051's half wing-body has the plane z = 0 as its symmetry plane. Its
faces on the plane come out of `premesh` on the plane to 4e-16, but as
`FREE` faces. `orthogonalize(0.5)` then moves them up to **0.026** off it,
in the O-grid layer blocks round the fuselage (b75–b83, `k_lo`), and
`wall_layer` keeps the error. ADR 0051's specs never checked those faces.
This was found while drawing the documentation figures (night of
2026-10-30). **The owner asked for the fix with an ADR and a frozen spec,
including a SYMMETRY notion in blocking (2026-10-30).**

The cause: a wall node on the plane is a boundary node of the wall's
surface grid.
- `np.gradient` takes a one-sided difference there.
- A half model has no mirror block to average with.

So the wall normal tilts out of the plane. The Hermite line leaves along
it, and the line, which lies in the plane face, leaves the plane.

## Decision

1. **Declare the plane.** `p.symmetry_plane(point, normal)`, the command
   `SymmetryPlane` (journaled, undoable), adds the plane (point, unit
   normal) to the current blocking: `Blocking.planes`, copied with it. A
   quarter model declares two.
   - Refused (`ParameterError`, nothing changed): a zero normal, a point or
     normal that is not three numbers, no blocking.
   - This mirrors the collars' `plane=` (ADR 0012, 0013): the plane is a
     user statement, not a guess from the geometry.
2. **Premesh** makes every boundary face whose four skeleton corners lie on
   a declared plane `SYMMETRY`, if it has no association (otherwise it is a
   `WALL`, as before).
   - Tolerance: 1e-9 of the skeleton's extent.
   - CGNS writes such faces as `BCSymmetryPlane`, `cc.par` as natural code
     2 (ADR 0035, 0036), the face graph as `SYMMETRY`.
3. **Orthogonalize** keeps the plane. The normal of every wall node lying on
   a declared plane is projected into the plane, before the lines are built.
   - The normal comes from the declaration, not from each block, and is
     keyed by the node's exact coordinates (ADR 0048). Blocks sharing the
     node use one normal, so seams stay closed bit for bit.
   - Each line's start tangent, end tangent and end nodes are then in the
     plane, so the whole line is.
   - `wall_layer` resamples each line on its own polyline, so it stays
     there too.
4. **Without a declared plane nothing changes**, bit for bit: the frozen
   ADR 0051 journal and grids are untouched. Projects that need the plane
   add one journal line after `create_blocking`.

## Alternatives considered

- **Detect planar FREE faces automatically:** any flat outer face of a box
  would become a symmetry plane, and the far-field faces of a box are
  planar too.
- **Snap the nodes back onto the plane after `orthogonalize`:** this hides
  the cause. The lines would still leave the wall tilted, and the snap
  would move nodes off their lines.
- **Fit each SYMMETRY face's plane from its nodes:** round-off differences
  between faces would break the bit-for-bit seams.

## Consequences

- The wing-body tutorial and the GUI tutorial declare the plane; the GUI has
  a *Symmetry plane* form.
- ADR 0051's grids without the declaration keep the drift: the spec of
  ADR 0051 is frozen, and the fix is opt-in by design.
- Not done: moving vertices onto the plane (`move_vertex` with a plane), a
  body that crosses the plane, and associated faces on the plane.

## Verification

`tests/spec/test_symmetry_plane.py` (CAD job), frozen before the
implementation, on ADR 0051's case with the plane declared:
- Y1: the faces on z = 0 are SYMMETRY, and only they.
- Y2: on the plane to 1e-12 after `orthogonalize` and `wall_layer`.
  Measured: 3.7e-16 (the premesh round-off, 28 SYMMETRY faces, min
  scaled Jacobian 0.173 as in ADR 0051), against 0.026 without the
  projection, which a mutation check confirmed.
- Y3: ADR 0051's quality holds.
- Y4: CGNS `BCSymmetryPlane` and the face graph.
- Y5: journal, replay bit for bit, undo.
- Y6: refusals.
- Y7: without the plane no SYMMETRY face (a guard).

At freeze: 1 passed, 6 xfailed.
