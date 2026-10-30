# 0049 — Feature curves and points of STL geometry

- **Status:** Proposed (spec frozen 2026-10-29)
- **Date:** 2026-10-29

## Context

ADR 0046 option 4: ICEM-style curve associations on STL. A CAD model has
edges and vertices to associate block edges and vertices with (ADR 0003,
0039). An STL had only the body, so an edge between two faces associated
with the body could not follow a sharp edge or a junction line.

**The user asked for all of ADR 0046's options (night of 2026-10-29).**

## Investigation (2026-10-29, scratchpad `wing46/feat.py`, `tecurve.py`; T1)

**Feature edges.**
- A mesh edge is a feature edge if its two facets' normals differ by more
  than the geometry's `crease_deg` (30°, the threshold of ADR 0014's crease
  normals), or if it is an open-boundary edge.
- Feature edges are chained into polylines (the CURVE entities), split at
  every vertex where other than two meet (junctions and open ends: the POINT
  entities). Closed loops without junctions start at their lowest vertex.
- Curves are ordered by their lowest vertex index, so the order is
  deterministic. After an edit, the fingerprints (centroid, length, box)
  re-find them anyway.

**Results:**
- a smooth sphere has no curves and no points, so ADR 0038's frozen display
  of an STL is unchanged;
- a box has 12 two-point curves and 8 corners;
- the blunt-TE wing on a plane (tip cap, no root cap) has 6 curves and 4
  points:
  - the two TE corner lines (21 points along the span);
  - the root and tip airfoil outlines (119 points; the root one is the open
    boundary);
  - the TE base at root and tip;
  - the 4 TE corners as points.

**Projection:**
- onto a curve: the closest point of the polyline, with NaN normals (the
  protocol's rule for curves);
- onto a point: the point.

The GUI draws and picks the curves like a CAD model's, with no GUI change.

**No quality gain on ADR 0046's wing.** Associating the hole's rear
vertical edges with the TE lines changes nothing (median wall angle 29.8°,
TE 46.7°, the same): closest-point projection already put those edges on
the TE corners. Curves matter where an edge must follow a line that
closest-point projection to the body would not find: a wing–body junction
in one STL, a plane cut, a sharp edge between faces associated with the
body.

**Spec pitfalls found while validating:**
- The two TE corner lines are 0.0025 apart, under a pixel in a whole-wing
  view. The GUI spec picks the tip outline instead.
- The out-of-range refusals (C5) already hold before the implementation,
  because the old STL refused every entity but the body. They are frozen as
  guards, not as pending cases.

## Decision

- **`StlGeometry`:**
  - `features()` (computed once);
  - `entities(CURVE | POINT)`;
  - `project` onto a curve or point;
  - `fingerprint` of a curve (`"polyline"`) or point;
  - `display()` with the curves and points.
- Out-of-range indices are refused (`ParameterError`), as before.
- `evaluate` stays refused: a polyline has no parametrisation.

## Alternatives considered

| Alternative | Why not (now) |
|---|---|
| Curves from plane sections (cut planes) | Needs the user's planes; a later command. Feature edges are intrinsic. |
| A separate crease threshold for curves | One threshold with ADR 0014's normals keeps "sharp" consistent. |
| Fit splines through feature polylines | The facets are the geometry of an STL; polylines are exact. |

## Consequences

- **STL curves and points** can be associated (`associate`,
  `associate_blocking`, the GUI), snapped and re-found after an edit, as
  CAD's are.
- **Not measured here:** a case where curves improve a blocking's quality.
  The wing–body junction (one STL with both parts) is the natural one.

## Verification

### Acceptance specs (frozen 2026-10-29, before implementation)

**Verdicts.**
- `tests/spec/test_stl_curves.py` (main job): required 7 passed, 0
  skipped, 0 xfailed. At freeze: 2 passed (C5 guards), 5 xfailed. With the
  prototype: 7 passed.
- `tests_gui/spec/test_gui_stl_curves.py` (GUI job): required 2 passed, 0
  skipped, 0 xfailed. At freeze: 2 xfailed. With the prototype: 2 passed.

| Case | Required |
|---|---|
| C1 | a smooth sphere: no curves, no points, display unchanged |
| C2 | a box: 12 curves between 8 corners, projection, NaN normals, display, distinct fingerprints |
| C3 | the wing: TE lines, root and tip outlines, 4 TE corners |
| C4 | a box edge associated with the TE line and snapped; replayed on a longer wing: on its TE line |
| C5 (× 2) | out-of-range curve and point refused (guards) |
| C6 | a vertex associated with a box corner and snapped |
| U1 | the curves drawn; a pick on the tip outline picks it |
| U2 | that curve associated with a box edge from the Association dock |
