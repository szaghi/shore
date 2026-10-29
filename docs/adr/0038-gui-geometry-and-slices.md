# 0038 — GUI slice 2: geometry display and slice planes

- **Status:** Accepted
- **Date:** 2026-10-29

## Context

GUI slice 1 (ADR 0037) shows the blocks' boundaries. Inspecting a grid also
needs:
- the **geometry** it was built on (STL bodies, and CAD models with their
  faces, edges and vertices). Association editing, the next slice, picks
  those entities.
- **interior planes** (ADR 0006: "interiors on demand", constant-i/j/k
  planes per block), coloured by quality where wanted.

User decisions of 2026-10-29:
- **This slice:** geometry display plus slices; inspection only.
- **View state** (what is visible, which slices) is not a command, and is
  not journaled.

## Investigation (2026-10-29, prototype in the session scratchpad `gui38/`; T1)

**What a geometry offers to draw.**
- `OcctGeometry` already keeps:
  - its tier-1 tessellation with the owning face of every triangle
    (`_tri_face`);
  - a polyline per edge (`_edge_polys`);
  - its vertices.
- `StlGeometry` is a single body.

A public `display()` on both backends returns them as one `Display`
(float64 points, triangles, each triangle's entity, curves, vertices), so the
GUI never reads backend internals.

**Drawing and picking.**
- Triangles are one pickable mesh per geometry. pygfx returns the triangle
  and its barycentric coordinates, so the pick reports the triangle's
  entity, and a point computed in float64.
- Curves are one line per geometry (polylines separated by NaN), points one
  point cloud. A pick's vertex index maps back to the curve or point.
- Edges and vertices lie on the faces' triangles, so with depth testing a
  pick there hits the face. They are drawn over the faces, as CAD viewers
  do, which keeps them pickable.
- A slice plane is a quad mesh of `nodes[index, :, :]` (or j, k). A quad
  maps to its cell on the plane's + side (the − side on the last plane).

**Spec pitfalls found while validating** (none was a prototype defect):
- picks on cells hidden behind other faces, or behind the curved slice
  itself, legitimately return the front object;
- straight CAD edges are sampled with two points, which are vertices;
- the default view direction (−1, −1, −1) sees axis-aligned faces at
  cos 0.577.

The spec now picks only points that a CPU ray cast (trimesh) finds
unoccluded, and uses polyline midpoints.

**Headless prototype run.** STL picks returned the body at the picked
point (7/7). A STEP box with a hole showed 7 faces, 15 edges and 10
vertices, and its face, edge and vertex picks were right. Slice picks
returned the slice, the cell and the exact centre. Quality coloured
slices. The journal was unchanged by every view action.

## Decision

- **Core:**
  - `shore.geometry.base.Display`;
  - `StlGeometry.display()`: one body;
  - `OcctGeometry.display()`: faces, edges, vertices.
- **Viewport:**
  - `set_geometries`, `set_slices` (`SliceKey(block, axis, index)`);
  - `show(kind, visible)` for blocks, geometries and slices;
  - `pick` returns `GeometryPick(geometry, EntityRef, xyz)` or a
    `PickResult` for a slice ("k=2");
  - quality (`set_face_scalars`) covers slices;
  - the float64 rebase covers everything drawn.
- **Window:**
  - a **Geometry** tab beside the model tree, one row per geometry (its
    triangles, or its faces, edges and vertices);
  - a **Display** dock: visibility of blocks, geometry and slices; slices
    added by block, axis and index, listed and removed;
  - the status bar names the picked entity;
  - slices of blocks that disappear (undo, a rebuild) are dropped.
- **View state is not a command.** It changes neither the journal nor the
  project, and it is not saved in this slice.
- **CI:** the `test-gui` job installs `.[dev,gui,cad]`. The CAD cases need
  build123d.

## Alternatives considered

| Alternative | Why not |
|---|---|
| The GUI reading `OcctGeometry` internals | Couples the GUI to the backend's private arrays. |
| Edges and vertices depth-tested | They would be unpickable where they lie on faces. |
| A geometry row in the model tree | ADR 0037's spec pins the tree's top level to the components. |
| Journaled view state | User decision: display settings are not project data. |

## Consequences

- **Association editing** (next) can pick CAD faces, edges and vertices,
  and block faces.
- **Hidden edges show through faces.** All of a CAD model's edges are
  visible, including the hidden ones.
- **Not in this slice:** saved view state, per-entity colours and
  visibility, transparency, a cell-quality threshold cloud.

## Verification

### Acceptance spec (frozen 2026-10-29, before implementation)

`tests_gui/spec/test_gui_inspect.py` runs in the `test-gui` job.

**Verdicts.**
- Required: 12 passed, 0 skipped, 0 xfailed.
- At freeze: 12 xfailed.
- With the prototype behind it: 12 passed.

| Case | Required |
|---|---|
| V1 | every project geometry shown and listed (`1280 triangles`; `7 faces, 15 edges, 10 vertices`); uploads = display points |
| V2 | STL picks (blocks hidden): the body, within 2 px of the picked point |
| V3 | CAD face picks: the owning face, within 2 px; at least 3 faces |
| V4 (× 2) | CAD edge and vertex picks: the right CURVE and POINT, at least 3 each |
| V5 (× 2) | slice picks: slice key, cell, exact centre; unknown block and out-of-range index refused |
| V6 | quality colours on a slice = colour map of the scaled Jacobian of its cells |
| V7 (× 2) | view state changes neither journal nor history; slices of removed blocks dropped |
| V8 (× 2) | `display()` of STL and OCCT: consistent arrays, one curve per edge, one point per vertex |

## Acceptance (2026-10-29): spec met

**Implementation.** As decided, plus:
- unit tests: the slice cell and value round trip, the last plane, the
  Display dock's buttons and toggles (`tests_gui/test_gui_units.py`), and
  `StlGeometry.display()` (`tests/test_geometry.py`, main suite);
- `docs/reference/gui.md`.

**Verdict** (T2 locally, offscreen Qt with llvmpipe; `scripts/check_spec.py`
reports all files ok): `tests_gui/spec/test_gui_inspect.py` gives **12
passed, 0 skipped, 0 xfailed**, and slice 1's spec still gives 16 passed.
