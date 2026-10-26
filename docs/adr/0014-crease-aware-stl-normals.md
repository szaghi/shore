# 0014 — Crease-aware STL surface normals

- **Status:** Accepted
- **Date:** 2026-10-26

## Context

`StlGeometry.project` returns, as the surface normal, the barycentric blend
of the closest triangle's trimesh vertex normals (ADR 0002, "as in
`volume.hyperbolic._stl_normals(interpolate=True)`"). A vertex normal
averages **every** face around the vertex, so it is wrong wherever faces
meet at a sharp edge, and the error is interpolated across whole triangles.
Measured on the 2,000-4,000 sampled surface points of this spec:

| Mesh | Worst normal error today | Mean |
|---|---|---|
| box 4 × 4 × 1 (flat faces) | 54.5° | |
| cylinder, 96 sections | 45.6° | 26.4° |
| icosphere, 3 subdivisions (no sharp edge) | 0.30° | 0.10° |

The defect reached ADR 0012: on a half model, the vertices of the cut see
faces on one side only, so the half mesh's normals are not the mirror
image of the whole mesh's. The STL half model could not be required to
march mirror-symmetrically.

Every protocol consumer reads these normals: `march_surface` (ADR 0008,
0012), the intersection tracer (ADR 0009, 0011), the collar blocks
(ADR 0010, 0013) and the OCCT tessellation proxy.

## Decision

User decisions of 2026-10-26:

- **Crease-angle corner normals.** Each corner of a triangle gets the
  corner-angle-weighted mean of the normals of the faces around its vertex
  that lie within `crease_deg` (default 30°) of the triangle's own normal.
  This is trimesh's own weighting (Thürrner & Wüthrich 1998), restricted to
  one side of every sharper edge. Flat faces bounded by sharp edges get
  their exact normal; curved regions stay C0-smooth across edges, which is
  why interpolation was introduced.
- **Planar open boundaries.** The corner normals of the vertices of an open
  boundary loop that lies in a plane (the cut of a half model) are
  projected into that plane, which restores the mirror normal. Other open
  boundaries are left alone. There is no new argument.
- **`StlGeometry` only.** `volume.hyperbolic._stl_normals`, used by the
  legacy `shore mesh` marching path, is not changed, so existing meshes are
  bit-identical. ADR 0002's identity claim narrows to points, distances and
  ray casts.
- **A frozen spec.**

Interface: `StlGeometry(mesh, crease_deg=30.0)`,
`StlGeometry.from_file(path, skip_hygiene=False, crease_deg=30.0)`.
`crease_deg = 180` on a closed mesh reproduces today's normal.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_stl_normals.py` runs in the main CI job.

**Verdict.** It must report 17 passed, 0 skipped, 0 xfailed. At freeze it
reports 2 passed, 15 xfailed. The 2 that pass are regression guards that
hold before and after: the icosphere keeps its normals, and an open
boundary that is not planar keeps them. The 15 fail on the missing
`crease_deg` argument or on today's normals.

| Case | Required (vectors, 1e-12) | Prototype | Today |
|---|---|---|---|
| N1 box | face normal at every sampled point | exact | 54.5° off |
| N2 cylinder; icosphere | exact radial / axial; icosphere unchanged | exact; unchanged | 45.6° off |
| N3 `crease_deg=180` (× 3 closed meshes) | today's normal | 2e-16 | |
| N4 half post, half box | the whole meshes' normals off crease edges; cut-vertex normals in the plane | round-off; 1.5e-33 | 1.9° off |
| N5 STL half-model march (× 2) | equal to the mirrored whole-mesh march, <= 1e-6 × first step | 2e-13 box, 7e-8 post | 0.47, 3.9e-3 |
| N6 icosphere with a curved open cut | unchanged | unchanged | |
| N7 interface (× 5) | `crease_deg` in (0, 180], checked; `from_file` passes it on | | |

The post's 7e-8 in N5 is the reference's, not the half model's. trimesh
splits every side quad of the cylinder along the same diagonal, so the
whole mesh is not exactly mirror-symmetric, and its own end node drifts
7e-10 off the plane by layer 12.

**Harness validation (T1).** The prototype was run behind the frozen
interface. It passes 12 of 17 cases and fails exactly the 5 interface
cases (no argument check, no `from_file` argument).

**Downstream (T1, the whole suite with the prototype patched into
`StlGeometry`):** 1392 passed, every frozen spec included. The STL collar
of ADR 0010 (box and post) gives min scaled Jacobian 0.517 (today 0.516);
the half-collar of ADR 0013, 0.517 (today 0.512). The CAD job, which builds
`OcctGeometry` on a tessellation proxy, was not run locally; its verdict
comes from CI.

**Out of scope:**
- points on a crease edge, where the closest triangle is either side's and
  so is the normal: a curve that runs along a crease;
- feature detection beyond the dihedral angle;
- the legacy `_stl_normals` path;
- run time.

**For review (your call before implementation starts):**
1. Default crease angle of 30°. Coarse curved meshes with facet angles
   above 30° (fewer than 12 sections round a cylinder) would be split into
   flat facets. Keep 30°, or choose 45°?
2. The legacy `shore mesh` path keeps the defective normals for bit
   identity. Should a follow-up ADR switch it too, with a re-baseline of
   the marching tests?
3. Corner normals are computed once per `StlGeometry` (memory: 9 floats
   per face). Fine?

## Review decisions (2026-10-26, user)

1. **Default crease angle 30°**, as frozen.
2. **The legacy `shore mesh` path gets a follow-up ADR** with a re-baseline
   of the marching tests; until then it keeps trimesh's vertex normals.
3. **Corner normals are cached**, computed when the `StlGeometry` is built.

## Acceptance (2026-10-26): crease-aware normals spec met

`src/shore/geometry/stl.py`: `StlGeometry(mesh, crease_deg=30.0)` checks
the angle and computes `_corner_normals` once; `project` blends them;
`from_file` passes `crease_deg` on. `_planar_boundaries` finds the planar
open boundaries (trimesh row grouping, SVD, deviation <= 1e-9 × diagonal).
The corner normals are bit-identical to the prototype's on the box, the
cylinder and the icosphere.

Cost of the accumulation, which runs one vertex-face slot at a time over
blocks of 65,536 faces: on a 1.31M-face icosphere, 1.4 s and 269 MB peak
beyond trimesh's own caches. The cached result itself is 94 MB
(72 bytes per face). The prototype's fully vectorised form peaked at
606 MB.

**Verdict for `tests/spec/test_stl_normals.py`: 17 passed, 0 skipped,
0 xfailed** (T2; `scripts/check_spec.py` reports all files ok). All 8
`@PENDING` lines are removed. Full suite: 1412 passed, 0 xfailed. The CAD
job's verdict comes from CI.

**Unit tests** (`tests/test_geometry.py`):
- chunking does not change the result;
- a half cylinder's cut is found as a plane, and a closed mesh has none;
- the legacy `_stl_normals` keeps its tilt on the box.

The existing check that `project` matches `_stl_normals` on the icosphere
is renamed: it holds because no edge there is sharper than 30°.
