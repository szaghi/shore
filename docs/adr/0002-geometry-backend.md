# 0002 — Geometry protocol with STL and OCCT backends

- **Status:** Proposed
- **Date:** 2026-10-25

## Context

SHORE's only geometry today is a `trimesh.Trimesh` loaded from STL, queried
through `Projector` (cached BVH ray-cast from an interior anchor). STEP/IGES
input needs an OpenCASCADE kernel. The candidate is **build123d** (0.13.0,
Apache-2.0), which depends on **cadquery-ocp** (8.0.1, i.e. OCCT 8.0;
LGPL-2.1 with exception; wheel ≈ 67 MB).

build123d is a CAD *authoring* API. Meshing needs a geometry *query* layer:

| Need | Where it lives |
|---|---|
| STEP import, names/colours via XCAF | build123d `import_step` (wraps `STEPCAFControl_Reader`) |
| IGES import | OCP `IGESControl_Reader`; build123d 0.13 has no IGES importer (to verify) |
| topology traversal, adjacency | build123d `Shape.faces()/edges()/vertices()` |
| surface evaluation at (u, v), normals | OCP `BRepAdaptor_Surface`, `BRepLProp_SLProps` |
| point → surface / curve projection | OCP `GeomAPI_ProjectPointOnSurf`, `GeomAPI_ProjectPointOnCurve`, `ShapeAnalysis_Surface` |
| closest point to a shape | OCP `BRepExtrema_DistShapeShape` |
| healing, sewing | OCP `ShapeFix_Shape`, `BRepBuilderAPI_Sewing` |
| tessellation with UV per node | OCP `BRepMesh_IncrementalMesh` + `BRep_Tool.Triangulation_s` |

**Performance constraint.** OCP calls go one point at a time through Python.
At a few µs per call, 10⁵ projections take about 1 s. Elliptic surface
smoothing re-projects every node on every sweep, so 100 sweeps × 10⁵ nodes
means about 100 s per surface. Interactive vertex dragging must answer in
under 16 ms. Per-point exact projection is fine for a final snap and too slow
for iterative or interactive work. Writing a compiled batch projector is
excluded (ADR 0006 sets Python + numpy only).

## Decision

1. **A `Geometry` protocol** in a new `shore.geometry` package is the only
   geometry API the mesh core sees:

   ```python
   class Geometry(Protocol):
       def entities(self, kind: EntityKind) -> Sequence[EntityRef]: ...
       def bbox(self) -> np.ndarray: ...                       # (2, 3)
       def ray_cast(self, origins, dirs) -> RayHits: ...       # existing topologies
       def project(self, points, entity: EntityRef | None, *, exact: bool = False
                   ) -> Projection: ...                        # xyz, uv|t, dist, status
       def evaluate(self, entity: EntityRef, params) -> np.ndarray: ...
       def normal(self, entity: EntityRef, uv) -> np.ndarray: ...
   ```

   All arguments and results are numpy arrays (batched). `entity=None`
   means "the whole body" (today's semantics).

2. **Two backends.**
   - `StlGeometry`: today's `trimesh` + `Projector`, unchanged behaviour. It
     has a single entity, the body.
   - `OcctGeometry`: STEP/BREP via build123d, IGES via OCP, healed with
     `ShapeFix` on load.

3. **Two-tier projection in `OcctGeometry`.** Each face and edge is tessellated
   once (`BRepMesh_IncrementalMesh`), keeping (u, v) or t at every
   tessellation node. Then:
   - **`exact=False` (default, vectorised):** BVH closest-point on the
     per-entity tessellation (trimesh, already a dependency). Parameters are
     interpolated barycentrically from the nodes' (u, v). Used inside
     smoothing loops and for interactive dragging.
   - **`exact=True`:** a per-point OCP Newton projection, seeded from the
     tier-1 (u, v). Used once at the end of each surface algorithm.

   The tessellation deflection is a geometry setting, recorded in the
   project. Tier-1 error is bounded by it, and tier 2 removes it.

4. **Optional dependency.** `shore[cad]` = `build123d`. The core imports
   `shore.geometry.occt` lazily. Without the extra, STEP input raises a
   `MeshInputError` that names the extra.

5. **`shapes/` migrates to build123d** when `shore[cad]` is installed.
   Analytic bodies become exact B-reps. The trimesh builders stay as the
   fallback that exports STL.

## Alternatives considered

- **Import build123d/OCP throughout the core.** Rejected: it makes a
  67 MB+ dependency mandatory, ties every test to OCCT, and lets the STL path
  rot.
- **Tessellate STEP to STL and keep the STL pipeline.** Rejected: SHORE would
  lose exact geometry, curve/surface identity (needed for association, ADR
  0003) and (u, v) parameters (needed for parametric-space surface meshing).
- **Vectorised numpy NURBS evaluator** (`BRepBuilderAPI_NurbsConvert`, then
  batched de Boor in numpy). Exact and fast for evaluation, but projection
  still needs Newton iterations and trimming loops. Kept as a later
  optimisation for tier 1 if the tessellation proxy proves too coarse.
- **Gmsh's OCC kernel as the backend.** Rejected: it is reachable only
  through Gmsh's model API and tags, and it adds a second large native
  dependency for no gain over OCP.
- **ESP/EGADS.** Its persistent attribution is the best precedent for ADR
  0003, but it is not on PyPI and needs a native build per platform.
  Studied, not adopted.

## Consequences

- `Projector`, `Body` and the surface builders take a `Geometry`, not a
  `trimesh.Trimesh`. `Body.from_stl` keeps its signature and builds a
  `StlGeometry`.
- The star-shaped restriction stays only on the centroid-ray-cast templates.
  Association-driven projection (ADR 0003) has no such restriction.
- CI gains a `cad` job that installs `shore[cad]`. The default job stays light.

## Verification

- **Spike:** `OcctGeometry` loads a STEP sphere and a STEP NACA 0012 wing. The
  existing cubed-sphere and C-grid templates run on them through `ray_cast`.
  Their k=0 surfaces match the STL runs to within the STL facet error.
- **Benchmark (recorded in the PR):**
  - tier-1 projection reaches at least 10⁶ points/s;
  - tier-2 throughput is measured;
  - tier-1 vs tier-2 distance error is no more than the tessellation
    deflection.
- **Check:** IGES support in build123d 0.13. If absent, wrap
  `IGESControl_Reader` directly.
