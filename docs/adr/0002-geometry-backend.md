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
   @runtime_checkable
   class Geometry(Protocol):
       @property
       def bounds(self) -> np.ndarray: ...                     # (2, 3)
       def entities(self, kind: EntityKind) -> tuple[EntityRef, ...]: ...
       def ray_cast(self, anchor, dirs, on_failure="raise"
                    ) -> ProjectionResult: ...                 # existing templates
       def project(self, points, entity: EntityRef | None = None, *,
                   exact: bool = False) -> Projection: ...     # points, dist, normals, params
       def evaluate(self, entity: EntityRef, params) -> np.ndarray: ...
   ```

   All arguments and results are numpy arrays (batched), and leading
   dimensions are preserved. `entity=None` means "the whole body" (today's
   semantics). Two choices made at implementation (2026-10-26), relative to
   the first draft of this record:
   - **`ray_cast` takes one interior anchor**, with `Projector` semantics
     and per-ray status. Every existing template casts that way. A
     many-origin cast is added only when a caller needs one.
   - **Normals come back in `Projection`** instead of from a separate
     `normal(entity, uv)` call. STL has no (u, v), and every caller wants
     the normal at a projected point.

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
   project. For queries near the surface, tier-1 error is a small multiple
   of it, and tier 2 removes it. See Results below for the measured bound.

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

  *As implemented (2026-10-26):*
  - **The k=0 builders keep their `projector` argument.** `cap`,
    `equator`, `remesh` and the topologies only call
    `.project(dirs)`/`.anchor`. For CAD, `Body.projector` is an
    `AnchoredCaster` that ray-casts on the tessellation and **snaps** hits
    to the exact surface. For STL it is still the legacy `Projector`, so
    results are bit-identical.
  - **Marcher normal blending** goes through `Body.normals` →
    `geometry.project(...).normals`.
  - **`Body.from_file`, `Mesh.from_file` and the `mesh`/`project` CLI
    commands** read CAD.
  - **Remaining trimesh consumers** (anchor resolution, C-grid slicing and
    TE detection) read `Body.mesh`, which for CAD is the merged,
    degenerate-free (watertight) tessellation. They run at tessellation
    accuracy until a CAD-native section replaces them.
- The star-shaped restriction stays only on the centroid-ray-cast templates.
  Association-driven projection (ADR 0003) has no such restriction.
- CI gains a `cad` job (`test-cad`) that installs `shore[cad]`. The default job stays light.

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

## Results (2026-10-26, `OcctGeometry` in `src/shore/geometry/occt.py`)

Measured on the development machine (WSL2, Python 3.12, OCCT 8.0 via
cadquery-ocp 8.0.1, build123d 0.13.0), default deflection 10⁻³ × the bbox
diagonal. Coverage is `tests/test_geometry_occt.py` (24 tests, run by the
`test-cad` CI job).

| Item | Result |
|---|---|
| IGES | build123d 0.13 has **no** IGES importer. IGES goes through OCP `IGESControl_Reader` and is then **sewn** (`BRepBuilderAPI_Sewing`), because IGES stores faces as independent trimmed surfaces: a box read back has 24 edges unsewn and 12 sewn |
| Tier-2 accuracy | Exact: sphere radius error 2·10⁻¹⁶; at a 10³ offset, 6·10⁻¹⁴; box faces exact. Outward normals hold on REVERSED faces. A foot outside a trimmed face falls back to `BRepExtrema_DistShapeShape` and lands on the face boundary |
| Tier-1 accuracy | Near the surface, error ≤ ~2× deflection. BRepMesh bounds edge sag by the deflection, but triangle interiors sag up to ~1.7× further (sphere: 0.006 at δ = 0.0035). For a query at distance h the error grows by about h × the facet normal error: 0.045 at h = 1 on the unit sphere. **Distant queries must use `exact=True`** |
| Tier-1 (u, v) | Interpolated per node. At singular parameters (sphere poles, cone apexes, detected as a vanishing ∂S/∂u or ∂S/∂v), the singular corners are excluded from that parameter's interpolation. Without that, pole errors reached 0.06 |
| **Tier-1 throughput** | **2.5·10⁴–1.3·10⁵ points/s** near the surface (trimesh `closest_point` with R-trees built once at load). **The 10⁶ points/s target is not met.** A vectorised custom index (KD-tree + certified candidate sets) was prototyped and dropped: with a warm R-tree trimesh was as fast or faster. An earlier "8·10³ points/s" figure was an artefact of timing cold R-tree builds |
| Tier-2 throughput | 1.4·10³–4·10⁴ points/s (per-point OCP Newton + face classifier; slowest on the filleted box, where distant queries hit the trimmed-face fallback) |
| Load | Box 0.04 s, sphere 0.05 s, box with all 12 edges filleted (26 faces) 0.25 s |

Consequence for phase 4. Elliptic surface smoothing that re-projects 10⁵
nodes per sweep for 100 sweeps (10⁷ queries) takes ~100 s at these rates.
Before phase 4, do one of the following and measure it here:
- smooth in (u, v) per face and evaluate only at the end;
- project only nodes that moved more than a tolerance;
- accept a compiled kernel in the core (the "Python + numpy only" rule of
  ADR 0006 binds the GUI, not the mesher).
