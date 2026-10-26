# Roadmap: from pipeline to mesher

The phases are ordered by dependency, not by visibility. The GUI comes last
because it has nothing to drive until the command layer exists, and nothing
worth driving until association and the surface algorithms exist. Each phase
ends on an acceptance test, never on "feature implemented".

| Phase | Delivers | ADR | Acceptance |
|---|---|---|---|
| **0. GUI spike** (parallel, throwaway) | PySide6 + pygfx viewport, 10⁷-node render, GPU pick, worker thread | 0006 | Budget table met on Win/macOS/Linux; WSL2 result recorded |
| **1. Geometry protocol** ✅ | `shore.geometry`, `StlGeometry` (refactor), `OcctGeometry` (STEP/IGES/BREP, two-tier projection), `shore[cad]`, consumers migrated (`Body.from_file`, `Mesh.from_file`, CLI) | 0002 | Met 2026-10-26: suite green on `StlGeometry`; STEP sphere k=0 exact and within the STL facet sag; STEP wing C-grid equals STL to 1e-8 (`tests/test_cad_pipeline.py`); throughput recorded (10⁶ target missed, see ADR 0002) |
| **2. Association** ✅ | `Association` on Vertex/Edge/Face, projection priority, fingerprint keys, `.adjacency.json` v3, propagation through `split_blocks` | 0003 | Met 2026-10-26: round-trip 100 %; +5 % span edit re-matches every face and edge (`changed`), +30 % reports losses; v2 sidecars load (`tests/test_association.py`). Open: C-grid TE edges → TE curve association |
| **3. Command layer** ✅ | `Project`, commands, Python journal, undo cache, `shore run`; CLI re-expressed over commands | 0004 | Met 2026-10-26: `shore mesh` 68/68 golden files byte-identical (all five topologies); `primitive *`, `cap-tip`, `split`, `grd`, `cc-par`, `proc-input` 66/66 byte-identical (the chimera-assembly example end to end, all five primitives, split, both tip caps); journal replay byte-identical for each, and one journal reproduces a whole Xall component (`tests/test_project.py`, `tests/test_project_pipeline.py`) |
| **4. Algorithms** (in progress) | surface TFI in (u, v) + elliptic smoothing with projection; true hyperbolic volume marching; hyperbolic surface (collar) marching; quality metrics module | 0005 | Frozen metric thresholds (orthogonality, concave-corner robustness). 2026-10-26: `shore.quality` done (hand-verified, `tests/test_quality.py`); hyperbolic-marching spec frozen (`tests/spec/`, 13 cases); hyperbolic scheme implemented; spec re-frozen (body-normal metric, per-case `@PENDING`): hyperbolic marching accepted: 13 passed, 0 skipped, 0 xfailed (seam-free explicit dissipation, bounded outer-face volume correction; ADR 0005). Surface TFI + elliptic spec frozen (ADR 0007, `tests/spec/test_surface_tfi*.py`, 23 + 4 cases); accepted: 23 + 4 passed, 0 skipped, 0 xfailed (`surface_grid` kernel, `Project.surface_grid` command). Collar grids: surface-marching kernel accepted, 23 passed, 0 skipped, 0 xfailed (ADR 0008, `shore.surface.march`); intersection loop accepted, 14 passed, 0 skipped, 0 xfailed after the review re-freeze (ADR 0009, `shore.surface.intersect`); review questions of ADR 0007-0009 decided 2026-10-26; collar block accepted, 26 passed, 0 skipped, 0 xfailed (ADR 0010, `shore.collar`, `Project.build_collar`: hyperbolic near the wall, trimmed-offset envelope, Winslow interior; open curves, end constraints, loop spacing, 3-D offset for thick collars and multi-component projects deferred); open intersection curves on a symmetry plane: spec frozen (ADR 0011, `tests/spec/test_open_curve.py`, 20 cases, awaiting review; marcher end constraint and half-collar block next). Order decided: hyperbolic volume marching, then surface TFI + elliptic smoothing in (u, v) per face (ADR 0002 throughput decision), then collar grids |
| **5. Exports** | CGNS (structured zones, 1-to-1 connectivity from the face graph), Plot3D | 0001 | Round-trip read by an independent reader (e.g. ParaView CGNS reader) |
| **6. GUI v1** | model tree, viewport, property panel, journal console; geometry + topology display; association editing; spacing; march with progress; quality view | 0006 | Every GUI action appears as a journal line and replays headless to identical output |
| **7. GUI v2** | interactive blocking: split, O-grid insertion, collapse, merge; vertex drag on associated entities | 0003, 0006 | A wing–body blocking built interactively, replayed from its journal after a CAD edit |

## Main risk

The main risk is not the GUI. It is phase 4, elliptic smoothing with
projection that stays robust across trimmed-surface seams and curve chains.
ICEM's reputation rests on it. Budget for it accordingly, and freeze its
acceptance metrics before any implementation starts.
