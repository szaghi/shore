# Roadmap: from pipeline to mesher

The phases are ordered by dependency, not by visibility. The GUI comes last
because it has nothing to drive until the command layer exists, and nothing
worth driving until association and the surface algorithms exist. Each phase
ends on an acceptance test, never on "feature implemented".

| Phase | Delivers | ADR | Acceptance |
|---|---|---|---|
| **0. GUI spike** (parallel, throwaway) | PySide6 + pygfx viewport, 10⁷-node render, GPU pick, worker thread | 0006 | Budget table met on Win/macOS/Linux; WSL2 result recorded |
| **1. Geometry protocol** | `shore.geometry`, `StlGeometry` (refactor), `OcctGeometry` (STEP/IGES/BREP, two-tier projection), `shore[cad]` | 0002 | Existing suite green on `StlGeometry`; STEP sphere and wing reproduce STL k=0 within facet error; projection throughput recorded |
| **2. Association** | `Association` on Vertex/Edge/Face, projection priority, fingerprint keys, `.adjacency.json` v2, propagation through `split_blocks` | 0003 | Round-trip 100 %; +5 % span edit re-matches; old sidecars load |
| **3. Command layer** | `Project`, commands, Python journal, undo cache, `shore run`; CLI re-expressed over commands | 0004 | Byte-identical outputs on all worked examples; journal replay bit-identical |
| **4. Algorithms** | surface TFI in (u, v) + elliptic smoothing with projection; true hyperbolic volume marching; hyperbolic surface (collar) marching; quality metrics module | 0005 | Frozen metric thresholds (orthogonality, concave-corner robustness) |
| **5. Exports** | CGNS (structured zones, 1-to-1 connectivity from the face graph), Plot3D | 0001 | Round-trip read by an independent reader (e.g. ParaView CGNS reader) |
| **6. GUI v1** | model tree, viewport, property panel, journal console; geometry + topology display; association editing; spacing; march with progress; quality view | 0006 | Every GUI action appears as a journal line and replays headless to identical output |
| **7. GUI v2** | interactive blocking: split, O-grid insertion, collapse, merge; vertex drag on associated entities | 0003, 0006 | A wing–body blocking built interactively, replayed from its journal after a CAD edit |

## Main risk

The main risk is not the GUI. It is phase 4, elliptic smoothing with
projection that stays robust across trimmed-surface seams and curve chains.
ICEM's reputation rests on it. Budget for it accordingly, and freeze its
acceptance metrics before any implementation starts.
