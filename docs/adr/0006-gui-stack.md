# 0006 — GUI: PySide6 + pygfx/wgpu, desktop, Python + numpy only

- **Status:** Proposed
- **Date:** 2026-10-25

## Context

The GUI must provide ICEM/Pointwise-class interaction:
- geometry display;
- block topology editing;
- association editing (ADR 0003);
- spacing control;
- marching with progress;
- quality inspection.

It is a front end over the command layer (ADR 0004), never a second owner of
state.

### Constraints

Agreed 2026-10-25, with untestable adjectives replaced by checkable
criteria:

| Requirement | Criterion |
|---|---|
| Portable | Runs on Windows 10+, macOS 13+ (Apple Silicon and x86), Linux (X11, Wayland) from `pip install shore[gui]`, with no system packages and no compiler |
| Minimal dependencies | GUI extra ≤ 100 MB of wheels; no dependency that is unmaintained or single-platform |
| Modern | Native GPU API per OS (Metal / D3D12 / Vulkan); HiDPI-correct; dark and light themes following the OS; dockable panels |
| Fast, responsive | See the performance budget below |
| Scale | Up to 10⁷ grid nodes interactive |
| Deployment | Desktop only; no remote/browser mode in v1 |
| Code | Python + numpy only; no compiled extensions of our own |

### Facts behind the choice

Measured 2026-10-25 (PyPI, linux x86_64 wheels):

| Package | Version | Wheel |
|---|---|---|
| cadquery-ocp (already required by `shore[cad]`) | 8.0.1.1 | 67 MB |
| vtk | 9.7.1 | 140 MB |
| PySide6-Essentials | 6.12.0 | 84 MB |
| imgui-bundle | 1.93.0 | 14 MB |
| wgpu + pygfx + rendercanvas | 0.32 / 0.17 / 2.7 | ≈ 5 MB |

macOS has deprecated OpenGL since 10.14 (frozen at 4.1). Every
OpenGL-rendering stack therefore rests on a platform API Apple is removing.
That includes VTK and pyvistaqt (today's `vis.py`), vispy and Dear PyGui.
wgpu (wgpu-native) renders through Metal, D3D12 and Vulkan.

## Decision

### Stack

- **UI toolkit: PySide6-Essentials** (LGPL, official Qt binding). It provides
  docking, tree and property widgets, native file dialogs, HiDPI support and
  OS colour-scheme tracking (`QStyleHints.colorScheme`). The style is Fusion
  with a SHORE palette for dark and light themes. Icons are a bundled SVG set
  under a permissive licence (e.g. Lucide, ISC). No QML, no Qt3D.
- **Renderer: pygfx on wgpu**, embedded through `rendercanvas`' Qt widget.
- **`shore[gui]`** = `PySide6-Essentials`, `pygfx`, `wgpu`, `rendercanvas`:
  about 90 MB.
- Entry points: `shore gui` (CLI) and a `gui-scripts` entry point
  `shore-gui`, so Windows launches without a console window.

### Isolation of the pre-1.0 renderer

pygfx is pre-1.0, and wgpu-py states its API will change with WebGPU. All
rendering goes through a SHORE-owned **`Viewport` interface** in
`shore.gui.viewport` with these operations:
- add, update and remove block surfaces, edges and slices;
- set scalars and colour map;
- camera control;
- pick → `PickResult`;
- screenshot.

Only `shore.gui.viewport.pygfx_backend` imports pygfx or wgpu. Upgrading
pygfx touches one module. A VTK backend remains possible behind the same
interface if pygfx is ever abandoned.

### Rendering architecture for 10⁷ nodes

1. **Never upload volume interiors by default.** Each block uploads its six
   boundary faces as indexed triangle meshes and its twelve edges as line
   segments. A 10⁷-node assembly in ~20 blocks has about 10⁶ surface nodes,
   which is cheap for the GPU.
2. **Interiors on demand**, as constant-i/j/k planes (O(ni·nj) upload per
   plane, natural for structured grids) and as cells filtered by a quality
   threshold (e.g. Jacobian < ε) drawn as a point or cell cloud.
3. **Float precision: camera-relative rebase.** The CPU keeps float64.
   With a 10³ m domain and 10⁻⁶ m wall spacing, float32 absolute resolution
   is about 6·10⁻⁵ m, so boundary-layer cells collapse on screen and picking
   is wrong. Spike results, at a fixed wall-cell aspect ratio of 10⁴:

   | Scheme | Result |
   |---|---|
   | Per-block centroid offset in the world transform | FAIL. pygfx applies `world_transform` and then `cam_transform` on the GPU in float32, so the large translation cancels in float32 |
   | Floating origin, one origin per block | Holds while `ds_wall / extent >~ 1e-8`. For y⁺ ≈ 1, `ds/L ≈ 1/(Re·sqrt(Cf/2))` is about 4·10⁻⁸ at Re = 10⁹ and 5·10⁻⁹ at Re = 10¹⁰, so large full-scale cases sit at or past the limit |
   | Per-tile local origins (quadtree) | Thousands of meshes and seconds of setup at `ds/extent` ≤ 10⁻⁹, and still misses picks. Rejected |
   | **Rebase** | Exact (18/18 picks) for `ds/extent` from 10⁻⁷ to 10⁻¹⁰, one mesh per block, 3.8 ms CPU for 10⁷ nodes |

   Decision:
   - Each block keeps its float64 surface positions on the CPU.
   - When the floating world origin moves, the GPU positions are re-uploaded
     as `float32(p − origin)`.
   - The origin follows the camera target. A rebase is triggered only when
     the target drifts more than ~10³ view widths from the origin, the bound
     beyond which visible geometry loses sub-pixel accuracy, or on explicit
     "look at" actions.
   - Reported coordinates always come from the float64 arrays, never from
     the GPU.

4. **Picking via the GPU pick buffer.** It returns (object, primitive index).
   SHORE maps that to (block, i, j, k) in O(1) from the structured index
   arithmetic, so no CPU BVH is needed for grid picking. Geometry picking
   uses the same buffer on the tessellation.
5. **Quality metrics** (Jacobian, skewness, aspect ratio, orthogonality,
   growth ratio) are computed vectorised in numpy on a worker thread and
   uploaded as per-vertex scalars.

### Threading

Qt's main thread owns widgets and rendering only. Commands run in a
`QThreadPool` worker, with the ADR 0004 progress and cancel hooks marshalled
back by signals. Numpy releases the GIL in its heavy kernels, so a march
keeps the UI responsive (spike: max UI lateness 1.8 ms). pybind11 bindings,
OCP included, do **not** release the GIL by default, and a single long call
freezes the UI for its whole duration (spike: 705 ms). Therefore:
- long OCP operations (STEP import, healing, `BRepMesh_IncrementalMesh`) run
  in a **subprocess**;
- a command that holds the GIL for more than 50 ms in one call is a defect.

### Performance budget

These are acceptance criteria, on a reference laptop with an integrated GPU:

| Interaction | Budget |
|---|---|
| orbit, pan, zoom on 10⁷-node assembly (surfaces + edges) | ≤ 16 ms/frame (60 fps) |
| pick → highlight + status bar (block, i, j, k, xyz) | ≤ 50 ms |
| drag an associated vertex (tier-1 projection, ADR 0002) | ≤ 16 ms per update |
| load project with 10⁷ nodes to first frame | ≤ 2 s |
| move a slice plane | ≤ 33 ms |
| UI input latency during a running march | ≤ 100 ms |

### Layout

The layout follows ICEM and Pointwise so that experienced users find their
way:
- **left:** a model tree (Geometry / Topology / Grids / Overset);
- **centre:** the viewport;
- **right:** a property panel for the selection;
- **bottom:** a tabbed console with the journal (ADR 0004), each GUI action
  shown as the Python line it executed, and the log.

The selection filter (geometry point/curve/surface, block
vertex/edge/face/block, grid node) is a toolbar toggle. Undo and redo are
bound to the command stack.

## Alternatives considered

| Option | Why rejected |
|---|---|
| PySide6 + pyvistaqt/VTK | ≈ 230 MB; OpenGL on macOS; reuses `vis.py` but fails the dependency and modernity criteria |
| imgui-bundle + pygfx | ≈ 20 MB and fastest, but forms, validation, tree views and file dialogs would be hand-built; SHORE's UI is parameter-heavy |
| trame (browser) | Remote use is not required in v1; carries VTK on the server; web-app feel |
| Dear PyGui | OpenGL; its 3D is not a scene graph; picking would be ours to write |
| Tauri/Electron + three.js | Second language and runtime; violates "Python only" |

## Consequences

- `vis.py` (PyVista) stays for scripted, offline and Jupyter visualisation
  and VTK export. The GUI does not depend on it.
- GUI logic is tested headless through the command layer (ADR 0004). GUI
  smoke tests run with `QT_QPA_PLATFORM=offscreen` and a wgpu software
  adapter (lavapipe) in Linux CI, plus a weekly manual check on macOS and
  Windows.
- **WSL2 is not a GUI development box** (spike, 2026-10-25). On Ubuntu 24.04
  WSLg, wgpu finds only llvmpipe (CPU). Mesa ships no `dzn` ICD, so the GPU
  is invisible, and orbit on 10⁷ nodes runs at about 160 ms/frame. GUI work
  uses native Windows Python on the same repo.
- The `Viewport` owns all camera state. pygfx's `OrbitController` is an input
  device only: its `target` setter calls `camera.look_at` in float32.

## Verification

A **spike** (≤ 1 week, throwaway code, results recorded in this ADR before it
moves to *Accepted*):

1. A PySide6 window with docked panels and a pygfx viewport via
   `rendercanvas` runs on Windows, macOS and Linux from a clean venv
   (`pip install PySide6-Essentials pygfx wgpu rendercanvas`).
2. A synthetic 10⁷-node, 20-block assembly is rendered as block surfaces and
   edges. Frame time is measured against the budget on all three OSes and on
   WSL2.
3. GPU pick → (block, i, j, k) is shown correct at a wall cell with 10⁻⁶
   spacing at 10³ offset (the float32 trap).
4. A 30 s march runs in a worker thread with the UI held at ≤ 100 ms input
   latency.

If (2) or (3) fails on pygfx, fall back to the imgui-bundle option before
considering VTK.

The spike code and per-machine results live in `spikes/gui/`; see its README.
Status on 2026-10-25:
- (3) and (4) pass;
- (2) is measured only on a CPU adapter (WSL2);
- (1) and (2) on native Windows, macOS and Linux GPUs are still open;
- the block-extent / `ds_wall` sweep for (3) is done; it selected rebase over
  per-tile origins;
- the rebase trigger (factor 100 view widths) is implemented and verified
  headless. The drift-error model `eps32 * D` matches the measured pick
  failures, and a user-style zoom stays exact with the trigger on and fails
  with it off. Checking it by eye on a native GPU is still open.
