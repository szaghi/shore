# 0037 — GUI slice 1: viewer and command shell

- **Status:** Proposed
- **Date:** 2026-10-29

## Context

Roadmap phase 6 is GUI v1: model tree, viewport, property panel and journal
console; geometry and topology display; association editing; spacing;
marching with progress; a quality view. Its acceptance is that every GUI
action appears as a journal line and replays headless to identical output.
ADR 0006 fixes the stack (PySide6 + pygfx/wgpu). It was accepted on
2026-10-29 on its native Windows spike; on that date the window was also
checked by eye.

User decisions of 2026-10-29:
- **The slice:** the first slice is a viewer and command shell:
  - open and save projects;
  - the model tree;
  - the viewport, with picks returning (block, face, i, j, k, xyz);
  - the journal console;
  - the existing commands run from forms, the march with progress and
    cancel;
  - quality colouring.

  Association editing, the spacing editor and slices come later.
- **Acceptance:** a headless spec (`QT_QPA_PLATFORM=offscreen`, software
  Vulkan). The performance budgets are measured on native Windows (T1), not
  in CI.
- **Dependencies:** a `[gui]` extra; the GUI tests run in their own CI job.

## Investigation (2026-10-29, prototype in the session scratchpad `gui37/`; T1)

**The session layer.** A Qt-free `Session` owns the `Project`. Every widget
action calls it, and each state change is a `Project` command, so the
journal is the record of the GUI session by construction. Forms are
specifications (`shore.gui.forms`): a form's `arguments(values)` builds the
method's arguments, so a form adds no logic between the user and the
command.

**The viewport.** The spike's renderer, with only the rebase mode kept
(ADR 0006).
- Each block is one pygfx mesh of its six boundary faces as quads (pygfx
  0.17 renders quads natively), plus its twelve edges.
- A pick returns the quad index, and the structured arithmetic
  (`quad_cell`) gives the face and the (i, j, k) cell. The reported
  centre comes from the float64 nodes.
- Quality colours are per quad (`color_mode="face"`): the scaled Jacobian of
  the cell behind each boundary face.

**Two defects found while prototyping.**
- Setting pygfx's orbit-controller target calls `camera.look_at`. When the
  camera sat on the target (a fresh window), that produced a NaN pose.
  Rebasing in `show_all` / `look_at_point` now leaves the camera alone and
  frames afterwards.
- Worker signals connected to plain closures run in the worker thread (Qt
  has no receiver to queue them to). Bound methods of the window are
  queued to the UI thread.

**Headless run of the prototype.** The cubed sphere, a flat-caps second
component, undo/redo, a bad form value and saving all work, and the saved
project's replay equals the session's grids bit for bit. The spec's G6 fails
when the rebase is disabled, so it does catch the float32 trap.

## Decision

- **Package `shore.gui`:**
  - `session` (Qt-free `Session`, `BlockKey`);
  - `forms` (`Form`, `Field`, `FORMS`: load geometry, use component, build
    topology, build primitive, march, snap, write Xall / CGNS / Plot3D);
  - `viewport` (interface and structured arithmetic);
  - `viewport.pygfx_backend` (the only import of pygfx, wgpu and
    rendercanvas);
  - `window` (`MainWindow`).
- **Window** (ADR 0006 layout): the model tree on the left (components →
  blocks → faces with BCs), the viewport in the centre, command forms on
  the right, the journal and log at the bottom.
  - Toolbar: New, Open, Save, Save as, Undo, Redo, Fit, Quality.
  - A double click picks; the status bar shows component/block, face, cell
    and xyz.
- **Threading:** long forms (builds, march, Xall) run on a one-thread
  `QThreadPool`. Progress and cancel use the ADR 0004 hooks (cancel applies
  to the march). Errors go to the log and the status bar, not to modal
  dialogs.
- **Entry points:** `shore gui [DIR] [--no-replay] [--smoke]` and the
  `gui-scripts` entry point `shore-gui`.
- **Packaging:** `shore[gui]` = PySide6-Essentials, pygfx 0.17, wgpu 0.32,
  rendercanvas 2.7.
- **Tests:** in `tests_gui/`, outside `testpaths`, run by the `test-gui` CI
  job (lavapipe, offscreen Qt).

## Alternatives considered

| Alternative | Why not |
|---|---|
| Widgets calling `Project` directly | The GUI logic would be testable only through Qt. |
| Forms generated from method signatures | `build_topology` / `build_primitive` take `**params`; signatures cannot describe them. |
| GUI tests under `tests/` with an import skip | The main job would report skips; the frozen `tests/conftest.py` cannot hide them. |
| Modal error dialogs | They block headless tests; the log and status bar suffice. |

## Consequences

- **The journal is the GUI's record.** A GUI session saved as a project
  replays headless.
- **Not in this slice:** association editing, the spacing editor, slices,
  geometry display, vertex drag, mouse-driven tests, macOS and native
  Linux.

## Verification

### Acceptance spec (frozen 2026-10-29, before implementation)

`tests_gui/spec/test_gui_shell.py` runs in the `test-gui` CI job, with
`tests_gui/conftest.py`.

**Verdicts.**
- Required: 16 passed, 0 skipped, 0 xfailed.
- At freeze: 1 passed (G11b), 15 xfailed.
- With the prototype behind it: 16 passed.

| Case | Required |
|---|---|
| G1 | a session driven through the forms has one journal line per command, the console shows the journal, and the project saved from the GUI replays (load and `run_journal`) to identical grids |
| G2 | undo and redo through the window |
| G3 | the model tree is the project's |
| G4 | the viewport holds every block of every component, its boundary nodes uploaded |
| G5 | picks on camera-facing cells: block, face, cell, float64 centre |
| G6 | the float32 trap: 1e-6 wall cells at 1e3, exact picks after `look_at_point` |
| G7 | quality colours = colour map of the scaled Jacobian per boundary cell |
| G8 (× 2) | a march on the worker keeps UI lateness ≤ 100 ms; cancel leaves mesh and journal unchanged |
| G9 (× 2) | a bad value and a failing command are logged; nothing changes; the window still works |
| G10 (× 2) | a saved two-component project opens, with and without replay; `shore gui DIR --smoke` |
| G11 (× 2) | import isolation; the core imports without the GUI stack |
| G12 | every form calls an existing Project method with arguments it accepts |

### Performance (to be measured on native Windows at acceptance; T1)

These are ADR 0006's budgets, measured on the real GUI:
- orbit at 10⁷ nodes;
- pick;
- load to first frame;
- UI latency during a march.
