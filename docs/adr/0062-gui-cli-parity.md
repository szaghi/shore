# 0062 — Every CLI capability in the GUI

- **Status:** Accepted
- **Date:** 2026-11-01

## Context

The GUI (ADR 0006, 0037–0045) had 20 hand-written forms over `Project`
methods. Much of what the CLI does was out of reach:

- **No usable primitive or topology form.** `build_primitive` and
  `build_topology` took their parameters as a raw Python dict literal, with
  no defaults: "I cannot create primitives in the GUI".
- **No generators:** `shore body` (sphere, cylinder, torus, box, capsule,
  extrude) and `shore profile naca`.
- **No Xall file tools:** `grd`, `cc-par`, `proc-input`, their merges,
  `balance`, `boxes-from-stl`, `boxes-export`.
- **No inspection tools:** `info`, `check`, `export`, `project`, `config`.
- **No forms for many Project commands:** collar, cap-tip, split,
  load-blocks, write-geo, the Xall writers, a surface grid on a CAD face.
- **No options on the writers:** `write_xall` had no options TOML, and
  CGNS / Plot3D had no component selection.

Every hand-written form also drifts from its CLI command: defaults, names
and new options. **The owner asked for every CLI feature in the GUI
(2026-11-01).**

## Decision

1. **The GUI's command forms are generated from the CLI.**
   - Every leaf command of the typer app becomes a form keyed
     `"cli:<command path>"`.
   - The fields are the command's options and arguments: named as the CLI's
     Python parameters, typed from their click types, with the CLI's
     defaults, choices and help.
   - The options the GUI owns are left out: `project`, `component`,
     `no_replay`, `journal`.
   - A new CLI option therefore appears in the GUI without GUI code.
2. **A CLI form runs the CLI in-process, on the session's project.**
   - `shore.cli.using_project(project)` makes `_open_project` and the
     commands' own `Project()` return the session's project, so a command
     applies its Project commands there: the same journal lines as on the
     command line, each one undoable.
   - Files are written as the CLI writes them, and the command's output
     goes to the log.
   - A failing command (exit code not 0) sets `last_error`, and its partial
     changes are undone, so the project is as it was.
3. **Mapped, not generated:**
   - `xall`, `cgns` and `plot3d` open a project directory, so they map to
     the session's `write_xall` (now with the options TOML),
     `write_cgns` and `write_plot3d` (with components);
   - `run` maps to `Session.run_journal` (the journal's project becomes the
     session's);
   - `view` maps to `Session.view_file`: a `.geo` loads as blocks (one
     journal line), and a `.npz` surface is displayed only;
   - `gui` is excluded.
4. **Every Project command has a form.** The hand-written forms the
   existing specs use stay; forms are added for the commands that had none.
   - `FORMS` keeps only the forms over Project methods: ADR 0037's G12
     requires every one to call `Project.<method>`.
   - `ALL_FORMS` adds the session's forms and the CLI's.
   - `form(key)` finds any of them.
5. **The forms are grouped.** A group combo in the Command dock and a
   Commands menu follow the CLI's structure: Project, mesh, primitive,
   body, profile, Xall files, inspect.

## Alternatives considered

- **Hand-written forms per command:** about 40 forms to keep in step with
  the CLI by hand, which is the drift this ADR removes.
- **Forms from the `Project` method signatures:** they miss the file tools
  and the CLI's defaults, and the CLI's mapping from options to commands
  (`shore mesh`'s 70 options) would be duplicated.
- **Running the CLI as a subprocess on a saved copy:** it loses undo and
  the live session, and costs a save and a load per command.

## Consequences

- Running `mesh` from the GUI writes its `.geo` files where the CLI would.
  The form's `output` field sets the location.
- The CLI's Rich output is captured, so commands run one at a time; the
  GUI already runs commands one at a time (ADR 0037).
- A CLI command that builds is two undo steps in the GUI (build, march),
  as its journal is two lines.

## Verification

`tests_gui/spec/test_gui_parity.py` (GUI job), frozen before the
implementation:
- P1: every CLI leaf command has a form, a mapping or an exclusion.
- P2: every Project command has a form.
- P3: the fields carry the CLI's defaults, without the project options.
- P4 (×5 primitives) and P5 (`mesh`): the CLI's journal lines and grids.
- P6 (generators) and P7 (`grd`, `proc-input`, `cc-par`): the CLI's bytes.
- P8: a failure sets `last_error` and changes nothing.
- P9: `write_xall` with options matches `shore xall --options`.
- P10: `run_journal` replays a journal.
- P11: `view_file` handles `.geo` and `.npz` files.
- P12: every form is in the Commands menu.

At freeze: 0 passed, 16 xfailed.
