# 0004 — Command layer and replayable journal

- **Status:** Proposed
- **Date:** 2026-10-25

## Context

Pointwise's Glyph (since Gridgen v14, 2002), ICEM's Tcl replay files and
Cubit's journals share one property. **Every GUI action is a script command,
recorded as it happens.** That is what makes those tools reproducible,
automatable and testable. A GUI that mutates model state directly yields a
tool whose results cannot be regenerated after a CAD change, and whose
behaviour can only be tested by clicking.

SHORE today has two drivers, the Typer CLI and the Python API, both stateless
per invocation. A GUI needs a third driver that is stateful, can undo, and
runs long operations without blocking the UI.

## Decision

1. **A `Project` document** (`shore.project`) owns the whole state:
   - geometries;
   - block topology with associations;
   - spacing and marching parameters;
   - generated grids;
   - overset metadata (levels, groups, boxes).

   On disk it is a directory:
   - `project.toml` holds parameters and references;
   - `adjacency.json` holds the topology and associations;
   - `grids/*.npz` holds the node arrays;
   - `journal.py` holds the command history.

   No new dependency: TOML, JSON and npz are already in use.

2. **Commands** are frozen dataclasses with:
   - `validate(project)`;
   - `apply(project) -> Undo`;
   - `to_source() -> str`, one line of public Python API.

   They are the *only* way to mutate a `Project`. The CLI, the Python API
   and the GUI all construct and execute commands.

3. **The journal is plain Python.** It calls the public `shore` API, in the
   way a Glyph script calls Glyph. `shore run journal.py` replays it
   headless. A journal is the reproducible recipe of a grid. Re-running it
   after a CAD edit regenerates the grid, with association re-matching from
   ADR 0003.

4. **Undo.**
   - Cheap commands (move vertex, set spacing, re-associate) store their
     inverse.
   - Commands that rewrite node arrays (march, smooth) snapshot the affected
     blocks' arrays to an on-disk undo cache under the project directory.
     They never hold them in RAM, because a 10⁷-node grid is 240 MB in
     float64.

5. **Long-running commands** take a `progress(fraction, message)` callback
   and a cancellation token, checked once per marching layer or smoothing
   sweep. The GUI runs them in a worker thread (ADR 0006). The CLI renders
   progress with `rich`.

6. **Determinism.** Replaying a journal on the same inputs gives
   bit-identical `.grd`. This is a tested invariant, since numpy
   reproducibility holds for a fixed operation order on one machine.

## Alternatives considered

- **GUI mutates the model directly and saves a binary state file.** Rejected:
  the result is not reproducible, cannot be tested headless, and cannot be
  regenerated after a CAD change.
- **A custom DSL or TOML command list instead of Python.** Rejected: it would
  need a second interpreter for loops, parameters and conditionals that
  Python already provides. TOML stays the format for *parameters*
  (`MeshConfig`), not for procedures.
- **Event sourcing with a persistent log database.** Rejected as
  over-engineering for a single-user desktop tool.

## Consequences

- The existing CLI commands are re-expressed as thin wrappers over commands.
  Their flags, config precedence and output stay unchanged. This is a
  refactor with the existing 999-test suite as the safety net.
- Every GUI feature is testable headless, by executing its commands and
  asserting on the `Project`.
- The journal is a user-facing artefact. Its API calls become a public
  compatibility surface, and breaking changes need deprecation.

## Verification

- `shore mesh` via the command layer produces byte-identical `.geo`, `.grd`
  and `cc.par` to the current implementation on every worked example.
- A journal recorded from a scripted session replays to bit-identical output.
- Undo followed by redo of a march restores the node arrays bit-identically.
