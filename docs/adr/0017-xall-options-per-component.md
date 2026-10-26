# 0017 — Per-component Xall options in `write_xall`

- **Status:** Accepted
- **Date:** 2026-10-26

## Context

`write_xall` (ADR 0016) writes one Xall run from every component of a
project, but only with the writers' defaults. A background component needs
natural boundary conditions on its outer faces, components need block
metadata and marker boxes, and `proc.input` needs body numbers. Without
these options, the "one project, one Xall run" promise of ADR 0016 holds
only for default setups.

`merge_cc_par` takes the `cc.par` header from the first component and
requires the others to agree (`_check_header_consistency`). The header
options are therefore one per run; only the block-level options can differ
per component.

## Decision

User decisions of 2026-10-26:

- **Keyword arguments on `write_xall`.** The command stays output-only,
  with no new project state; the journal line carries the options.
  ```
  write_xall(directory, *, grd="mesh.grd", np=1, components=None,
             output_basename="cc.", multigrid_levels=(4, 1),
             boundary_layer_thickness=-0.1, numerical_beach=-1.0)
  ```
  - The header options apply to every component's `write_cc_par`.
  - `components = {name: {...}}` gives each component its `cc.par`
    options, with `write_cc_par`'s names: `free_face_bc`, `meta` (block
    label → level / group / priority / comment / free_face_bc) and
    `boxes`.
- **`proc.input` overrides per component.** `body` sets the Xall body
  number of every block of the component. `proc_meta` (block label → proc
  / group / body) overrides single blocks, and its `body` wins over the
  component's. The overrides are passed to `write_proc_input` at the merged
  1-based block indices; the run is still balanced once (ADR 0016, review
  decision 1).
- **Checked before writing.** An unknown or empty component, an unknown
  key, or a block label the component does not have raises
  `ParameterError`, and nothing is written.
- **A new frozen spec module**; the accepted ADR 0016 spec is untouched.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_xall_options.py` runs in the main CI job.

**Verdict.** It must report 9 passed, 0 skipped, 0 xfailed. At freeze it
reports 0 passed, 9 xfailed, every case failing only with the `TypeError`
of the missing keyword arguments.

| Case | Required |
|---|---|
| X1 per-component cc.par options | merged and per-component `cc.par` byte-identical to the accepted writers with the same options, merged |
| X2 run-wide header options | the same, with all four header options changed |
| X3 proc.input | component bodies and per-block overrides, byte-identical to `write_proc_input` at the merged indices |
| X4 journal | the `write_xall` line with options replays to the same files |
| X5 errors (× 5) | unknown component, unknown key, unknown label in `meta` / `proc_meta`, a component without a mesh; nothing written |

**Harness validation (T1).** A stand-in `write_xall` with the options,
calling the existing writers, was run behind the frozen interface.
- It passes 8 of 9 cases. It fails exactly X4: as a stand-in it is not a
  command, so it leaves no journal line.
- With the merged block indices off by one, it also fails X3.

**Out of scope:**
- options for `.geo` / `.grd` writing;
- per-component rank sets;
- CLI flags.

## Acceptance (2026-10-26): Xall options spec met

`WriteXall` gains `components`, `output_basename`, `multigrid_levels`,
`boundary_layer_thickness` and `numerical_beach`; `Project.write_xall`
passes them on. Per component it builds `WriteCcPar` with the run-wide
header and the component's `free_face_bc` (default `WriteCcPar`'s own),
`meta` and `boxes`. It collects `body` and the `proc_meta` overrides at
the merged 1-based block indices for one `write_proc_input`. `_check`
validates every option before anything is written; besides the frozen
cases it also rejects unknown keys inside `proc_meta`.

**Verdict for `tests/spec/test_xall_options.py`: 9 passed, 0 skipped,
0 xfailed** (T2; `scripts/check_spec.py` reports all files ok). All 5
`@PENDING` lines are removed. `test_project_components.py` still reports
14 passed. Full suite: 1441 passed, 0 xfailed.

**Unit test** (`tests/test_project.py`): an unknown `proc_meta` key raises
before anything is written.
