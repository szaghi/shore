# 0016 — Multi-component, portable projects

- **Status:** Proposed
- **Date:** 2026-10-26

## Context

A `Project` (ADR 0004) holds one body and one mesh, and its saved form is
not portable:
- `load_geometry` journals the geometry path as given, and `save` does not
  copy the file. **No saved project survives moving it or its inputs**,
  not even with a single body.
- A collar (ADR 0010, 0013) reads its second body from a file path inside
  the command, and it **replaces** the mesh. A body-fitted mesh and its
  collar cannot live in one project; they are separate Chimera
  components, combined afterwards with `grd-merge`, `cc-par-merge` and
  `proc-input-merge`.

## Decision

User decisions of 2026-10-26:

- **Named geometries and named components.**
  - `load_geometry(path, ..., name="body")` registers a named geometry.
    Reloading a name replaces it and discards the components built on it.
  - `use_component(name)` makes a component current, creating it empty if
    new; it is journaled and undoable. Every mesh command acts on the
    current component, and `mesh` / `build` are the current component's.
    The default component is `"main"`.
  - The commands that read a geometry take `geometry="body"`. For
    `build_collar`, *other* is a geometry name when one is loaded, else a
    file path.
  - One project is one Xall run: `write_xall(directory, *, grd, np)`
    writes each component's files and then the merged `.grd`, `cc.par`
    and `proc.input`. These are byte-identical to the accepted per-component
    writers followed by `merge_grd` and `merge_cc_par`. `proc.input` is
    balanced once over the merged `.grd`, not merged per component.
- **Inputs copied into the project.** `save` copies every file the
  journal reads into `inputs/`, writes the journal against those copies,
  and records their SHA-256 in `project.toml` (format 2). `load` resolves
  them against the project directory wherever it now is, and fails with
  `MeshInputError` on a changed copy.
- **Old journals and projects keep working.** No `name` means `"body"`; a
  file path still works for `build_collar`; there is one implicit
  component `"main"`; format-1 projects load. Saving writes format 2.
- **A frozen acceptance spec**, although this is architecture, not
  numerics: what users rely on is behaviour (a moved project loads, the
  merged Xall files are right), and that can be frozen.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_project_components.py` runs in the main CI job.

**Verdict.** It must report 14 passed, 0 skipped, 0 xfailed. At freeze it
reports 1 passed and 13 xfailed:
- P3 passes: format-1 projects load today, and it guards that they still
  do;
- 12 cases fail on the missing interface (`name=`, `use_component`,
  `geometries`);
- P1 fails because today's save is not portable: the moved project looks
  for the deleted STL.

| Case | Required |
|---|---|
| P1 | a saved project, moved, with its original STL deleted, loads bit-identically |
| P2 | a changed input copy raises `MeshInputError` on load |
| P3 | a format-1 project (absolute STL path) still loads (regression guard) |
| P4 | `build_collar` by geometry name equals `build_collar` by file path, bit for bit |
| P5 | an O-grid and a collar coexist as two components; the first is unchanged |
| P6 | `write_xall`: merged `.grd`, `cc.par`, `proc.input` and the per-component `cc.par` byte-identical to the accepted tools |
| P7 | undo / redo across components |
| P8 | save / load of two geometries and two components, inputs deleted |
| P9 | reloading a geometry discards only the components built on it |
| P10 (× 4) | unknown geometry, empty names, `write_xall` without a mesh; reporting |

**Harness validation (T1).** There is no prototype: for an architecture
change the implementation is the first attempt. Every *reference* the spec
compares with was run with today's API and works: the O-grid on the
icosphere STL, the collar by file path (11 × 16 × 9, 0 folds, 0.2 s), and
the per-component writers followed by `merge_grd`, `merge_cc_par` and
`write_proc_input`. At freeze every pending case fails for the reasons
listed above, and only those.

**Out of scope:**
- per-component `cc.par` options in `write_xall` (block metadata, marker
  boxes, natural BCs);
- the *output* paths of journaled write commands, which replay writes
  where the journal says;
- CLI flags;
- the GUI.

**For review (your call before implementation starts):**
1. `proc.input` is balanced once over the whole run, not concatenated
   from per-component balances (as `proc-input-merge` does). Right for
   Xall?
2. `write_xall` takes no per-component `cc.par` options yet, so a
   background component with natural BCs still needs the file tools. Add
   them now (a re-freeze), or later?
3. Replaying a saved journal re-runs its output commands, which write to
   the paths recorded in the journal, relative to the current directory.
   Should `save` drop output commands from the saved journal (they are
   never undone and do not change state)?
