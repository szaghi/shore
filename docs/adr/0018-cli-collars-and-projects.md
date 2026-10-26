# 0018 — CLI for collars and multi-component projects

- **Status:** Proposed
- **Date:** 2026-10-26

## Context

Collars (ADR 0010-0013), named geometries and components (ADR 0016) and
`write_xall` with options (ADR 0017) can be reached only from Python or a
journal. `shore mesh` and `shore primitive *` build one throwaway project
each; `shore run` replays journals. (`shore project` is the existing
surface-projection command, so the name is taken.)

A project directory *is* its journal: `Project.load` replays everything.

## Decision

User decisions of 2026-10-26:

- **Scope:** `shore collar`, a `--project DIR [--component C]` option on
  `shore mesh`, `shore collar` and `shore primitive *`, and
  `shore xall PROJECT_DIR --out RUN [--options xall.toml]`.
  - `shore collar BODY_A BODY_B ...` loads each file under a geometry name
    (default: its stem), selects the component (default `"collar"`) and
    runs `build_collar`. `--plane` builds a half-collar. `-o` writes the
    `.geo` files as `write_geo` does.
  - `shore mesh` loads its input under `--geometry-name` (default: the
    stem); its default component is `"main"`.
- **Editing a project: load, apply, save.** Each command loads the project
  (a full replay and bit-identity check), applies itself and saves. That is
  simple and always consistent. Each call re-meshes the existing components
  once: seconds for collars, longer for large meshes.
- **Geometries are reused.** A geometry name already loaded from a file with
  the same SHA-256 is not reloaded. Reloading would discard the components
  built on it (ADR 0016), so two collars on the same files could not share
  a project.
- **Options in a TOML file.** A `[run]` table holds the header options and
  `np` / `grd`. `[components.<name>]` tables hold `free_face_bc` and `body`,
  plus `meta.<label>` and `proc_meta.<label>` subtables. `boxes` is the
  path of a `shore boxes-from-stl` JSON file, relative to the options file.
- **A frozen spec.**

Found while validating the spec: **naming an open half-model STL failed.**
`Body` resolves its ray-casting anchor when it loads, and an open half
body has none, so `load_geometry(half_post, name=...)` raised. ADR 0013's
command never hit this, because it read body B without a `Body`. The
anchor is needed only for ray-cast topologies, and collars never use it.
So `Body` keeps the resolution error and raises it from `projector`, the
only place the anchor is used. A body that is not star-shaped then fails
where it is first ray-cast (`build_topology`), not where it loads. It must
not fall back silently to the centroid, which `Body.__post_init__` does
for a missing anchor.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_cli_projects.py` runs in the main CI job.

**Verdict.** It must report 13 passed, 0 skipped, 0 xfailed. At freeze it
reports 0 passed, 13 xfailed, every case failing on the missing command or
option (Typer's usage error, exit code 2). The error cases require exit
code 1, so a usage error does not satisfy them.

| Case | Required |
|---|---|
| L1 `shore collar --project` | equals the Python API's collar, bit for bit |
| L2 half collar `--plane` | equals the Python API's; SYMMETRY j-faces (needs the deferred anchor) |
| L3 `shore collar -o` | `.geo` files byte-identical to `write_geo` of the API's collar |
| L4 `shore mesh --project` | equals the replay of the same run's `--journal`; geometry named by the stem |
| L5 mesh, collar and primitive into one project | three components, the earlier ones unchanged |
| L6 two collars on the same files | both survive; geometries reused |
| L7 `shore xall --options` | `.grd`, `cc.par`, `proc.input` byte-identical to `write_xall` with the options the file describes |
| L8 errors (× 6) | unknown key or component in the options file, not a project, `--component` without `--project` (collar, mesh), collar with nothing to write: exit 1 |

**Harness validation (T1).** Every reference the spec compares with was
run through the Python API, with the anchor deferred: the full and the
half collar (11 × 16 × 9, 0 folds, faces exactly on the plane), and the
two-component `write_xall` with the options of L7. Default STL hygiene
leaves the open half meshes open.

**Out of scope:**
- a CLI for every `Project` method;
- progress output;
- CAD inputs (the same code path; checked by the CAD job);
- run time.

**For review (your call before implementation starts):**
1. With load, apply, save, adding a collar to a project with a large
   body-fitted mesh re-marches that mesh first. Acceptable for now, or
   should projects cache their grids and replay lazily (a format change)?
2. Geometry names default to file stems. Two different files with the
   same stem, for example `wing/body.stl` and `fuselage/body.stl`, would
   collide: the second load replaces the first and discards its
   components. Should a stem collision with different content be an error
   asking for `--name-a` / `--geometry-name`?
