# 0033 — Loading a saved project without replay

- **Status:** Proposed
- **Date:** 2026-10-28

## Context

`Project.load` replays the whole journal and checks every grid bit for bit
(ADR 0004, 0016), which makes a load also a determinism check. It also
costs as much as rebuilding the project, which can be minutes since ADR
0027's sub-stepping at clustered corners. `shore xall` reads a saved
project only to write its Xall run, yet replays it all.

User decisions of 2026-10-28:
- an opt-in fast load;
- the default load unchanged;
- a `--no-replay` flag on the commands that open a project.

## Investigation (2026-10-28, prototype in the session scratchpad `fl/`; T1)

**What a saved project held:**
- the journal;
- every component's grids;
- the current component's face graph only;
- the input hashes.

So a load without replay could not rewire the other components' blocks.
Saving each component's face graph next to its grids
(`grids/<component>/adjacency.json`) is enough.

**Rebuilding the state without the kernels.** The journal is executed in a
recording mode:
- `LoadGeometry` and `UseComponent` are applied: they are cheap, and later
  commands and saves need them;
- every other command is only recorded.

Each component's blocks are then read from its grids and rewired from its
face graph (`wire_blocks_from_adjacency`, as `load_blocks` does).

**Result** on a primitive plus a collar on STL bodies:
- the same components, current component, geometries, blocks (nodes bit
  for bit, face BCs, partners, axis maps) and journal as a replayed load;
- saving again, or writing the Xall run, gives identical files;
- the load takes 0.01 s, against 0.35 s replayed. The gap grows with the
  meshes.

**Format.** The face graphs are additive, so the project format stays 3:
ADR 0026's spec pins it. A project saved before this ADR lacks them; it
loads normally and must be saved once before it can load without replay.

**An older behaviour noticed:** each re-save prefixes the input copies
again (`00_plate.stl` → `00_00_plate.stl`). It is harmless and out of
scope here.

## Decision

- **Saving.** `Project.save` writes `grids/<component>/adjacency.json` for
  every non-empty component.
- **Loading without replay.** `Project.load(directory, replay=False)`:
  1. checks the input hashes;
  2. records the journal (the context variable `recording`;
     `Project.execute` applies only `LoadGeometry` and `UseComponent`);
  3. reads each component's blocks from its grids and rewires them from
     its face graph.

  A loaded component's `BuildState` is `"blocks"`, not marchable, built on
  the geometries its last recorded command names.
- **Undo** across the recorded history raises `ParameterError`.
- **CLI.** `--no-replay` on `shore xall` and on the mesh, collar and
  primitive commands with `--project`.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Memoised replay (heavy commands restore cached results) | Same state, far more machinery: topology objects, per-command caches. |
| Bump the project format to 4 | ADR 0026's spec pins format 3; the new files are additive. |
| Rebuild cheap commands (primitives), record heavy ones | Mixed states (a marchable topology over loaded nodes) invite wrong marches. |

## Consequences

- **Fast:** loading a saved project for its outputs (`shore xall`) or to
  add a component takes a fraction of a second.
- **No determinism check:** a load without replay does not re-check the
  grids. The default load still does.
- **Limits:** loaded meshes cannot be marched again, and the loaded history
  cannot be undone.

## Verification

### Acceptance spec (frozen 2026-10-28, before implementation)

**Spec.** `tests/spec/test_fast_load.py` runs in the main CI job.

**Verdict.** It must report 10 passed, 0 skipped, 0 xfailed. At freeze it
reports 1 passed and 9 xfailed; the 1 is a regression guard (the replayed
load still rejects a changed grid). Validated with the prototype behind
the interface: 10 passed.

| Case | Required |
|---|---|
| L1 | face graphs saved per component; format still 3 |
| L2 | same components, current, geometries, blocks and face graphs as a replayed load |
| L3 | no kernel runs (collar and primitive builders patched to fail) |
| L4 | same journal; re-saved directories identical |
| L5 | identical Xall runs |
| L6 (× 3) | undo and march refused; a project without face graphs refused, naming the replay |
| L7 | `shore xall --no-replay` writes the same files |
| L8 | the replayed load still rejects a changed grid (guard) |
