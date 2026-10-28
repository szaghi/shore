# 0034 — Remove the marcher's normal blending (`blend_normals_k`)

- **Status:** Proposed
- **Date:** 2026-10-28
- **Supersedes:** the deprecation interface of [ADR 0019](./0019-deprecate-normal-blending.md)

## Context

ADR 0019 deprecated `blend_normals_k` in release 0.2.0. Blending body
normals into the first layers' marching directions showed no benefit over
the grid normals, and folded every C-grid it was tried on. Its removal was
left for a later release. That release is now 0.5.0, a breaking minor
release on the 0.x line.

The option reaches five entry points:
- `Mesh.march`;
- `Project.march` (the `March` command and its journal line);
- `shore mesh --blend-normals-k`;
- the `[march] blend_normals_k` config key;
- saved projects whose journal contains `p.march(blend_normals_k=K)`.

User decisions of 2026-10-28:
- **K > 0 is refused everywhere,** with an error that names the option and
  the removal. A project saved with blending still loads without replay
  (ADR 0033), because its grids are saved and the march is only recorded.
- **An explicit 0 is still accepted,** with a notice. Old configs and
  scripts that only spell out the default keep running.
- **ADR 0019's spec is replaced** by a removal spec. Deleting it is the
  re-freeze.

## Investigation (2026-10-28, prototype in the session scratchpad `br/`; T1)

**What the blend touched.**
- Each marcher (cubed sphere and OH wrap, C-grid and CH wrap, single-block
  C-grid, O-grid) computed `alpha * body normal + (1 - alpha) * grid
  normal` for k < K and passed it down as `normals_override`.
- Nothing else passes `normals_override`, so the hook of `_march_one_layer`,
  `_march_cap_layer` and `hyperbolic_layer` goes with it.
- `Body.normals` stays: ADR 0005's spec measures k-line angles against it.

**Unblended marches are unchanged.** Golden marches taken before the removal
are reproduced bit for bit on the prototype:
- O-grid and cubed sphere, algebraic and hyperbolic schemes;
- a 3-block C-grid.

**Old projects.** A journal line `p.march(blend_normals_k=3)` must still
parse and record, so that ADR 0033's load without replay keeps working. The
`March` command therefore keeps the field, and refuses it only when it
runs. The field is a tombstone, not a feature.

## Decision

- **`Mesh.march(..., blend_normals_k=None)`.**
  - `None` (the new default) does nothing.
  - 0 emits a `DeprecationWarning` saying it was removed.
  - Anything else raises `ParameterError` naming `blend_normals_k` and the
    removal, before any layer is marched.
  - The keyword keeps its position, so positional calls still bind
    `smooth_mode` and `march_scheme` correctly.
- **`Project.march(..., blend_normals_k=None)`.**
  - 0 warns and journals nothing extra.
  - K > 0 is recorded in the `March` command, which raises the same
    `ParameterError` when it runs, before the undo snapshot: the mesh and
    the journal stay unchanged.
  - A replayed load of a blended project therefore fails, and a load
    without replay succeeds.
- **CLI and config.**
  - `--blend-normals-k` is hidden from `--help`.
  - The flag or the `[march] blend_normals_k` key set to K > 0 exits 1 with
    an error naming `--blend-normals-k` and the removal.
  - Set to 0, it prints a notice and runs.
  - The config template no longer lists the key.
- **Code removed:** the blend in every marcher, and the `normals_override`
  hook of the layer helpers.
- **Documentation and examples** drop the option. Examples 01 and 05 stop
  setting it, and the fuselage config used K = 5.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Ignore K > 0 with a warning | A blended project's replay would fail the grid check anyway; it silently changes what scripts produce. |
| Remove the name completely (unknown key, `TypeError`) | Breaks configs and scripts that only spell out the default 0. |
| Refuse blended journals even without replay | Their saved grids are valid; ADR 0033 can still serve them. |
| Edit ADR 0019's spec in place | The removal is a new decision; a new spec keeps 0019's record intact. |

## Consequences

- **Breaking change (0.5.0).** Any `blend_normals_k` other than 0 fails.
  This includes old journals on a replayed load; re-mesh, or load without
  replay.
- **The marchers lose a code path** that every layer branched on; there is
  no change to unblended grids.
- **The tombstones** (the keyword, the hidden flag, the config key and the
  `March` field) can go in a later release, once old projects are re-saved.
  Deleting them would then break only explicit zeros.

## Verification

### Acceptance spec (frozen 2026-10-28, before implementation)

`tests/spec/test_blend_removal.py` runs in the main CI job, with its data
`data/blend_removal_default.npz` (the golden marches) and
`data/blended_project.zip` (an O-grid sphere marched with
`blend_normals_k=3` and saved, both made at the commit before the removal).
It replaces `tests/spec/test_blend_deprecation.py` (ADR 0019), which is
deleted together with its hash.

**Verdicts.**
- Required: 18 passed, 0 skipped, 0 xfailed.
- At freeze: 7 passed (R7, R8 × 5, R9a: regression guards), 11 xfailed.
- With the prototype behind the interface: 18 passed, and the five golden
  marches bit-identical.

| Case | Required |
|---|---|
| R1 (× 3) | `Mesh.march` refuses K > 0 (O-grid, cubed sphere, C-grid); mesh unchanged |
| R2 | `Project.march` and a journal replay refuse K > 0; mesh and journal unchanged |
| R3, R4 | `--blend-normals-k 3` and `[march] blend_normals_k = 3` exit 1, naming the removal |
| R5 (× 2) | an explicit 0 in `Mesh.march` / `Project.march` warns; same grids and journal |
| R6 (× 2) | an explicit 0 as CLI flag / config key prints a notice; same files |
| R7 | no warning without the option |
| R8 (× 5) | unblended marches equal the golden (atol 1e-10) |
| R9a, R9b | the blended project loads without replay; its replayed load is refused |
