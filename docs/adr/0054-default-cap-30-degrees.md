# 0054 — The cubed sphere's default cap: 30 degrees

- **Status:** Accepted
- **Date:** 2026-10-31

## Context

Since ADR 0052 a cubed-sphere march that folds a cell is refused. Its
default cap angle, `theta_cap = 0.05π` (9°), folds the cap corners in every
algebraic case of ADR 0052's sweep: sphere and capsule at three sizes,
4–104 folded cells. It also trips the cap/equator C1-pinning fallback
warning. So `shore mesh body.stl`, with the defaults, now fails on most
bodies.

At 30°, the same sweep has no fold and no warning: both bodies, the three
sizes, both schemes, minimum scaled Jacobian 0.57 or more. ADR 0052 left
the change to the owner, because the default backs ADR 0034's frozen R8
goldens. **The owner chose the per-topology default with a legacy guard
(option A, 2026-10-31).**

## Decision

1. **The cubed sphere's default cap is π/6 (30°).** It applies in
   `CubedSphereTopology.build`, `Mesh.from_stl` / `from_file`,
   `Project.build_topology` and `shore mesh`.
   - The **O-grid keeps 0.05π**. There `theta_cap` is the polar opening:
     the O-grid covers colatitudes [θ, π − θ], so 30° would leave 30° holes.
   - `topology.CUBED_SPHERE_THETA_CAP` and `topology.LEGACY_THETA_CAP` name
     the two values.
   - OH still pins 45°.
2. **Journals carry the value.** `Project.build_topology("cubed_sphere",
   …)` without `theta_cap` adds the default to its parameters, so the
   journal line says `theta_cap=0.5235987755982988`. No journal depends on
   a default again. `shore mesh` already wrote the value it used.
3. **Projects become format 4.**
   - **Legacy guard:** while a project saved before format 4 replays, a
     cubed sphere built without `theta_cap` gets the old 0.05π
     (`commands.legacy_theta_cap`, like ADR 0026's `legacy_collar_volume`).
     Its saved grids match bit for bit, and a re-save journals 0.05π
     explicitly.
   - **Loading without replay** (ADR 0033) accepts formats 3 and 4, both of
     which save face graphs.
4. **Re-freezes, approved by the owner 2026-10-31.**
   - ADR 0034 R8's cubed-sphere goldens are regenerated at the new default.
     The other topologies' entries are unchanged, checked byte for byte.
   - `test_cap_folds` F1, F2 and F4 give 9° explicitly.
   - `test_fast_load` L1 and `test_collar_default` D1 check `format = 4`.
   - The golden must change with the default. A separate re-freeze commit
     would fail its own gate, so the golden is in the implementation commit.

## Alternatives considered

- **Keep 9° and document it:** the default command fails on most bodies.
- **Change the default only in the CLI:** the CLI and Python defaults would
  disagree.
- **Treat a missing `theta_cap` in any replayed journal as 9°:** a
  hand-written script run through `shore run` would silently get 9°.

## Consequences

- `shore mesh body.stl` (cubed sphere) marches the sweep's cases cleanly,
  without the C1 warning.
- **Bare journal scripts** written before ADR 0054 that build a cubed
  sphere without `theta_cap` replay at 30°. They carry no format, so
  nothing can tell them apart. Add `theta_cap=0.15707963267948966` to
  replay the old grid. Saved projects are protected by the guard.
- The docs no longer need `--theta-cap-deg 35` to avoid folds; the
  tutorials keep their explicit value.

## Verification

`tests/spec/test_default_cap.py` (main job), frozen before the
implementation:
- C1–C3: the Python, CLI and `Mesh.from_stl` defaults give the same 30°
  grid, fold-free, without the warning, journaled explicitly.
- C4: the O-grid's default is unchanged (a guard).
- C5: a format-3 project made at the base commit
  (`data/cap9_project.zip`) loads bit for bit (a guard).
- C6: saved again, it is format 4 and journals 0.05π.
- C7: new projects save format 4.

At freeze: 2 passed, 5 xfailed.

A mutation check confirmed the guard is live: with a 30° legacy value, C5
fails, because the replay no longer matches the saved grids.
