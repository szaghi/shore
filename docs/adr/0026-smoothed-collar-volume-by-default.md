# 0026 — The smoothed collar volume by default

- **Status:** Accepted
- **Date:** 2026-10-28

## Context

ADR 0025 added `CollarParams.volume = "smoothed"` as an opt-in. On 49 test
junctions it cut the folded cells outside a folded wall from 2367 to 72.
The user decided (2026-10-28) to make it the default. That changes every
collar, so it re-freezes the specs that pin the old default bit for bit.

## Investigation (2026-10-28; T1, then T2 on the frozen specs)

**With `"smoothed"` as the default, twelve frozen cases fail:**
- 9 pin the old default bit for bit:
  - ADR 0024's S6: 5 cases;
  - ADR 0025's V3: 4 cases.
- **ADR 0020's S5:** `Project.build_collar` equals the kernel. The
  project's default was still `"stations"`.
- **ADR 0010's K4 and ADR 0013's H3** (the outer boundary within ±1 % of
  the thickness from the bodies) on cylinder × cylinder: 0.986. Averaging
  offsets whose directions rotate along a curved loop shortens them.

**Two remedies tried:**

| Remedy | Result |
|---|---|
| Restore each offset's length after smoothing | K4 and H3 still fail; ADR 0025's V1/V2 fail (3 cases) |
| Smooth only near corners, detected on the near-wall layer at 30° per node | K4 and H3 pass; the blunt wing at 11.8 × h with curvature folds 11 cells |

The blunt wing still folds because of the leading edge:
- with curvature its nodes turn less than 30° each, and on the near-wall
  layer the corners are already rounded;
- yet its radius (0.016) is far below the collar thickness (0.059), so
  smoothing is needed there.

The trigger that works is a corner at the collar's scale: the loop turns
by 60° or more within ± one thickness of the station.
- **Cylinder × cylinder:** radius 0.5, thickness 0.118, turns 27-32°, so it
  is not smoothed.
- **Leading edge:** about 180°, so it is smoothed.
- **Corners:** 90°, so they are smoothed.

The weight is 1 within one thickness of such stations and falls to 0 at
two thicknesses.

**Results with that trigger:**
- every frozen spec passes except the 9 bit-identity cases;
- the 49-case bench is unchanged from ADR 0025 (1843 folded cells; 72
  outside a folded wall).

**Saved projects.** `Project.load` replays the journal and requires every
grid to match bit for bit. A collar journalled without a volume would
replay with the new default, and older projects with collars would stop
loading.

## Decision

User decision of 2026-10-28: `"smoothed"` is the default.

- `CollarParams.volume` defaults to `"smoothed"`. `_smoothed_outer`
  smooths the offsets only near corners at the collar's scale
  (`_near_corners`: 60° within ± the thickness, measured on the wall's
  junction row). The local repair where the boundary steps back and the
  second trim stay as in ADR 0025.
- `BuildCollar.volume`, `Project.build_collar(volume=)` and
  `shore collar --collar-volume` default to `None`: the current default.
- **Projects are format 3.** While a format 1 or 2 project replays
  (`legacy_collar_volume`), a `BuildCollar` without a volume uses
  `"stations"`, so the saved grids still match.
- **Re-freezes (this ADR's spec commit):**
  - ADR 0024's S6 asks for `"stations"` by name. The golden collars are
    that construction, unchanged.
  - ADR 0025's V3 states the new default.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Keep `"stations"` the default | The user decided otherwise; it folds thick collars at corners. |
| Restore offset lengths after smoothing | Breaks ADR 0025's validity cases. |
| Smooth everywhere | Shrinks the envelope by 1.4 % on smooth curved loops (ADR 0010's K4). |
| Corners by node turn ≥ 30° on the near-wall layer | Misses bends tighter than the collar (the blunt leading edge). |
| Project default `"smoothed"` without a format bump | Older projects with collars fail their bit-identity check on load. |
| Loosen K4/H3 to ±2 % | Weakens accepted criteria to fit the change. |

## Consequences

- **New collars change.** At corners and tight bends the new default is
  valid where the old one folded (ADR 0025's numbers). Elsewhere it equals
  the 3-D offset-line construction, which is the same as `"stations"` on
  planar stations.
- **Old projects keep their collars.** Re-saving such a project writes
  format 3, with its collars still journalled as before.
- **To reproduce an old grid outside a project,** use
  `CollarParams(volume="stations")` or `--collar-volume stations`.

## Verification

### Acceptance spec (frozen 2026-10-28, before implementation)

**Spec.** `tests/spec/test_collar_default.py` runs in the main CI job. The
re-frozen `test_split_collar.py` (S6) and `test_collar_volume.py` (V3) go
with it.

**Verdict.**
- `test_collar_default.py`: 3 passed, 0 skipped, 0 xfailed. At freeze: 1
  passed (old projects load) and 2 xfailed.
- `test_collar_volume.py`: 25 passed. At the re-freeze: 21 passed, 4
  xfailed.
- `test_split_collar.py`: 33 passed throughout.

Validated with the implementation behind the interface: 61 passed.

| Case | Required |
|---|---|
| D1 | saved projects are format 3 |
| D2 | a format-2 project with a collar journalled without a volume loads; its collar is the `"stations"` block |
| D3 | a new project's collar without a volume is the `"smoothed"` block and round-trips |

### Review (2026-10-28, user decisions)

1. The trigger stays at 60° within ± one thickness, as an internal
   constant.
2. Project format 3 with legacy replay of format 1/2 collars is accepted.
3. Next step: release 0.4.0.

## Acceptance (2026-10-28): smoothed-default spec met

**Kernel.** As decided. The changes:
- `CollarParams.volume = "smoothed"` is the default;
- `_near_corners` weights the smoothing near loop corners at the collar's
  scale;
- `BuildCollar.volume` is `None`, resolved via `legacy_collar_volume`;
- `PROJECT_FORMAT = 3`;
- `--collar-volume` defaults to `None`.

**Verdicts** (T2; `scripts/check_spec.py` reports all files ok):
- `tests/spec/test_collar_default.py`: **3 passed, 0 skipped, 0 xfailed**;
- `tests/spec/test_collar_volume.py`: **25 passed, 0 skipped, 0 xfailed**;
- `tests/spec/test_split_collar.py`: **33 passed, 0 skipped, 0 xfailed**.

All 6 `@PENDING` lines are removed.

**Unit test** (`tests/test_collar.py`): the trigger.
- A loop of radius 0.5 is not smoothed at thickness 0.118.
- A loop of radius 0.05 is smoothed everywhere.
- A square's corners are smoothed and its mid-sides are not.

**49-case bench:** 1841 folded cells, against 1843 for ADR 0025's global
smoothing.
