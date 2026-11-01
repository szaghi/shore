# 0058 — Refuse folded wing-tip caps

- **Status:** Accepted
- **Date:** 2026-10-31

## Context

The C-grid's wing-tip caps are built, not marched:
- `cap_tips` on `CGridTopology.build`, which is `shore mesh --cap-tips`;
- `CGridTipCap`, `Project.build_tip_cap` and `shore cap-tip`.

So none of the marchers' guards (ADR 0052, 0056) sees them. ADR 0057
measured them on a NACA 0012 wing. Every tip has mirrored blocks (every
cell at about −1) and a degenerate or twisted LE block. Its TE blocks have
one node along i, so they have no cells. All of it is written without a
word. Since ADR 0056, "SHORE never writes a folded grid" holds for the
marchers but not for the caps. **The owner asked for the guard (2026-10-31),
before the redesign (ADR 0059).**

## Decision

1. **A cap build checks its blocks.** Once the cap blocks exist, every one
   must have at least one cell along each axis, and every cell a minimum
   corner scaled Jacobian above 0 (`shore.quality`, the frozen metric). The
   check runs in `CGridTopology.build`'s cap path and in
   `CGridTipCap.build`.
2. **Otherwise the build is refused** with a `ParameterError`.
   - The message names the block and either its worst cell or its missing
     cells.
   - The refused command is not recorded: the project keeps its mesh, and
     the CLI exits 1 without writing files.
3. **The caps of the base commit fail on every tested wing,** so until ADR
   0059 `--cap-tips` and `shore cap-tip` are refused. This is deliberate:
   they produce grids no solver should get. The documentation says so.
4. **Builds without caps are unchanged,** bit for bit.

## Alternatives considered

- **A warning:** the grids would still be written, and a warning in a long
  log is the silent failure ADR 0052 removed for the marchers.
- **Deprecating the feature:** the redesign follows the same night
  (ADR 0059), so a deprecation would be withdrawn at once.

## Consequences

- `shore mesh --cap-tips` and `shore cap-tip` fail until ADR 0059, with a
  message naming the block.
- The guard stays after the redesign: a cap that folds on unexpected input
  is refused, not written.
- The cap unit tests that built the old caps and checked their structure
  now check the refusal; the wing tutorial drops its cap sections until
  ADR 0059.

## Verification

`tests/spec/test_tip_cap_guard.py` (main job), frozen before the
implementation:
- P1–P4: every path either returns only valid cap blocks or refuses the
  build, changing nothing. This holds for any construction:
  - P1 (×2): the conformal caps on a sharp and a blunt wing;
  - P2: `Project.build_tip_cap`;
  - P3: `shore cap-tip`;
  - P4: `shore mesh --cap-tips both`.
- P5: a C-grid without caps equals its grid at the base commit
  (`data/tip_cap_guard_clean.npz`).

At freeze: 1 passed, 5 xfailed.
