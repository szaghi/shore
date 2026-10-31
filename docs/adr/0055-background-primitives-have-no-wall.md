# 0055 — Background primitives have no wall

- **Status:** Accepted
- **Date:** 2026-10-31

## Context

`HexBlock` gives every new block `k_lo = WALL`. That default is right for
marched topologies, whose body is at k = 0. The background primitives (box,
annular cylinder and flat-caps cylinder) set none of their face conditions,
so each kept a wall at `k_lo`, which is a far-field face for them:
- `cc.par` wrote BC 1 (`parete`, a viscous wall) on the end of the
  background cylinder;
- CGNS wrote `BCWall`;
- `free_face_bc`, which sets the FREE faces, could not change it.

The documentation had to work round it: the Xall-bundle tutorial patched the
sidecar by hand, and the multi-component recipe warned about it. Found while
writing the documentation (2026-10-30). **The owner asked for the fix
(2026-10-31).**

## Decision

1. **`box`, `annulus` and `flat_caps` make every WALL face FREE**
   (`primitive.base.free_boundaries`, in their constructors). Their faces are
   FREE, or SHARED where the flat-caps blocks meet. `free_face_bc` then sets
   all their boundary codes: 9 for a far field, 40 for an overset fringe.
2. **The duct primitives are unchanged.** `h_grid_rect` keeps its j and k
   walls, and `h_grid_tube` its outer `k_hi` wall: their walls are real.
   `h_grid_tube` marks its centreline singular axis (`k_lo`) WALL by design;
   that is not changed here.
3. **Re-freeze, approved with the owner's request.** The flat-caps entries of
   `tests/spec/data/xall_before.zip` (ADR 0036's "the writer before" archive)
   carried the wall.
   - They are rewritten as the old files with only this change: the five
     `k_lo` patches 1 → 40 in `flat_caps.par` and WALL → FREE in
     `flat_caps.json`.
   - ADR 0036's lo-face connection codes stay as they were, so the archive
     still means "before ADR 0036".
   - The archive changes with the code, so both are in one commit.

## Alternatives considered

- **Change `HexBlock`'s default to FREE:** every marched topology relies on
  it, so the risk is out of proportion.
- **A per-call option on the primitives:** a far field never wants a wall,
  so it would be an option that is always set.

## Consequences

- The Xall-bundle tutorial no longer edits the sidecar; the recipes no
  longer warn.
- A user who wants a wall on a primitive face sets it on the block
  (`_set_face_bc`), as the split unit test now does.
- A primitive's BCs in the face graph of a project saved before ADR 0055
  change on replay. A load without replay keeps the saved face graph.

## Verification

`tests/spec/test_background_bcs.py` (main job), frozen before the
implementation:
- B1 (×3): box, annulus and flat-caps have no WALL face.
- B2: flat-caps through `write_xall` with code 9: no BC 1, all 9.
- B3: a box in CGNS has no `BCWall`.
- B4: the duct walls are kept (a guard).
- B5: a marched cubed sphere keeps its body wall (a guard).
- B6: `shore primitive box --project` saves no WALL.

At freeze: 2 passed, 6 xfailed.
