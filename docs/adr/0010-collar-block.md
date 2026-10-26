# 0010 — The collar block of two intersecting bodies

- **Status:** Accepted
- **Date:** 2026-10-26

## Context

A Chimera collar is a volume block wrapped around the junction of two
components, such as a wing root on a fuselage or a post on a plate. ADRs
0008 and 0009 provide its wall:
- `intersection_loop` gives the junction curve;
- `march_surface` grows a surface grid on each body away from the other.

This ADR joins the two halves and extrudes them into a volume block.

The hard part is the **concave junction**. The wall is creased along the
junction row, at an angle of 63–117 deg in the spec's cases. Extruding it
with the accepted hyperbolic volume kernel (ADR 0005) folds every case.

## Diagnosis (2026-10-26, prototypes in the session scratchpad)

**Plain hyperbolic extrusion folds.** On all four analytic junctions it
inverts 384–1,408 cells:
- only in the three rows at the junction, at every loop station;
- starting exactly when the layer thickness reaches the junction's surface
  spacing, independently of the first step (1e-4 to 5e-3);
- with the cells recovering in later layers.

**The mechanism is geometric.** Orthogonal offsetting at a concave crease
collapses the cells beside it. On a 90 deg junction, the cell between the
junction node (offset to (d, d)) and its neighbour on the floor (offset to
(h, d)) has width h − d, which vanishes at d = h. Any scheme that keeps
k-lines orthogonal at the crease folds there. Grids survive concave corners
by *splaying*, not by following offsets.

H3 (ADR 0005) never showed this: its grid crosses the concave edges
obliquely. A collar's junction is a grid line.

**Hypotheses tested and falsified** (all reverted, H spec 13/13
throughout):

| Hypothesis | Result |
|---|---|
| A first step too large for the crease cells | Wall-resolved steps (1e-4) fold the same way |
| An isotropic cell-size measure in the kernel's ρ | Measuring by the smallest spacing changes nothing |
| The junction node lagging behind the exact offset (√2) | Correcting its local solution: the node does move √2 × step, but the solve still explodes once the junction cells close, at k = 17 |

**Treatments prototyped**, on a collar 12 times the junction spacing thick:

| Treatment | Junctions fold-free |
|---|---|
| Arc-length splay of every layer | 1 of 4 |
| Hyperbolic, then elliptic, with an offset envelope | 0 of 4; folds only in the outermost layer |
| … with the envelope smoothed | 2 of 4 at 20 iterations; 100 pulls it onto the walls |
| … with an exact distance-level-set envelope | 0 of 4; nodes out of order at the envelope's own crease |
| … with an offset whose swallowtail is trimmed | 3 of 4; the tilted junction folds 19 cells at its sharpest (63-74 deg) stations |
| **… trimmed, with the crease mapped to the junction row** | **4 of 4**, min scaled Jacobian 0.08-0.62 |

## Decision

`shore.collar.build_collar(body_a, body_b, seed, *, n_loop, steps_a,
steps_b, steps_k, label="collar") -> HexBlock`:

1. **Wall.** The composition of the accepted kernels, frozen exactly in the
   spec. The loop runs along j (periodic, the O-grid convention). i runs
   from B's far edge, through the junction row, to A's far edge. j is
   reversed so that k leaves both walls.
2. **Near-wall layers.** The hyperbolic kernel of ADR 0005, while the
   layers are thinner than half the junction spacing, where it cannot fold.
3. **Outer boundary.** Per loop station, the near-wall layer's i-line is
   offset by the remaining thickness. The swallowtail that the concave
   junction makes is trimmed at its self-intersection, and the nodes are
   placed by arc length with the trimmed crease mapped to the junction
   row.
4. **Interior.** Per loop station, TFI between the near-wall line and the
   outer boundary, using the k schedule, then Winslow smoothing with
   Thomas–Middlecoff control functions. This is the method of ADR 0007,
   which keeps the wall clustering. Elliptic grids cannot fold in the
   continuum.

`Project.build_collar(other, *, seed, n_loop, ds_surface, growth_surface,
m_a, m_b, ds_wall, growth_wall, nk, label="collar")`:
- body A is the project's geometry; body B is loaded from the file *other*,
  whose path the journal records;
- the collar replaces the mesh, like a primitive, and becomes a separate
  Chimera component for the existing `*-merge` tools.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_collar.py` runs in the main CI job.

**Verdict.** It must report 26 passed, 0 skipped, 0 xfailed. At freeze it
reports 0 passed, 26 xfailed: the 24 kernel cases fail only with
`ModuleNotFoundError` and the 2 command cases only with the missing
`Project.build_collar`.

**Cases.** The kernel cases use the four analytic junctions of ADR 0009,
with `n_loop` = 64, 15 surface layers `0.01 × 1.12^j` on each side, and 30
wall-normal layers `1e-4 × 1.2^k`. The collar is 0.118 thick, 12 times the
junction spacing.

| Case | Required | Prototype |
|---|---|---|
| K1 wall (× 4) | the composition above (<= 1e-12); rows on their bodies | exact |
| K2 validity (× 4) | 0 folded cells; min scaled Jacobian >= 0.05 | 0; 0.081 (tilted), 0.48-0.62 (others) |
| K3 first cell (× 4) | 0.95-1.05 × `steps_k[0]` | 0.984-1.001 |
| K4 envelope (× 4) | outer layer within 1 % of the thickness | 0.997-1.002 |
| K5 topology (× 4) | shape, periodicity, label, face BCs | |
| K6 interface | determinism; `ParameterError` for bad steps and `n_loop` < 3 | |
| C1, C2 command (STL box + post) | one valid block; journal, replay, undo | prototype on the same STL pair: 0 folded, min sJ 0.52, first cell 0.988-1.000, envelope 1.000 |

**Harness validation (T1).** The prototype was run behind the frozen kernel
interface. It passes 23 of the 26 cases. It fails exactly three:
- `steps_k` validation, which the prototype does not do;
- the two command cases, because `Project.build_collar` does not exist.

**Out of scope:**
- collars much thicker than 12 times the junction spacing at junctions
  sharper than 90 deg (the prototype folds the tilted junction at 32 times,
  1,203 cells);
- open intersection curves and symmetry planes;
- spacing laws along the loop;
- run time (the prototype takes about 12 s per collar with plain Jacobi
  iterations; the implementation is expected to use Picard + SuperLU, as in
  ADR 0007).

**For review (your call before implementation starts):**
1. The quality floor (min scaled Jacobian 0.05) has a 1.6× margin on the
   tilted junction and about 10× on the others. Tighten it per junction,
   or keep one floor?
2. Collars much thicker than 12 times the junction spacing at sharp
   junctions are out of scope. The likely fix is to offset in 3-D instead
   of per loop-station plane. Leave that for a follow-up?
3. Body B enters the project by file path inside the command, and the
   project still holds one body. Multi-component projects (named
   geometries) would be the general answer. Defer?

## Review decisions (2026-10-26, user)

1. **One quality floor** (min scaled Jacobian >= 0.05). Per-junction floors
   would freeze the prototype's numbers rather than a requirement; K2's
   fold count catches regressions.
2. **Thick collars at sharp junctions stay out of scope.** The likely fix,
   offsetting in 3-D instead of per loop-station plane, is a follow-up with
   its own prototype.
3. **Multi-component projects deferred.** Known limitation: the journal
   records body B by path, so a saved project replays only while that file
   is still there.

## Acceptance (2026-10-26): collar-block spec met

`src/shore/collar.py` implements the decision; `Project.build_collar`
(`BuildCollar`) wraps it. Differences from the prototype:
- Winslow solved by Picard iteration with one SuperLU factorisation per
  step for all loop-station sheets and the three coordinates at once,
  converged to 1e-12 of the collar thickness (58-85 steps), instead of
  4,000 Jacobi sweeps;
- the swallowtail search is vectorised;
- the junction spacing *h* is the smaller first surface step on either
  side; the out-of-plane part of the offset follows the B / A split also
  when `len(steps_a) != len(steps_b)`;
- argument checks (`n_loop`, `steps_*`; `m_a`, `m_b`, `nk` and a loaded
  body in the command). The collar mesh carries no body association: its
  wall lies on two bodies.

**Verdict for `tests/spec/test_collar.py`: 26 passed, 0 skipped,
0 xfailed** (T2; `scripts/check_spec.py` reports all files ok). All 10
`@PENDING` lines are removed. Full suite: 1324 passed, 0 xfailed.

Measured (T1, the spec's settings):

| Pair | Time | Folded | Min scaled Jacobian | First cell / `steps_k[0]` | Envelope / thickness |
|---|---|---|---|---|---|
| plane × cylinder | 1.7 s | 0 | 0.491 | 0.992-1.000 | 1.000-1.000 |
| tilted plane × cylinder | 2.5 s | 0 | 0.081 | 0.984-1.000 | 0.998-1.000 |
| cylinder × cylinder | 1.8 s | 0 | 0.532 | 0.985-1.000 | 0.998-1.002 |
| sphere × offset cylinder | 1.6 s | 0 | 0.660 | 0.990-1.001 | 0.997-1.001 |
| required | | 0 | >= 0.05 | 0.95-1.05 | 0.99-1.01 |

The fully converged solve gives the tilted junction the prototype's
0.081: the floor's 1.6× margin did not depend on stopping early.

**Unit tests** (`tests/test_collar.py`): the swallowtail cut and the
junction's image; a perturbed tensor-product sheet returns to the exact
TFI solution under Winslow with Thomas-Middlecoff (1e-9); the command
without a body.
