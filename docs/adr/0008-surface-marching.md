# 0008 — Hyperbolic surface marching: the collar-grid kernel

- **Status:** Proposed
- **Date:** 2026-10-26

## Context

Chimera assemblies of several components need **collar grids** where the
components meet: a wing–body junction, an appendage on a hull. A collar is
built from the intersection curve. A surface grid is marched away from
that curve over each component, and the result is extruded into a volume
block. ADR 0005 (Decision 3) planned this as hyperbolic surface marching
(SURGRD type; Chan & Steger 1992; Chan & Buning 1995), reusing the volume
scheme's ideas in two dimensions on the surface.

This ADR covers the core of that pipeline only: marching a structured grid
over a surface from a given curve. The intersection curve, two-sided
marching and the volume collar block come later, with their own specs
(user decision, 2026-10-26).

## Decision

`shore.surface.march.march_surface(geometry, curve, steps, *,
periodic=False, entity=BODY) -> (n, m + 1, 3)`:
- Layer j + 1 is marched from layer j to the `n × t` side, solving the
  linearised surface system as one sparse solve per layer:
  - `r_ξ·Δr = 0` (k-lines orthogonal to the layer);
  - `n·(r_ξ × Δr) = ΔA` (prescribed cell area);
  - `n·Δr = 0` (stay in the tangent plane).
- Each new layer is projected onto the surface. Only `project` (with
  normals) is needed, so STL bodies work as well as CAD faces.
- Stabilisation follows the volume scheme: smoothing of the area target,
  implicit and explicit dissipation ramped by `ρ = min(1, step / h)`, and a
  front-convergence sensor. The implementer chooses the parameters; no case
  passes them.
- Open curves march their ends freely. Periodic curves close on themselves.

## Alternatives considered

- **Pure normal marching (step along `n × t`, then project).** Rejected:
  fronts collide at concave corners. On the spec's concave L it folds 2
  cells.
- **Elliptic surface grids (ADR 0007) for collars.** Rejected as the
  primary tool: a collar needs orthogonal, wall-clustered grid lines
  leaving the intersection curve, with no opposite boundary to interpolate
  towards.

## Verification

### Acceptance spec (frozen 2026-10-26, before implementation)

**Spec.** `tests/spec/test_surface_march.py` runs in the main CI job and
needs no CAD extra.

**Verdict.** It must report 23 passed, 0 skipped, 0 xfailed. At freeze it
reports 0 passed, 23 xfailed, every case failing only with
`ModuleNotFoundError`.

**Rules.** As in ADR 0005: per-case `@PENDING` markers, with the hashes
in `tests/spec/SHA256SUMS`. The cases run on exact analytic surfaces
defined in the spec (plane, unit cylinder, unit sphere). Every case uses
the same schedule: 25 layers, `steps = 0.01 × 1.08^j`, total 0.731.

| Case | Required | Prototype |
|---|---|---|
| M1 line on a plane | exact parallel lines (1e-12) | 6e-16 |
| M2 ring on a cylinder, along the axis | exact rings (1e-12) | 2e-16 |
| M3 generator on a cylinder, around it | geodesic distance within 0.5 % | 0.056 % |
| M4 equator on a sphere, north (converging) | longitudes exact; latitude within 5 % | exact; 3.5 % |
| M5 helix on a cylinder | orthogonality mean <= 0.5 deg, max <= 1 deg; 0 folded | 0.02 / 0.04 deg |
| M6 concave L (90 deg) | pure normal folds; 0 folded, min angle >= 5 deg | pure 2 folded; 17 deg |
| M7 convex L | 0 folded, min angle >= 30 deg | 45 deg |
| M8 all 7 cases | shape, row 0 = curve, on the surface (1e-12), `n × t` side, <= 4 projected points per node and layer | 2 per node and layer |
| M9 interface | defaults, determinism, `ParameterError` for bad steps, bad curves and off-surface curves | |

**Harness validation (T1).** The prototype (`area_smooth` 0.5, `eps_i` 1,
`eps_e` 0.5, convergence gain 16) was run behind the frozen interface. It
passes all 15 exact, quality and invariant cases. It fails exactly the 8
interface-error cases, because it has no input validation. This shows that
the bounds are attainable and that the error cases discriminate.

**The robustness-versus-accuracy trade-off (measured).**
- Without area smoothing the prototype folds the concave L worse than pure
  normal marching: 14 cells against 2. Smoothing plus stronger dissipation
  removes the folds, reaching 17-44 deg.
- That same stabilisation lengthens the steps on converging fronts. On the
  sphere the latitude overshoots by 3.5 % (5.2 % with the stronger
  dissipation).
- The outer-edge excess correction of ADR 0005, applied with both signs,
  levels off at 1.6 % on the sphere, and the concave L folds again for
  bounds of 0.1 and above. The remaining error comes from the convergence
  sensor's dissipation, not from the area target.
- The spec therefore freezes 5 % on the sphere. This is the H2/H3 tension
  of ADR 0005 again.

**For review (your call before implementation starts):**
1. Is 5 % on a converging front acceptable for collars? It decides how far
   accuracy may trade against concave robustness. The alternative is a
   tighter bound with a smaller concave corner.
2. Open-curve ends march freely. A wing-root collar needs an end held on a
   symmetry plane or on another surface. Leave that for the collar-block
   spec?
3. The concave L requires a minimum angle of >= 5 deg; the prototype
   reaches 17 deg. Tighten it?
4. Is the work budget of 4 projected points per node and layer right? The
   need is 2: the normals and the final projection.

## Acceptance (2026-10-26): surface-marching spec met

`src/shore/surface/march.py` implements the decision. It uses the
prototype's defaults (`MarchParams`: area smoothing 0.5, `eps_implicit` 1,
`eps_explicit` 0.5, convergence gain 16) and no excess correction.

**Verdict for `tests/spec/test_surface_march.py`: 23 passed, 0 skipped,
0 xfailed** (T2; `scripts/check_spec.py` reports all files ok). All 12
`@PENDING` lines are removed; the parametrised cases share them. Full
suite: 1269 passed, 0 xfailed.

Measured (T1, the same defaults for every case):

| Case | Measured | Required | Prototype |
|---|---|---|---|
| M1 / M2 | 5.6e-16 / 1.7e-16 | <= 1e-12 | the same |
| M3 generator | 0.056 % | <= 0.5 % | 0.056 % |
| M4 sphere equator | 3.50 % | <= 5 % | 3.5 % |
| M5 helix | 0.018 / 0.036 deg | <= 0.5 / 1 | 0.02 / 0.04 |
| M6 concave L | 0 folded, 23.9 deg | >= 5 deg | 17.0 deg |
| M7 convex L | 0 folded, 45.0 deg | >= 30 deg | 45.0 deg |
| M8 work | 1.00 projected point per node and layer | <= 4 | 2 |

Two differences from the prototype:
- **Work.** Each layer's final projection also returns the normals of the
  next layer, so the geometry is queried once per node and layer
  (`tests/test_surface_march.py` pins this exactly).
- **Concave corner.** The contraction sensor reads the area-smoothed trial
  layer; the prototype read the unsmoothed one. This raises the concave
  corner's minimum angle from 17 to 23.9 deg.

The kernel also marches on an STL body (`StlGeometry` icosphere, unit
test).

The review questions above are still open, the 5 % bound on converging
fronts first among them. Changing any answer is a re-freeze. Next for
collars: the intersection curve between components and the collar block.

## Review decisions and re-freeze (2026-10-26, user)

1. **The 5 % bound on converging fronts is accepted**, as a documented
   robustness trade-off. Collar thickness tolerates a few percent of Chimera
   overlap.
2. **Open-curve end constraints move to the collar-block spec**, together
   with open intersection curves (ADR 0009, decision 3).
3. **Re-freeze:** the M6 concave-corner floor rises from 5 deg to
   **15 deg**. That is a real regression floor below the accepted 23.9 deg;
   5 deg would let a near-degenerate corner pass.
4. **Re-freeze:** the M8 work budget falls from 4 to **2** projected points
   per node and layer. That leaves room for one extra normal query but rules
   out iterative re-projection; the kernel uses 1.

Verdict after re-freeze: 23 passed, 0 skipped, 0 xfailed.
