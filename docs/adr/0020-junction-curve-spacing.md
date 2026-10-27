# 0020 — Spacing laws along the junction curve

- **Status:** Proposed
- **Date:** 2026-10-27

## Context

A collar's j direction runs along the junction curve. ADR 0009 resamples
that curve with nodes equally spaced in arc length; ADR 0009's review moved
spacing laws to the collar step, and ADR 0010 and ADR 0011 left them out of
scope.

Equal arc-length spacing fails where the junction bends sharply. A
thick-wing-like junction, a plate pierced by an elliptic cylinder
`x^2 + (y / 0.2)^2 = 1` (curvature ratio 125 between the ends and the
flanks), gives the following with ADR 0010's collar steps and 128 nodes:
- a uniform loop turns 42.6° per cell at the leading and trailing ends;
- the collar's minimum scaled Jacobian is 0.008.

On an STL post of the same section, a uniform loop folds 4 collar cells.

## Investigation (2026-10-27, prototypes in the session scratchpad; T1)

**Combining the two controls.** The two controls compose in one map:
- the user's `Spacing1D` law acts on a computational coordinate;
- curvature clustering maps that coordinate to arc length.

So the law is matched exactly by arc length when the curvature weight is 0.

**Curvature distribution (equidistribution).** Nodes equidistribute the
monitor `m = 1 + alpha * kappa * L / (2 pi)` along the curve, where `L` is
the curve length:
- at `alpha = 0` the nodes are equally spaced;
- on a circle the monitor is constant, so the loop stays uniform for any
  `alpha`;
- the weight is dimensionless.

**Curvature estimate.** Three steps:
1. The traced polyline (step `4e-3` x the model diagonal) is too coarse at
   a leading edge of radius 0.04. It is first refined: segments next to a
   vertex that turns more than 2° are bisected, and the new points are
   polished onto both bodies. A minimum length of `1e-5` x the diagonal
   stops the refinement at real corners. Refining took the ellipse
   polyline from 151 to 348 points in 0.01 s.
2. Curvature at single vertices is useless on a faceted STL, where every
   facet edge is a spike. With `alpha = 3`, vertex curvature folded 7 cells
   on the STL post. So the curvature is averaged over half a uniform cell,
   with each vertex's turning spread over its half-segments and integrated
   exactly. An earlier count of vertices per window was noisy when the
   window was shorter than a segment, and the frozen spec's
   constant-curvature case caught it.
3. On an open curve the curvature is measured on the curve together with
   its mirror image. The ends on the symmetry plane therefore see the
   bend that the plane cuts through.

**Cell-ratio limit.** Plain equidistribution gave adjacent cells up to
1.9x apart. A limiter keeps them within 1.2x: forward and backward
min-sweeps in log space, then a rescale to the curve length.

Measured results, all with ADR 0010's collar steps:

| Case | Uniform | `alpha = 1` | `alpha = 3` | `alpha = 10` |
|---|---|---|---|---|
| analytic ellipse loop, n = 128: max turn / cell | 42.6° | 5.5° | 4.4° | 4.0° |
| analytic ellipse collar, n = 128: min scaled Jacobian | 0.008 | 0.23 | 0.25 | 0.27 |
| STL elliptic post collar, n = 96, half the steps: folded cells / min sJ | 4 / −0.017 | 0 / 0.15 | 0 / 0.18 | 0 / 0.05 |
| half ellipse (plane y = 0), n = 65: min sJ | 0.012 | 0.16 | — | — |
| ADR 0009's four pairs: cell spacing max / min | 1.00 | ≤ 1.04 | ≤ 1.04 | ≤ 1.04 |

**More notes:**
- **Accuracy of the law.** On the exact circle, the law is matched to
  2e-7 of the loop length.
- **Too few nodes.** With n = 64 on the analytic ellipse, the limit of 1.2x
  caps how much clustering the nodes allow. The ends then still turn about
  20° per cell, and one collar cell folds. The uniform loop folds 8 there.
- **Corners.** Junction curves with corners are out of scope: a sharp or a
  blunt trailing edge through a plate. The ADR 0009 tracer stalls at the
  corner and never closes the loop. Fixing that is a separate step.

## Decision

User decisions of 2026-10-27: `Spacing1D` laws plus curvature clustering;
the law lives in the curve kernel and the collar forwards it; exposed in the
API, the project and the CLI; the four acceptance criteria below.

- `intersection_loop(..., spacing=None, curvature=0.0)` and
  `intersection_curve(..., spacing=None, curvature=0.0)`:
  - `spacing` is a `Spacing1D`:
    - loops: `uniform` or `tanh2` (with an optional `cluster_at` in
      [0, 1]); `tanh` on a loop would jump in spacing at the seed;
    - open curves: `uniform`, `tanh` or `tanh2`.
  - `curvature` is the weight `alpha >= 0` above.
  - The defaults keep the ADR 0009 / 0011 resampling bit for bit (saved
    projects replay unchanged). `Spacing1D("uniform")` with `curvature = 0`
    is the same.
  - The new settings go in `TraceParams`: `refine_turn_deg = 2`,
    `refine_min = 1e-5`, `curvature_window = 0.5` and `max_cell_ratio = 1.2`.
- `build_collar(..., loop_spacing=None, loop_curvature=0.0)` forwards
  them.
- `Project.build_collar(..., loop_law="uniform", loop_beta=3.0,
  loop_cluster_at=None, loop_curvature=0.0)`. The fields are on
  `BuildCollar` and journaled when they are not the defaults.
- `shore collar` gains `--loop-law`, `--loop-beta`, `--loop-cluster-at` and
  `--loop-curvature`.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Collar-only redistribution of the uniform curve | Couples spacing to the collar. The kernel is the natural owner and stays reusable. |
| Clustering at user-given 3-D points | Larger interface. Curvature finds the ends of a section by itself. Deferred. |
| Curvature at polyline vertices | Spikes at STL facet edges fold collars at `alpha >= 3`. |
| No cell-ratio limit | Adjacent cells up to 1.9x apart. |
| Golden digests of the default outputs | One prototype value lies 1.6e-14 from a rounding boundary. A golden `.npz` compared at 1e-12 is robust across machines. |

## Consequences

- Collars around thick, rounded sections (struts, wing-body junctions with
  a rounded trailing edge) become valid at modest loop resolution:
  `--loop-curvature 1`.
- Combining a law and curvature is allowed. The law then acts on the
  curvature-equidistributed coordinate, so only `curvature = 0` matches the
  law by arc length.
- Junctions with corners still fail in the tracer (ADR 0009). That is the
  natural next step.

## Verification

### Acceptance spec (frozen 2026-10-27, before implementation)

**Spec.** `tests/spec/test_loop_spacing.py` runs in the main CI job. The
golden file `tests/spec/data/loop_default.npz` holds today's default
curves and is hashed in `SHA256SUMS`.

**Verdict.** It must report 29 passed, 0 skipped, 0 xfailed. At freeze it
reports 5 passed (the S2 golden comparisons: regression guards) and
24 xfailed.

**Harness validation.** With the prototype behind the interface, the spec
reported 29 passed after the curvature-window fix above. Before that fix,
S3 failed: the constant-curvature loop was off uniform by 2e-3, against a
2e-6 bound.

| Case | Required |
|---|---|
| S1 (× 4) | law matched by arc length on the exact circle, loop and open half (≤ 2e-6; prototype 2e-7); nodes on both bodies; node 0 unchanged; ends on the plane |
| S2 (× 6) | defaults equal the golden curves within 1e-12 (ADR 0009's four loops, the open half circle); `Spacing1D("uniform")` is the default, bit for bit |
| S3 (× 3) | ellipse loop `alpha = 1`: max turn ≤ 0.25 × uniform (0.13), cells ≤ 1.21 apart; half ellipse: end cells ≤ 0.25 × uniform (0.17); circle stays uniform (≤ 2e-6) |
| S4 (× 5) | the wall composes the spaced curve (full and half, ≤ 1e-12); ellipse collars valid with min scaled Jacobian ≥ 0.1: `alpha = 1` (0.23), `tanh2` + `alpha` (0.17), half (0.16) |
| S5 (× 4) | STL elliptic post: valid collar (min sJ ≥ 0.05, prototype 0.11; uniform folds 4), equal to the kernel, journal + bit-exact replay; `shore collar --loop-*` equals the API |
| S6 (× 7) | `ParameterError`: `tanh` on a loop, `geometric`, `pinned`, `cluster_at` outside [0, 1], negative or NaN curvature, unknown `loop_law` (mesh unchanged) |

**Out of scope:**
- junction curves with corners;
- clustering at 3-D points;
- run time.

### Review questions

1. **The default `alpha`.** It stays 0, so existing projects and journals
   are unchanged. Should `shore collar` instead recommend
   `--loop-curvature 1` in its help, or default to it in the CLI only?
2. **The cell-ratio limit** is fixed at 1.2 in `TraceParams`. Should it be
   exposed in the CLI and the project?
3. **Next step:** tracing junctions with corners (sharp and blunt trailing
   edges), which the ADR 0009 tracer cannot close.
