# 0005 — Name the current marcher "algebraic"; add a true hyperbolic scheme

- **Status:** Proposed
- **Date:** 2026-10-25

## Context

`volume/hyperbolic.py` advances each layer as

$$
P_{i,j,k+1} = P_{i,j,k} + \Delta s_k\,\hat n_{i,j,k},
$$

with $\hat n$ from the cross product of i- and j-tangents, followed by
Laplacian blending. This is **algebraic normal extrusion with smoothing**. It
does not solve a hyperbolic system. In the literature and in the tools SHORE
users will compare against (HYPGEN, pyHyp, Pointwise), "hyperbolic grid
generation" means the Steger–Chaussee formulation. It marches the linearised
system

$$
\mathbf{x}_\xi\cdot\mathbf{x}_\zeta = 0,\quad
\mathbf{x}_\eta\cdot\mathbf{x}_\zeta = 0,\quad
\det\!\left[\frac{\partial(x,y,z)}{\partial(\xi,\eta,\zeta)}\right] = \Delta V,
$$

which enforces orthogonality and a prescribed cell volume. It is solved
implicitly in the marching direction with ADI, with explicit and implicit
dissipation and a volume-averaging that spreads fronts in concave regions
(Steger & Chaussee 1980; Chan & Steger 1992).

The difference matters exactly where CAD geometry takes SHORE: concave
corners, wing–body junctions and colliding fronts. There, normal extrusion
crosses itself and SHORE aborts on a negative Jacobian.

## Decision

1. **Rename, do not remove.** The current scheme becomes
   `march_scheme="algebraic"` and stays the default until the new scheme
   passes the verification below. Module and docs text stop calling it
   hyperbolic. `volume/hyperbolic.py` becomes `volume/algebraic.py`, with a
   re-export shim for one minor release.

2. **New `march_scheme="hyperbolic"`** in `volume/hyperbolic.py`:
   Steger–Chaussee with Chan–Steger enhancements:
   - cell-volume smoothing;
   - spatially variable dissipation;
   - grid-convergence sensors at concave corners;
   - boundary conditions for the seam ghost rows SHORE already builds
     (`ghost_rows`), constant-plane and symmetry planes.

   The implementation stays pure numpy + scipy: block-tridiagonal 3×3 solves
   batched over lines, with ADI in the two surface directions.

3. **Hyperbolic surface marching** (SURGRD-style), marching a curve over
   `Geometry` surfaces (ADR 0002), is planned as a separate feature for
   **collar grids** at component intersections. Multi-component Chimera
   assemblies (wing–body, appendages) need it. It reuses the same solver
   kernel in two dimensions on the surface.

## Alternatives considered

- **Keep the name, improve the algebraic smoother.** Rejected: the name
  misleads users coming from HYPGEN or pyHyp, and smoothing alone cannot
  enforce the volume constraint that controls concave-region behaviour.
- **Depend on pyHyp.** Rejected: it is Fortran + PETSc + CGNS and
  single-block-surface oriented. It would break the pip-installable
  constraint and does not know SHORE's seam ghost rows.

## Consequences

- Every topology marcher (`_march_caps`, `_march_cgrid_three`,
  `_march_ogrid`, the OH/CH wraps) dispatches on `march_scheme`. The seam
  machinery (ghost rows, Dirichlet pins, restitching) is shared.
- The `Spacing1D` k-schedule now prescribes the marching distance. The cell
  volume per node also derives from the surface cell area.

## Verification

- **Orthogonality:** on a sphere, the hyperbolic scheme's mean
  |cos(ξ, ζ)| at k=1 is no worse than the algebraic scheme.
- **Concave robustness:** on a 90° internal-corner body (L-section extrusion)
  where the algebraic scheme fails the Jacobian check, the hyperbolic scheme
  reaches the target thickness with positive Jacobian everywhere.
- **Reference:** grids from HYPGEN or pyHyp are not available as golden data,
  so acceptance rests on the metric tests above, with fixed inputs and
  thresholds frozen *before* implementation (ai-research-verification
  discipline).

## Finding (2026-10-26): smoothing folds thin first layers

While capturing phase-3 golden outputs, the current (algebraic) marcher
failed with a negative Jacobian at layer k=1 on a 6-block cubed sphere
(icosphere STL, `--ni 10 --nj 12`) whenever the first step was small:
`ds = 1e-3` geometric, or a tanh law whose first step is 5.6e-3. The
failure persists at `--smooth 0.05 --smooth-iters 1` and disappears with
`--smooth 0`.

**Mechanism (confirmed by instrumenting layer k=1 of sub-block 0, h = 0.32):**
two sweeps at strength 0.2 move nodes *along the normal* by up to
-3.7e-2, independent of the step. That is the Laplacian of positions on a
convex surface, curvature shrinkage of about `strength * h**2 * kappa`.
Their tangential motion is up to 1.8e-2. For any step below ~0.04, 16 nodes
end behind the wall. A first guess that the tangential motion was
responsible was wrong: the normal component dominates.

**Fix (smoothing limiter, `smooth_mode='limited'`, the default):**
`volume.hyperbolic._limit_smoothing` splits each node's smoothing motion
into normal and tangential parts:
- the normal part may pull the node back by at most 0.25 step;
- the tangential part is capped at 0.25 step;
- nodes within both limits keep the plain Laplacian result bit for bit.

Two rejected alternatives:
- Removing the normal component everywhere fixed the folds, but changed
  every grid. On the 60x60x10 cubed-sphere regression grid it raised the
  cap–equator C1 seam drift from 0.452 to 0.670 degrees.
- Limits (0.5, 1.0) left the regression grid untouched, but tilted coarse
  wall cells by 64–68 degrees.

The chosen (0.25, 0.25) is the tightest pair that leaves the regression
grid bit-identical.

| Case | `legacy` | `limited` |
|---|---|---|
| 10x12 sphere: cubed `ds` 1e-3, 1e-5; tanh first step 5.6e-3; O-grid `ds` 1e-3 | negative Jacobian at k=1 | marches; worst wall-cell tilt 19–24 deg (8 deg unsmoothed) |
| 60x60x10 cubed-sphere regression grid | — | bit-identical to `legacy` |
| `shore mesh` goldens with smoothing on (cubed, O-grid, OH at `ds` 0.05) | — | change: the limiter binds on those coarse grids |
| all 68 goldens with `--smooth-mode legacy` | byte-identical to the pre-fix CLI | — |

Regression tests: `tests/test_smoothing.py`. The limiter treats the
symptom inside the algebraic scheme. The hyperbolic scheme of this ADR
remains the structural answer for near-wall orthogonality.

## Acceptance spec (frozen 2026-10-26, before implementation)

`tests/spec/test_hyperbolic_marching.py` is the acceptance spec for the
hyperbolic scheme. It was written, and its helpers validated, before any
implementation.

**Verdict.** `pytest tests/spec` must report **13 passed, 0 skipped,
0 xfailed**. Until then the module is `xfail(strict=True)`. Deleting that
`pytestmark` line is the only permitted edit, and strict mode turns a
passing-but-still-marked suite into a failure.

**Frozen files.** Implementing work may not edit:
- the spec file and `tests/spec/__init__.py`;
- `tests/conftest.py`;
- the pytest configuration in `pyproject.toml`;
- `src/shore/quality.py`, which defines the metrics the thresholds are
  written in.

`tests/spec/SHA256SUMS` holds their hashes. The spec file has two
accepted hashes: frozen, and with the `pytestmark` line removed.
`scripts/check_spec.py` verifies the manifest and runs in CI's lint job.

**Cases and thresholds.** The spec's docstring carries the full table.
Every threshold sits next to the algebraic (limited) reference measured on
the same case with `shore.quality`:

| Case | Required |
|---|---|
| H1 coarse sphere, `ds` 1e-3 | no folds; wall orthogonality mean <= 3 deg, max <= 10 deg (algebraic 12.15 / 18.4); first cell 0.95–1.05 `ds` (algebraic down to 0.751); outer-layer distance 0.90–1.10 of the schedule; seams to round-off; bit-identical rerun |
| H2 60x60x10 regression grid | not worse than algebraic: mean <= 2.0, max <= 5.9 deg, min scaled Jacobian >= 0.70, thickness 0.95–1.05 |
| H3 plus prism, 4 concave 90-deg edges, cubed-sphere and O-grid | completes nk = 15 with 0 folded cells. Algebraic fails at k = 13 (legacy at k = 1); a case keeps checking that it still does |
| H4 3-block C-grid NACA 0012 | no folds; wall orthogonality mean <= 1.0, max <= 20 deg (TE); seams to round-off |
| H5 interface | `march_scheme` (default `"algebraic"`, unknown values raise), journaled by `Project.march`, CLI `--march-scheme` |

**Harness validation (T1, before freezing):**
- The helpers reproduce the algebraic baselines exactly on H1.
- The seam check passes on algebraic cubed-sphere and C-grid meshes. It
  detects a 1e-9 corruption on the C-grid wake cut, the one seam that is
  not a shared numpy view.

Validating the harness changed the spec once. The seam requirement was
first bitwise, but the existing algebraic C-grid wake cut already differs
at 1e-17. It is now "to round-off" (1e-12 × model diagonal).

**For review (your call before implementation starts):**
1. Is H1's max <= 10 deg right at the cubed-sphere 3-block corners, where
   the surface normal is ill-defined? Unsmoothed algebraic reaches 7.7 deg
   there.
2. H2 requires "no worse than algebraic" on a grid where algebraic is
   already good. Should hyperbolic be required to be *better* there?
3. H3 accepts outer-layer distances of 0.5–1.5 × the schedule. Fronts merge
   at concave edges, so the distance there is genuinely shorter. Is the
   band too generous?
4. H4 allows 20 deg at the trailing edge, a sharp edge with no unique
   normal. Should the TE columns be excluded and the rest held tighter?

Out of scope (stated in the spec): blunt-TE C-grids, CH/OH channel blocks
(they do not march), non-star-shaped CAD bodies, finer grids, run time.

## Implementation status (2026-10-26): spec NOT yet accepted

`volume/hyperbolic_scheme.py` implements the linearised Steger–Chaussee
step. Per block layer it solves

  `(I + C⁻¹A δξ + C⁻¹B δη − ε_i ∇²) Δr = r_ζ⁰ + ε_e T(∇² r)`

as one sparse system (SuperLU, deterministic). Its parts:
- volume smoothing, with the marching distance clipped to [½, 2] × step;
- a front-convergence sensor;
- implicit dissipation, and explicit dissipation restricted to the tangent
  plane;
- ramps of ρ (implicit) and ρ² (explicit), with ρ = min(1, step/h).

Edge rows with a neighbour block (ghost rows) are unknowns; unconnected and
pinned (cap) edges are Dirichlet. It is reached through
`Mesh.march(march_scheme="hyperbolic")`, `Project.march`,
`--march-scheme` and `[march] march_scheme`. The default stays `algebraic`,
and the 66 pipeline goldens are byte-identical. Defaults (chosen by a sweep
over all spec cases, none per case): `eps_implicit` 0.5, `eps_explicit`
0.25, `volume_smoothing` 0.2, `convergence_gain` 16, `distance_clip` 2.

Spec verdict: **10 of 13 pass, 3 fail** (T1: reproduced by scripts;
`pytest tests/spec` reports the 10 as strict XPASS failures).

| Case | Result |
|---|---|
| H3 plus prism, cubed-sphere and O-grid | pass: completes where algebraic fails |
| H4 C-grid | pass |
| H5 interface | pass |
| H1 validity, first cell, thickness, seams, determinism | pass |
| H1 orthogonality | fail: mean 5.9 / max 13.8 (required 3 / 10) |
| H2 orthogonality and validity | fail: max 9.6 (<= 5.9), min scaled Jacobian 0.69 (>= 0.70); interior max 2.4 vs algebraic 3.7 |
| H2 thickness | fail: 0.937 .. 1.075 (0.95 .. 1.05) |

Spec defects found during implementation, for a deliberate re-freeze:

1. **The orthogonality metric is wrong at block edges.**
   `quality.k_orthogonality_deg` takes the layer normal from one-sided
   differences on block-edge nodes, where the surface normal spans the
   seam. Unsmoothed algebraic marching is exactly orthogonal in the
   interior (0.00 deg), yet scores mean 9.4 / max 13.9 deg on edge nodes.
   On H1's coarse grid edge nodes dominate, so **H1's mean <= 3 / max
   <= 10 is unattainable for any seam-consistent scheme**. The angle to the
   exact sphere normal is <= 7.8 deg. H2's max <= 5.9 is the same artefact.
   Fix: measure against the body normal at the wall
   (`body.geometry.project(...).normals`), or interior nodes only.
2. **The H2 thickness floor is below what algebraic achieves.** The 0.95
   floor (algebraic 0.932) is missed by both schemes away from block
   corners. The cause is not the marching scheme and is not yet
   identified. Investigate before fixing the band.
3. **Module-level strict xfail is all-or-nothing.** Partial progress shows
   as failures and turns CI red. Mark the cases individually.

H2's min scaled Jacobian (0.69 < 0.70) and its 7 % thickness overshoot at
block corners (volume smoothing) are genuine shortfalls of the scheme, not
spec defects.

## Re-freeze (2026-10-26, user decision)

The spec was re-frozen as a deliberate, separate change to fix the three
defects listed above.

1. **Wall orthogonality is measured against the body's surface normal.** A
   new frozen metric, `quality.k_angle_deg`, measures the angle between the
   first k-segment and given directions. The spec takes the body's normals
   at the wall, orientation-free. Nodes with an undefined normal (the
   C-grid TE columns and the wing-tip stations) are excluded in H4. Before
   the change, H4 read 17.6 / 88.7 deg for every scheme alike.
2. **Per-case `@PENDING` strict-xfail markers** replace the module-level
   `pytestmark`. `scripts/check_spec.py` hashes spec modules with `@PENDING`
   lines removed: deleting one is the only permitted edit, and any other
   change is reported as CHANGED (checked both ways).
3. **H1's mean bound is anchored to pure normal extrusion**, computed live
   in the test: <= 1.15 × the mean of unsmoothed algebraic marching on the
   same grid. **This is a threshold change, not only a metric fix.** Under
   the corrected metric, pure normal extrusion scores 3.49 deg on this
   coarse faceted grid, because the discrete grid normal differs from the
   STL surface normal. The original "<= 3" would have excluded even an
   exact-normal march.

H2's references are re-measured under the new metric: algebraic 2.05 /
4.6 deg, so the bound is <= 2.0 / <= 4.6. H2's thickness band
(0.95–1.05) is deliberately unchanged until the 0.93 floor is explained.

Harness validation (T1): the helper reproduces the pure-normal reference
(3.49 / 7.8) exactly, and reads a forced 20-degree tilt as 20.00 / 20.0.

**Verdict at re-freeze: 11 passed, 2 xfailed.** Both H2 cases are pending:
- orthogonality and validity (max 8.1 vs <= 4.6; scaled Jacobian 0.69 vs
  >= 0.70);
- thickness (0.937 .. 1.075 vs 0.95 .. 1.05).

Acceptance remains 13 passed, 0 skipped, 0 xfailed.
