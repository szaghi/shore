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
