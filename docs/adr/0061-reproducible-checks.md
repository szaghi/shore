# 0061 — Checks that hold across machines and dependency versions

- **Status:** Accepted
- **Date:** 2026-11-01

## Context

The main CI job (`check`) has failed on every push since ADR 0053, so no
release since 0.5.0 reached PyPI. Investigated on 2026-11-01 (T1).

1. **F3 failed on every run.** `test_cap_folds.py` F3 (a capsule with a
   45° cap, hyperbolic scheme, must be refused) never raised on CI. CI
   installs the latest packages: numpy 2.5.4, scipy 1.18.1, trimesh 5.2.0.
   - trimesh 5 tessellates `trimesh.creation.capsule` differently (226 →
     258 vertices). So `shore.shapes.capsule` builds a different body, and
     the case's worst cell goes from −0.0056 to +0.0040.
   - Swapping trimesh alone flips it, in a venv with CI's stack.
   - The case was marginal: a fold of −0.006 is within reach of any
     upstream change.
2. **The golden checks failed on two runs.** F6, H3, G5, ADR 0054's C5 and
   C6 failed on two runs with the same package versions as passing runs,
   and pass locally with CI's stack.
   - Locally, `OPENBLAS_CORETYPE=PRESCOTT` reproduces F6. The grid then
     differs from the golden by 1.45e-15 of its size.
   - numpy's SIMD groups (`NPY_DISABLE_CPU_FEATURES=X86_V3`) change
     nothing.
   - So the BLAS kernel OpenBLAS picks for the runner's CPU changes the
     last bits. The failing runs were also twice as fast, which fits
     different hardware.
3. **`Project.load` compares bit for bit.** It replays the journal and
   requires every grid to be **bit-identical** to the saved one (ADR 0004).
   So a project saved on one machine can fail to load on another. That is
   C5 and C6, and also the workstation-to-cluster case.

**The owner asked for all three fixes, re-freezes approved (2026-11-01).**

## Decision

1. **A: the capsule is frozen.** `test_cap_folds.py` reads the capsule from
   `data/cap_folds_capsule.stl` (`shapes.capsule` with trimesh 4.11.5,
   the body the spec was validated on). F3 moves to a 50° cap at the same
   size: its worst cell is −0.32, and −0.30 on trimesh 5's tessellation
   too.
2. **B: grids are compared to 1e-10.**
   - The golden checks (F6, G5, H3 and ADR 0058's P5) use
     `assert_allclose(rtol=0, atol=1e-10)`, the convention ADR 0034's R8
     already used. The bodies are about unit size, and the measured
     machine difference is about 1e-15.
   - `Project.load` allows each node to differ from its saved grid by up to
     1e-10 of the block's largest bounding-box extent
     (`project.GRID_TOLERANCE`). The refusal message gives the largest
     deviation.
   - A replay is still a determinism check: any change to an algorithm
     moves nodes by far more than 1e-10.
3. **C: CI is reproducible, and the latest packages are watched.**
   - The blocking jobs install with `ci/constraints.txt`, which pins the
     packages the results depend on to the versions the goldens were made
     with.
   - The `check` job prints numpy's runtime and BLAS configuration, so a
     failing run shows its kernels.
   - A separate `test-latest` job runs the suite on the newest packages
     with `continue-on-error`. An upstream change then shows up as a
     warning, not as a blocked release. Bumping the constraints is a
     deliberate commit.

## Alternatives considered

- **Pin `trimesh<5` for users:** it would hold every user back, and the
  marginal case would still be marginal.
- **Disable BLAS kernel dispatch in CI** (`OPENBLAS_CORETYPE`): fixes CI
  but not `Project.load` across machines.
- **Re-make the goldens on CI:** they would then fail on the developers'
  machines.

## Consequences

- CI's `check` job passes again, and the next release (0.7.1) reaches
  PyPI.
- `shore body capsule` still depends on the installed trimesh version.
  Saved projects are safe: they copy their input files in.
- Specs written from now on compare grids made on one machine with grids
  made on another to a tolerance.

## Verification

- `tests/spec/test_cap_folds.py`, `test_marcher_folds.py`,
  `test_ch_handedness.py` and `test_tip_cap_guard.py`, re-frozen, pass in
  three settings: the local stack, the local stack with
  `OPENBLAS_CORETYPE=PRESCOTT`, and CI's stack (numpy 2.5.4, scipy 1.18.1,
  trimesh 5.2.0).
- `tests/spec/test_load_tolerance.py` (new, main job), frozen before the
  implementation:
  - L1: grids moved by 1e-13 of their size load;
  - L2: moved by 1e-6, they are refused (a guard);
  - L3: the refusal names the block and the largest deviation;
  - L4: an unchanged project loads (a guard).

  At freeze: 2 passed, 2 xfailed.
