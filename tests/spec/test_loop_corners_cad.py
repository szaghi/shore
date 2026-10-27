"""FROZEN acceptance spec: junction curves with corners on CAD bodies (ADR 0021).

Needs the ``cad`` extra; the CI CAD job runs it.

Frozen 2026-10-27, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (CAD job): ``pytest tests/spec/test_loop_corners_cad.py``
reports 2 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 2 xfailed.

Interface under test: as ``test_loop_corners.py`` (``intersection_loop``
through corners), on ``OcctGeometry`` bodies.

Cases
-----
- K6 a build123d box plate (top face at z = 0) pierced by a NACA 0012 wing
  extruded from a section of two splines (upper and lower surface, 20
  points each from ``shore.profiles.naca``): sharp TE, the splines meet at
  the TE; blunt TE, a line closes the base. The loop (n = 64) closes within
  10 s (prototype 0.6 s and 0.8 s) and every TE point is a node within
  1e-7 x the plate's diagonal (the tessellation-based projection of
  ADR 0002; prototype 2e-11 sharp, exact blunt).

Out of scope of this spec: as ``test_loop_corners.py``; exact (tier-2)
projection.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

bd = pytest.importorskip("build123d")

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-27; implementation pending")


def cad_wing(closed_te: bool):
    from shore.geometry.occt import OcctGeometry
    from shore.profiles import naca

    section = naca("0012", n_pts=40, closed_te=closed_te)  # upper TE -> LE -> lower TE
    pts = [bd.Vector(float(x), float(y), -0.3) for x, y in section]
    half = len(pts) // 2
    edges = [bd.Spline(*pts[: half + 1]), bd.Spline(*pts[half:])]
    if not closed_te:
        edges.append(bd.Line(pts[-1], pts[0]))
    wing = bd.extrude(bd.Face(bd.Wire(edges)), amount=1.0)
    plate = bd.Box(4, 4, 1, align=(bd.Align.CENTER, bd.Align.CENTER, bd.Align.MAX))
    plate = plate.moved(bd.Location((0.5, 0.0, 0.0)))
    te = np.unique(np.c_[section[[0, -1]], [0.0, 0.0]], axis=0)
    return OcctGeometry(plate.wrapped), OcctGeometry(wing.wrapped), te


class TestCadWings:
    @pytest.mark.parametrize("closed_te", [True, False], ids=["sharp", "blunt"])
    def test_k6_cad_wing_through_a_plate(self, closed_te):
        from shore.surface.intersect import TraceParams, intersection_loop

        plate, wing, te = cad_wing(closed_te)
        start = time.perf_counter()
        fast = TraceParams(max_steps=2000)  # a stalled trace fails fast
        loop = intersection_loop(plate, wing, np.array([0.5, 0.06, 0.0]), 64, params=fast)
        assert time.perf_counter() - start <= 10.0
        diag = float(np.linalg.norm(np.diff(plate.bounds, axis=0)))
        dist = np.linalg.norm(loop[:, None] - te[None], axis=-1).min(0)
        assert dist.max() <= 1e-7 * diag
