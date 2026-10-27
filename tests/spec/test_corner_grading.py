"""FROZEN acceptance spec: graded spacing across corners, collars at blunt corners (ADR 0022).

Frozen 2026-10-27, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_corner_grading.py``
reports 13 passed, 0 skipped, 0 xfailed. At freeze: 5 passed (the box-edge
loop, graded already, and the four collar cases: regression guards that pin
what the single block does at blunt corners), 8 xfailed.

Interface under test
--------------------
``intersection_loop`` / ``intersection_curve`` with corners (ADR 0021) and
``curvature > 0`` (ADR 0020), and ``build_collar`` on such loops:

- with ``curvature > 0`` the nodes follow one graded spacing over the whole
  curve: adjacent cells differ by at most ``TraceParams.max_cell_ratio``
  (1.2) everywhere, across the corners and node 0 included (plus the
  discretisation: bound 1.22 on chords); every corner is still a node;
- ``curvature = 0`` keeps ADR 0021's uniform cells per segment (its spec);
  curves without corners keep ADR 0020's distribution (ADR 0021's golden);
- the single collar block is valid at convex corners of up to 90 deg within
  the thickness measured in ADR 0022.

Cases
-----
- G1 grading, ``curvature`` 1 and 3: the NACA 0012 STL wing (ADR 0021's
  K5) through a plate, blunt and sharp TE, n = 96 (prototype 1.198; before:
  8.7 and 6.6 on the blunt wing); the analytic lens, rectangle, L-shape and
  box edge of ADR 0021 (n = 64, prototype <= 1.207); the open rectangle cut
  at x = 0 (n = 33, prototype 1.207): adjacent chords within 1.22, corners
  on nodes within 1e-12 x diagonal (STL: 1e-9).
- G2 collars at convex 90 deg corners (ADR 0010's surface steps
  ``0.005 * 1.12**j``, 15 layers; wall-normal ``5e-5 * 1.2**k``): the plane x
  rectangle prism and the box [-2, 1] x [-2, 2] x [-1, 0] x a sphere of radius
  0.5 on its edge (corners on body A), n_loop 96; uniform loop, 25 layers
  (4.7 x the junction spacing): 0 folded cells, min scaled Jacobian >= 0.025
  (prototype 0.053 and 0.148); curvature 1, 20 layers (1.9 x): >= 0.1
  (prototype 0.28 and 0.22).

Out of scope of this spec: thicker collars at corners (the rectangle folds
13 cells at 11.8 x, ADR 0022), concave corners (the L-shape folds at
1.9 x), sharp corners, the STL wing collar (the surface marcher of ADR 0008
flips at a clustered leading edge: a later ADR), ADR 0020's distribution on
smooth curves (pinned by ADR 0021's golden file).
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.quality import scaled_jacobian
from tests.spec.test_intersection import AnalyticBody
from tests.spec.test_loop_corners import (
    CORNER_PAIRS,
    DIAG,
    FAST,
    Box,
    corner_nodes,
    stl_wing,
)

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-27; implementation pending")

STEPS_S = 0.005 * 1.12 ** np.arange(15)


def kernels():
    from shore.surface.intersect import intersection_curve, intersection_loop

    return intersection_loop, intersection_curve


def max_ratio(nodes: np.ndarray, closed: bool) -> float:
    pts = np.vstack([nodes, nodes[:1]]) if closed else nodes
    h = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    h = np.r_[h, h[:1]] if closed else h
    return float(np.max(np.maximum(h[1:] / h[:-1], h[:-1] / h[1:])))


class Ball(AnalyticBody):
    """Sphere of radius *r* about *centre* (exact)."""

    def __init__(self, centre=(1.0, 0.0, 0.0), r: float = 0.5) -> None:
        self.c, self.r = np.asarray(centre, dtype=np.float64), r

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        q = p - self.c
        n = q / np.linalg.norm(q, axis=-1, keepdims=True)
        return self.c + self.r * n, n


def check_analytic(name: str) -> None:
    loop_kernel, _ = kernels()
    make_a, make_b, seed, corners = CORNER_PAIRS[name]
    loop = loop_kernel(make_a(), make_b(), np.asarray(seed), 64, curvature=1.0, **FAST)
    assert max_ratio(loop, closed=True) <= 1.22
    assert corner_nodes(loop, corners)[1].max() <= 1e-12 * DIAG


# ---------------------------------------------------------------------------
# G1: graded across the corners
# ---------------------------------------------------------------------------


class TestGrading:
    @pytest.mark.parametrize("curvature", [1.0, 3.0])
    @pytest.mark.parametrize("closed_te", [False, True], ids=["blunt", "sharp"])
    def test_g1_stl_wing(self, closed_te, curvature):
        loop_kernel, _ = kernels()
        plate, wing, te = stl_wing(closed_te)
        loop = loop_kernel(plate, wing, np.array([0.5, 0.06, 0.0]), 96, curvature=curvature, **FAST)
        assert max_ratio(loop, closed=True) <= 1.22
        diag = float(np.linalg.norm(np.diff(plate.bounds, axis=0)))
        assert corner_nodes(loop, te)[1].max() <= 1e-9 * diag

    @pytest.mark.parametrize("name", ["lens", "rectangle", "l_shape"])
    def test_g1_analytic_loops(self, name):
        check_analytic(name)

    def test_g1_box_edge_loop(self):
        check_analytic("box_edge")  # graded already: no segment much shorter than its cells

    def test_g1_open_curve(self):
        _, curve_kernel = kernels()
        make_a, make_b, _, _ = CORNER_PAIRS["rectangle"]
        plane = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
        curve = curve_kernel(
            make_a(), make_b(), np.array([0.2, 0.06, 0.0]), 33, plane=plane, curvature=1.0, **FAST
        )
        assert max_ratio(curve, closed=False) <= 1.22
        corners = [[0.5, -0.05, 0.0], [0.5, 0.05, 0.0]]
        assert corner_nodes(curve, corners)[1].max() <= 1e-12 * DIAG


# ---------------------------------------------------------------------------
# G2: collars at convex 90 deg corners
# ---------------------------------------------------------------------------

BLUNT = {
    "rectangle": (CORNER_PAIRS["rectangle"][0], CORNER_PAIRS["rectangle"][1], [0.0, 0.06, 0.0]),
    "box_edge": (Box, Ball, [0.7, 0.0, 0.05]),
}


class TestBluntCollars:
    @pytest.mark.parametrize(
        "nk, curvature, bound", [(25, 0.0, 0.025), (20, 1.0, 0.1)], ids=["uniform", "curvature"]
    )
    @pytest.mark.parametrize("name", list(BLUNT))
    def test_g2_valid_collar(self, name, nk, curvature, bound):
        from shore.collar import build_collar

        make_a, make_b, seed = BLUNT[name]
        steps_k = 5e-5 * 1.2 ** np.arange(nk)
        block = build_collar(
            make_a(), make_b(), np.asarray(seed), n_loop=96, steps_a=STEPS_S, steps_b=STEPS_S,
            steps_k=steps_k, loop_curvature=curvature,
        )  # fmt: skip
        jac = scaled_jacobian(np.concatenate([block.nodes, block.nodes[:, :1]], axis=1))
        assert (jac <= 0.0).sum() == 0
        assert jac.min() >= bound
