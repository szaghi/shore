"""FROZEN acceptance spec: junction curves with corners (ADR 0021).

Frozen 2026-10-27, before any implementation. Implementing work may not edit
this file, ``tests/spec/data/loop_bends_default.npz``,
``tests/spec/data/loop_default.npz``, ``tests/spec/__init__.py``,
``tests/conftest.py``, the pytest configuration in ``pyproject.toml``, or
``src/shore/quality.py``. The ONE permitted edit is deleting an ``@PENDING``
line once that case passes (strict xfail turns a passing pending case into a
failure). ``tests/spec/SHA256SUMS`` holds the hashes, taken with ``@PENDING``
lines removed; ``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_loop_corners.py``
reports 22 passed, 0 skipped, 0 xfailed. At freeze: 10 passed (K3, the
golden smooth curves: regression guards), 12 xfailed. The CAD wings are in
``test_loop_corners_cad.py`` (CAD job).

Interface under test
--------------------
``intersection_loop`` and ``intersection_curve`` (ADR 0009, 0011, 0020),
unchanged in signature, now also on curves with corners: points where the
curve's tangent jumps because it crosses an edge (a crease) of either body,
convex or concave, at any angle the two branches still differ by at least
30 deg (``TraceParams.corner_deg``):

- the trace passes every corner and closes (loop) or reaches the plane
  (curve), within the run times below;
- every corner is a node, on the exact corner (to the trace tolerance);
- the nodes between fixed nodes (node 0 / the plane ends, and the corners)
  form segments; each segment gets cells in proportion to its length (its
  ADR 0020 monitor integral with ``curvature``), at least one, and the
  ``spacing`` law of ADR 0020 by itself: by default uniform per segment;
- curves without corners, tight smooth bends included (a leading edge that
  turns more than 30 deg per trace step), are unchanged bit for bit;
- ``ParameterError`` when *n* leaves a segment without a cell.

The corner cases pass ``params=TraceParams(max_steps=2000)`` (their traces
need fewer than 300 steps), so that a tracer that stalls on a corner fails
fast instead of spinning through the default 200 000 steps.

Cases
-----
Exact analytic bodies (frozen below): a plane through z-extruded prisms
(biconvex lens with two 22.8 deg tips; rectangle 1 x 0.1, the blunt base;
an L-shape with five convex corners and one concave), and a sphere centred
on an edge of a box (the curve crosses the edge of body A).

- K1 loops, n = 64: corners on nodes within 1e-12 x diagonal (prototype
  2e-15), each on its own node, every node on both bodies within 1e-9 x
  diagonal.
- K2 open curves: the rectangle cut by x = 0 (two corners), the lens cut
  by x = 0.25 (one tip): corners on nodes, ends on the plane.
- K3 unchanged (golden files written by the code before this ADR, within
  1e-12): ADR 0020's ``loop_default.npz`` (ADR 0009's four loops, the open
  half circle) and ``loop_bends_default.npz`` (tight bends: the elliptic
  junction of ADR 0020 as a loop, with curvature and as a half curve; its
  STL post; a faceted sphere x cylinder).
- K4 per-segment spacing on the rectangle: uniform and ``tanh2`` laws hold
  on every straight segment within 1e-9 (prototype 1e-15); cell counts
  within 1 of the segment's proportional share; curvature 1 on the lens
  shrinks the cells next to the corners to at most 0.5 x the uniform ones
  (prototype 0.35) with adjacent cells within 1.21 inside each segment;
  ``n`` too small raises, naming the corners.
- K5 STL wings (NACA 0012 from ``shore.profiles``, extruded, through a box
  plate): sharp TE, the TE a node within 1e-9 x diagonal (prototype exact);
  blunt TE, both base corners nodes; within 5 s each (prototype 0.5 s and
  0.8 s).

Out of scope of this spec: the collar at corners (single-block collars
fold there: measured in ADR 0021, a later ADR), a seed on a corner, a
curve meeting the plane at a corner, corners of less than 30 deg (traced
as smooth bends), cusps (branches 180 deg apart).
"""

from __future__ import annotations

import itertools
import time

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.mesh.spacing import Spacing1D
from shore.surface.intersect import TraceParams
from tests.spec.test_intersection import PAIRS as SMOOTH_PAIRS
from tests.spec.test_intersection import AnalyticBody, CylinderZ, Plane, Sphere
from tests.spec.test_loop_spacing import (
    ELLIPSE_SEED,
    GOLDEN,
    HALF_CIRCLE_SEED,
    HALF_ELLIPSE_SEED,
    SYM,
    EllipticCylinder,
)

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-27; implementation pending")

BENDS = GOLDEN.parent / "loop_bends_default.npz"
# corner cases: a trace needs < 300 steps; without corner handling it stalls on
# the corner and would spin through the default 200 000 before failing
FAST = {"params": TraceParams(max_steps=2000)}
DIAG = float(np.linalg.norm([4.0, 4.0, 4.0]))  # the analytic bodies' bounds


def loop_kernel():
    from shore.surface.intersect import intersection_loop

    return intersection_loop


def curve_kernel():
    from shore.surface.intersect import intersection_curve

    return intersection_curve


# ---------------------------------------------------------------------------
# Exact bodies with edges (frozen with the spec)
# ---------------------------------------------------------------------------


class PolygonPrism(AnalyticBody):
    """z-extruded simple polygon (CCW vertices): exact foot onto its sides."""

    def __init__(self, vertices) -> None:
        self.v = np.asarray(vertices, dtype=np.float64)
        self.e0, self.e1 = self.v, np.roll(self.v, -1, axis=0)
        d = self.e1 - self.e0
        self.len = np.linalg.norm(d, axis=1)
        self.dir = d / self.len[:, None]
        self.nrm = np.c_[self.dir[:, 1], -self.dir[:, 0]]  # outward for CCW

    def inside(self, q: np.ndarray) -> np.ndarray:
        x, y = q[..., 0, None], q[..., 1, None]
        x0, y0, x1, y1 = self.e0[:, 0], self.e0[:, 1], self.e1[:, 0], self.e1[:, 1]
        straddle = (y0 > y) != (y1 > y)
        with np.errstate(divide="ignore", invalid="ignore"):
            xc = (x1 - x0) * (y - y0) / (y1 - y0) + x0
        return (straddle & (x < xc)).sum(-1) % 2 == 1

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        q = p[..., :2]
        rel = q[..., None, :] - self.e0
        t = np.clip(np.einsum("...ej,ej->...e", rel, self.dir), 0.0, self.len)
        f = self.e0 + t[..., None] * self.dir
        k = np.linalg.norm(q[..., None, :] - f, axis=-1).argmin(-1)
        foot2 = np.take_along_axis(f, k[..., None, None], -2)[..., 0, :]
        tk = np.take_along_axis(t, k[..., None], -1)[..., 0]
        off = q - foot2
        r = np.linalg.norm(off, axis=-1, keepdims=True)
        sign = np.where(self.inside(q), -1.0, 1.0)[..., None]
        at_vertex = (r[..., 0] > 1e-12) & ((tk <= 0.0) | (tk >= self.len[k]))
        n2 = np.where(at_vertex[..., None], sign * off / np.maximum(r, 1e-300), self.nrm[k])
        n = np.concatenate([n2, np.zeros_like(n2[..., :1])], axis=-1)
        return np.concatenate([foot2, p[..., 2:]], axis=-1), n


class LensPrism(AnalyticBody):
    """z-extruded biconvex lens (two disks intersected): sharp tips at (+-a, 0)."""

    def __init__(self, a: float = 0.5, t: float = 0.05) -> None:
        self.R = (a * a + t * t) / (2.0 * t)
        self.d = self.R - t
        self.a = a

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        q = p[..., :2]
        R, d, a = self.R, self.d, self.a
        tip = np.where(q[..., :1] >= 0.0, [a, 0.0], [-a, 0.0])
        cands, norms = [], []
        for cy, side in ((-d, 1.0), (d, -1.0)):  # the upper arc is on the circle at (0, -d)
            v = q - [0.0, cy]
            u = v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-300)
            f = [0.0, cy] + R * u
            on_arc = side * f[..., 1:2] >= 0.0
            cands.append(np.where(on_arc, f, tip))
            norms.append(np.where(on_arc, u, np.nan))
        cands, norms = np.stack(cands, -2), np.stack(norms, -2)
        k = np.linalg.norm(q[..., None, :] - cands, axis=-1).argmin(-1)
        foot2 = np.take_along_axis(cands, k[..., None, None], -2)[..., 0, :]
        n2 = np.take_along_axis(norms, k[..., None, None], -2)[..., 0, :]
        off = q - foot2
        r = np.linalg.norm(off, axis=-1, keepdims=True)
        inside = (np.linalg.norm(q - [0.0, -d], axis=-1) < R) & (
            np.linalg.norm(q - [0.0, d], axis=-1) < R
        )
        sign = np.where(inside, -1.0, 1.0)[..., None]
        arc = np.where(q[..., 1:2] >= 0.0, (foot2 - [0.0, -d]) / R, (foot2 - [0.0, d]) / R)
        at_tip = np.where(r > 1e-12, sign * off / np.maximum(r, 1e-300), arc)
        n2 = np.where(np.isnan(n2), at_tip, n2)
        n = np.concatenate([n2, np.zeros_like(n2[..., :1])], axis=-1)
        return np.concatenate([foot2, p[..., 2:]], axis=-1), n


class Box(AnalyticBody):
    """The box [-2, 1] x [-2, 2] x [-1, 0]: exact foot onto its surface."""

    lo = np.array([-2.0, -2.0, -1.0])
    hi = np.array([1.0, 2.0, 0.0])

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        lo, hi = self.lo, self.hi
        inside = np.all((p > lo) & (p < hi), axis=-1)
        k = np.concatenate([p - lo, hi - p], axis=-1).argmin(-1)
        axis, upper = k % 3, k >= 3
        idx = np.indices(axis.shape)
        pushed = p.copy()
        pushed[(*idx, axis)] = np.where(upper, hi[axis], lo[axis])
        n_in = np.zeros_like(p)
        n_in[(*idx, axis)] = np.where(upper, 1.0, -1.0)
        foot = np.where(inside[..., None], pushed, np.clip(p, lo, hi))
        off = p - foot
        r = np.linalg.norm(off, axis=-1, keepdims=True)
        face = (np.abs(foot - hi) <= 1e-15).astype(float) - (np.abs(foot - lo) <= 1e-15)
        first = np.argmax(face != 0.0, axis=-1)
        n_on = np.zeros_like(p)
        n_on[(*idx, first)] = np.take_along_axis(face, first[..., None], -1)[..., 0]
        n = np.where(
            inside[..., None], n_in, np.where(r > 1e-12, off / np.maximum(r, 1e-300), n_on)
        )
        return foot, n


RECT = [[-0.5, -0.05], [0.5, -0.05], [0.5, 0.05], [-0.5, 0.05]]
L_SHAPE = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.0], [0.0, 0.0], [0.0, 0.5], [-0.5, 0.5]]
CORNER_PAIRS = {
    # name: (body A, body B, seed, exact corners)
    "lens": (Plane, LensPrism, [0.0, 0.06, 0.0], [[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]]),
    "rectangle": (Plane, lambda: PolygonPrism(RECT), [0.0, 0.06, 0.0], np.c_[RECT, [0.0] * 4]),
    "l_shape": (
        Plane, lambda: PolygonPrism(L_SHAPE), [-0.2, -0.55, 0.0], np.c_[L_SHAPE, [0.0] * 6]
    ),
    "box_edge": (
        Box, lambda: Sphere((1.0, 0.0, 0.0)), [0.6, 0.0, 0.05], [[1.0, 1.0, 0.0], [1.0, -1.0, 0.0]]
    ),
}  # fmt: skip


def corner_nodes(nodes: np.ndarray, corners) -> tuple[np.ndarray, np.ndarray]:
    """Index of the nearest node to each exact corner, and its distance."""
    d = np.linalg.norm(nodes[:, None] - np.asarray(corners, dtype=np.float64)[None], axis=-1)
    return d.argmin(0), d.min(0)


# ---------------------------------------------------------------------------
# K1-K2: the corners are traced and become nodes
# ---------------------------------------------------------------------------


class TestCorners:
    @pytest.mark.parametrize("name", list(CORNER_PAIRS))
    def test_k1_loop_through_corners(self, name):
        make_a, make_b, seed, corners = CORNER_PAIRS[name]
        a, b = make_a(), make_b()
        loop = loop_kernel()(a, b, np.asarray(seed, dtype=np.float64), 64, **FAST)
        assert loop.shape == (64, 3)
        index, dist = corner_nodes(loop, corners)
        assert dist.max() <= 1e-12 * DIAG
        assert len(set(index.tolist())) == len(index)
        assert a.distance(loop).max() <= 1e-9 * DIAG
        assert b.distance(loop).max() <= 1e-9 * DIAG

    @pytest.mark.parametrize(
        "name, seed, x0, corners",
        [
            ("rectangle", [0.2, 0.06, 0.0], 0.0, [[0.5, -0.05, 0.0], [0.5, 0.05, 0.0]]),
            ("lens", [0.3, 0.03, 0.0], 0.25, [[0.5, 0.0, 0.0]]),
        ],
    )
    def test_k2_open_curve_through_corners(self, name, seed, x0, corners):
        make_a, make_b, _, _ = CORNER_PAIRS[name]
        a, b = make_a(), make_b()
        plane = ((x0, 0.0, 0.0), (1.0, 0.0, 0.0))
        curve = curve_kernel()(a, b, np.asarray(seed, dtype=np.float64), 33, plane=plane, **FAST)
        assert curve.shape == (33, 3)
        _, dist = corner_nodes(curve, corners)
        assert dist.max() <= 1e-12 * DIAG
        assert np.abs(curve[[0, -1], 0] - x0).max() <= 1e-9 * DIAG
        assert b.distance(curve).max() <= 1e-9 * DIAG


# ---------------------------------------------------------------------------
# K3: curves without corners are unchanged
# ---------------------------------------------------------------------------


def _stl_bodies(name: str):
    from shore.geometry import StlGeometry

    if name == "stl_elliptic_post":
        box = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
        box.apply_translation((0.0, 0.0, -0.5))
        post = trimesh.creation.cylinder(radius=1.0, height=1.6, sections=256)
        post.apply_scale((0.5, 0.1, 1.0))
        post.apply_translation((0.0, 0.0, 0.2))
        return StlGeometry(box), StlGeometry(post), np.array([0.0, 0.11, 0.02]), 96
    sphere = trimesh.creation.icosphere(subdivisions=4, radius=1.0)
    cylinder = trimesh.creation.cylinder(radius=0.4, height=4.0, sections=48)
    return StlGeometry(sphere), StlGeometry(cylinder), np.array([0.4, 0.0, 0.9]), 48


class TestUnchanged:
    @pytest.mark.parametrize("name", list(SMOOTH_PAIRS))
    def test_k3_smooth_loops(self, name):
        make_a, make_b, seed, _, _ = SMOOTH_PAIRS[name]
        loop = loop_kernel()(make_a(), make_b(), np.asarray(seed, dtype=np.float64), 64)
        np.testing.assert_allclose(loop, np.load(GOLDEN)[name], rtol=0.0, atol=1e-12)

    def test_k3_smooth_open_curve(self):
        curve = curve_kernel()(Plane(), CylinderZ(), HALF_CIRCLE_SEED, 33, plane=SYM)
        golden = np.load(GOLDEN)["plane_cylinder_curve"]
        np.testing.assert_allclose(curve, golden, rtol=0.0, atol=1e-12)

    @pytest.mark.parametrize("name", ["ellipse", "ellipse_curvature", "half_ellipse"])
    def test_k3_tight_bends(self, name):
        a, b = Plane(), EllipticCylinder()
        if name == "half_ellipse":
            out = curve_kernel()(a, b, HALF_ELLIPSE_SEED, 65, plane=SYM)
        else:
            out = loop_kernel()(a, b, ELLIPSE_SEED, 128, curvature=float(name != "ellipse"))
        np.testing.assert_allclose(out, np.load(BENDS)[name], rtol=0.0, atol=1e-12)

    @pytest.mark.parametrize("name", ["stl_elliptic_post", "stl_sphere_cylinder"])
    def test_k3_faceted_bends(self, name):
        a, b, seed, n = _stl_bodies(name)
        loop = loop_kernel()(a, b, seed, n)
        np.testing.assert_allclose(loop, np.load(BENDS)[name], rtol=0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# K4: spacing per segment
# ---------------------------------------------------------------------------


def segments(nodes: np.ndarray, corners) -> list[np.ndarray]:
    """The loop's node rows between node 0 and the corners, ends included."""
    index, _ = corner_nodes(nodes, corners)
    fixed = [0, *sorted(index.tolist()), len(nodes)]
    closed = np.vstack([nodes, nodes[:1]])
    return [closed[i : j + 1] for i, j in itertools.pairwise(fixed)]


class TestSegments:
    @pytest.mark.parametrize(
        "law", [Spacing1D("uniform"), Spacing1D("tanh2", beta=2.0)], ids=["uniform", "tanh2"]
    )
    def test_k4_law_per_segment(self, law):
        make_a, make_b, seed, corners = CORNER_PAIRS["rectangle"]
        spacing = None if law.law == "uniform" else law
        loop = loop_kernel()(make_a(), make_b(), np.asarray(seed), 64, spacing=spacing, **FAST)
        rows = segments(loop, corners)
        lengths = np.array([np.linalg.norm(r[-1] - r[0]) for r in rows])  # straight segments
        share = lengths / lengths.sum() * 64
        for row, part in zip(rows, share, strict=True):
            cells = len(row) - 1
            assert cells >= 1 and abs(cells - part) < 1.0
            frac = np.linalg.norm(row - row[0], axis=1) / np.linalg.norm(row[-1] - row[0])
            assert np.abs(frac - law.build(cells + 1, 0.0, 1.0)).max() <= 1e-9

    def test_k4_curvature_clusters_at_corners(self):
        make_a, make_b, seed, corners = CORNER_PAIRS["lens"]
        a, b = make_a(), make_b()
        uniform = loop_kernel()(a, b, np.asarray(seed), 64, **FAST)
        loop = loop_kernel()(a, b, np.asarray(seed), 64, curvature=1.0, **FAST)

        def next_to_corners(nodes):
            h = np.linalg.norm(np.diff(np.vstack([nodes, nodes[:1]]), axis=0), axis=1)
            index, _ = corner_nodes(nodes, corners)
            return np.r_[h[index], h[index - 1]]

        assert next_to_corners(loop).max() <= 0.5 * next_to_corners(uniform).min()
        for row in segments(loop, corners):
            h = np.linalg.norm(np.diff(row, axis=0), axis=1)
            if len(h) > 1:
                assert np.max(np.maximum(h[1:] / h[:-1], h[:-1] / h[1:])) <= 1.21

    def test_k4_too_few_nodes(self):
        make_a, make_b, seed, _ = CORNER_PAIRS["rectangle"]
        with pytest.raises(ParameterError, match="corner"):  # the reason, not a failed trace
            loop_kernel()(make_a(), make_b(), np.asarray(seed), 4, **FAST)  # 5 segments


# ---------------------------------------------------------------------------
# K5: STL wings through a plate
# ---------------------------------------------------------------------------


def stl_wing(closed_te: bool):
    from shore.geometry import StlGeometry
    from shore.profiles import naca
    from shore.shapes import extrude_section

    section = naca("0012", n_pts=80, closed_te=closed_te)
    wing = extrude_section(section, span=1.0, n_stations=2, span_axis="z", end_caps=True)
    wing.apply_translation((0.0, 0.0, -0.3))  # chord along x from 0 to 1
    plate = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
    plate.apply_translation((0.5, 0.0, -0.5))  # top face at z = 0
    te = np.unique(np.c_[section[[0, -1]], [0.0, 0.0]], axis=0)
    return StlGeometry(plate), StlGeometry(wing), te


class TestStlWings:
    @pytest.mark.parametrize("closed_te", [True, False], ids=["sharp", "blunt"])
    def test_k5_wing_through_a_plate(self, closed_te):
        plate, wing, te = stl_wing(closed_te)
        start = time.perf_counter()
        loop = loop_kernel()(plate, wing, np.array([0.5, 0.06, 0.0]), 96, **FAST)
        assert time.perf_counter() - start <= 5.0
        diag = float(np.linalg.norm(np.diff(plate.bounds, axis=0)))
        _, dist = corner_nodes(loop, te)
        assert dist.max() <= 1e-9 * diag
