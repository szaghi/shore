"""FROZEN acceptance spec: intersection loop of two bodies, for collar grids (ADR 0009).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_intersection.py``
reports 13 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 13 xfailed (the
kernel module does not exist yet).

Interface under test
--------------------
``shore.surface.intersect.intersection_loop(body_a, body_b, seed, n, *,
entity_a=BODY, entity_b=BODY)`` returns ``(n, 3)``:

- the closed intersection loop of the two bodies through the junction
  nearest *seed* (a point near the curve), resampled to *n* nodes, n >= 3,
  equally spaced in arc length; the last node is not repeated;
- every node on both bodies;
- ordered along ``t = n_b x n_a`` (the bodies' outward normals), so that
  ``march_surface(body_a, loop, ...)`` marches on A away from B and
  ``march_surface(body_b, loop[::-1], ...)`` on B away from A: the two halves
  of a collar;
- ``loop[0]`` within one node of the loop's point nearest *seed*;
- computed with ``Geometry.project`` (with normals) only, so it works for STL,
  CAD and analytic bodies alike;
- ``ParameterError`` for n < 3, a *seed* that is not a 3-vector, and no
  intersection near *seed*.

The bodies are exact analytic ``Geometry`` implementations defined here
(numpy only). Every loop case uses n = 64.

Thresholds and references (prototype measurements, ADR 0009):

- X1 plane x cylinder (circle), X2 tilted plane x cylinder (ellipse),
  X3 cylinder x cylinder (two loops; the seed's), X4 sphere x offset cylinder
  (two loops; the seed's): every node on both bodies within 1e-9 x diagonal
  (prototype 1e-15); chord lengths max / min <= 1.01 (prototype 1.0003);
  polygon length within 0.5 % of the exact curve length; every node on the
  seed's loop; ``loop[0]`` at the seed.
- X5 each pair: one ``march_surface`` layer on A from the loop lies outside B,
  and on B from the reversed loop outside A.
- X6 interface: determinism, errors.

Out of scope of this spec: tangential (grazing) intersections, curves that
end on a boundary (open intersection curves), more than two bodies, spacing
laws other than uniform, run time.
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.geometry import BODY, EntityKind, EntityRef, Projection

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

N = 64


def kernel():
    """The function under test (imported late: absent until implemented)."""
    from shore.surface.intersect import intersection_loop

    return intersection_loop


# ---------------------------------------------------------------------------
# Exact analytic bodies (frozen with the spec)
# ---------------------------------------------------------------------------


class AnalyticBody:
    """An exact surface behind the ``Geometry`` protocol (BODY only)."""

    @property
    def bounds(self) -> np.ndarray:
        return np.array([[-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]])

    def entities(self, kind: EntityKind) -> tuple[EntityRef, ...]:
        return (BODY,) if kind is EntityKind.BODY else ()

    def ray_cast(self, *args: object, **kwargs: object) -> object:
        raise NotImplementedError("not used")

    def evaluate(self, *args: object, **kwargs: object) -> object:
        raise NotImplementedError("not used")

    def fingerprint(self, *args: object) -> object:
        raise NotImplementedError("not used")

    def project(
        self, points: np.ndarray, entity: EntityRef | None = None, *, exact: bool = False
    ) -> Projection:
        pts = np.asarray(points, dtype=np.float64)
        foot, normal = self.foot(pts)
        return Projection(foot, np.linalg.norm(pts - foot, axis=-1), normal, None, None)

    def distance(self, p: np.ndarray) -> np.ndarray:
        return np.linalg.norm(p - self.foot(p)[0], axis=-1)

    def side(self, p: np.ndarray) -> np.ndarray:
        """Signed distance, > 0 outside (along the outward normal)."""
        foot, normal = self.foot(p)
        return np.einsum("...j,...j->...", p - foot, normal)


class Plane(AnalyticBody):
    """The plane through the origin with unit normal *normal*."""

    def __init__(self, normal=(0.0, 0.0, 1.0)) -> None:
        self.n = np.asarray(normal, dtype=np.float64) / np.linalg.norm(normal)

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        d = p @ self.n
        return p - d[..., None] * self.n, np.broadcast_to(self.n, p.shape).copy()


class CylinderZ(AnalyticBody):
    """Cylinder of radius *r* about the z-parallel axis through (cx, 0)."""

    def __init__(self, r: float = 1.0, cx: float = 0.0) -> None:
        self.r, self.c = r, np.array([cx, 0.0, 0.0])

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        q = p - self.c
        rr = np.linalg.norm(q[..., :2], axis=-1, keepdims=True)
        n = np.concatenate([q[..., :2] / rr, np.zeros_like(rr)], axis=-1)
        return self.c + self.r * n + q * [0.0, 0.0, 1.0], n


class CylinderX(AnalyticBody):
    """Cylinder of radius *r* about the x axis."""

    def __init__(self, r: float = 0.5) -> None:
        self.r = r

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        rr = np.linalg.norm(p[..., 1:], axis=-1, keepdims=True)
        n = np.concatenate([np.zeros_like(rr), p[..., 1:] / rr], axis=-1)
        return self.r * n + p * [1.0, 0.0, 0.0], n


class Sphere(AnalyticBody):
    """Sphere of radius 1 about *centre*."""

    def __init__(self, centre=(0.0, 0.0, 0.0)) -> None:
        self.c = np.asarray(centre, dtype=np.float64)

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        q = p - self.c
        n = q / np.linalg.norm(q, axis=-1, keepdims=True)
        return self.c + n, n


# ---------------------------------------------------------------------------
# The pairs, their exact curves and the checks (frozen with the spec)
# ---------------------------------------------------------------------------


def _theta(m: int = 200001) -> np.ndarray:
    return np.linspace(0.0, 2.0 * np.pi, m)


def _length(curve: np.ndarray) -> float:
    return float(np.linalg.norm(np.diff(curve, axis=0), axis=1).sum())


def exact_circle() -> np.ndarray:
    t = _theta()
    return np.c_[np.cos(t), np.sin(t), 0 * t]


def exact_ellipse() -> np.ndarray:
    t = _theta()
    return np.c_[np.cos(t), np.sin(t), 0.5 * np.cos(t)]


def exact_steinmetz() -> np.ndarray:  # x > 0 branch of x^2 + y^2 = 1, y^2 + z^2 = 1/4
    t = _theta()
    y, z = 0.5 * np.cos(t), 0.5 * np.sin(t)
    return np.c_[np.sqrt(1.0 - y * y), y, z]


def exact_sphere_cylinder() -> np.ndarray:  # upper branch
    t = _theta()
    x, y = 0.3 + 0.5 * np.cos(t), 0.5 * np.sin(t)
    return np.c_[x, y, np.sqrt(1.0 - x * x - y * y)]


PAIRS = {
    "plane_cylinder": (
        Plane, CylinderZ, [1.1, 0.1, 0.05], exact_circle, lambda p: np.ones(len(p), bool)
    ),
    "tilted_plane_cylinder": (
        lambda: Plane((-0.5, 0.0, 1.0)), CylinderZ, [1.0, 0.1, 0.6], exact_ellipse,
        lambda p: np.ones(len(p), bool),
    ),
    "cylinder_cylinder": (
        CylinderZ, CylinderX, [1.05, 0.05, 0.45], exact_steinmetz, lambda p: p[:, 0] > 0.0
    ),
    "sphere_cylinder": (
        Sphere, lambda: CylinderZ(0.5, 0.3), [0.8, 0.0, 0.6], exact_sphere_cylinder,
        lambda p: p[:, 2] > 0.0,
    ),
}  # fmt: skip


def loop_of(name: str) -> tuple[AnalyticBody, AnalyticBody, np.ndarray, np.ndarray]:
    make_a, make_b, seed, _, _ = PAIRS[name]
    a, b = make_a(), make_b()
    seed = np.asarray(seed, dtype=np.float64)
    return a, b, seed, kernel()(a, b, seed, N)


# ---------------------------------------------------------------------------
# X1-X4: the loop
# ---------------------------------------------------------------------------


class TestLoop:
    @pytest.mark.parametrize("name", list(PAIRS))
    def test_x1_x4_loop(self, name):
        a, b, seed, loop = loop_of(name)
        _, _, _, exact, on_branch = PAIRS[name]
        assert loop.shape == (N, 3) and loop.dtype == np.float64
        diag = float(np.linalg.norm(np.ptp(a.bounds, axis=0)))
        assert a.distance(loop).max() <= 1e-9 * diag
        assert b.distance(loop).max() <= 1e-9 * diag
        chords = np.linalg.norm(np.roll(loop, -1, axis=0) - loop, axis=1)
        assert chords.max() / chords.min() <= 1.01
        assert abs(chords.sum() / _length(exact()) - 1.0) <= 0.005
        assert on_branch(loop).all()
        nearest = int(np.argmin(np.linalg.norm(loop - seed, axis=1)))
        assert nearest in (0, 1, N - 1)


# ---------------------------------------------------------------------------
# X5: orientation, composed with the surface marcher (ADR 0008)
# ---------------------------------------------------------------------------


class TestCollarSides:
    @pytest.mark.parametrize("name", list(PAIRS))
    def test_x5_marching_leaves_the_other_body(self, name):
        from shore.surface.march import march_surface

        a, b, _, loop = loop_of(name)
        steps = np.full(2, 0.02)
        on_a = march_surface(a, loop, steps, periodic=True)
        on_b = march_surface(b, loop[::-1], steps, periodic=True)
        assert (b.side(on_a[:, 1]) > 0.0).all()
        assert (a.side(on_b[:, 1]) > 0.0).all()


# ---------------------------------------------------------------------------
# X6: interface
# ---------------------------------------------------------------------------


class TestInterface:
    def test_x6_deterministic(self):
        a, b, seed = Sphere(), CylinderZ(0.5, 0.3), np.array([0.8, 0.0, 0.6])
        first = kernel()(a, b, seed, N)
        np.testing.assert_array_equal(first, kernel()(a, b, seed, N, entity_a=BODY, entity_b=BODY))

    def test_x6_too_few_nodes(self):
        with pytest.raises(ParameterError):
            kernel()(Plane(), CylinderZ(), np.array([1.1, 0.1, 0.05]), 2)

    @pytest.mark.parametrize("seed", [np.zeros(2), np.zeros((2, 3))])
    def test_x6_bad_seed(self, seed):
        with pytest.raises(ParameterError):
            kernel()(Plane(), CylinderZ(), seed, N)

    def test_x6_no_intersection(self):
        with pytest.raises(ParameterError):
            kernel()(Plane(), Sphere((0.0, 0.0, 5.0)), np.array([0.0, 0.0, 2.0]), N)
