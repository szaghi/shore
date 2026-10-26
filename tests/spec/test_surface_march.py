"""FROZEN acceptance spec: hyperbolic surface marching, the collar-grid kernel (ADR 0008).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_surface_march.py``
reports 23 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 23 xfailed (the
kernel module does not exist yet).

Interface under test
--------------------
``shore.surface.march.march_surface(geometry, curve, steps, *, periodic=False,
entity=BODY)`` returns the surface grid as an array ``(n, m + 1, 3)``:

- *curve* ``(n, 3)``, n >= 3, lies on the surface (on *entity* of *geometry*;
  only ``project`` with its normals is needed, so an STL body works); with
  *periodic* it is a closed loop (the last node is not repeated).
- *steps* ``(m,)``, m >= 1, positive: the marching distance of each layer.
- Row ``[:, 0]`` is the curve; every node lies on the surface; layer j + 1 is
  marched from layer j to the side ``n x t`` (``n`` the surface normal, ``t``
  the curve direction), so that ``(e_i x e_j) . n > 0``.
- Hyperbolic (Steger-Chaussee type on the surface): k-lines orthogonal to the
  layers, prescribed cell areas, with smoothing and dissipation of the
  implementer's choice, not passed per case.
- ``ParameterError`` for a malformed *curve* or *steps*, and for a curve
  farther than 1e-6 x the model diagonal from the surface.
- At most 4 points projected per node per layer (no iterative re-projection).

The surfaces are exact analytic ``Geometry`` implementations defined here
(numpy only). Every case marches the same schedule: 25 layers,
``steps = 0.01 * 1.08**j``, total 0.731.

Thresholds and references (prototype measurements, ADR 0008; "pure normal" =
each node stepped along ``n x t`` and projected, computed in the cases):

- M1 straight line on a plane: exact parallel lines (<= 1e-12).
- M2 ring on a cylinder, marching along the axis: exact rings (<= 1e-12).
- M3 generator of a cylinder, marching around it: geodesic distance within
  0.5 % of the schedule (prototype 0.056 %), z unchanged (<= 1e-12).
- M4 equator of a sphere, marching north (fronts converge): longitudes kept
  (<= 1e-12), latitude within 5 % of the schedule (prototype 3.5 %).
- M5 helix on a cylinder (no symmetry): k-lines orthogonal to the layers,
  |90 deg - angle| mean <= 0.5, max <= 1 (prototype 0.02 / 0.04); 0 folded.
- M6 concave L on a plane (fronts collide at 90 deg): pure normal folds
  (prototype 2 cells); 0 folded, min corner angle >= 5 deg (prototype 17).
- M7 convex L on a plane (fronts fan out): 0 folded, min corner angle
  >= 30 deg (prototype 45).
- M8 every case: shape, row 0 = curve, nodes on the surface (<= 1e-12 diag),
  orientation, work budget.
- M9 interface: defaults, determinism, errors.

Out of scope of this spec: intersection curves between components, marching
on both sides of a curve, symmetry-plane or other end constraints for open
curves (ends march freely), multi-face CAD surfaces, the volume collar block,
run time.
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.geometry import BODY, EntityKind, EntityRef, Projection

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

STEPS = 0.01 * 1.08 ** np.arange(25)
DIST = np.r_[0.0, np.cumsum(STEPS)]  # distance of each layer from the curve
TOTAL = float(STEPS.sum())


def kernel():
    """The function under test (imported late: absent until implemented)."""
    from shore.surface.march import march_surface

    return march_surface


# ---------------------------------------------------------------------------
# Exact analytic surfaces (frozen with the spec)
# ---------------------------------------------------------------------------


class AnalyticBody:
    """An exact surface behind the ``Geometry`` protocol (BODY only).

    Counts the points it is asked to project (``work``).
    """

    def __init__(self) -> None:
        self.work = 0

    @property
    def bounds(self) -> np.ndarray:
        return np.array([[-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]])

    def entities(self, kind: EntityKind) -> tuple[EntityRef, ...]:
        return (BODY,) if kind is EntityKind.BODY else ()

    def ray_cast(self, *args: object, **kwargs: object) -> object:
        raise NotImplementedError("not used by surface marching")

    def evaluate(self, *args: object, **kwargs: object) -> object:
        raise NotImplementedError("not used by surface marching")

    def fingerprint(self, *args: object) -> object:
        raise NotImplementedError("not used by surface marching")

    def project(
        self, points: np.ndarray, entity: EntityRef | None = None, *, exact: bool = False
    ) -> Projection:
        pts = np.asarray(points, dtype=np.float64)
        self.work += pts[..., 0].size
        foot, normal = self.foot(pts)
        return Projection(foot, np.linalg.norm(pts - foot, axis=-1), normal, None, None)

    def distance(self, pts: np.ndarray) -> np.ndarray:  # exact, not counted
        return np.linalg.norm(pts - self.foot(pts)[0], axis=-1)

    def normal(self, pts: np.ndarray) -> np.ndarray:  # exact, not counted
        return self.foot(pts)[1]


class Plane(AnalyticBody):
    """z = 0, normal +z."""

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        f = p.copy()
        f[..., 2] = 0.0
        return f, np.broadcast_to([0.0, 0.0, 1.0], p.shape).copy()


class Cylinder(AnalyticBody):
    """Unit cylinder about z, outward normal."""

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        r = np.linalg.norm(p[..., :2], axis=-1, keepdims=True)
        n = np.concatenate([p[..., :2] / r, np.zeros_like(r)], axis=-1)
        return n + p * [0.0, 0.0, 1.0], n


class Sphere(AnalyticBody):
    """Unit sphere, outward normal."""

    def foot(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        n = p / np.linalg.norm(p, axis=-1, keepdims=True)
        return n, n


# ---------------------------------------------------------------------------
# Curves, references and checks (frozen with the spec)
# ---------------------------------------------------------------------------


def wrap(a: np.ndarray) -> np.ndarray:
    return (a + np.pi) % (2.0 * np.pi) - np.pi


def ring(n: int = 48) -> np.ndarray:
    phi = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return np.c_[np.cos(phi), np.sin(phi), np.zeros(n)]


def l_curve(turn: int) -> np.ndarray:
    """41 nodes, legs of length 1 meeting at the origin; turn +1 left, -1 right."""
    a = np.linspace(-1.0, 0.0, 21)[:-1]
    b = np.linspace(0.0, 1.0, 21)
    return np.vstack([np.c_[a, 0 * a, 0 * a], np.c_[0 * b, turn * b, 0 * b]])


def helix() -> np.ndarray:
    t = np.linspace(0.0, 1.2, 41)
    return np.c_[np.cos(t), np.sin(t), 0.6 * t]


def generator() -> np.ndarray:
    z = np.linspace(0.0, 1.0, 41)
    return np.c_[np.ones_like(z), 0 * z, z]


def line() -> np.ndarray:
    x = np.linspace(0.0, 1.0, 41)
    return np.c_[x, 0 * x, 0 * x]


def pure_normal(surface: AnalyticBody, curve: np.ndarray) -> np.ndarray:
    """Reference: each node stepped along n x t (one-sided t at the ends), then projected."""
    out = [curve]
    for s in STEPS:
        r = out[-1]
        t = np.empty_like(r)
        t[1:-1] = 0.5 * (r[2:] - r[:-2])
        t[0], t[-1] = r[1] - r[0], r[-1] - r[-2]
        b = np.cross(surface.normal(r), t)
        b /= np.linalg.norm(b, axis=-1, keepdims=True)
        out.append(surface.foot(r + s * b)[0])
    return np.stack(out, axis=1)


def corner_sines(grid: np.ndarray, normals: np.ndarray) -> np.ndarray:
    """Signed sine of each cell's smallest corner angle, oriented by the normal; <= 0 folded."""
    p = [grid[:-1, :-1], grid[1:, :-1], grid[1:, 1:], grid[:-1, 1:]]
    n = [normals[:-1, :-1], normals[1:, :-1], normals[1:, 1:], normals[:-1, 1:]]
    out = []
    for c in range(4):
        a = p[(c + 1) % 4] - p[c]
        b = p[(c - 1) % 4] - p[c]
        s = np.einsum("...j,...j->...", np.cross(a, b), n[c])
        out.append(s / (np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1)))
    return np.minimum.reduce(out)


def orthogonality_error(grid: np.ndarray) -> np.ndarray:
    """|90 - angle| (deg) between each marching segment and its layer's tangent, interior."""
    seg = grid[1:-1, 1:] - grid[1:-1, :-1]
    tan = 0.5 * (grid[2:, :-1] - grid[:-2, :-1])
    c = np.einsum("...j,...j->...", seg, tan)
    c = c / (np.linalg.norm(seg, axis=-1) * np.linalg.norm(tan, axis=-1))
    return np.abs(90.0 - np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def min_angle_deg(grid: np.ndarray, surface: AnalyticBody) -> float:
    s = corner_sines(grid, surface.normal(grid))
    return float(np.degrees(np.arcsin(np.clip(s.min(), -1.0, 1.0))))


CASES = {
    "line_plane": (Plane, line, False),
    "ring_cylinder": (Cylinder, ring, True),
    "generator_cylinder": (Cylinder, generator, False),
    "equator_sphere": (Sphere, ring, True),
    "helix_cylinder": (Cylinder, helix, False),
    "concave_l_plane": (Plane, lambda: l_curve(+1), False),
    "convex_l_plane": (Plane, lambda: l_curve(-1), False),
}


def run(name: str) -> tuple[AnalyticBody, np.ndarray, np.ndarray]:
    make_surface, make_curve, periodic = CASES[name]
    surface, curve = make_surface(), make_curve()
    return surface, curve, kernel()(surface, curve, STEPS, periodic=periodic)


# ---------------------------------------------------------------------------
# M1-M4: exact references
# ---------------------------------------------------------------------------


class TestExact:
    @PENDING
    def test_m1_line_on_plane(self):
        _, curve, g = run("line_plane")
        np.testing.assert_allclose(
            g[..., 0], np.broadcast_to(curve[:, :1], g.shape[:2]), atol=1e-12
        )
        np.testing.assert_allclose(g[..., 1], np.broadcast_to(DIST, g.shape[:2]), atol=1e-12)

    @PENDING
    def test_m2_ring_on_cylinder(self):
        _, curve, g = run("ring_cylinder")
        phi0 = np.arctan2(curve[:, 1], curve[:, 0])[:, None]
        assert np.abs(wrap(np.arctan2(g[..., 1], g[..., 0]) - phi0)).max() <= 1e-12
        np.testing.assert_allclose(g[..., 2], np.broadcast_to(DIST, g.shape[:2]), atol=1e-12)

    @PENDING
    def test_m3_generator_on_cylinder(self):
        _, curve, g = run("generator_cylinder")
        geodesic = np.abs(np.arctan2(g[..., 1], g[..., 0]))  # from the generator at phi = 0
        assert np.abs(geodesic - DIST)[1:-1].max() <= 0.005 * TOTAL
        np.testing.assert_allclose(
            g[..., 2], np.broadcast_to(curve[:, 2:], g.shape[:2]), atol=1e-12
        )

    @PENDING
    def test_m4_equator_on_sphere(self):
        _, curve, g = run("equator_sphere")
        lon0 = np.arctan2(curve[:, 1], curve[:, 0])[:, None]
        assert np.abs(wrap(np.arctan2(g[..., 1], g[..., 0]) - lon0)).max() <= 1e-12
        lat = np.arcsin(np.clip(g[..., 2], -1.0, 1.0))
        assert np.abs(lat - DIST).max() <= 0.05 * TOTAL


# ---------------------------------------------------------------------------
# M5-M7: orthogonality and robustness
# ---------------------------------------------------------------------------


class TestQuality:
    @PENDING
    def test_m5_helix_orthogonal(self):
        surface, _, g = run("helix_cylinder")
        err = orthogonality_error(g)
        assert err.mean() <= 0.5 and err.max() <= 1.0
        assert (corner_sines(g, surface.normal(g)) <= 0).sum() == 0

    @PENDING
    def test_m6_concave_corner(self):
        surface, curve, g = run("concave_l_plane")
        ref = pure_normal(surface, curve)
        assert (corner_sines(ref, surface.normal(ref)) <= 0).sum() >= 1  # the case discriminates
        assert (corner_sines(g, surface.normal(g)) <= 0).sum() == 0
        assert min_angle_deg(g, surface) >= 5.0

    @PENDING
    def test_m7_convex_corner(self):
        surface, _, g = run("convex_l_plane")
        assert (corner_sines(g, surface.normal(g)) <= 0).sum() == 0
        assert min_angle_deg(g, surface) >= 30.0


# ---------------------------------------------------------------------------
# M8: invariants on every case
# ---------------------------------------------------------------------------


class TestInvariants:
    @PENDING
    @pytest.mark.parametrize("name", list(CASES))
    def test_m8_shape_curve_surface_orientation_work(self, name):
        make_surface, make_curve, periodic = CASES[name]
        surface, curve = make_surface(), make_curve()
        g = kernel()(surface, curve, STEPS, periodic=periodic)
        assert g.shape == (len(curve), len(STEPS) + 1, 3) and g.dtype == np.float64
        np.testing.assert_array_equal(g[:, 0], curve)
        assert surface.distance(g).max() <= 1e-12 * float(np.linalg.norm(np.ptp(surface.bounds, 0)))
        # orientation: the first marching step goes to the n x t side
        t = np.gradient(curve, axis=0)
        side = np.einsum("ij,ij->i", g[:, 1] - g[:, 0], np.cross(surface.normal(curve), t))
        assert (side > 0).all()
        assert surface.work <= 4 * g.shape[0] * g.shape[1]


# ---------------------------------------------------------------------------
# M9: interface
# ---------------------------------------------------------------------------


class TestInterface:
    @PENDING
    def test_m9_defaults_and_determinism(self):
        a = kernel()(Cylinder(), helix(), STEPS)
        b = kernel()(Cylinder(), helix(), STEPS, periodic=False, entity=BODY)
        np.testing.assert_array_equal(a, b)

    @PENDING
    @pytest.mark.parametrize(
        "steps", [np.array([]), np.array([0.01, -0.01]), np.zeros((2, 2)) + 0.01, [0.01, np.nan]]
    )
    def test_m9_bad_steps(self, steps):
        with pytest.raises(ParameterError):
            kernel()(Plane(), line(), steps)

    @PENDING
    @pytest.mark.parametrize("curve", [np.zeros((2, 3)), np.zeros((5, 2)), np.zeros(5)])
    def test_m9_bad_curve(self, curve):
        with pytest.raises(ParameterError):
            kernel()(Plane(), curve, STEPS)

    @PENDING
    def test_m9_curve_off_surface(self):
        curve = line()
        curve[7, 2] = 1e-3
        with pytest.raises(ParameterError):
            kernel()(Plane(), curve, STEPS)
