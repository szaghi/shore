"""Unit checks of the surface TFI kernel internals (the frozen spec is tests/spec/)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.geometry import EntityKind, EntityRef, Projection
from shore.surface.tfi import EllipticParams, _control, _Metric, _tfi, surface_grid

FACE = EntityRef(EntityKind.SURFACE, 0)


class Sphere:
    """Unit sphere, (u, v) = (longitude, latitude); counts evaluated points."""

    bounds = np.array([[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]])

    def __init__(self) -> None:
        self.work = 0

    def evaluate(self, entity: EntityRef, params: np.ndarray) -> np.ndarray:
        self.work += params[..., 0].size
        u, v = params[..., 0], params[..., 1]
        return np.stack([np.cos(v) * np.cos(u), np.cos(v) * np.sin(u), np.sin(v)], axis=-1)

    def project(self, points, entity=None, *, exact=False) -> Projection:
        self.work += points[..., 0].size
        n = points / np.linalg.norm(points, axis=-1, keepdims=True)
        uv = np.stack([np.arctan2(n[..., 1], n[..., 0]), np.arcsin(n[..., 2])], axis=-1)
        return Projection(n, np.linalg.norm(points - n, axis=-1), n, uv)


class TestTfi:
    def test_bilinear_field_exact_and_boundary_kept(self):
        s, t = np.meshgrid(np.linspace(0, 1, 7), np.linspace(0, 1, 5), indexing="ij")
        f = np.stack([2 + 3 * s - t + 0.5 * s * t, -s + 4 * t], axis=-1)
        g = f.copy()
        g[1:-1, 1:-1] = np.nan
        np.testing.assert_allclose(_tfi(g), f, atol=1e-14)


class TestMetric:
    def test_lattice_derivatives_match_sphere(self):
        lo, hi = np.array([0.1, 0.0]), np.array([1.4, 1.1])
        metric = _Metric(Sphere(), FACE, lo, hi, (41, 33))
        uv = np.stack(np.meshgrid(np.linspace(0.2, 1.3, 9), np.linspace(0.1, 1.0, 7)), -1)
        su, sv, suu, suv, svv = metric.at(uv)
        u, v = uv[..., 0], uv[..., 1]
        cu, su_, cv, sv_ = np.cos(u), np.sin(u), np.cos(v), np.sin(v)
        exact = {
            "S_u": (su, np.stack([-cv * su_, cv * cu, 0 * u], -1)),
            "S_v": (sv, np.stack([-sv_ * cu, -sv_ * su_, cv], -1)),
            "S_uu": (suu, np.stack([-cv * cu, -cv * su_, 0 * u], -1)),
            "S_uv": (suv, np.stack([sv_ * su_, -sv_ * cu, 0 * u], -1)),
            "S_vv": (svv, np.stack([-cv * cu, -cv * su_, -sv_], -1)),
        }
        for name, (got, want) in exact.items():
            np.testing.assert_allclose(got, want, atol=2e-3, err_msg=name)


class TestControl:
    def test_tanh_row_and_clamp(self):
        x = 1.0 - np.tanh(3.0 * (1.0 - np.linspace(0, 1, 21))) / np.tanh(3.0)
        edge = np.stack([x, 0 * x, 0 * x], axis=-1)
        p = _control(edge, limit=10.0)
        d1, d2 = 0.5 * (x[2:] - x[:-2]), x[2:] + x[:-2] - 2 * x[1:-1]
        np.testing.assert_allclose(p[1:-1], -d2 / d1, rtol=1e-12)
        assert p[0] == p[-1] == 0.0
        assert np.abs(_control(edge, limit=0.05)).max() == 0.05


class TestWork:
    def test_geometry_work_independent_of_iterations(self):
        u, v = np.meshgrid(np.linspace(0.2, 1.2, 15), np.linspace(0.1, 0.9, 11), indexing="ij")
        grid = Sphere().evaluate(FACE, np.stack([u, v], axis=-1))
        work = []
        for iters in (1, 400):
            sphere = Sphere()
            surface_grid(sphere, FACE, grid, params=EllipticParams(max_iterations=iters))
            work.append(sphere.work)
        assert work[0] == work[1]


class PeriodicCylinder(Sphere):
    """Unit cylinder about z, u = atan2 in (-pi, pi]: its seam is at u = +-pi."""

    bounds = np.array([[-1.0, -1.0, 0.0], [1.0, 1.0, 1.0]])

    def evaluate(self, entity, params):
        self.work += params[..., 0].size
        return np.stack([np.cos(params[..., 0]), np.sin(params[..., 0]), params[..., 1]], -1)

    def project(self, points, entity=None, *, exact=False):
        u = np.arctan2(points[..., 1], points[..., 0])
        foot = np.stack([np.cos(u), np.sin(u), points[..., 2]], -1)
        n = foot * [1.0, 1.0, 0.0]
        return Projection(
            foot, np.linalg.norm(points - foot, axis=-1), n, np.stack([u, points[..., 2]], -1)
        )


class TestOutOfScopeFaces:
    def test_seam_crossing_rejected(self):
        u, z = np.meshgrid(np.linspace(2.8, 3.5, 9), np.linspace(0.0, 1.0, 7), indexing="ij")
        grid = np.stack([np.cos(u), np.sin(u), z], -1)  # straddles u = pi
        with pytest.raises(ParameterError, match="seam"):
            surface_grid(PeriodicCylinder(), FACE, grid)
        ok = np.stack([np.cos(u - 1.0), np.sin(u - 1.0), z], -1)  # the same patch off the seam
        surface_grid(PeriodicCylinder(), FACE, ok)

    def test_collapsed_edge_rejected(self):
        u, v = np.meshgrid(np.linspace(0.2, 1.2, 9), np.linspace(0.1, 1.0, 7), indexing="ij")
        grid = Sphere().evaluate(FACE, np.stack([u, v], -1))
        grid[:, -1] = grid[0, -1]  # the j_hi edge collapsed to one point, as at a pole
        with pytest.raises(ParameterError, match="collapsed"):
            surface_grid(Sphere(), FACE, grid)
