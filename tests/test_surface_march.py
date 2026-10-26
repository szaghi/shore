"""Unit checks of the surface marcher outside the frozen spec (tests/spec/)."""

from __future__ import annotations

import numpy as np
import trimesh

from shore.geometry import StlGeometry
from shore.surface.march import march_surface


class TestStlBody:
    def test_marches_on_a_triangulated_sphere(self):
        geometry = StlGeometry(trimesh.creation.icosphere(subdivisions=4, radius=1.0))
        phi = np.linspace(0.0, 2.0 * np.pi, 40, endpoint=False)
        curve = geometry.project(np.c_[np.cos(phi), np.sin(phi), 0 * phi]).points
        steps = np.full(10, 0.03)
        grid = march_surface(geometry, curve, steps, periodic=True)
        assert grid.shape == (40, 11, 3)
        assert geometry.project(grid).distances.max() < 1e-9
        assert np.all(np.diff(grid[..., 2], axis=1) > 0.0)  # every layer moves north


class TestWork:
    def test_one_projection_per_node_and_layer(self):
        calls = []

        class Plane:
            bounds = np.array([[-1.0, -1.0, 0.0], [2.0, 2.0, 0.0]])

            def project(self, points, entity=None, *, exact=False):
                from shore.geometry import Projection

                calls.append(points.shape[0])
                foot = points * [1.0, 1.0, 0.0]
                normal = np.broadcast_to([0.0, 0.0, 1.0], points.shape).copy()
                return Projection(foot, np.abs(points[:, 2]), normal)

        x = np.linspace(0.0, 1.0, 21)
        march_surface(Plane(), np.c_[x, 0 * x, 0 * x], np.full(7, 0.05))
        assert calls == [21] * 8  # the curve once, then each new layer once
