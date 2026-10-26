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


class TestMirror:
    """The symmetry-plane helper of the end constraint (ADR 0012)."""

    def test_ghosts_are_mirror_images(self):
        from shore.surface.march import _Mirror

        mirror = _Mirror.of(((0.0, 1.0, 0.0), (0.0, 2.0, 0.0)))  # y = 1, normal unnormalised
        f = np.array([[0.0, 1.0, 0.0], [1.0, 1.5, 2.0], [3.0, 2.0, 1.0], [4.0, 1.0, 5.0]])
        ext = mirror.extend(f)
        np.testing.assert_allclose(ext[0], [1.0, 0.5, 2.0])
        np.testing.assert_allclose(ext[-1], [3.0, 0.0, 1.0])
        np.testing.assert_allclose(mirror.reflect(mirror.reflect(f)), f)
        np.testing.assert_allclose(mirror.matrix @ [0.0, 1.0, 0.0], [0.0, -1.0, 0.0])

    def test_end_tangent_is_the_plane_normal(self):
        from shore.surface.march import _central, _Mirror

        mirror = _Mirror.of(((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)))
        theta = np.linspace(0.0, np.pi, 9)
        arc = np.c_[np.cos(theta), np.sin(theta), 0.0 * theta]  # ends on y = 0
        t = _central(arc, False, mirror)
        np.testing.assert_allclose(t[[0, -1], [0, 0]], 0.0, atol=1e-15)
