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


class TestCrease:
    """The marching frame at a corner of the front on a crease (ADR 0023)."""

    def test_square_prism_marches_nearly_as_a_translation(self):
        # a square around a z-extruded square prism: every corner sits on a crease, where
        # the node normal (one face's) leans 45 deg along the front tangent
        side = np.linspace(-0.5, 0.5, 9)[:-1]
        front = np.vstack(
            [np.c_[side, -0.5 + 0 * side], np.c_[0.5 + 0 * side, side],
             np.c_[-side, 0.5 + 0 * side], np.c_[-0.5 + 0 * side, -side]]
        )  # fmt: skip
        curve = np.c_[front, np.zeros(len(front))]
        prism = StlGeometry(trimesh.creation.box(extents=(1.0, 1.0, 4.0)))
        grid = march_surface(prism, curve, np.full(10, 0.05), periodic=True)
        # nearly a translation (before the fix: 0.071 drift, top layer 0.507-0.583)
        assert np.abs(grid[:, :, :2] - front[:, None]).max() <= 0.02
        assert 0.48 <= grid[:, -1, 2].min() and grid[:, -1, 2].max() <= 0.53


class TestSharpCrease:
    """A front folded back at a sharp edge (ADR 0024)."""

    @staticmethod
    def folded_front():
        a, h = np.radians(8.0), 0.001  # a 16 deg wedge, the front around its edge
        k = np.arange(1, 8)[::-1]
        upper = np.c_[-k * h * np.cos(a), k * h * np.sin(a), 0.0 * k]
        lower = np.c_[-k[::-1] * h * np.cos(a), -k[::-1] * h * np.sin(a), 0.0 * k]
        r = np.vstack([upper, [[0.0, 0.0, 0.0]], lower])
        n_up, n_lo = np.array([np.sin(a), np.cos(a), 0.0]), np.array([np.sin(a), -np.cos(a), 0.0])
        normals = np.vstack([np.tile(n_up, (8, 1)), np.tile(n_lo, (7, 1))])
        return r, normals

    def test_edge_node_marches_the_step(self):
        from shore.surface.march import MarchParams, _layer

        r, normals = self.folded_front()
        fixed = np.linalg.norm(_layer(r, normals, 0.01, False, MarchParams()) - r, axis=1)
        assert np.abs(fixed - 0.01).max() <= 1e-4
        # without the fix the edge node's area turns into a longer step
        old = np.linalg.norm(
            _layer(r, normals, 0.01, False, MarchParams(sharp_crease=0.0)) - r, axis=1
        )
        assert old.max() >= 0.0103

    def test_arc_spacing(self):
        from shore.surface.march import _arc_spacing

        r = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 2.0, 0.0], [0.0, 2.0, 0.0]])
        np.testing.assert_allclose(_arc_spacing(r, False, None), [1.0, 1.5, 1.5, 1.0])
        np.testing.assert_allclose(_arc_spacing(r, True, None), [1.5, 1.5, 1.5, 1.5])


class TestEnds:
    """Two end planes (ADR 0024): the same plane twice is the plane mode."""

    def test_same_plane_twice_is_bit_identical(self):
        from shore.surface.march import _Ends, _march_between, _Mirror

        cylinder = StlGeometry(trimesh.creation.cylinder(radius=1.0, height=4.0, sections=64))
        theta = np.linspace(0.0, np.pi, 17)
        arc = np.c_[np.cos(theta), np.sin(theta), 0.0 * theta]
        arc = np.asarray(cylinder.project(arc).points)
        plane = ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))
        steps = np.full(6, 0.05)
        one = march_surface(cylinder, arc, steps, plane=plane)
        mirror = _Mirror.of(plane)
        two = _march_between(cylinder, arc, steps, _Ends(mirror, mirror))
        assert np.array_equal(one, two)

    def test_each_end_on_its_plane(self):
        from shore.surface.march import _Ends, _march_between, _Mirror

        cylinder = StlGeometry(trimesh.creation.cylinder(radius=1.0, height=4.0, sections=64))
        theta = np.linspace(0.0, 0.5 * np.pi, 17)  # a quarter: ends on y = 0 and x = 0
        arc = np.asarray(cylinder.project(np.c_[np.cos(theta), np.sin(theta), 0.0 * theta]).points)
        first, last = _Mirror.of(((0, 0, 0), (0, 1, 0))), _Mirror.of(((0, 0, 0), (1, 0, 0)))
        grid = _march_between(cylinder, arc, np.full(6, 0.05), _Ends(first, last))
        assert np.abs(grid[0, :, 1]).max() <= 1e-9
        assert np.abs(grid[-1, :, 0]).max() <= 1e-9


class TestClusteredCorners:
    """Sub-stepping and untangling (ADR 0027)."""

    def test_substeps_ignore_converging_nodes(self):
        from shore.surface.march import DEFAULT, _substeps

        x = np.r_[0.0, 0.1, 0.2, 0.3, 0.301, 0.302, 0.303, 0.4, 0.5]
        clustered = np.c_[x, 0.0 * x, 0.0 * x]  # a run of 1e-3 cells, as at a clustered corner
        uniform = np.c_[np.linspace(0.0, 0.5, 9), 0.0 * x, 0.0 * x]
        assert _substeps(uniform, uniform, 0.05, False, None, DEFAULT) == 1
        # clustered from the start: 0.1 / 1e-3 = 100 > 24 -> ceil(100 / 2) sub-steps
        assert _substeps(clustered, clustered, 0.1, False, None, DEFAULT) == 50
        # the same cells after a front converged there from uniform: left alone
        assert _substeps(clustered, uniform, 0.1, False, None, DEFAULT) == 1

    def test_inverted_finds_a_reversed_segment(self):
        from shore.surface.march import _inverted

        old = np.c_[np.linspace(0.0, 1.0, 6), np.zeros(6), np.zeros(6)]
        new = old + np.array([0.0, 0.1, 0.0])
        normals = np.tile([0.0, 0.0, 1.0], (6, 1))
        assert not _inverted(old, new, normals, False).any()
        new[[2, 3]] = new[[3, 2]]  # two nodes overtake each other
        assert _inverted(old, new, normals, False)[2]
