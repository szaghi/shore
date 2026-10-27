"""Unit checks of the intersection tracer outside the frozen spec (tests/spec/)."""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.geometry import StlGeometry
from shore.surface.intersect import intersection_curve, intersection_loop
from shore.surface.march import march_surface


@pytest.fixture(scope="module")
def stl_pair() -> tuple[StlGeometry, StlGeometry]:
    """A faceted sphere pierced by a faceted vertical cylinder (two loops)."""
    sphere = StlGeometry(trimesh.creation.icosphere(subdivisions=4, radius=1.0))
    cylinder = StlGeometry(trimesh.creation.cylinder(radius=0.4, height=4.0, sections=48))
    return sphere, cylinder


@pytest.fixture(scope="module")
def stl_loop(stl_pair) -> np.ndarray:
    sphere, cylinder = stl_pair
    return intersection_loop(sphere, cylinder, np.array([0.4, 0.0, 0.9]), 48)


class TestStlBodies:
    def test_loop_on_two_triangulated_bodies(self, stl_pair, stl_loop):
        sphere, cylinder = stl_pair
        loop = stl_loop
        assert loop.shape == (48, 3)
        assert sphere.project(loop).distances.max() < 1e-8
        assert cylinder.project(loop).distances.max() < 1e-8
        assert (loop[:, 2] > 0.0).all()  # the seed's (upper) loop only
        chords = np.linalg.norm(np.roll(loop, -1, axis=0) - loop, axis=1)
        assert chords.max() / chords.min() < 1.05

    def test_collar_halves_leave_the_other_body(self, stl_pair, stl_loop):
        sphere, _ = stl_pair
        on_sphere = march_surface(sphere, stl_loop, np.full(2, 0.03), periodic=True)
        assert (np.linalg.norm(on_sphere[:, 1, :2], axis=1) > 0.4).all()  # outside the cylinder


class TestErrors:
    def test_integer_n_required(self):
        with pytest.raises(ParameterError):
            intersection_loop(None, None, np.zeros(3), 3.5)  # type: ignore[arg-type]


class TestOpenCurve:
    """intersection_curve (ADR 0011) beyond the frozen spec."""

    def test_stl_sphere_cylinder_cut_by_a_plane(self, stl_pair):
        sphere, cylinder = stl_pair
        plane = ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))
        curve = intersection_curve(sphere, cylinder, np.array([0.3, 0.25, 0.9]), 24, plane=plane)
        assert sphere.project(curve).distances.max() < 1e-8
        assert cylinder.project(curve).distances.max() < 1e-8
        assert np.abs(curve[[0, -1], 1]).max() < 1e-12
        assert curve[:, 1].min() >= 0.0 and (curve[:, 2] > 0.0).all()

    def test_tangent_plane_raises(self):
        from tests.spec.test_intersection import CylinderZ, Plane

        # the circle x^2 + y^2 = 1 touches the plane x = 1 without crossing it
        with pytest.raises(ParameterError):
            intersection_curve(
                Plane(), CylinderZ(), np.array([0.0, 1.05, 0.05]), 16,
                plane=((1.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
            )  # fmt: skip

    def test_ends_polished_to_output_tolerance(self):
        # on a half model the trace lands on the plane within the on-curve tolerance;
        # its end point must still be polished to the output tolerance (1e-13), not
        # left at the trace tolerance (1e-10: ends off the sphere by 7e-10)
        from tests.spec.test_open_curve import HALF, SYMMETRY

        make_a, make_b, seed, _, _ = HALF["sphere_cylinder"]
        a, b = make_a(), make_b()
        curve = intersection_curve(a, b, np.array(seed), 33, plane=SYMMETRY)
        assert a.distance(curve[[0, -1]]).max() < 1e-11


class TestDistribution:
    """Helpers of the node distribution along the curve (ADR 0020)."""

    def test_windowed_curvature_of_a_polygon_circle(self):
        from shore.surface.intersect import DEFAULT, _windowed_curvature

        t = np.linspace(0.0, 2.0 * np.pi, 400, endpoint=False)
        loop = 2.0 * np.c_[np.cos(t), np.sin(t), 0.0 * t]
        closed = np.vstack([loop, loop[:1]])
        s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))]
        # windows shorter and longer than a segment: smooth either way
        for n in (40, 2000):
            k = _windowed_curvature(loop, s, n, DEFAULT)
            np.testing.assert_allclose(k, 0.5, rtol=1e-4)

    @pytest.mark.parametrize("closed", [True, False])
    def test_limit_ratio(self, closed):
        from shore.surface.intersect import _limit_ratio

        h = np.array([1.0, 1.0, 0.01, 1.0, 1.0, 1.0, 0.5, 1.0])
        out = _limit_ratio(h, 1.2, closed)
        assert np.all(out <= h)
        ring = np.r_[out, out[:1]] if closed else out
        assert np.max(np.maximum(ring[1:] / ring[:-1], ring[:-1] / ring[1:])) <= 1.2 + 1e-12
        assert out[2] == pytest.approx(0.01, rel=1e-14)  # the smallest cell is kept

    def test_refine_stops_at_a_corner(self):
        from shore.surface.intersect import DEFAULT, _Tracer

        class Still:  # projection is the identity: refinement only bisects
            bounds = np.array([[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]])

            def project(self, points, entity=None):
                from shore.geometry import Projection

                p = np.asarray(points, dtype=np.float64)
                return Projection(p, np.zeros(len(p)), None, None, None)

        square = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
        tracer = _Tracer(Still(), Still(), None, None, 2.0 * np.sqrt(3.0), DEFAULT)
        fine = tracer.refine(square, closed=True)
        seg = np.linalg.norm(np.diff(np.vstack([fine, fine[:1]]), axis=0), axis=1)
        assert seg.min() >= DEFAULT.refine_min * tracer.diag / 2.0
        assert len(fine) < 4 * 2**20  # it terminates near the corners only

    def test_tanh_is_allowed_on_an_open_curve_only(self):
        from shore.mesh.spacing import Spacing1D
        from shore.surface.intersect import _check_distribution

        law = Spacing1D("tanh", beta=2.0)
        assert _check_distribution(law, 0.0, ("uniform", "tanh", "tanh2")) is False
        with pytest.raises(ParameterError):
            _check_distribution(law, 0.0, ("uniform", "tanh2"))


class TestCollarJournal:
    def test_default_loop_keywords_stay_out_of_the_journal(self):
        from shore.project.commands import BuildCollar

        cmd = BuildCollar("post.stl", (0.5, 0.0, 0.0), 16, 0.02, 1.15, 5, 5, 1e-3, 1.25, 8)
        assert "loop_" not in cmd.to_source()
        spaced = BuildCollar(
            "post.stl", (0.5, 0.0, 0.0), 16, 0.02, 1.15, 5, 5, 1e-3, 1.25, 8, loop_curvature=1.0
        )
        assert "loop_curvature=1.0" in spaced.to_source()


class TestCornerHelpers:
    """Helpers of the corner-aware tracer (ADR 0021)."""

    def test_crossing_of_two_lines(self):
        from shore.surface.intersect import _crossing

        c = _crossing(
            np.array([0.0, 0, 0]),
            np.array([1.0, 0, 0]),
            np.array([1.0, 1, 0]),
            np.array([0, 1.0, 0]),
        )
        np.testing.assert_allclose(c, [1.0, 0.0, 0.0], atol=1e-15)
        # parallel lines: the midpoint of the two points
        p = _crossing(
            np.array([0.0, 0, 0]),
            np.array([1.0, 0, 0]),
            np.array([0.0, 2, 0]),
            np.array([2.0, 0, 0]),
        )
        np.testing.assert_allclose(p, [0.0, 1.0, 0.0])

    def test_allocate_keeps_one_cell_per_segment(self):
        from shore.surface.intersect import _allocate

        cells = _allocate(np.array([10.0, 0.01, 0.01, 10.0]), 20)
        assert cells.sum() == 20 and cells.min() >= 1
        assert list(cells) == [9, 1, 1, 9]
        np.testing.assert_array_equal(_allocate(np.array([1.0, 2.0, 1.0]), 8), [2, 4, 2])

    def test_add_corner_replaces_a_point_on_it(self):
        from shore.surface.intersect import _add_corner

        points, corners = [np.zeros(3), np.array([1.0, 0, 0])], []
        _add_corner(points, corners, np.array([1.0, 0, 1e-12]), 1e-9)
        assert len(points) == 2 and corners == [1] and points[1][2] == 1e-12
        _add_corner(points, corners, np.array([2.0, 0, 0]), 1e-9)
        assert len(points) == 3 and corners == [1, 2]

    def test_refine_renumbers_fixed_vertices(self):
        from shore.surface.intersect import DEFAULT, _Tracer

        class Still:
            bounds = np.array([[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]])

            def project(self, points, entity=None):
                from shore.geometry import Projection

                p = np.asarray(points, dtype=np.float64)
                return Projection(p, np.zeros(len(p)), None, None, None)

        square = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
        tracer = _Tracer(Still(), Still(), None, None, 2.0 * np.sqrt(3.0), DEFAULT)
        fine, fixed = tracer.refine(square, closed=True, fixed=[0, 1, 2, 3])
        np.testing.assert_array_equal(fine[fixed], square)

    def test_ring_finds_both_branches(self):
        from shore.surface.intersect import DEFAULT, _Tracer
        from tests.spec.test_intersection import CylinderZ, Plane

        a, b = Plane(), CylinderZ()
        tracer = _Tracer(a, b, None, None, float(np.linalg.norm([4.0, 4, 4])), DEFAULT)
        ring = tracer._ring(np.array([1.0, 0, 0]), 0.1, ())
        assert len(ring) == 2  # the unit circle crosses the sphere twice
        for q in ring:
            assert abs(np.linalg.norm(q - [1.0, 0, 0]) - 0.1) <= 1e-9
            assert abs(np.hypot(q[0], q[1]) - 1.0) <= 1e-9 and abs(q[2]) <= 1e-12
