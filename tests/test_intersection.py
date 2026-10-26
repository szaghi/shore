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
