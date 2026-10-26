"""Unit checks of the intersection tracer outside the frozen spec (tests/spec/)."""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.geometry import StlGeometry
from shore.surface.intersect import intersection_loop
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
