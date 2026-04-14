"""Tests for the jacobian module."""

from __future__ import annotations

import numpy as np
import pytest

from shore.volume.jacobian import jacobian_field, jacobian_per_layer, min_jacobian_two_layers
from tests.conftest import sphere_layer


def _box_grid(ni: int, nj: int, nk: int) -> np.ndarray:
    """Uniform Cartesian grid — all Jacobians should be identical and positive."""
    i = np.linspace(0.0, 1.0, ni)
    j = np.linspace(0.0, 1.0, nj)
    k = np.linspace(0.0, 1.0, nk)
    ig, jg, kg = np.meshgrid(i, j, k, indexing="ij")
    return np.stack([ig, jg, kg], axis=-1)


class TestJacobianField:
    def test_shape(self):
        grid = _box_grid(5, 6, 4)
        jac = jacobian_field(grid)
        assert jac.shape == (4, 6, 3)

    def test_ogrid_all_positive(self):
        # Use a sphere O-grid: j-direction is genuinely periodic.
        # A Cartesian box grid is NOT a valid O-grid topology —
        # the periodic j-wrap would create a fold.
        surface = sphere_layer(8, 10, radius=1.0)
        from shore.volume.hyperbolic import march  # noqa: PLC0415
        grid = march(surface, nk=5, ds=0.05, growth=1.1, smooth_strength=0.0)
        assert jacobian_field(grid).min() > 0.0

    def test_sphere_ogrid_all_positive(self):
        surface = sphere_layer(16, 20, radius=1.0)
        from shore.volume.hyperbolic import march
        grid = march(surface, nk=8, ds=0.05, growth=1.1, smooth_strength=0.0)
        assert jacobian_field(grid).min() > 0.0

    def test_inverted_layer_negative(self):
        """Moving the outer layer inside the inner one must give negative Jacobian."""
        grid = _box_grid(5, 6, 3)
        # Reverse the k-direction: outer shell moved to k=0 position
        grid_inv = grid[:, :, ::-1, :]
        assert jacobian_field(grid_inv).min() < 0.0

    def test_j_periodicity_accounted(self):
        """The cell at j=nj-1 must be present in the output (nj cells, not nj-1)."""
        grid = _box_grid(4, 8, 3)
        jac = jacobian_field(grid)
        assert jac.shape[1] == 8   # nj cells (periodic)


class TestJacobianPerLayer:
    def test_shape(self):
        grid = _box_grid(5, 6, 4)
        jpl = jacobian_per_layer(grid)
        assert jpl.shape == (3,)   # nk-1

    def test_values_match_field_min(self):
        grid = _box_grid(5, 6, 4)
        jac = jacobian_field(grid)
        jpl = jacobian_per_layer(grid)
        for k in range(3):
            assert np.isclose(jpl[k], jac[:, :, k].min())


class TestMinJacobianTwoLayers:
    def test_positive_for_outward_extrusion(self):
        surface = sphere_layer(10, 12, radius=1.0)
        outer = sphere_layer(10, 12, radius=1.5)
        assert min_jacobian_two_layers(surface, outer) > 0.0

    def test_negative_for_inward_extrusion(self):
        surface = sphere_layer(10, 12, radius=1.0)
        inner = sphere_layer(10, 12, radius=0.5)
        assert min_jacobian_two_layers(surface, inner) < 0.0

    def test_consistent_with_jacobian_field(self):
        """Must match jacobian_field on a 2-layer stack."""
        l0 = sphere_layer(10, 12, radius=1.0)
        l1 = sphere_layer(10, 12, radius=1.5)
        two = np.stack([l0, l1], axis=2)
        expected = float(jacobian_field(two).min())
        assert np.isclose(min_jacobian_two_layers(l0, l1), expected)
