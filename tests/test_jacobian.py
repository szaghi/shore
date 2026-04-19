"""Tests for the jacobian module."""

from __future__ import annotations

import numpy as np

from shore.volume.jacobian import _jacobian_field_open, jacobian_field
from tests.conftest import sphere_layer


def _box_grid(ni: int, nj: int, nk: int) -> np.ndarray:
    """Uniform Cartesian grid — all Jacobians should be identical and positive."""
    i = np.linspace(0.0, 1.0, ni)
    j = np.linspace(0.0, 1.0, nj)
    k = np.linspace(0.0, 1.0, nk)
    ig, jg, kg = np.meshgrid(i, j, k, indexing="ij")
    return np.stack([ig, jg, kg], axis=-1)


def _sphere_shell(ni: int, nj: int, radii: list[float]) -> np.ndarray:
    """Stack sphere O-grid layers at the given radii into a (ni, nj, nk, 3) volume."""
    layers = [sphere_layer(ni, nj, radius=r) for r in radii]
    return np.stack(layers, axis=2)


class TestJacobianField:
    def test_shape(self):
        grid = _box_grid(5, 6, 4)
        jac = jacobian_field(grid)
        assert jac.shape == (4, 6, 3)

    def test_sphere_ogrid_all_positive(self):
        grid = _sphere_shell(16, 20, [1.0, 1.1, 1.25, 1.45, 1.7])
        assert jacobian_field(grid).min() > 0.0

    def test_inverted_layer_negative(self):
        grid = _box_grid(5, 6, 3)
        grid_inv = grid[:, :, ::-1, :]
        assert jacobian_field(grid_inv).min() < 0.0

    def test_j_periodicity_accounted(self):
        """The cell at j=nj-1 must be present in the output (nj cells, not nj-1)."""
        grid = _box_grid(4, 8, 3)
        jac = jacobian_field(grid)
        assert jac.shape[1] == 8   # nj cells (periodic)


class TestJacobianFieldOpen:
    def test_shape(self):
        grid = _box_grid(5, 6, 4)
        jac = _jacobian_field_open(grid)
        assert jac.shape == (4, 5, 3)   # non-periodic j → nj-1

    def test_sphere_shell_positive(self):
        grid = _sphere_shell(10, 12, [1.0, 1.2, 1.5])
        assert _jacobian_field_open(grid).min() > 0.0
