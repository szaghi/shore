"""Unit checks of the hyperbolic marching kernel (exact cases, not the frozen spec)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh import CubedSphereTopology, Mesh
from shore.quality import scaled_jacobian
from shore.volume.hyperbolic_scheme import hyperbolic_layer


def _plane(ni: int = 7, nj: int = 6, h: float = 0.1) -> np.ndarray:
    x, y = np.meshgrid(np.arange(ni) * h, np.arange(nj) * h, indexing="ij")
    return np.stack([x, y, np.zeros_like(x)], axis=-1)


class TestExactCases:
    @pytest.mark.parametrize("step", [1e-6, 1e-3, 0.05, 0.3])
    def test_uniform_plane_is_pure_normal_extrusion(self, step):
        # uniform flat grid: coupling and dissipation terms vanish identically
        layer = _plane()
        new = hyperbolic_layer(layer, step, (False, False))
        np.testing.assert_allclose(new[..., :2], layer[..., :2], atol=1e-14)
        np.testing.assert_allclose(new[..., 2], step, rtol=1e-12)

    def test_periodic_cylinder_ring_moves_radially(self):
        # an i-periodic band on a cylinder, (theta, z) so that t_i x t_j = theta x z
        # is the outward radial direction: every node moves outward by ~step
        theta = np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False)
        z = np.linspace(0.0, 0.5, 6)
        t, zz = np.meshgrid(theta, z, indexing="ij")
        layer = np.stack([np.cos(t), np.sin(t), zz], axis=-1)
        new = hyperbolic_layer(layer, 1e-3, (True, False))
        dr = np.linalg.norm(new[..., :2], axis=-1) - 1.0
        np.testing.assert_allclose(dr, 1e-3, rtol=2e-2)
        np.testing.assert_allclose(new[..., 2], layer[..., 2], atol=1e-6)

    def test_boundary_target_is_exact(self):
        layer = _plane()
        target = layer.copy()
        target[..., 2] = 0.02
        target[0, :, 0] -= 0.003  # pinned boundary row moves sideways too
        new = hyperbolic_layer(layer, 0.02, (False, False), boundary_target=target)
        for edge in (np.s_[0], np.s_[-1], np.s_[:, 0], np.s_[:, -1]):
            np.testing.assert_array_equal(new[edge], target[edge])

    def test_positive_volume_and_deterministic(self):
        rng = np.random.default_rng(0)
        layer = _plane()
        layer[1:-1, 1:-1, :2] += 0.01 * rng.standard_normal(layer[1:-1, 1:-1, :2].shape)
        a = hyperbolic_layer(layer, 0.02, (False, False))
        b = hyperbolic_layer(layer, 0.02, (False, False))
        np.testing.assert_array_equal(a, b)
        assert scaled_jacobian(np.stack([layer, a], axis=2)).min() > 0.0


class TestMeshInterface:
    def test_default_is_algebraic_and_unknown_rejected(self, sphere_stl):
        mesh = Mesh.from_file(sphere_stl, CubedSphereTopology(), ni=10, nj=12, nk=4)
        with pytest.raises(ParameterError, match="march_scheme"):
            mesh.march(ds=0.05, march_scheme="elliptic")

    def test_hyperbolic_marches_sphere(self, sphere_stl):
        mesh = Mesh.from_file(sphere_stl, CubedSphereTopology(), ni=10, nj=12, nk=6)
        mesh.march(ds=1e-3, growth=1.2, march_scheme="hyperbolic")
        for block in mesh.blocks:
            assert scaled_jacobian(block.nodes).min() > 0.0
