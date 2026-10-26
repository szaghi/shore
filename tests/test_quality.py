"""Hand-computed checks of shore.quality (each metric on cells with known answers)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.quality import (
    aspect_ratio,
    k_angle_deg,
    k_growth,
    k_orthogonality_deg,
    k_spacing,
    scaled_jacobian,
    summary,
)


def _grid(xs, ys, zs) -> np.ndarray:
    x, y, z = np.meshgrid(xs, ys, zs, indexing="ij")
    return np.stack([x, y, z], axis=-1).astype(float)


def _rotation(seed: int = 0) -> np.ndarray:
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q *= np.sign(np.diag(r))
    return q * np.sign(np.linalg.det(q))  # proper rotation (det = +1)


class TestScaledJacobian:
    def test_unit_cube(self):
        np.testing.assert_allclose(scaled_jacobian(_grid([0, 1], [0, 1], [0, 1])), 1.0)

    def test_rectangular_box_is_one(self):
        nodes = _grid([0, 2], [0, 1], [0, 0.5])
        np.testing.assert_allclose(scaled_jacobian(nodes), 1.0)

    @pytest.mark.parametrize("theta_deg", [10.0, 30.0, 60.0])
    def test_shear_gives_cosine(self, theta_deg):
        # k-edges tilted by theta toward +i: every corner's scaled J is cos(theta)
        nodes = _grid([0, 1], [0, 1], [0, 1])
        nodes[:, :, 1, 0] += np.tan(np.radians(theta_deg))
        np.testing.assert_allclose(scaled_jacobian(nodes), np.cos(np.radians(theta_deg)))

    def test_inverted_cell_is_negative(self):
        nodes = _grid([0, 1], [0, 1], [0, -1])  # k points the wrong way
        np.testing.assert_allclose(scaled_jacobian(nodes), -1.0)

    def test_degenerate_corner_flagged(self):
        nodes = _grid([0, 1], [0, 1], [0, 1])
        nodes[0, 0, 1] = nodes[0, 0, 0]  # zero-length k-edge at one corner
        assert scaled_jacobian(nodes).min() == pytest.approx(-1.0)  # degenerate flagged

    def test_rotation_invariant(self):
        nodes = _grid([0, 1, 2], [0, 0.5, 1.5], [0, 0.1, 0.3])
        nodes[1, 1, 1] += [0.05, -0.02, 0.01]
        rotated = nodes @ _rotation().T + [3.0, -1.0, 7.0]
        np.testing.assert_allclose(scaled_jacobian(rotated), scaled_jacobian(nodes), atol=1e-12)

    def test_shape(self):
        assert scaled_jacobian(_grid(range(4), range(5), range(6))).shape == (3, 4, 5)


class TestOrthogonality:
    def test_orthogonal_stack(self):
        np.testing.assert_allclose(k_orthogonality_deg(_grid(range(3), range(3), [0, 1])), 0.0)

    @pytest.mark.parametrize("theta_deg", [5.0, 45.0])
    def test_tilt_angle(self, theta_deg):
        nodes = _grid(range(3), range(3), [0, 1])
        nodes[:, :, 1, 0] += np.tan(np.radians(theta_deg))
        np.testing.assert_allclose(k_orthogonality_deg(nodes), theta_deg, atol=1e-9)

    def test_backward_segment_over_90(self):
        nodes = _grid(range(3), range(3), [0, -1])
        np.testing.assert_allclose(k_orthogonality_deg(nodes), 180.0)

    def test_other_layer(self):
        nodes = _grid(range(3), range(3), [0, 1, 2])
        nodes[:, :, 2, 1] += 1.0  # tilt 45 deg between layers 1 and 2 only
        np.testing.assert_allclose(k_orthogonality_deg(nodes, k=0), 0.0)
        np.testing.assert_allclose(k_orthogonality_deg(nodes, k=1), 45.0, atol=1e-9)


class TestAngleToDirections:
    def test_against_given_normals(self):
        nodes = _grid(range(3), range(3), [0, 1])
        up = np.broadcast_to([0.0, 0.0, 1.0], (3, 3, 3))
        np.testing.assert_allclose(k_angle_deg(nodes, up), 0.0)
        tilted = np.broadcast_to([0.0, np.sin(0.3), np.cos(0.3)], (3, 3, 3))
        np.testing.assert_allclose(k_angle_deg(nodes, 7.0 * tilted), np.degrees(0.3))

    def test_edge_nodes_unbiased(self):
        # a curved layer: k_orthogonality_deg is one-sided on edges, k_angle_deg is exact
        x, y = np.meshgrid(np.linspace(-0.5, 0.5, 5), np.linspace(-0.5, 0.5, 5), indexing="ij")
        r = 1.0
        base = np.stack([x, y, np.sqrt(r**2 - x**2 - y**2)], axis=-1)
        normal = base / np.linalg.norm(base, axis=-1, keepdims=True)
        nodes = np.stack([base, base + 0.01 * normal], axis=2)
        np.testing.assert_allclose(k_angle_deg(nodes, normal), 0.0, atol=1e-6)
        assert k_orthogonality_deg(nodes)[0, 0] > 1.0  # the biased one-sided value


class TestSpacing:
    def test_geometric_stack(self):
        z = np.concatenate([[0.0], np.cumsum(0.01 * 1.2 ** np.arange(5))])
        nodes = _grid([0, 1], [0, 1], z)
        np.testing.assert_allclose(k_spacing(nodes)[0, 0], 0.01 * 1.2 ** np.arange(5))
        np.testing.assert_allclose(k_growth(nodes), 1.2)

    def test_aspect_ratio(self):
        np.testing.assert_allclose(aspect_ratio(_grid([0, 2], [0, 1], [0, 0.5])), 4.0)
        np.testing.assert_allclose(aspect_ratio(_grid([0, 1], [0, 1], [0, 1])), 1.0)


class TestSummary:
    def test_keys_and_values(self):
        nodes = _grid([0, 2], [0, 1], np.concatenate([[0.0], np.cumsum([0.1, 0.12, 0.144])]))
        s = summary([nodes, nodes])
        assert s["scaled_jacobian_min"] == pytest.approx(1.0)
        assert s["folded_cells"] == 0
        assert s["wall_orthogonality_max_deg"] == pytest.approx(0.0, abs=1e-9)
        assert s["k_growth_max"] == pytest.approx(1.2)
        assert s["aspect_ratio_max"] == pytest.approx(2.0 / 0.1)
