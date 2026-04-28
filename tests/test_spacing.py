"""Tests for the 1D spacing laws in shore.volume.spacing."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.volume.spacing import (
    build_i_coords,
    build_j_coords,
    build_k_steps,
    tanh_one_sided,
    tanh_two_sided,
    uniform,
)

# ---------------------------------------------------------------------------
# Low-level primitives
# ---------------------------------------------------------------------------


class TestUniform:
    def test_shape(self):
        assert uniform(10, 0.0, 1.0).shape == (10,)

    def test_endpoints(self):
        c = uniform(10, 0.5, 2.5)
        assert c[0] == pytest.approx(0.5)
        assert c[-1] == pytest.approx(2.5)

    def test_equal_spacing(self):
        c = uniform(5, 0.0, 4.0)
        diffs = np.diff(c)
        assert np.allclose(diffs, diffs[0])


class TestTanhOneSided:
    def test_shape(self):
        assert tanh_one_sided(10, 0.0, 1.0, beta=3.0).shape == (10,)

    def test_endpoints(self):
        c = tanh_one_sided(10, 0.1, 0.9, beta=3.0)
        assert c[0] == pytest.approx(0.1)
        assert c[-1] == pytest.approx(0.9)

    def test_clusters_at_start(self):
        c = tanh_one_sided(20, 0.0, 1.0, beta=5.0)
        diffs = np.diff(c)
        # First interval must be smaller than last
        assert diffs[0] < diffs[-1]

    def test_monotone(self):
        c = tanh_one_sided(20, 0.0, 1.0, beta=3.0)
        assert np.all(np.diff(c) > 0)


class TestTanhTwoSided:
    def test_shape(self):
        assert tanh_two_sided(10, 0.0, 1.0, beta=3.0).shape == (10,)

    def test_endpoints(self):
        c = tanh_two_sided(10, 0.2, 0.8, beta=3.0)
        assert c[0] == pytest.approx(0.2)
        assert c[-1] == pytest.approx(0.8)

    def test_symmetric_clusters_at_both_ends(self):
        c = tanh_two_sided(20, 0.0, 1.0, beta=5.0)
        diffs = np.diff(c)
        mid = len(diffs) // 2
        # Interior spacing must be larger than boundary spacing
        assert diffs[mid] > diffs[0]
        assert diffs[mid] > diffs[-1]

    def test_monotone(self):
        c = tanh_two_sided(20, 0.0, 1.0, beta=3.0)
        assert np.all(np.diff(c) > 0)

    def test_cluster_at_shifts_density(self):
        # cluster toward 0.2 (near start): spacing before 0.2 should be finer
        c = tanh_two_sided(40, 0.0, 1.0, beta=4.0, cluster_at=0.2)
        diffs = np.diff(c)
        # index corresponding to 0.2 (approx)
        idx = int(np.searchsorted(c, 0.2))
        assert diffs[max(0, idx - 2)] < diffs[idx + 5]


# ---------------------------------------------------------------------------
# Domain builders
# ---------------------------------------------------------------------------


class TestBuildICoords:
    def test_uniform_default(self):
        c = build_i_coords(10, 0.1, 3.0)
        assert np.allclose(c, np.linspace(0.1, 3.0, 10))

    def test_tanh_clusters_near_start(self):
        c = build_i_coords(20, 0.1, np.pi - 0.1, spacing="tanh", beta=4.0)
        diffs = np.diff(c)
        assert diffs[0] < diffs[-1]

    def test_tanh2_clusters_at_both_ends(self):
        c = build_i_coords(20, 0.1, np.pi - 0.1, spacing="tanh2", beta=4.0)
        diffs = np.diff(c)
        assert diffs[len(diffs) // 2] > diffs[0]

    def test_unknown_spacing_raises(self):
        with pytest.raises(ParameterError, match="Unknown i-spacing"):
            build_i_coords(10, 0.1, 1.0, spacing="bad")


class TestBuildJCoords:
    def test_uniform_default(self):
        c = build_j_coords(12)
        expected = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
        assert np.allclose(c, expected)

    def test_tanh2_shape_and_range(self):
        c = build_j_coords(16, spacing="tanh2", beta=3.0)
        assert c.shape == (16,)
        assert c[0] >= 0.0
        assert c[-1] < 2.0 * np.pi

    def test_tanh2_monotone(self):
        c = build_j_coords(16, spacing="tanh2", beta=3.0)
        assert np.all(np.diff(c) > 0)

    def test_tanh2_cluster_lon(self):
        target = np.pi / 2
        c = build_j_coords(32, spacing="tanh2", beta=5.0, cluster_lon=target)
        diffs = np.diff(c)
        idx = int(np.searchsorted(c, target))
        # Spacing near target should be finer than far from it
        assert diffs[max(0, idx - 1)] < diffs[0]

    def test_unknown_spacing_raises(self):
        with pytest.raises(ParameterError, match="Unknown j-spacing"):
            build_j_coords(12, spacing="tanh")


class TestBuildKSteps:
    def test_geometric_default(self):
        steps = build_k_steps(6, ds=0.01, growth=1.1)
        assert steps.shape == (5,)
        assert steps[0] == pytest.approx(0.01)
        assert steps[1] == pytest.approx(0.01 * 1.1)

    def test_tanh_shape_and_positive(self):
        steps = build_k_steps(10, ds=0.01, spacing="tanh", beta=3.0, total_thickness=1.0)
        assert steps.shape == (9,)
        assert np.all(steps > 0)

    def test_tanh_total_thickness(self):
        steps = build_k_steps(10, ds=0.01, spacing="tanh", beta=3.0, total_thickness=1.0)
        assert steps.sum() == pytest.approx(1.0, rel=1e-10)

    def test_tanh2_total_thickness(self):
        steps = build_k_steps(10, ds=0.01, spacing="tanh2", beta=3.0, total_thickness=2.0)
        assert steps.sum() == pytest.approx(2.0, rel=1e-10)

    def test_tanh_requires_thickness(self):
        with pytest.raises(ParameterError, match="total_thickness"):
            build_k_steps(10, ds=0.01, spacing="tanh")

    def test_unknown_spacing_raises(self):
        with pytest.raises(ParameterError, match="Unknown k-spacing"):
            build_k_steps(10, ds=0.01, spacing="linear", total_thickness=1.0)
