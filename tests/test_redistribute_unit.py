"""Unit tests of the edge redistribution helpers (ADR 0041)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh.redistribute import edge_axis, fractions, resample
from shore.mesh.spacing import Spacing1D


def test_edge_axis():
    assert edge_axis("j_lo_k_lo") == "i"
    assert edge_axis("i_hi_k_lo") == "j"
    assert edge_axis("i_lo_j_hi") == "k"


def test_fractions_span_zero_to_one():
    nodes = np.stack([np.linspace(0, 3, 7), np.zeros(7), np.zeros(7)], axis=1)
    t = fractions(Spacing1D("tanh", beta=2.0), nodes)
    assert t[0] == 0.0 and t[-1] == 1.0 and np.all(np.diff(t) > 0)
    with pytest.raises(ParameterError):
        fractions(Spacing1D("geometric"), nodes)
    with pytest.raises(ParameterError):
        fractions(Spacing1D("uniform"), np.zeros((7, 3)))


def test_resample_keeps_ends_and_collapsed_lines():
    t = np.array([0.0, 0.1, 0.5, 1.0])
    line = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 1, 0]])
    point = np.zeros((4, 3))  # a collapsed (polar) line
    out = resample(np.stack([line, point]), t)
    assert np.array_equal(out[0, 0], line[0]) and np.array_equal(out[0, -1], line[-1])
    np.testing.assert_allclose(out[0, 1], [0.3, 0, 0])  # 10 % of the arc length 3
    np.testing.assert_allclose(out[0, 2], [1, 0.5, 0])  # half-way
    assert np.array_equal(out[1], point)
