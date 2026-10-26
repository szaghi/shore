"""Unit checks of the collar-block internals (the frozen spec is tests/spec/)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.collar import DEFAULT, _trim, _winslow
from shore.exceptions import ParameterError


class TestTrim:
    # segment (2, 1)-(1, -1) crosses segment (0, 0)-(2, 0) at (1.5, 0)
    POLY = np.array([[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, -1.0], [1.0, -2.0]])

    def test_swallowtail_cut_at_crossing(self):
        pts, _ = _trim(self.POLY, 0)
        np.testing.assert_allclose(pts, [[0, 0], [1.5, 0], [1, -1], [1, -2]], atol=1e-15)

    @pytest.mark.parametrize(("junction", "image"), [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3)])
    def test_junction_image(self, junction, image):
        # vertices 1 and 2 are cut out: both map to the crossing point
        assert _trim(self.POLY, junction)[1] == image

    def test_simple_polyline_unchanged(self):
        poly = np.c_[np.linspace(0, 1, 7), np.linspace(0, 1, 7) ** 2]
        pts, image = _trim(poly, 3)
        np.testing.assert_array_equal(pts, poly)
        assert image == 3


class TestWinslow:
    def test_tensor_product_grid_kept(self):
        # uniform in i, geometric in k: Thomas-Middlecoff makes it an exact solution
        s = np.linspace(0.0, 2.0, 9)
        t = np.r_[0.0, np.cumsum(1e-3 * 1.3 ** np.arange(12))]
        grid = np.stack(np.meshgrid(s, t, indexing="ij"), axis=-1)
        sheets = np.concatenate([grid, np.zeros((*grid.shape[:2], 1))], axis=-1)[None]
        perturbed = sheets.copy()
        perturbed[0, 1:-1, 1:-1, 1] += 0.2 * (t[-1] - t[0]) * np.sin(np.pi * s[1:-1, None] / 2)
        out = _winslow(perturbed, float(t[-1]), DEFAULT)
        np.testing.assert_allclose(out, sheets, atol=1e-9)


class TestCommandErrors:
    def test_needs_a_body(self):
        from shore.project import Project

        with pytest.raises(ParameterError, match="load_geometry"):
            Project().build_collar(
                "other.stl", seed=(0, 0, 0), n_loop=8, ds_surface=0.1, growth_surface=1.0,
                m_a=2, m_b=2, ds_wall=0.01, growth_wall=1.0, nk=3,
            )  # fmt: skip
