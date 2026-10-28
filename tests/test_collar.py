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


class TestExtrudeCollar:
    def test_bad_wall_and_junction(self):
        from shore.collar import extrude_collar

        wall = np.zeros((5, 4, 3))
        with pytest.raises(ParameterError):
            extrude_collar(np.zeros((5, 4)), np.array([0.1]), 2)
        with pytest.raises(ParameterError):
            extrude_collar(wall, np.array([0.1]), 0)
        with pytest.raises(ParameterError):
            extrude_collar(wall, np.array([-0.1]), 2)

    def test_end_angle(self):
        from shore.collar import _end_angle_deg

        theta = np.linspace(0.0, np.pi, 17)
        arc = np.c_[np.cos(theta), np.sin(theta), 0.0 * theta]  # meets y = 0 at right angles
        assert _end_angle_deg(arc, np.array([0.0, 1.0, 0.0])) < 0.5
        assert _end_angle_deg(arc, np.array([1.0, 0.0, 0.0])) > 80.0


class TestHalfCollarCommand:
    def test_malformed_plane(self):
        from shore.project import Project

        with pytest.raises(ParameterError, match="plane"):
            Project().build_collar(
                "other.stl", seed=(0, 0, 0), n_loop=8, ds_surface=0.1, growth_surface=1.0,
                m_a=2, m_b=2, ds_wall=0.01, growth_wall=1.0, nk=3, plane=(1.0, 2.0),
            )  # fmt: skip


class TestSplitCollarHelpers:
    """Internals of split_collar (ADR 0024)."""

    def test_wake_cells(self):
        from shore.collar import _wake_cells

        cells = _wake_cells(0.001, 0.2, 1.2)
        assert cells[0] == 0.001
        assert abs(cells.sum() - 0.2) <= 1e-12
        ratio = cells[1:] / cells[:-1]
        assert 1.0 <= ratio.min() and ratio.max() <= 1.2
        np.testing.assert_array_equal(_wake_cells(0.3, 0.2, 1.2), [0.2])

    def test_loop_corner_nodes(self):
        from shore.geometry.base import BODY
        from shore.surface.intersect import DEFAULT as TRACE
        from shore.surface.intersect import _loop
        from tests.spec.test_intersection import Plane
        from tests.spec.test_loop_corners import RECT, PolygonPrism

        args = (Plane(), PolygonPrism(RECT), np.array([0.0, 0.06, 0.0]), 64, None, 0.0)
        nodes, corners = _loop(*args, BODY, BODY, TRACE)
        assert len(corners) == 4
        corner_xy = np.sort(np.round(nodes[corners, :2], 12), axis=0)
        np.testing.assert_allclose(corner_xy, np.sort(np.asarray(RECT), axis=0), atol=1e-12)

    def test_bad_wake_length(self):
        from shore.collar import split_collar
        from tests.spec.test_intersection import Plane
        from tests.spec.test_loop_corners import LensPrism

        steps = np.full(3, 0.01)
        with pytest.raises(ParameterError, match="wake_length"):
            split_collar(
                Plane(), LensPrism(), np.array([0.0, 0.06, 0.0]), n_loop=32,
                steps_a=steps, steps_b=steps, steps_k=steps, wake_length="long",
            )  # fmt: skip


class TestSmoothedVolume:
    """Internals of the ``"smoothed"`` outer boundary (ADR 0025)."""

    @staticmethod
    def station(n_b=8, n_a=10):
        """A planar L-section: up a wall (x = 0) and along a floor (z = 0), y = 0."""
        z = np.cumsum(0.01 * 1.15 ** np.arange(n_b))[::-1]
        x = np.cumsum(0.01 * 1.15 ** np.arange(n_a))
        line = np.vstack([np.c_[0 * z, 0 * z, z], [[0.0, 0.0, 0.0]], np.c_[x, 0 * x, 0 * x]])
        corner = [[np.sqrt(0.5), 0.0, np.sqrt(0.5)]]
        normals = np.vstack(
            [np.tile([1.0, 0.0, 0.0], (n_b, 1)), corner, np.tile([0.0, 0.0, 1.0], (n_a, 1))]
        )
        return line, normals, n_b

    def test_planar_station_matches_adr0010(self):
        from shore.collar import _offset_line, _outer_line

        line, normals, junction = self.station()
        old = _outer_line(line, normals, 0.03, junction)
        new = _offset_line(line, line + 0.03 * normals, junction)
        np.testing.assert_allclose(new, old, atol=1e-12)

    def test_trim_origins(self):
        from shore.collar import _trim, _trim_origins

        line, normals, junction = self.station()
        offset = (line + 0.03 * normals)[:, [0, 2]]
        pts, _ = _trim(offset, junction)
        first, second = _trim_origins(offset)
        assert len(first) == len(pts)
        cut = np.flatnonzero(first != second)
        assert cut.size == 1  # one swallowtail, at the crease
        assert first[cut[0]] < junction < second[cut[0]]

    def test_loop_kernel(self):
        from shore.collar import _loop_kernel

        theta = np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False)
        ring = np.c_[np.cos(theta), np.sin(theta), 0.0 * theta]
        w = _loop_kernel(ring, 0.3, periodic=True)
        np.testing.assert_allclose(w.sum(axis=1), 1.0)
        np.testing.assert_allclose(w, np.roll(np.roll(w, 3, axis=0), 3, axis=1))  # rotation
        open_w = _loop_kernel(ring[:12], 0.3, periodic=False)
        assert open_w[0, -1] < open_w[-1, -1]

    def test_unknown_volume(self):
        from shore.collar import CollarParams, extrude_collar

        wall = np.zeros((4, 5, 3))
        wall[..., 0] = np.arange(4)[:, None]
        wall[..., 1] = np.arange(5)[None, :]
        with pytest.raises(ParameterError, match="volume"):
            extrude_collar(wall, np.full(3, 0.1), 1, params=CollarParams(volume="x"))
