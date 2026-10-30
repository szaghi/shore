"""Unit tests of the wall-normal grid lines (ADR 0048)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh.orthogonal import grid_normals, hermite_lines, orthogonal_walls


def test_grid_normals_of_a_plane():
    a, b = np.meshgrid(np.linspace(0, 1, 4), np.linspace(0, 2, 5), indexing="ij")
    plane = np.stack([a, b, np.zeros_like(a)], axis=-1)
    np.testing.assert_allclose(grid_normals(plane), np.broadcast_to([0, 0, 1], (4, 5, 3)))


def test_hermite_lines_keep_ends_and_fractions():
    t = np.linspace(0.0, 1.0, 5) ** 2
    line = np.stack([t, t, np.zeros_like(t)], axis=-1)[None]  # a diagonal line, clustered
    new = hermite_lines(line, np.array([[0.0, 1.0, 0.0]]), 1.0, samples=400)
    np.testing.assert_array_equal(new[0, [0, -1]], line[0, [0, -1]])
    first = new[0, 1] - new[0, 0]
    assert first[1] / np.linalg.norm(first) > 0.99  # leaves along the normal
    arc = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(new[0], axis=0), axis=1))])
    assert np.all(np.diff(arc) > 0)  # nodes keep their order


def test_orthogonal_walls_on_a_box_changes_nothing():
    from shore.mesh import blocking as bk
    from shore.mesh.association import Association

    class Flat:  # projection onto the plane z = 0
        def project(self, assoc, pts):
            pts = np.array(pts, dtype=float)
            pts[..., 2] = 0.0
            return pts

    b = bk.associate(bk.box((0, 0, 0), (1, 1, 1), 4, 4, 4), "b0", "k_lo", Association.body())
    blocks = bk.premesh(b, Flat())
    before = blocks[0].nodes.copy()
    assert orthogonal_walls(blocks) == ["b0"]
    np.testing.assert_allclose(blocks[0].nodes, before, atol=1e-12)  # already orthogonal
    with pytest.raises(ParameterError):
        orthogonal_walls(blocks, 0.0)


def test_layer_fractions():
    from shore.mesh.orthogonal import layer_fractions

    t = layer_fractions(10.0, 6, 0.5)
    cells = np.diff(t) * 10.0
    assert t[0] == 0.0 and t[-1] == 1.0
    np.testing.assert_allclose(cells[0], 0.5, rtol=1e-12)
    ratios = cells[1:] / cells[:-1]
    np.testing.assert_allclose(ratios, ratios[0], rtol=1e-9)  # one geometric ratio
    assert ratios[0] > 1.0
    with pytest.raises(ParameterError):
        layer_fractions(1.0, 6, 0.2)  # = the uniform spacing


def test_wall_layer_on_a_box():
    from shore.mesh import blocking as bk
    from shore.mesh.association import Association
    from shore.mesh.orthogonal import wall_layer

    class Flat:
        def project(self, assoc, pts):
            pts = np.array(pts, dtype=float)
            pts[..., 2] = 0.0
            return pts

    b = bk.associate(bk.box((0, 0, 0), (1, 1, 1), 3, 3, 9), "b0", "k_lo", Association.body())
    blocks = bk.premesh(b, Flat())
    assert wall_layer(blocks, 0.01) == ["b0"]
    z = blocks[0].nodes[1, 1, :, 2]
    np.testing.assert_allclose(z[1], 0.01, rtol=1e-12)
    assert z[-1] == 1.0 and np.all(np.diff(np.diff(z)) > 0)
    with pytest.raises(ParameterError):
        wall_layer(blocks, -1.0)
