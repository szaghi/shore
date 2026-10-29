"""Unit tests of the blocking skeleton and its premesh (ADR 0042)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh import blocking as bk
from shore.mesh.association import Association
from shore.quality import scaled_jacobian


def test_edge_and_face_tables():
    assert bk.EDGES["j_lo_k_lo"] == ("i", 0, 4)
    assert bk.EDGES["i_hi_j_hi"] == ("k", 6, 7)
    assert bk.face_corners("i_hi").tolist() == [[4, 5], [6, 7]]
    assert bk.face_corners("j_lo").tolist() == [[0, 1], [4, 5]]


def test_box_and_entities():
    b = bk.box((0, 0, 0), (1, 2, 3), 3, 4, 5)
    assert b.shape("b0") == (3, 4, 5)
    np.testing.assert_array_equal(b.points[7], [1, 2, 3])
    assert b.entity("b0", "v5") == (("v", 5), "vertex")
    assert b.entity("b0", "i_lo_k_hi") == (("e", 1, 3), "edge")
    assert b.entity("b0", "k_lo") == (("f", 0, 2, 4, 6), "face")
    with pytest.raises(ParameterError):
        b.entity("b0", "v8")
    with pytest.raises(ParameterError):
        bk.box((0, 0, 0), (1, 0, 1), 3, 3, 3)
    with pytest.raises(ParameterError):
        bk.box((0, 0, 0), (1, 1, 1), 1, 3, 3)


def test_split_keeps_counts_and_associations():
    b = bk.box((0, 0, 0), (2, 1, 1), 9, 5, 5)
    b = bk.associate(b, "b0", "j_lo_k_lo", Association.body())
    b = bk.associate(b, "b0", "k_lo", Association.body("other"))
    s = bk.split(b, "b0", "i", 0.25)
    assert s.shape("b0") == (3, 5, 5) and s.shape("b1") == (7, 5, 5)
    assert ("e", 0, 4) not in s.associations  # the cut edge's halves carry it
    halves = [k for k, a in s.associations.items() if k[0] == "e"]
    assert len(halves) == 2
    faces = [a.geometry for k, a in s.associations.items() if k[0] == "f"]
    assert faces == ["other", "other"]
    assert b.labels() == ["b0"]  # the input is not modified


def test_split_refusals():
    b = bk.box((0, 0, 0), (1, 1, 1), 3, 3, 3)
    for args in (("b0", "q", 0.5), ("b0", "i", 0.0), ("b0", "i", 0.1), ("bx", "i", 0.5)):
        with pytest.raises(ParameterError):
            bk.split(b, *args)
    with pytest.raises(ParameterError):
        bk.delete(b, ["b0"])


def test_tfi_face_bilinear():
    a = np.linspace(0, 1, 4)[:, None]
    b = np.linspace(0, 2, 3)[:, None]
    i_lo = np.hstack([np.zeros_like(b), b, np.zeros_like(b)])
    i_hi = np.hstack([np.ones_like(b), b, np.zeros_like(b)])
    j_lo = np.hstack([a, np.zeros_like(a), np.zeros_like(a)])
    j_hi = np.hstack([a, 2 * np.ones_like(a), np.zeros_like(a)])
    grid = bk.tfi_face(i_lo, i_hi, j_lo, j_hi)
    assert grid.shape == (4, 3, 3)
    np.testing.assert_allclose(grid[2, 1], [2 / 3, 1.0, 0.0], atol=1e-15)


def test_tfi_volume_reproduces_trilinear():
    x = np.stack(np.meshgrid(*(np.linspace(0, 1, n) for n in (4, 5, 6)), indexing="ij"), -1)
    x = x @ np.array([[1.0, 0.2, 0.0], [0.0, 1.0, 0.3], [0.1, 0.0, 1.0]])
    y = x.copy()
    y[1:-1, 1:-1, 1:-1] = 0.0
    bk.tfi_volume(y)
    np.testing.assert_allclose(y, x, atol=1e-14)


def rotated_pair():
    """Two unit cubes sharing x = 1; the second's axes are (z, -y, x)."""
    pts = np.array([[x, y, z] for x in (0, 1, 2) for y in (0, 1) for z in (0, 1)], float)
    idx = {tuple(p): n for n, p in enumerate(pts)}
    a = tuple(idx[(i, j, k)] for i in (0, 1) for j in (0, 1) for k in (0, 1))
    b = tuple(idx[(1 + k, 1 - j, i)] for i in (0, 1) for j in (0, 1) for k in (0, 1))
    blocking = bk.Blocking(pts, [("A", a), ("B", b)], {})
    for edge in blocking.edges():
        blocking.nodes[edge] = 4
    return blocking


def test_rotated_neighbour_axis_map_matches_nodes():
    from shore.mesh.face import _face_slice

    blocks = bk.premesh(rotated_pair(), None)
    a, b = blocks
    assert min(float(scaled_jacobian(x.nodes).min()) for x in blocks) == pytest.approx(1.0)
    face = a.face("i_hi")
    assert face.partner is b.face("k_lo")
    assert face.axis_map.along == (("j", True), ("i", False))
    mine, theirs = _face_slice(a, "i_hi"), _face_slice(b, "k_lo")  # (j, k) and (i, j)
    np.testing.assert_array_equal(mine, theirs.transpose(1, 0, 2)[::-1])


def test_effective_associations():
    b = bk.box((0, 0, 0), (1, 1, 1), 3, 3, 3)
    b = bk.associate(b, "b0", "i_lo", Association.body("a"))
    b = bk.associate(b, "b0", "j_lo", Association.body("b"))
    eff = bk.effective_associations(b)
    assert ("e", 0, 1) not in eff and ("v", 0) not in eff
    assert eff[("e", 2, 3)].geometry == "a"  # i_lo_j_hi: i_lo only
    assert eff[("v", 3)].geometry == "a"
    b = bk.associate(b, "b0", "i_lo_j_lo", Association.body("c"))
    eff = bk.effective_associations(b)
    assert eff[("e", 0, 1)].geometry == "c" and eff[("v", 0)].geometry == "c"


def test_premesh_needs_a_resolver_for_associations():
    b = bk.associate(bk.box((0, 0, 0), (1, 1, 1), 3, 3, 3), "b0", "k_lo", Association.body())
    with pytest.raises(ParameterError):
        bk.premesh(b, None)


def test_polylines():
    b = bk.split(bk.box((0, 0, 0), (2, 1, 1), 9, 5, 5), "b0", "i", 0.5)
    lines = b.polylines()
    assert len(lines) == 20
    *_, pts = next(x for x in lines if x[:2] == ("b1", "j_lo_k_lo"))
    np.testing.assert_allclose(pts[:, 0], np.linspace(1, 2, 5))
