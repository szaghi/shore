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


def test_ogrid_one_block_layers():
    b = bk.box((-1, -1, -1), (1, 1, 1), 5, 5, 5)
    o = bk.ogrid(b, ["b0"], factor=0.25, nodes=3)
    np.testing.assert_allclose(o.points[list(o.block("b0"))], 0.25 * b.points[:8], atol=1e-15)
    assert o.labels() == ["b0", "b1", "b2", "b3", "b4", "b5", "b6"]
    assert all(o.shape(label)[0] == 3 for label in ("b1", "b2"))
    blocks = bk.premesh(o, None)
    assert min(float(scaled_jacobian(x.nodes).min()) for x in blocks) > 0.5
    assert b.labels() == ["b0"]  # the input is not modified


def test_ogrid_layer_handedness_and_seams():
    from shore.mesh.seams import seams

    o = bk.ogrid(bk.box((0, 0, 0), (2, 1, 1), 5, 4, 3), factor=0.6, nodes=3)
    blocks = bk.premesh(o, None)
    assert min(float(scaled_jacobian(x.nodes).min()) for x in blocks) > 0.0
    assert sum(len(v) for v in seams(blocks, 1e-12).values()) == 36


def test_ogrid_open_face_keeps_plane():
    b = bk.box((0, 0, 0), (1, 1, 1), 5, 5, 5)
    o = bk.ogrid(b, exclude=[("b0", "k_lo")], factor=0.5, nodes=3)
    core = o.points[list(o.block("b0"))]
    np.testing.assert_array_equal(core[0::2, 2], 0.0)
    assert len(o.labels()) == 6


def _cut(b, axis, x):
    """Split the first block spanning coordinate *x* along *axis*."""
    d = "ijk".index(axis)
    for label, v in b.blocks:
        lo, hi = b.points[v[0]][d], b.points[v[bk._BIT[axis]]][d]
        if lo < x < hi:
            return bk.split(b, label, axis, (x - lo) / (hi - lo))
    raise AssertionError(x)


def test_ogrid_hole_moves_away_from_it():
    b = bk.box((-3, -3, -3), (3, 3, 3), 7, 7, 7)
    for axis in "ijk":
        b = _cut(_cut(b, axis, -1.0), axis, 1.0)
    centre = next(lab for lab, v in b.blocks if np.allclose(b.points[list(v)].mean(axis=0), 0))
    hole = {b.entity(centre, f)[0] for f in bk.FACES}
    b = bk.delete(b, [centre])
    faces = [(lab, f) for lab, _ in b.blocks for f in bk.FACES if b.entity(lab, f)[0] in hole]
    assert len(faces) == 6
    o = bk.ogrid(b, faces=faces, factor=0.5, nodes=3)
    assert len(o.labels()) == 32
    moved = o.points[len(b.points) :]
    assert len(moved) == 8
    # (1, 1, 1) moves by (1 - 0.5) * ((2, 0, 0) + (0, 2, 0) + (0, 0, 2)) / 2
    assert any(np.allclose(p, [1.5, 1.5, 1.5]) for p in moved)


def test_ogrid_refusals():
    b = bk.split(bk.box((0, 0, 0), (2, 1, 1), 5, 3, 3), "b0", "i", 0.5)
    for kwargs in (
        {"factor": 0.0},
        {"nodes": 1},
        {"labels": ["b0"], "exclude": [("b0", "i_hi")]},
        {"labels": ["b0", "b1"], "faces": [("b0", "i_hi")]},
        {"labels": ["b0"], "exclude": [("b0", f) for f in bk.FACES]},
    ):
        with pytest.raises(ParameterError):
            bk.ogrid(b, **kwargs)


def test_centres():
    b = bk.split(bk.box((0, 0, 0), (2, 1, 1), 5, 3, 3), "b0", "i", 0.5)
    blocks = b.centres("block")
    assert [(lab, f) for lab, f, _ in blocks] == [("b0", ""), ("b1", "")]
    np.testing.assert_allclose(blocks[1][2], [1.5, 0.5, 0.5])
    faces = b.centres("face")
    assert len(faces) == 11
    shared = [c for lab, f, c in faces if (lab, f) == ("b0", "i_hi")]
    np.testing.assert_allclose(shared[0], [1.0, 0.5, 0.5])
    assert ("b1", "i_lo") not in {(lab, f) for lab, f, _ in faces}  # once, by its first block
