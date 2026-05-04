"""Tests for shore.topology.oh_channel.build_channel_block — frustum builder.

The frustum builder takes:
    - inner_face: the wrap block's k_hi outer envelope, shape (ni, nj, 3)
    - outer_quad: 4 corners of the box face this channel block "owns"
    - nk_channel + spacing law
    - optional c1_first_step (pin to wrap's last step for C1 continuity)

and returns a (ni, nj, nk_channel, 3) frustum whose:
    - k=0 slice == inner_face exactly
    - k=nk-1 slice == bilinear (a, b) parameterisation of outer_quad
    - interior is ruled along k.

These tests are pure-numpy — no STL, no Projector, no marching.
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.topology.oh_channel import build_channel_block


def _fake_inner_face(ni: int, nj: int, *, scale: float = 1.0) -> np.ndarray:
    """Synthetic 'wrap envelope' patch on a sphere (ni, nj, 3).

    A small lat-lon patch centred on +z, parametrically (i, j) ↔ (lat, lon).
    """
    lats = np.linspace(0.3, 0.6, ni)  # colatitude band
    lons = np.linspace(-0.2, 0.2, nj)
    la, lo = np.meshgrid(lats, lons, indexing="ij")
    x = scale * np.sin(la) * np.cos(lo)
    y = scale * np.sin(la) * np.sin(lo)
    z = scale * np.cos(la)
    return np.stack([x, y, z], axis=-1)


def _box_face_quad_zp(half: float = 1.0) -> np.ndarray:
    """4 corners of the +z face of a centred [-half, +half]^3 box.

    Corners ordered (a=0,b=0), (a=1,b=0), (a=1,b=1), (a=0,b=1) — i.e.
    the same bilinear (a, b) parameterisation as build_cap_k0.
    """
    return np.array(
        [
            [-half, -half, half],
            [+half, -half, half],
            [+half, +half, half],
            [-half, +half, half],
        ],
        dtype=np.float64,
    )


class TestChannelShape:
    def test_returns_correct_shape(self):
        inner = _fake_inner_face(5, 7)
        outer = _box_face_quad_zp()
        block = build_channel_block(inner, outer, nk_channel=4)
        assert block.shape == (5, 7, 4, 3)

    def test_dtype_is_float64(self):
        inner = _fake_inner_face(4, 4)
        outer = _box_face_quad_zp()
        block = build_channel_block(inner, outer, nk_channel=3)
        assert block.dtype == np.float64


class TestChannelInnerFaceBitExact:
    """k=0 slice must equal inner_face exactly — bit-exact, no smoothing."""

    def test_k0_is_inner_face(self):
        inner = _fake_inner_face(6, 8)
        outer = _box_face_quad_zp()
        block = build_channel_block(inner, outer, nk_channel=5)
        np.testing.assert_array_equal(block[:, :, 0, :], inner)

    def test_k0_is_inner_face_with_c1_pin(self):
        """C1 pin only affects interior spacing; k=0 stays bit-exact."""
        inner = _fake_inner_face(6, 8)
        outer = _box_face_quad_zp()
        block = build_channel_block(inner, outer, nk_channel=5, c1_first_step=0.05)
        np.testing.assert_array_equal(block[:, :, 0, :], inner)


class TestChannelOuterFaceBilinear:
    """k=nk-1 slice must be a bilinear (a, b) sweep of the outer quad."""

    def test_outer_face_corners_match_quad(self):
        ni, nj = 5, 7
        inner = _fake_inner_face(ni, nj)
        outer = _box_face_quad_zp(half=2.5)
        block = build_channel_block(inner, outer, nk_channel=3)
        # (a=0, b=0) → outer[0]; (a=1, b=0) → outer[1];
        # (a=1, b=1) → outer[2]; (a=0, b=1) → outer[3].
        # Parametric (a, b) = (i / (ni-1), j / (nj-1)).
        np.testing.assert_allclose(block[0, 0, -1, :], outer[0], atol=1e-12)
        np.testing.assert_allclose(block[-1, 0, -1, :], outer[1], atol=1e-12)
        np.testing.assert_allclose(block[-1, -1, -1, :], outer[2], atol=1e-12)
        np.testing.assert_allclose(block[0, -1, -1, :], outer[3], atol=1e-12)

    def test_outer_face_lies_on_box_face_plane(self):
        """All outer-face nodes have z = +half on the +z box face."""
        inner = _fake_inner_face(5, 5)
        half = 1.5
        outer = _box_face_quad_zp(half=half)
        block = build_channel_block(inner, outer, nk_channel=4)
        outer_layer = block[:, :, -1, :]
        np.testing.assert_allclose(outer_layer[..., 2], half, atol=1e-12)

    def test_outer_face_uniform_param(self):
        """Default param: (a, b) is uniform → outer-face is a uniform grid."""
        inner = _fake_inner_face(4, 4)
        outer = _box_face_quad_zp()
        block = build_channel_block(inner, outer, nk_channel=3)
        # Outer face: x along i, y along j, z = const.  Spacings are uniform.
        x = block[:, 0, -1, 0]
        dx = np.diff(x)
        assert np.allclose(dx, dx[0], atol=1e-12)


class TestChannelInteriorRuled:
    """Interior k-layers are ruled-surface interpolation between inner / outer."""

    def test_uniform_law_midlayer_is_midpoint(self):
        """With nk_channel=3 and uniform law, layer k=1 = 0.5 * (inner + outer)."""
        ni, nj = 4, 4
        inner = _fake_inner_face(ni, nj)
        outer = _box_face_quad_zp()
        block = build_channel_block(inner, outer, nk_channel=3, k_law="uniform")
        # The bilinear sweep of outer at the same (i, j) uses the same (a, b)
        # grid as the channel, so we can compute it explicitly:
        a = np.linspace(0.0, 1.0, ni)[:, None]
        b = np.linspace(0.0, 1.0, nj)[None, :]
        outer_layer = (
            (1 - a)[..., None] * (1 - b)[..., None] * outer[0]
            + a[..., None] * (1 - b)[..., None] * outer[1]
            + a[..., None] * b[..., None] * outer[2]
            + (1 - a)[..., None] * b[..., None] * outer[3]
        )
        expected_mid = 0.5 * (inner + outer_layer)
        np.testing.assert_allclose(block[:, :, 1, :], expected_mid, atol=1e-12)

    def test_monotonic_distance_from_inner(self):
        """At any (i, j), node-to-inner distance grows monotonically in k."""
        inner = _fake_inner_face(5, 5)
        outer = _box_face_quad_zp(half=2.0)
        block = build_channel_block(inner, outer, nk_channel=6)
        # Distance from each node to its k=0 anchor.
        dist = np.linalg.norm(block - block[:, :, :1, :], axis=-1)  # (ni, nj, nk)
        # Should be monotonically nondecreasing along k for every (i, j).
        for k in range(dist.shape[2] - 1):
            assert np.all(dist[:, :, k] <= dist[:, :, k + 1] + 1e-12)


class TestChannelSpacingLaws:
    def test_tanh_clusters_near_inner(self):
        inner = _fake_inner_face(4, 4)
        outer = _box_face_quad_zp(half=2.0)
        block = build_channel_block(inner, outer, nk_channel=20, k_law="tanh", k_beta=4.0)
        # Distance-from-inner along k at the centre of the (i, j) patch.
        ic, jc = 2, 2
        dist = np.linalg.norm(block[ic, jc, :, :] - block[ic, jc, 0, :], axis=-1)
        ds = np.diff(dist)
        # tanh (one-sided, cluster near k=0) → ds[0] significantly < ds[-1].
        assert ds[0] < 0.5 * ds[-1]

    def test_c1_pin_matches_first_step(self):
        """With c1_first_step set, ||x[k=1] - x[k=0]|| ≈ first_step at every (i,j)."""
        inner = _fake_inner_face(4, 4)
        outer = _box_face_quad_zp(half=2.0)
        first_step = 0.07
        block = build_channel_block(inner, outer, nk_channel=10, c1_first_step=first_step)
        # Distance from k=0 to k=1 should equal first_step (within FP noise).
        ds_0 = np.linalg.norm(block[:, :, 1, :] - block[:, :, 0, :], axis=-1)
        np.testing.assert_allclose(ds_0, first_step, atol=1e-10)


class TestChannelErrors:
    def test_invalid_nk_rejected(self):
        inner = _fake_inner_face(4, 4)
        outer = _box_face_quad_zp()
        with pytest.raises(ParameterError, match="nk_channel"):
            build_channel_block(inner, outer, nk_channel=1)

    def test_outer_quad_wrong_shape_rejected(self):
        inner = _fake_inner_face(4, 4)
        with pytest.raises(ParameterError, match="outer_quad"):
            build_channel_block(inner, np.zeros((3, 3)), nk_channel=3)

    def test_inner_face_wrong_rank_rejected(self):
        outer = _box_face_quad_zp()
        with pytest.raises(ParameterError, match="inner_face"):
            build_channel_block(np.zeros((4, 4)), outer, nk_channel=3)

    def test_unknown_law_rejected(self):
        inner = _fake_inner_face(4, 4)
        outer = _box_face_quad_zp()
        with pytest.raises(ParameterError, match="k_law"):
            build_channel_block(inner, outer, nk_channel=3, k_law="bogus")

    def test_negative_c1_first_step_rejected(self):
        inner = _fake_inner_face(4, 4)
        outer = _box_face_quad_zp()
        with pytest.raises(ParameterError, match="c1_first_step"):
            build_channel_block(inner, outer, nk_channel=4, c1_first_step=-0.01)
