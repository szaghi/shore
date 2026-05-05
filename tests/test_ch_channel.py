"""Unit tests for the CH channel-block helpers (resolve_frame + outer_quads_for_cgrid)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.topology.ch_channel import (
    CHFrame,
    outer_quads_for_cgrid,
    resolve_frame,
)

# ---------------------------------------------------------------------------
# resolve_frame
# ---------------------------------------------------------------------------


class TestResolveFrame:
    def test_default_freestream_x_chord(self):
        frame = resolve_frame(
            freestream_world=None,
            span_axis="z",
            chord_min=-1.0,
            chord_max=2.0,
            normal_min=-1.0,
            normal_max=1.0,
        )
        assert frame.chord_axis == 0
        assert frame.normal_axis == 1
        assert frame.span_axis_idx == 2
        assert frame.chord_sign == 1

    def test_negative_freestream_sets_chord_sign(self):
        frame = resolve_frame(
            freestream_world=np.array([-1.0, 0.0, 0.0]),
            span_axis="z",
            chord_min=-1.0,
            chord_max=2.0,
            normal_min=-1.0,
            normal_max=1.0,
        )
        assert frame.chord_sign == -1

    def test_y_chord_with_z_span(self):
        frame = resolve_frame(
            freestream_world=np.array([0.0, 1.0, 0.0]),
            span_axis="z",
            chord_min=-1.0,
            chord_max=2.0,
            normal_min=-1.0,
            normal_max=1.0,
        )
        assert frame.chord_axis == 1
        assert frame.normal_axis == 0

    def test_x_chord_with_y_span(self):
        frame = resolve_frame(
            freestream_world=np.array([1.0, 0.0, 0.0]),
            span_axis="y",
            chord_min=-1.0,
            chord_max=2.0,
            normal_min=-1.0,
            normal_max=1.0,
        )
        assert frame.chord_axis == 0
        assert frame.normal_axis == 2

    def test_invalid_span_axis(self):
        with pytest.raises(ParameterError, match="span_axis must be"):
            resolve_frame(
                freestream_world=None,
                span_axis="w",
                chord_min=-1.0,
                chord_max=1.0,
                normal_min=-1.0,
                normal_max=1.0,
            )

    def test_inverted_extents(self):
        with pytest.raises(ParameterError, match="chord_max must exceed"):
            resolve_frame(
                freestream_world=None,
                span_axis="z",
                chord_min=2.0,
                chord_max=1.0,
                normal_min=-1.0,
                normal_max=1.0,
            )


# ---------------------------------------------------------------------------
# outer_quads_for_cgrid
# ---------------------------------------------------------------------------


def _fake_envelope(ni: int, nj: int, *, chord_range, normal, span_lo, span_hi) -> np.ndarray:
    """Build a fake (ni, nj, 3) envelope with chord/normal/span = (x, y, z)."""
    chord = np.linspace(chord_range[0], chord_range[1], ni)
    span = np.linspace(span_lo, span_hi, nj)
    nodes = np.zeros((ni, nj, 3))
    nodes[:, :, 0] = chord[:, None]  # chord = x
    nodes[:, :, 1] = normal  # constant normal
    nodes[:, :, 2] = span[None, :]  # span = z
    return nodes


def _fake_three_envelopes() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build (upper, body, lower) envelope triplets for outer-quad tests."""
    e_kw = {"span_lo": 0.0, "span_hi": 1.0}
    upper = _fake_envelope(20, 5, chord_range=(15.0, 1.0), normal=0.5, **e_kw)
    body = _fake_envelope(80, 5, chord_range=(1.0, 0.0), normal=0.5, **e_kw)
    body[:80, :, 1] = np.linspace(0.5, -0.5, 80)[:, None]
    lower = _fake_envelope(20, 5, chord_range=(0.0, 15.0), normal=-0.5, **e_kw)
    return upper, body, lower


class TestOuterQuads:
    @pytest.fixture
    def frame(self):
        return resolve_frame(
            freestream_world=None,
            span_axis="z",
            chord_min=-2.0,
            chord_max=15.0,
            normal_min=-3.0,
            normal_max=3.0,
        )

    def test_ch3_yields_7_quads(self, frame):
        upper, body, lower = _fake_three_envelopes()
        quads = outer_quads_for_cgrid(
            frame=frame,
            upper_envelope=upper,
            body_envelope=body,
            lower_envelope=lower,
            i_split_body_upper=20,
            i_split_body_lower=60,
            span_lo=0.0,
            span_hi=1.0,
            is_blunt=False,
        )
        # Public quad labels for outer-face anchors of each channel.
        # (Sentinel "__outflow_*_inner__" entries are documentation
        # carriers used by fill_channel_blocks; tested separately.)
        public = {k for k in quads if not k.startswith("__")}
        assert public == {
            "upper",
            "body_top",
            "body_front",
            "body_bot",
            "lower",
            "outflow_upper",
            "outflow_lower",
        }

    def test_quad_corners_have_shape_4_3(self, frame):
        upper, body, lower = _fake_three_envelopes()
        quads = outer_quads_for_cgrid(
            frame=frame,
            upper_envelope=upper,
            body_envelope=body,
            lower_envelope=lower,
            i_split_body_upper=20,
            i_split_body_lower=60,
            span_lo=0.0,
            span_hi=1.0,
            is_blunt=False,
        )
        for label, q in quads.items():
            assert q.shape == (4, 3), f"{label}: shape {q.shape}"

    def test_chan_upper_quad_lies_on_top_frame(self, frame):
        upper, body, lower = _fake_three_envelopes()
        quads = outer_quads_for_cgrid(
            frame=frame,
            upper_envelope=upper,
            body_envelope=body,
            lower_envelope=lower,
            i_split_body_upper=20,
            i_split_body_lower=60,
            span_lo=0.0,
            span_hi=1.0,
            is_blunt=False,
        )
        # All 4 corners of the chan_upper quad should sit at normal_max (y).
        np.testing.assert_allclose(quads["upper"][:, 1], frame.normal_max)

    def test_chan_body_front_quad_lies_on_front_frame(self, frame):
        upper, body, lower = _fake_three_envelopes()
        quads = outer_quads_for_cgrid(
            frame=frame,
            upper_envelope=upper,
            body_envelope=body,
            lower_envelope=lower,
            i_split_body_upper=20,
            i_split_body_lower=60,
            span_lo=0.0,
            span_hi=1.0,
            is_blunt=False,
        )
        # All 4 corners of chan_body_front sit at chord_min (x).
        np.testing.assert_allclose(quads["body_front"][:, 0], frame.chord_min)

    def test_invalid_split_indices_rejected(self, frame):
        upper, body, lower = _fake_three_envelopes()
        with pytest.raises(ParameterError, match="invalid body splits"):
            outer_quads_for_cgrid(
                frame=frame,
                upper_envelope=upper,
                body_envelope=body,
                lower_envelope=lower,
                i_split_body_upper=60,  # > i_split_body_lower
                i_split_body_lower=20,
                span_lo=0.0,
                span_hi=1.0,
                is_blunt=False,
            )

    def test_blunt_requires_wake_base_envelope(self, frame):
        upper, body, lower = _fake_three_envelopes()
        with pytest.raises(ParameterError, match="wake_base_envelope is required"):
            outer_quads_for_cgrid(
                frame=frame,
                upper_envelope=upper,
                body_envelope=body,
                lower_envelope=lower,
                i_split_body_upper=20,
                i_split_body_lower=60,
                span_lo=0.0,
                span_hi=1.0,
                is_blunt=True,
            )


# ---------------------------------------------------------------------------
# CHFrame is a dataclass-like immutable carrier
# ---------------------------------------------------------------------------


class TestCHFrameDataclass:
    def test_frame_attributes(self):
        frame = CHFrame(
            chord_axis=0,
            normal_axis=1,
            span_axis_idx=2,
            chord_sign=1,
            chord_min=-1.0,
            chord_max=2.0,
            normal_min=-1.0,
            normal_max=1.0,
        )
        assert frame.chord_axis == 0
        assert frame.chord_sign == 1
