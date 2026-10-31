"""Tests for the CHTopology (C-grid wrap + H-channel far-field)."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh import Body, CGridTopology, CHTopology, FaceBC, Mesh

_CH_KW = {
    "chord_min": -2.0,
    "chord_max": 20.0,  # wake_length=15 + body chord ~ 1, +pad
    "normal_min": -2.5,
    "normal_max": 2.5,
    "nk_channel": 5,
}


def _build_ch_mesh(
    stl_path,
    *,
    nk=10,
    ni_body=80,
    ni_wake=20,
    n_stations=5,
    wake_length=15.0,
    nk_channel=5,
    k_law_channel="uniform",
    c1_seam=True,
    chord_min=-2.0,
    chord_max=20.0,
    normal_min=-2.5,
    normal_max=2.5,
    march=True,
):
    """Build (and optionally march) a CH-topology mesh from an STL."""
    body = Body.from_stl(stl_path, skip_hygiene=True)
    topo = CHTopology(
        chord_min=chord_min,
        chord_max=chord_max,
        normal_min=normal_min,
        normal_max=normal_max,
        nk_channel=nk_channel,
        k_law_channel=k_law_channel,
        c1_seam=c1_seam,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        blocks = topo.build(
            surface=None,
            nk=nk,
            mesh=body.mesh,
            ni_body=ni_body,
            ni_wake=ni_wake,
            n_stations=n_stations,
            wake_length=wake_length,
        )
        m = Mesh(body=body, blocks=blocks, topology=topo)
        if march:
            m.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)
    return m, topo


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


class TestCHConstruction:
    """Block count, labels, and shapes."""

    def test_sharp_te_yields_12_blocks(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, march=False)
        assert len(m.blocks) == 12
        assert [b.label for b in m.blocks] == [
            "upper",
            "body",
            "lower",
            "chan_upper",
            "chan_body_top",
            "chan_body_front",
            "chan_body_bot",
            "chan_lower",
            "chan_outflow_top_upper",
            "chan_outflow_wake_upper",
            "chan_outflow_wake_lower",
            "chan_outflow_top_lower",
        ]

    def test_blunt_te_yields_14_blocks(self, naca0012_blunt_stl):
        m, _ = _build_ch_mesh(naca0012_blunt_stl, march=False)
        assert len(m.blocks) == 14
        labels = [b.label for b in m.blocks]
        assert labels[:4] == ["upper", "body", "lower", "wake_base"]
        assert "chan_outflow_wake_base" in labels
        assert "chan_outflow_top_upper" in labels
        assert "chan_outflow_wake_upper" in labels

    def test_wrap_block_shapes(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, nk=10, ni_body=80, ni_wake=20, march=False)
        assert m.blocks[0].nodes.shape == (20, 5, 10, 3)  # upper
        assert m.blocks[1].nodes.shape == (80, 5, 10, 3)  # body
        assert m.blocks[2].nodes.shape == (20, 5, 10, 3)  # lower

    def test_channel_block_shapes(self, naca0012_stl):
        m, _ = _build_ch_mesh(
            naca0012_stl,
            nk=10,
            ni_body=80,
            ni_wake=20,
            nk_channel=5,
            march=False,
        )
        # k_hi-side channel depth = nk_channel; i-extent matches the wrap slice.
        i_split_upper = 80 // 3  # 26
        i_split_lower = (2 * 80) // 3  # 53
        assert m.blocks[3].nodes.shape == (20, 5, 5, 3)  # chan_upper
        assert m.blocks[4].nodes.shape == (i_split_upper + 1, 5, 5, 3)  # chan_body_top
        assert m.blocks[5].nodes.shape == (
            i_split_lower - i_split_upper + 1,
            5,
            5,
            3,
        )  # chan_body_front
        assert m.blocks[6].nodes.shape == (80 - i_split_lower, 5, 5, 3)  # chan_body_bot
        assert m.blocks[7].nodes.shape == (20, 5, 5, 3)  # chan_lower
        # Outflow channels (4-block conformal layout):
        #   top_*  : (ni_outflow, n_stations, nk_channel, 3)
        #   wake_* : (ni_outflow, n_stations, nk_wrap,    3)
        # ni_outflow defaults to 8.  Order: top_upper, wake_upper,
        # wake_lower, top_lower.
        assert m.blocks[8].nodes.shape == (8, 5, 5, 3)  # chan_outflow_top_upper
        assert m.blocks[9].nodes.shape == (8, 5, 10, 3)  # chan_outflow_wake_upper
        assert m.blocks[10].nodes.shape == (8, 5, 10, 3)  # chan_outflow_wake_lower
        assert m.blocks[11].nodes.shape == (8, 5, 5, 3)  # chan_outflow_top_lower

    def test_blunt_channel_shapes(self, naca0012_blunt_stl):
        m, _ = _build_ch_mesh(
            naca0012_blunt_stl,
            nk=8,
            ni_body=60,
            ni_wake=15,
            nk_channel=4,
            march=False,
        )
        # wake_base_outflow has i = ni_outflow (default 8), j = n_stations,
        # k = nk_w (= 3 by default, the wake_base wrap k-extent).
        wb_chan = next(b for b in m.blocks if b.label == "chan_outflow_wake_base")
        assert wb_chan.nodes.shape == (8, 5, 3, 3)
        # Wake-side outflow channels conformal to wrap: nk = nk_wrap (= 8).
        wake_up = next(b for b in m.blocks if b.label == "chan_outflow_wake_upper")
        wake_lo = next(b for b in m.blocks if b.label == "chan_outflow_wake_lower")
        assert wake_up.nodes.shape == (8, 5, 8, 3)
        assert wake_lo.nodes.shape == (8, 5, 8, 3)


# ---------------------------------------------------------------------------
# View aliasing (k_hi-side channels share memory with wrap)
# ---------------------------------------------------------------------------


class TestSharedViewAliasing:
    """Wrap.k_hi and chan.k_lo are bit-identical via numpy view aliasing."""

    def test_chan_upper_aliased_to_upper_khi(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, march=True)
        upper = m.blocks[0]
        chan_upper = m.blocks[3]
        # Modify a node through chan_upper.k=0; check upper.k=nk-1 sees it.
        original = chan_upper.nodes[5, 2, 0, :].copy()
        chan_upper.nodes[5, 2, 0, :] = -999.0
        assert upper.nodes[5, 2, -1, 0] == -999.0
        chan_upper.nodes[5, 2, 0, :] = original  # restore

    def test_chan_body_top_aliased_to_body_khi_slice(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, march=True)
        body = m.blocks[1]
        chan_body_top = m.blocks[4]
        # chan_body_top covers body i = 0..i_split_upper.
        assert np.shares_memory(chan_body_top.nodes, body.nodes)
        # Bit-exact equality on the shared slice.
        np.testing.assert_array_equal(
            chan_body_top.nodes[:, :, 0, :],
            body.nodes[: chan_body_top.nodes.shape[0], :, -1, :],
        )

    def test_body_split_overlap_is_aliased(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, march=True)
        chan_body_top = m.blocks[4]
        chan_body_front = m.blocks[5]
        # Last column of chan_body_top must equal first column of chan_body_front.
        np.testing.assert_array_equal(
            chan_body_top.nodes[-1, :, :, :],
            chan_body_front.nodes[0, :, :, :],
        )

    def test_outflow_wake_seam_geometric_match(self, naca0012_stl):
        # Phase-1.2 layout: upper and lower outflow stacks live in
        # separate buffers, joined at the wake-cut by SHARED face-graph.
        # The wake-cut row data must be bit-identical between
        # wake_upper.k=0 and wake_lower.k=0 (geometric coincidence
        # comes from the wrap supplying both sides identical wake-cut
        # data, which fill_channel_blocks copies from upper.i_lo / lower.i_hi).
        m, _ = _build_ch_mesh(naca0012_stl, march=True)
        labels = [b.label for b in m.blocks]
        wake_upper = m.blocks[labels.index("chan_outflow_wake_upper")]
        wake_lower = m.blocks[labels.index("chan_outflow_wake_lower")]
        # the upper side's j runs against the lower side's (ADR 0057: right-handed)
        np.testing.assert_allclose(
            wake_upper.nodes[:, ::-1, 0, :],
            wake_lower.nodes[:, :, 0, :],
            atol=1e-9,
        )


# ---------------------------------------------------------------------------
# Face BCs / wiring
# ---------------------------------------------------------------------------


class TestCHBCs:
    """Face boundary conditions are correctly wired."""

    def test_wrap_k_lo_wall_on_body(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, march=False)
        body = m.blocks[1]
        assert body.face("k_lo").bc == FaceBC.WALL

    def test_wrap_khi_shared_with_chan_klo(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, march=False)
        upper = m.blocks[0]
        chan_upper = m.blocks[3]
        assert upper.face("k_hi").bc == FaceBC.SHARED
        assert chan_upper.face("k_lo").bc == FaceBC.SHARED

    def test_chan_khi_is_free(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, march=False)
        # k_hi-side channels (indices 3..7) all have k_hi = FREE
        # (top / bottom / front far-field frame edges).
        for chan in m.blocks[3:8]:
            assert chan.face("k_hi").bc == FaceBC.FREE, (
                f"{chan.label}.k_hi should be FREE, got {chan.face('k_hi').bc}"
            )
        # Outflow channels (4-block layout, channel-k pointing outward
        # from wake-cut on each side):
        #   top_upper:  k=0 envelope-top (SHARED wake_upper.k_hi),
        #               k=last top frame (FREE)
        #   top_lower:  k=0 envelope-bot (SHARED wake_lower.k_hi),
        #               k=last bottom frame (FREE)
        #   wake_upper: k=0 wake-cut (SHARED to wake_lower),
        #               k=last envelope-top (SHARED top_upper.k_lo)
        #   wake_lower: k=0 wake-cut (SHARED to wake_upper),
        #               k=last envelope-bot (SHARED top_lower.k_lo)
        labels = [b.label for b in m.blocks]
        top_upper = m.blocks[labels.index("chan_outflow_top_upper")]
        top_lower = m.blocks[labels.index("chan_outflow_top_lower")]
        wake_upper = m.blocks[labels.index("chan_outflow_wake_upper")]
        wake_lower = m.blocks[labels.index("chan_outflow_wake_lower")]
        assert top_upper.face("k_hi").bc == FaceBC.FREE
        assert top_lower.face("k_hi").bc == FaceBC.FREE
        assert top_upper.face("k_lo").bc == FaceBC.SHARED
        assert top_lower.face("k_lo").bc == FaceBC.SHARED
        assert wake_upper.face("k_lo").bc == FaceBC.SHARED
        assert wake_lower.face("k_lo").bc == FaceBC.SHARED

    def test_chan_jlo_jhi_are_free(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, march=False)
        for chan in m.blocks[3:]:
            assert chan.face("j_lo").bc == FaceBC.FREE
            assert chan.face("j_hi").bc == FaceBC.FREE

    def test_chan_chan_lateral_seams_shared(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, march=False)
        chan_upper = m.blocks[3]
        chan_body_top = m.blocks[4]
        chan_body_front = m.blocks[5]
        chan_body_bot = m.blocks[6]
        chan_lower = m.blocks[7]
        assert chan_upper.face("i_hi").bc == FaceBC.SHARED
        assert chan_body_top.face("i_lo").bc == FaceBC.SHARED
        assert chan_body_top.face("i_hi").bc == FaceBC.SHARED
        assert chan_body_front.face("i_lo").bc == FaceBC.SHARED
        assert chan_body_front.face("i_hi").bc == FaceBC.SHARED
        assert chan_body_bot.face("i_lo").bc == FaceBC.SHARED
        assert chan_body_bot.face("i_hi").bc == FaceBC.SHARED
        assert chan_lower.face("i_lo").bc == FaceBC.SHARED

    def test_outflow_channels_conformal_to_wrap(self, naca0012_stl):
        # Phase-1.2 redesign: outflow split into top + wake sub-blocks
        # per side.  Wake sub-blocks are conformal SHARED to the wrap
        # upper.i_lo / lower.i_hi faces (same nk = nk_wrap, identity
        # AxisMap).  Top sub-blocks are conformal SHARED to the
        # k_hi-side chan_upper.i_lo / chan_lower.i_hi faces.
        m, _ = _build_ch_mesh(naca0012_stl, march=False)
        labels = [b.label for b in m.blocks]
        upper = m.blocks[0]
        lower = m.blocks[2]
        chan_upper = m.blocks[labels.index("chan_upper")]
        chan_lower = m.blocks[labels.index("chan_lower")]
        wake_upper = m.blocks[labels.index("chan_outflow_wake_upper")]
        wake_lower = m.blocks[labels.index("chan_outflow_wake_lower")]
        top_upper = m.blocks[labels.index("chan_outflow_top_upper")]
        top_lower = m.blocks[labels.index("chan_outflow_top_lower")]

        # Wrap-side i_lo/i_hi now SHARED, not FREE.
        assert upper.face("i_lo").bc == FaceBC.SHARED
        assert lower.face("i_hi").bc == FaceBC.SHARED
        # k_hi-side channel i_lo/i_hi at the wake far-tip end SHARED to
        # the top outflow sub-blocks.
        assert chan_upper.face("i_lo").bc == FaceBC.SHARED
        assert chan_lower.face("i_hi").bc == FaceBC.SHARED
        # Outflow sub-blocks' i_lo SHARED back.
        assert wake_upper.face("i_lo").bc == FaceBC.SHARED
        assert wake_lower.face("i_lo").bc == FaceBC.SHARED
        assert top_upper.face("i_lo").bc == FaceBC.SHARED
        assert top_lower.face("i_lo").bc == FaceBC.SHARED


# ---------------------------------------------------------------------------
# Frame resolution (chord/normal axes from freestream)
# ---------------------------------------------------------------------------


class TestCHFrame:
    """Freestream / span axis selection."""

    def test_default_frame_x_chord_z_span(self, naca0012_stl):
        _, topo = _build_ch_mesh(naca0012_stl, march=False)
        assert topo.frame.chord_axis == 0  # x
        assert topo.frame.normal_axis == 1  # y
        assert topo.frame.span_axis_idx == 2  # z

    def test_explicit_freestream_y_picks_y_chord(self):
        topo = CHTopology(
            chord_min=-1.0,
            chord_max=2.0,
            normal_min=-1.0,
            normal_max=1.0,
            nk_channel=3,
            freestream_world=np.array([0.0, 1.0, 0.0]),
            span_axis="z",
        )
        assert topo.frame.chord_axis == 1
        assert topo.frame.normal_axis == 0
        assert topo.frame.span_axis_idx == 2


# ---------------------------------------------------------------------------
# C1 seam pinning
# ---------------------------------------------------------------------------


class TestC1Seam:
    """First-channel-step matches wrap last-step at every (i, j)."""

    def test_c1_pin_matches_wrap_last_step(self, naca0012_stl):
        m, _ = _build_ch_mesh(naca0012_stl, c1_seam=True, march=True)
        upper = m.blocks[0]
        chan_upper = m.blocks[3]
        wrap_last_step = np.linalg.norm(
            upper.nodes[:, :, -1, :] - upper.nodes[:, :, -2, :], axis=-1
        )
        chan_first_step = np.linalg.norm(
            chan_upper.nodes[:, :, 1, :] - chan_upper.nodes[:, :, 0, :], axis=-1
        )
        np.testing.assert_allclose(chan_first_step, wrap_last_step, atol=1e-9)

    def test_no_c1_pin_falls_back_to_law(self, naca0012_stl):
        m_pin, _ = _build_ch_mesh(naca0012_stl, c1_seam=True, march=True)
        m_nopin, _ = _build_ch_mesh(naca0012_stl, c1_seam=False, march=True)
        # The two channels should have different geometry away from the seam
        # (since the c1 path uses geometric tail vs the law's distribution).
        a = m_pin.blocks[3].nodes
        b = m_nopin.blocks[3].nodes
        assert not np.allclose(a, b)


# ---------------------------------------------------------------------------
# Wrap geometry unaffected by CH wrapper (vs CGridTopology)
# ---------------------------------------------------------------------------


class TestMarchingUnaffectsWrap:
    """CH wrap geometry matches a standalone CGridTopology run, bit-exact."""

    def test_wrap_matches_cgrid_standalone(self, naca0012_stl):
        kw = {"ni_body": 80, "ni_wake": 20, "n_stations": 5, "wake_length": 15.0}
        body = Body.from_stl(naca0012_stl, skip_hygiene=True)

        # CH path
        topo_ch = CHTopology(**_CH_KW)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks_ch = topo_ch.build(surface=None, nk=10, mesh=body.mesh, **kw)
            m_ch = Mesh(body=body, blocks=blocks_ch, topology=topo_ch)
            m_ch.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)

        # Standalone CGrid path
        topo_cg = CGridTopology()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks_cg = topo_cg.build(surface=None, nk=10, mesh=body.mesh, **kw)
            m_cg = Mesh(body=body, blocks=blocks_cg, topology=topo_cg)
            m_cg.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)

        # Compare the 3 wrap blocks (upper / body / lower).  Channel views
        # share memory with the wrap, so they are inside the same buffer —
        # but the wrap-only k-range is what we compare here.
        for label in ("upper", "body", "lower"):
            ch_b = next(b for b in m_ch.blocks if b.label == label)
            cg_b = next(b for b in m_cg.blocks if b.label == label)
            np.testing.assert_array_equal(ch_b.nodes, cg_b.nodes)


# ---------------------------------------------------------------------------
# Channel filling produces sensible geometry
# ---------------------------------------------------------------------------


class TestChannelFilling:
    """Channels reach the frame and stay within it."""

    def test_chan_upper_outer_layer_on_top_frame(self, naca0012_stl):
        m, topo = _build_ch_mesh(naca0012_stl, march=True)
        chan_upper = m.blocks[3]
        # Outer face k=nk_channel-1 should sit at normal=normal_max.
        nax = topo.frame.normal_axis
        outer = chan_upper.nodes[:, :, -1, nax]
        np.testing.assert_allclose(outer, topo.frame.normal_max, atol=1e-9)

    def test_chan_body_front_outer_layer_on_front_frame(self, naca0012_stl):
        m, topo = _build_ch_mesh(naca0012_stl, march=True)
        chan_body_front = m.blocks[5]
        cax = topo.frame.chord_axis
        outer = chan_body_front.nodes[:, :, -1, cax]
        np.testing.assert_allclose(outer, topo.frame.chord_min, atol=1e-9)

    def test_chan_outflow_outer_layer_on_rear_frame(self, naca0012_stl):
        # In the rectangular outflow design, channel-i runs along chord,
        # so the rear frame edge is at i = -1.  Check all 4 outflow
        # sub-blocks land their i_hi face at chord = chord_max.
        m, topo = _build_ch_mesh(naca0012_stl, march=True)
        labels = [b.label for b in m.blocks]
        cax = topo.frame.chord_axis
        for blk_label in (
            "chan_outflow_top_upper",
            "chan_outflow_wake_upper",
            "chan_outflow_wake_lower",
            "chan_outflow_top_lower",
        ):
            blk = m.blocks[labels.index(blk_label)]
            outer = blk.nodes[-1, :, :, cax]
            np.testing.assert_allclose(
                outer, topo.frame.chord_max, atol=1e-9, err_msg=f"{blk_label} rear frame"
            )

    def test_envelope_inside_frame(self, naca0012_stl):
        m, topo = _build_ch_mesh(naca0012_stl, march=True)
        upper, body, lower = m.blocks[:3]
        cax = topo.frame.chord_axis
        nax = topo.frame.normal_axis
        for blk in (upper, body, lower):
            env = blk.nodes[:, :, -1, :]
            assert env[..., cax].min() >= topo.frame.chord_min - 1e-9
            assert env[..., cax].max() <= topo.frame.chord_max + 1e-9
            assert env[..., nax].min() >= topo.frame.normal_min - 1e-9
            assert env[..., nax].max() <= topo.frame.normal_max + 1e-9


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


class TestCHParameterValidation:
    """Constructor / build() parameter checks."""

    def test_nk_channel_below_2_rejected(self):
        with pytest.raises(ParameterError, match="nk_channel must be >= 2"):
            CHTopology(
                chord_min=-1.0,
                chord_max=1.0,
                normal_min=-1.0,
                normal_max=1.0,
                nk_channel=1,
            )

    def test_invalid_k_law_rejected(self):
        with pytest.raises(ParameterError, match="k_law_channel must be"):
            CHTopology(
                chord_min=-1.0,
                chord_max=1.0,
                normal_min=-1.0,
                normal_max=1.0,
                nk_channel=3,
                k_law_channel="bad",
            )

    def test_inverted_chord_extents_rejected(self):
        with pytest.raises(ParameterError, match="chord_max must exceed chord_min"):
            CHTopology(
                chord_min=2.0,
                chord_max=1.0,
                normal_min=-1.0,
                normal_max=1.0,
                nk_channel=3,
            )

    def test_inverted_normal_extents_rejected(self):
        with pytest.raises(ParameterError, match="normal_max must exceed normal_min"):
            CHTopology(
                chord_min=-1.0,
                chord_max=1.0,
                normal_min=1.0,
                normal_max=-1.0,
                nk_channel=3,
            )

    def test_zero_freestream_rejected(self):
        with pytest.raises(ParameterError, match="zero magnitude"):
            CHTopology(
                chord_min=-1.0,
                chord_max=1.0,
                normal_min=-1.0,
                normal_max=1.0,
                nk_channel=3,
                freestream_world=np.zeros(3),
            )

    def test_freestream_parallel_to_span_rejected(self):
        with pytest.raises(ParameterError, match="parallel to the span axis"):
            CHTopology(
                chord_min=-1.0,
                chord_max=1.0,
                normal_min=-1.0,
                normal_max=1.0,
                nk_channel=3,
                freestream_world=np.array([0.0, 0.0, 1.0]),
                span_axis="z",
            )

    def test_envelope_outside_frame_rejected(self, naca0012_stl):
        # Tiny frame: envelope can't fit.
        body = Body.from_stl(naca0012_stl, skip_hygiene=True)
        topo = CHTopology(
            chord_min=-0.05,
            chord_max=0.05,
            normal_min=-0.05,
            normal_max=0.05,
            nk_channel=3,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks = topo.build(
                surface=None,
                nk=4,
                mesh=body.mesh,
                ni_body=40,
                ni_wake=10,
                n_stations=4,
                wake_length=2.0,
            )
            m = Mesh(body=body, blocks=blocks, topology=topo)
            with pytest.raises(ParameterError, match="outside the frame"):
                m.march(ds=5e-4, growth=1.15, smooth_strength=0.0, smooth_iters=0)
