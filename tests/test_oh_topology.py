"""Tests for shore.mesh.OHTopology — Style-2 conforming hybrid (cubed-sphere + H-channel).

The OH topology produces 12 blocks:
    - 6 wrap blocks (cubed-sphere, k=0 on STL, marched outward to nk_wrap-1)
    - 6 channel blocks (frustum, k=0 ≡ wrap.k_hi, k=nk_c-1 on the box)

Wrap-channel SHARED-by-view: the wrap and its mirror channel share one
underlying buffer of depth (nk_wrap + nk_channel - 1).  The channel's
k=0 row IS the wrap's k=nk_wrap-1 row by numpy view aliasing.
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh.face import FaceBC
from shore.mesh.mesh import Mesh
from shore.mesh.spacing import Spacing1D
from shore.mesh.topology import OHTopology

# Common build parameters used across the test classes.
NK_WRAP = 5
NK_CHANNEL = 4
NI = 9
NJ = 16  # divisible by 4 — cubed-sphere requirement
BOX_VMIN = (-3.0, -3.0, -3.0)
BOX_VMAX = (3.0, 3.0, 3.0)


# ---------------------------------------------------------------------------
# Construction smoke tests (no marching)
# ---------------------------------------------------------------------------


class TestOHConstruction:
    def test_returns_12_blocks(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(
                box_vmin=BOX_VMIN,
                box_vmax=BOX_VMAX,
                nk_channel=NK_CHANNEL,
            ),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        assert len(m.blocks) == 12

    def test_block_labels(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        labels = [b.label for b in m.blocks]
        assert labels == [
            "sub0",
            "sub1",
            "sub2",
            "sub3",
            "cap_north",
            "cap_south",
            "chan_sub0",
            "chan_sub1",
            "chan_sub2",
            "chan_sub3",
            "chan_north",
            "chan_south",
        ]

    def test_block_shapes_match_pairs(self, sphere_stl):
        """Each channel block has (ni, nj) of its paired wrap block."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        for q in range(6):
            wrap = m.blocks[q]
            chan = m.blocks[6 + q]
            assert wrap.nodes.shape[:2] == chan.nodes.shape[:2]
            assert wrap.nodes.shape[2] == NK_WRAP
            assert chan.nodes.shape[2] == NK_CHANNEL


# ---------------------------------------------------------------------------
# View-aliasing: wrap.k_hi shares memory with channel.k_lo
# ---------------------------------------------------------------------------


class TestSharedViewAliasing:
    def test_buffers_share_memory(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        for q in range(6):
            wrap = m.blocks[q]
            chan = m.blocks[6 + q]
            # Mutating the wrap's k_hi slice must show up in the channel's k_lo.
            wrap.nodes[:, :, -1, 0] = 99.0
            np.testing.assert_array_equal(chan.nodes[:, :, 0, 0], 99.0)
            # And vice versa.
            chan.nodes[:, :, 0, 1] = 42.0
            np.testing.assert_array_equal(wrap.nodes[:, :, -1, 1], 42.0)

    def test_post_march_seam_is_bit_exact(self, sphere_stl):
        """After marching, channel.k_lo == wrap.k_hi node-for-node."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        m.march(ds=0.05, growth=1.1)
        for q in range(6):
            wrap = m.blocks[q]
            chan = m.blocks[6 + q]
            np.testing.assert_array_equal(wrap.nodes[:, :, -1, :], chan.nodes[:, :, 0, :])


# ---------------------------------------------------------------------------
# BC layout
# ---------------------------------------------------------------------------


class TestOHBCs:
    def test_wrap_k_hi_is_shared_with_channel(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        for q in range(6):
            wrap = m.blocks[q]
            chan = m.blocks[6 + q]
            f_w = wrap.face("k_hi")
            f_c = chan.face("k_lo")
            assert f_w.bc == FaceBC.SHARED
            assert f_c.bc == FaceBC.SHARED
            assert f_w.partner is f_c
            assert f_c.partner is f_w

    def test_channel_k_hi_is_free(self, sphere_stl):
        """Outer channel face = box → FREE (far-field)."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        for q in range(6, 12):
            assert m.blocks[q].face("k_hi").bc == FaceBC.FREE

    def test_wrap_k_lo_still_wall(self, sphere_stl):
        """Wrap k_lo (STL surface) BC is unchanged from CubedSphereTopology."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        for q in range(6):
            assert m.blocks[q].face("k_lo").bc == FaceBC.WALL

    def test_channel_lateral_seams_match_wrap_topology(self, sphere_stl):
        """Channel mirrors wrap's lateral seam structure: chan_sub_q.j_hi ↔ chan_sub_(q+1).j_lo."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        chans = m.blocks[6:10]  # 4 equatorial channel sub-blocks
        for q in range(4):
            q_next = (q + 1) % 4
            f_hi = chans[q].face("j_hi")
            f_lo = chans[q_next].face("j_lo")
            assert f_hi.bc == FaceBC.SHARED
            assert f_lo.bc == FaceBC.SHARED
            assert f_hi.partner is f_lo
            assert f_lo.partner is f_hi


# ---------------------------------------------------------------------------
# Channel construction (post-march)
# ---------------------------------------------------------------------------


class TestChannelConstruction:
    def test_channel_outer_face_lies_on_box(self, sphere_stl):
        """After build, channel.k_hi nodes lie on the corresponding box face."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        m.march(ds=0.05, growth=1.1)

        # Each channel block's outer face lies on exactly one of the 6 box
        # planes — i.e. one of (x, y, z) is constant at one of (vmin, vmax).
        x_planes = (BOX_VMIN[0], BOX_VMAX[0])
        y_planes = (BOX_VMIN[1], BOX_VMAX[1])
        z_planes = (BOX_VMIN[2], BOX_VMAX[2])
        tol = 1e-9
        for q in range(6, 12):
            outer = m.blocks[q].nodes[:, :, -1, :]  # (ni, nj, 3)
            on_x = np.any([np.allclose(outer[..., 0], v, atol=tol) for v in x_planes])
            on_y = np.any([np.allclose(outer[..., 1], v, atol=tol) for v in y_planes])
            on_z = np.any([np.allclose(outer[..., 2], v, atol=tol) for v in z_planes])
            assert on_x or on_y or on_z, (
                f"Channel block {m.blocks[q].label} outer face not on any "
                f"box face (extents x={outer[..., 0].min():.3f}/{outer[..., 0].max():.3f}, "
                f"y={outer[..., 1].min():.3f}/{outer[..., 1].max():.3f}, "
                f"z={outer[..., 2].min():.3f}/{outer[..., 2].max():.3f})"
            )

    def test_six_distinct_box_faces_covered(self, sphere_stl):
        """Each of the 6 box faces is covered by exactly one channel block."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        m.march(ds=0.05, growth=1.1)

        # For each channel block, identify which of the 6 box faces it lies on.
        faces_hit = []
        tol = 1e-9
        for q in range(6, 12):
            outer = m.blocks[q].nodes[:, :, -1, :]
            cx = outer[..., 0].mean()
            cy = outer[..., 1].mean()
            cz = outer[..., 2].mean()
            # Find the dominant axis.
            for axis, (ax_min, ax_max), _c in (
                ("x", (BOX_VMIN[0], BOX_VMAX[0]), cx),
                ("y", (BOX_VMIN[1], BOX_VMAX[1]), cy),
                ("z", (BOX_VMIN[2], BOX_VMAX[2]), cz),
            ):
                if np.allclose(outer[..., {"x": 0, "y": 1, "z": 2}[axis]], ax_min, atol=tol):
                    faces_hit.append(f"{axis}-")
                    break
                if np.allclose(outer[..., {"x": 0, "y": 1, "z": 2}[axis]], ax_max, atol=tol):
                    faces_hit.append(f"{axis}+")
                    break
        assert sorted(faces_hit) == ["x+", "x-", "y+", "y-", "z+", "z-"]

    def test_north_cap_owns_zminus_face(self, sphere_stl):
        """cap_north (z<0 side per shore convention) → -z box face."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        m.march(ds=0.05, growth=1.1)

        chan_north = m.blocks[10]  # paired with cap_north
        outer = chan_north.nodes[:, :, -1, :]
        np.testing.assert_allclose(outer[..., 2], BOX_VMIN[2], atol=1e-9)

    def test_south_cap_owns_zplus_face(self, sphere_stl):
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
        )
        m.march(ds=0.05, growth=1.1)

        chan_south = m.blocks[11]  # paired with cap_south
        outer = chan_south.nodes[:, :, -1, :]
        np.testing.assert_allclose(outer[..., 2], BOX_VMAX[2], atol=1e-9)


# ---------------------------------------------------------------------------
# Marching does not regress the wrap (the channel must not interfere)
# ---------------------------------------------------------------------------


class TestMarchingUnaffectsWrap:
    def test_wrap_radius_grows(self, sphere_stl):
        """Wrap blocks still expand outward like vanilla CubedSphereTopology."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        m.march(ds=0.05, growth=1.1)
        for q in range(6):
            blk = m.blocks[q]
            r0 = np.linalg.norm(blk.nodes[:, :, 0, :], axis=-1).mean()
            rN = np.linalg.norm(blk.nodes[:, :, -1, :], axis=-1).mean()
            assert rN > r0 * 1.05


# ---------------------------------------------------------------------------
# C1 seam: wrap's last marching step matches channel's first step
# ---------------------------------------------------------------------------


class TestC1Seam:
    def test_seam_step_match_default(self, sphere_stl):
        """Default c1_seam=True: ||chan[k=1]-chan[k=0]|| ≈ ||wrap[-1]-wrap[-2]|| at every (i, j)."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=NK_CHANNEL),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        m.march(ds=0.05, growth=1.1)
        for q in range(6):
            wrap = m.blocks[q]
            chan = m.blocks[6 + q]
            wrap_last = np.linalg.norm(wrap.nodes[:, :, -1, :] - wrap.nodes[:, :, -2, :], axis=-1)
            chan_first = np.linalg.norm(chan.nodes[:, :, 1, :] - chan.nodes[:, :, 0, :], axis=-1)
            np.testing.assert_allclose(chan_first, wrap_last, rtol=1e-8)

    def test_channel_geometric_growth_smooth(self, sphere_stl):
        """With c1_seam, channel cells should grow geometrically (constant ratio)."""
        # Use a denser channel so the ratio is well-defined.
        nk_channel = 8
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=nk_channel),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        m.march(ds=0.05, growth=1.1)
        # On any wrap-channel pair, take the channel column at the centre
        # (i, j) and verify cell-to-cell ratio is constant — i.e. a single
        # smooth geometric law from the seam out to the box face.
        for q in range(6):
            chan = m.blocks[6 + q]
            ic, jc = chan.nodes.shape[0] // 2, chan.nodes.shape[1] // 2
            steps = np.linalg.norm(np.diff(chan.nodes[ic, jc, :, :], axis=0), axis=-1)
            ratios = steps[1:] / steps[:-1]
            # All ratios within 0.5% of each other — geometric growth.
            assert ratios.std() / ratios.mean() < 5e-3, (
                f"chan_{chan.label}: cell ratios {ratios} not constant — "
                f"std/mean = {ratios.std() / ratios.mean():.4f}"
            )

    def test_c1_seam_off_makes_jump(self, sphere_stl):
        """c1_seam=False: spacing law dominates, no per-(i, j) match."""
        m = Mesh.from_stl(
            sphere_stl,
            topology=OHTopology(
                box_vmin=BOX_VMIN,
                box_vmax=BOX_VMAX,
                nk_channel=NK_CHANNEL,
                c1_seam=False,
            ),
            ni=NI,
            nj=NJ,
            nk=NK_WRAP,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        m.march(ds=0.05, growth=1.1)
        # On the equatorial sub-blocks, the unpinned channel divides the
        # remaining gap into nk_channel-1 ≈ uniform cells, which is much
        # bigger than the wrap's last step (~0.07).  The seam must show a
        # visible jump — assert at least 2x mismatch on sub0.
        wrap = m.blocks[0]
        chan = m.blocks[6]
        wrap_last = np.linalg.norm(
            wrap.nodes[:, :, -1, :] - wrap.nodes[:, :, -2, :], axis=-1
        ).mean()
        chan_first = np.linalg.norm(chan.nodes[:, :, 1, :] - chan.nodes[:, :, 0, :], axis=-1).mean()
        assert chan_first > 2.0 * wrap_last


# ---------------------------------------------------------------------------
# Parameter validation
# ---------------------------------------------------------------------------


class TestOHParameterValidation:
    def test_box_vmin_ge_vmax_rejected(self):
        with pytest.raises(ParameterError, match="box_vmin"):
            OHTopology(
                box_vmin=(1.0, -3.0, -3.0),
                box_vmax=(1.0, 3.0, 3.0),
                nk_channel=4,
            )

    def test_nk_channel_too_small_rejected(self):
        with pytest.raises(ParameterError, match="nk_channel"):
            OHTopology(box_vmin=BOX_VMIN, box_vmax=BOX_VMAX, nk_channel=1)

    def test_invalid_k_law_rejected(self):
        with pytest.raises(ParameterError, match="k_law"):
            OHTopology(
                box_vmin=BOX_VMIN,
                box_vmax=BOX_VMAX,
                nk_channel=4,
                k_law_channel="bogus",
            )

    def test_box_must_contain_sphere(self, sphere_stl):
        """A box too small to contain the wrap envelope should be rejected."""
        with pytest.raises(ParameterError, match="box"):
            Mesh.from_stl(
                sphere_stl,
                topology=OHTopology(
                    # Sphere has radius 1; this box is too small to contain
                    # even k=0, let alone the marched envelope.
                    box_vmin=(-0.5, -0.5, -0.5),
                    box_vmax=(0.5, 0.5, 0.5),
                    nk_channel=NK_CHANNEL,
                ),
                ni=NI,
                nj=NJ,
                nk=NK_WRAP,
            ).march(ds=0.05, growth=1.1)
