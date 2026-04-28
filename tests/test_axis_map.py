"""Tests for the AxisMap data model and its enforcement on Face / HexBlock.

The AxisMap field on every connection face (SHARED / DIRICHLET /
PERIODIC) is what lets cc.par-style consumers compute orientation
codes without per-topology lookup tables.  These tests cover:

- AxisMap validation (axis labels, flip flags, partner-axis consistency)
- AxisMap.identity / from_pairs / inverse helpers
- _set_face_bc requiring axis_map for connection BCs and rejecting
  axis_map for non-connection BCs
"""

from __future__ import annotations

import numpy as np
import pytest

from shore.mesh.block import HexBlock
from shore.mesh.face import FACE_ALONG_AXES, AxisMap, FaceBC


def _two_blocks() -> tuple[HexBlock, HexBlock]:
    """Two adjacent unit-cube blocks, geometry irrelevant for face-graph tests."""
    nodes_a = np.zeros((3, 3, 3, 3))
    nodes_b = np.zeros((3, 3, 3, 3))
    return (
        HexBlock(label="a", nodes=nodes_a, periodic_i=False, periodic_j=False),
        HexBlock(label="b", nodes=nodes_b, periodic_i=False, periodic_j=False),
    )


class TestAxisMapValidation:
    def test_along_axes_must_differ(self):
        with pytest.raises(ValueError, match="must differ"):
            AxisMap.from_pairs(("i", False), ("i", True))

    def test_axis_label_must_be_ijk(self):
        with pytest.raises(ValueError, match="invalid axis label"):
            AxisMap.from_pairs(("x", False), ("k", False))

    def test_flip_flag_must_be_bool(self):
        with pytest.raises(TypeError, match="must be bool"):
            AxisMap.from_pairs(("i", 1), ("k", False))  # type: ignore[arg-type]

    def test_frozen(self):
        amap = AxisMap.from_pairs(("i", False), ("k", False))
        # Frozen dataclass rejects mutation; the concrete exception is
        # FrozenInstanceError but that's a dataclasses internal — match
        # on either by accepting AttributeError too.
        with pytest.raises((AttributeError, TypeError)):
            amap.along = ()  # type: ignore[misc]


class TestAxisMapIdentity:
    @pytest.mark.parametrize(
        ("face_name", "expected"),
        [
            ("i_lo", (("j", False), ("k", False))),
            ("i_hi", (("j", False), ("k", False))),
            ("j_lo", (("i", False), ("k", False))),
            ("j_hi", (("i", False), ("k", False))),
            ("k_lo", (("i", False), ("j", False))),
            ("k_hi", (("i", False), ("j", False))),
        ],
    )
    def test_identity_matches_face_along_axes(self, face_name, expected):
        amap = AxisMap.identity(face_name)
        assert amap.along == expected

    def test_identity_partner_axes_match_canonical(self):
        for face_name in FACE_ALONG_AXES:
            amap = AxisMap.identity(face_name)
            assert amap.partner_axes() == FACE_ALONG_AXES[face_name]


class TestAxisMapInverse:
    def test_same_axis_no_flip_self_inverse(self):
        # j_hi <-> j_lo, identity on each side.
        amap = AxisMap.identity("j_hi")
        inv = amap.inverse("j_hi", "j_lo")
        # Both sides describe (i, k) -> (i, k) with no flips: same payload.
        assert inv.along == amap.along

    def test_double_inverse_is_identity(self):
        amap = AxisMap.from_pairs(("i", True), ("k", False))
        # Roundtrip: inverse then inverse-back returns the original.
        inv = amap.inverse("i_hi", "j_hi")
        back = inv.inverse("j_hi", "i_hi")
        assert back.along == amap.along

    def test_cross_axis_inverse_swaps_labels(self):
        # sub.i_hi -> cap.j_hi, sub.j -> cap.i, sub.k -> cap.k.
        amap = AxisMap.from_pairs(("i", False), ("k", False))
        inv = amap.inverse("i_hi", "j_hi")
        # cap.i (cap_along[0]) -> sub.j (no flip).
        # cap.k (cap_along[1]) -> sub.k.
        assert inv.along == (("j", False), ("k", False))

    def test_flip_preserved_through_inverse(self):
        amap = AxisMap.from_pairs(("j", True), ("k", False))
        inv = amap.inverse("i_hi", "i_hi")
        # Same-axis with j-flip: inverse is also j-flipped.
        assert inv.along == (("j", True), ("k", False))


class TestAxisMapValidateAgainstPartner:
    def test_accepts_consistent_partner(self):
        amap = AxisMap.from_pairs(("i", False), ("k", False))
        amap.validate_against_partner("i_hi", "j_hi")  # i_hi along=(j,k); cap=(i,k)

    def test_rejects_partner_axis_not_along(self):
        # Claim sub.j -> cap.k as second along — but cap.j_hi's along is (i, k),
        # so saying both pairs use ("j", ...) is invalid (j is the perpendicular
        # axis for j_hi).
        amap = AxisMap.from_pairs(("j", False), ("i", False))
        with pytest.raises(ValueError, match="partner axes"):
            amap.validate_against_partner("i_hi", "j_hi")


class TestSetFaceBCEnforcement:
    def test_shared_requires_partner(self):
        a, _b = _two_blocks()
        with pytest.raises(ValueError, match="requires a partner"):
            a._set_face_bc("j_hi", FaceBC.SHARED, axis_map=AxisMap.identity("j_hi"))

    def test_shared_requires_axis_map(self):
        a, b = _two_blocks()
        with pytest.raises(ValueError, match="requires an axis_map"):
            a._set_face_bc("j_hi", FaceBC.SHARED, partner=b.face("j_lo"))

    def test_dirichlet_requires_axis_map(self):
        a, b = _two_blocks()
        with pytest.raises(ValueError, match="requires an axis_map"):
            a._set_face_bc("i_hi", FaceBC.DIRICHLET, partner=b.face("k_lo"))

    def test_periodic_requires_axis_map(self):
        a, _b = _two_blocks()
        with pytest.raises(ValueError, match="requires an axis_map"):
            a._set_face_bc("j_hi", FaceBC.PERIODIC, partner=a.face("j_lo"))

    def test_wall_must_not_have_partner(self):
        a, b = _two_blocks()
        with pytest.raises(ValueError, match="must not have a partner"):
            a._set_face_bc("k_lo", FaceBC.WALL, partner=b.face("k_lo"))

    def test_free_must_not_have_axis_map(self):
        a, _b = _two_blocks()
        with pytest.raises(ValueError, match="must not have an axis_map"):
            a._set_face_bc("i_lo", FaceBC.FREE, axis_map=AxisMap.identity("i_lo"))

    def test_axis_map_validated_against_partner_face(self):
        a, b = _two_blocks()
        # sub.i_hi -> partner b.j_hi: along (j, k) on sub side.  Claiming
        # both along axes are partner-i is invalid (only one slot can map
        # to partner-i; b.j_hi's along axes are (i, k), so the AxisMap
        # must use one i and one k slot).
        bad_axis_map = AxisMap.from_pairs(("i", False), ("j", False))
        with pytest.raises(ValueError, match="partner axes"):
            a._set_face_bc(
                "i_hi",
                FaceBC.SHARED,
                partner=b.face("j_hi"),
                axis_map=bad_axis_map,
            )

    def test_set_face_bc_succeeds_with_consistent_axis_map(self):
        a, b = _two_blocks()
        a._set_face_bc(
            "j_hi",
            FaceBC.SHARED,
            partner=b.face("j_lo"),
            axis_map=AxisMap.identity("j_hi"),
        )
        face = a.face("j_hi")
        assert face.bc == FaceBC.SHARED
        assert face.partner is b.face("j_lo")
        assert face.axis_map == AxisMap.identity("j_hi")
