"""Per-edge spacing data model tests.

Covers:
- public ``pinned_distribution`` (lifted from surface/equator.py)
- ``Spacing1D(law='pinned')``
- ``Edge.first_step`` / ``Edge.last_step`` fields
- per-block k-spacing override through ``Mesh.march()``
- per-edge i/j on BoxMesh and AnnularCylinderMesh via TFI
- FlatCapsCylinder and OGridTopology rejection / no-op of per-edge requests
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.mesh import (
    _block_user_k_spacing,
    _resolve_k_steps,
    _spacings_equal,
)
from shore.mesh.spacing import Spacing1D
from shore.mesh.topology import CubedSphereTopology
from shore.volume.spacing import pinned_distribution

# ---------------------------------------------------------------------------
# pinned_distribution (public)
# ---------------------------------------------------------------------------


class TestPinnedDistribution:
    def test_uniform_interior_pins_endpoints_exactly(self):
        out = pinned_distribution(6, 0.0, 1.0, first_step=0.1, last_step=0.05)
        assert out[0] == 0.0
        assert out[-1] == 1.0
        assert out[1] - out[0] == pytest.approx(0.1, abs=1e-15)
        assert out[-1] - out[-2] == pytest.approx(0.05, abs=1e-15)

    def test_tanh_interior_pins_endpoints_exactly(self):
        out = pinned_distribution(
            8, 0.0, 1.0, first_step=0.1, last_step=0.05, spacing="tanh", beta=4.0
        )
        assert out[1] - out[0] == pytest.approx(0.1, abs=1e-15)
        assert out[-1] - out[-2] == pytest.approx(0.05, abs=1e-15)
        # Monotone non-decreasing
        diffs = np.diff(out)
        assert np.all(diffs > 0)

    def test_n_equals_2_endpoints_only(self):
        out = pinned_distribution(2, 0.0, 1.0, first_step=1.0, last_step=1.0)
        np.testing.assert_array_equal(out, [0.0, 1.0])

    def test_n_equals_3_collapses_to_pin_only(self):
        out = pinned_distribution(3, 0.0, 1.0, first_step=0.4, last_step=0.6)
        np.testing.assert_allclose(out, [0.0, 0.4, 1.0], atol=1e-15)

    def test_rejects_non_positive_steps(self):
        with pytest.raises(ParameterError, match="positive"):
            pinned_distribution(5, 0.0, 1.0, first_step=0.0, last_step=0.1)
        with pytest.raises(ParameterError, match="positive"):
            pinned_distribution(5, 0.0, 1.0, first_step=0.1, last_step=-0.05)

    def test_rejects_oversized_pins(self):
        with pytest.raises(ParameterError, match="exceeds"):
            pinned_distribution(5, 0.0, 1.0, first_step=0.6, last_step=0.5)

    def test_rejects_unknown_interior_law(self):
        with pytest.raises(ParameterError, match="Unknown interior"):
            pinned_distribution(5, 0.0, 1.0, first_step=0.1, last_step=0.1, spacing="bogus")

    def test_rejects_n_below_2(self):
        with pytest.raises(ParameterError, match=">= 2"):
            pinned_distribution(1, 0.0, 1.0, first_step=0.5, last_step=0.5)


# ---------------------------------------------------------------------------
# Spacing1D pinned law
# ---------------------------------------------------------------------------


class TestSpacing1DPinnedLaw:
    def test_build_delegates_to_pinned_distribution(self):
        sp = Spacing1D(law="pinned", first_step=0.1, last_step=0.05)
        out = sp.build(6, 0.0, 1.0)
        ref = pinned_distribution(6, 0.0, 1.0, first_step=0.1, last_step=0.05)
        np.testing.assert_array_equal(out, ref)

    def test_interior_law_propagates(self):
        sp = Spacing1D(law="pinned", first_step=0.1, last_step=0.05, interior="tanh", beta=4.0)
        out = sp.build(8, 0.0, 1.0)
        ref = pinned_distribution(
            8, 0.0, 1.0, first_step=0.1, last_step=0.05, spacing="tanh", beta=4.0
        )
        np.testing.assert_array_equal(out, ref)

    def test_missing_steps_raises(self):
        sp = Spacing1D(law="pinned")
        with pytest.raises(ParameterError, match="first_step and last_step"):
            sp.build(6, 0.0, 1.0)

    def test_existing_laws_unchanged(self):
        # Backwards-compat: the new fields are optional and the existing laws
        # ignore them.
        sp = Spacing1D(law="uniform", first_step=99.0, last_step=99.0)
        out = sp.build(5, 0.0, 1.0)
        np.testing.assert_allclose(out, [0.0, 0.25, 0.5, 0.75, 1.0])


# ---------------------------------------------------------------------------
# Edge dataclass: new fields
# ---------------------------------------------------------------------------


class TestEdgeFields:
    def test_default_steps_are_none(self):
        nodes = np.zeros((4, 4, 3, 3))
        b = HexBlock("t", nodes)
        for e in b.edges():
            assert e.first_step is None
            assert e.last_step is None

    def test_step_assignment_persists(self):
        nodes = np.zeros((4, 4, 3, 3))
        b = HexBlock("t", nodes)
        e = b.edge("j_lo_k_lo")
        e.first_step = 0.05
        e.last_step = 0.04
        assert b.edge("j_lo_k_lo").first_step == 0.05
        assert b.edge("j_lo_k_lo").last_step == 0.04

    def test_block_records_default_spacings(self):
        nodes = np.zeros((4, 4, 3, 3))
        sp_k = Spacing1D(law="geometric", ds=2e-3)
        b = HexBlock("t", nodes, spacing_k=sp_k)
        assert b._default_spacing_k is sp_k
        # All 4 k-edges share the default instance (identity, not equality).
        for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
            assert b.edge(name).spacing is sp_k


# ---------------------------------------------------------------------------
# Per-block k-spacing override resolver
# ---------------------------------------------------------------------------


def _make_block(nk: int = 5, sp_k: Spacing1D | None = None) -> HexBlock:
    nodes = np.zeros((4, 4, nk, 3))
    return HexBlock("t", nodes, spacing_k=sp_k)


class TestKSpacingResolver:
    def test_no_override_returns_none(self):
        blk = _make_block()
        assert _block_user_k_spacing(blk) is None

    def test_full_override_returns_user_spacing(self):
        blk = _make_block()
        new_sp = Spacing1D(law="geometric", ds=5e-4, growth=1.2)
        for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
            blk.edge(name).spacing = new_sp
        assert _block_user_k_spacing(blk) is new_sp

    def test_partial_override_warns_and_returns_none(self):
        blk = _make_block()
        new_sp = Spacing1D(law="geometric", ds=5e-4)
        # Override only 2 of 4 k-edges.
        blk.edge("i_lo_j_lo").spacing = new_sp
        blk.edge("i_hi_j_hi").spacing = new_sp
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = _block_user_k_spacing(blk)
        assert result is None
        assert any("not others" in str(item.message) for item in w)

    def test_inconsistent_full_override_warns_and_returns_none(self):
        blk = _make_block()
        sp_a = Spacing1D(law="geometric", ds=5e-4)
        sp_b = Spacing1D(law="geometric", ds=1e-3)
        blk.edge("i_lo_j_lo").spacing = sp_a
        blk.edge("i_hi_j_lo").spacing = sp_a
        blk.edge("i_lo_j_hi").spacing = sp_b
        blk.edge("i_hi_j_hi").spacing = sp_b
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = _block_user_k_spacing(blk)
        assert result is None
        assert any("4 k-edges carry different" in str(item.message) for item in w)

    def test_resolve_falls_back_to_ds_growth_without_override(self):
        b1 = _make_block()
        b2 = _make_block()
        steps = _resolve_k_steps([b1, b2], nk=5, ds=0.1, growth=1.0, steps=None)
        # geometric law with growth=1.0 → constant ds across all (nk-1) layers
        np.testing.assert_allclose(steps, [0.1, 0.1, 0.1, 0.1])

    def test_resolve_uses_user_spacing_when_overridden(self):
        b = _make_block()
        new_sp = Spacing1D(law="geometric", ds=0.01, growth=1.5)
        for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
            b.edge(name).spacing = new_sp
        steps = _resolve_k_steps([b], nk=5, ds=99.0, growth=99.0, steps=None)
        np.testing.assert_allclose(steps, [0.01, 0.015, 0.0225, 0.03375])

    def test_resolve_explicit_steps_override_everything(self):
        b = _make_block()
        explicit = np.array([0.5, 0.5, 0.5, 0.5])
        out = _resolve_k_steps([b], nk=5, ds=0.001, growth=1.0, steps=explicit)
        np.testing.assert_array_equal(out, explicit)

    def test_resolve_warns_on_cross_block_disagreement(self):
        b1 = _make_block()
        b2 = _make_block()
        sp_a = Spacing1D(law="geometric", ds=0.01, growth=1.0)
        sp_b = Spacing1D(law="geometric", ds=0.02, growth=1.0)
        for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
            b1.edge(name).spacing = sp_a
            b2.edge(name).spacing = sp_b
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            steps = _resolve_k_steps([b1, b2], nk=5, ds=99.0, growth=99.0, steps=None)
        # Block 0 wins.
        np.testing.assert_allclose(steps, [0.01, 0.01, 0.01, 0.01])
        assert any("different spacing" in str(item.message) for item in w)


class TestSpacingsEqual:
    def test_none_handling(self):
        assert _spacings_equal(None, None) is True
        assert _spacings_equal(Spacing1D(law="uniform"), None) is False

    def test_value_equality(self):
        a = Spacing1D(law="geometric", ds=0.01, growth=1.1)
        b = Spacing1D(law="geometric", ds=0.01, growth=1.1)
        assert a is not b
        assert _spacings_equal(a, b) is True

    def test_value_difference(self):
        a = Spacing1D(law="geometric", ds=0.01)
        b = Spacing1D(law="geometric", ds=0.02)
        assert _spacings_equal(a, b) is False


# ---------------------------------------------------------------------------
# End-to-end: per-edge k-spacing flows through Mesh.march()
# ---------------------------------------------------------------------------


class TestMarchHonoursPerEdgeKSpacing:
    """Override k-spacing on block 0's k-edges and verify the marcher uses
    it instead of the (ds, growth) defaults from .march()."""

    def test_override_changes_layer_thickness(self, sphere_stl):
        from shore.mesh.mesh import Mesh

        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=9,
            nj=16,
            nk=5,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )

        # Override every block's k-edges with a *larger* ds — chosen so the
        # marcher's stability is independent of the call-site ds/growth.
        override = Spacing1D(law="geometric", ds=0.08, growth=1.1)
        for blk in m.blocks:
            for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
                blk.edge(name).spacing = override

        # Pass a different ds/growth that, if (incorrectly) consumed, would
        # produce a clearly different first layer thickness.
        m.march(ds=0.02, growth=1.1, smooth_strength=0.2, smooth_iters=2)

        blk = m.blocks[0]
        # First-layer thickness (radial) on the wall — should match
        # override.ds (~0.08) and NOT the .march() ds (0.02).
        first_layer_thickness = float(np.linalg.norm(blk.nodes[0, 0, 1, :] - blk.nodes[0, 0, 0, :]))
        assert 0.06 < first_layer_thickness < 0.10

    def test_unset_falls_back_to_ds_growth(self, sphere_stl):
        from shore.mesh.mesh import Mesh

        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=9,
            nj=16,
            nk=5,
        )
        # No edge overrides; .march() supplies ds=0.05.
        m.march(ds=0.05, growth=1.1, smooth_strength=0.2, smooth_iters=2)

        blk = m.blocks[0]
        first_layer_thickness = float(np.linalg.norm(blk.nodes[0, 0, 1, :] - blk.nodes[0, 0, 0, :]))
        # ds=0.05 expected; smoothing perturbs ~few percent.
        assert 0.03 < first_layer_thickness < 0.07


# ---------------------------------------------------------------------------
# Per-edge i/j on BoxMesh via transfinite interpolation
# ---------------------------------------------------------------------------


class TestBoxMeshPerEdge:
    def test_no_pins_reduces_to_tensor_product(self):
        from shore.primitive import BoxMesh

        b1 = BoxMesh.build(vmin=(0, 0, 0), vmax=(1, 1, 1), ni=11, nj=5, nk=5)
        b2 = BoxMesh.build(
            vmin=(0, 0, 0),
            vmax=(1, 1, 1),
            ni=11,
            nj=5,
            nk=5,
            edge_pins={},  # empty pin dict — same as None
        )
        np.testing.assert_array_equal(b1.blocks[0].nodes, b2.blocks[0].nodes)

    def test_single_i_edge_pin_blends_through_interior(self):
        from shore.primitive import BoxMesh

        # Pin j_lo_k_lo (i-edge at j=0, k=0) to tanh; other 3 i-edges remain
        # at the default uniform i_law.
        box = BoxMesh.build(
            vmin=(0, 0, 0),
            vmax=(1, 1, 1),
            ni=11,
            nj=5,
            nk=5,
            edge_pins={"j_lo_k_lo": Spacing1D(law="tanh", beta=4.0)},
        )
        nodes = box.blocks[0].nodes
        # Edge (j=0, k=0) — exact tanh.
        x_pinned = nodes[:, 0, 0, 0]
        # Opposite corner (j=nj-1, k=nk-1) — exact uniform.
        x_uniform = nodes[:, -1, -1, 0]
        # Endpoints match in both:
        assert x_pinned[0] == 0.0 and x_pinned[-1] == 1.0
        assert x_uniform[0] == 0.0 and x_uniform[-1] == 1.0
        # Tanh clusters near 0: x_pinned[1] should be much smaller than 0.1.
        assert x_pinned[1] < 0.05
        # Uniform: x_uniform[1] == 0.1 within FP noise.
        assert abs(x_uniform[1] - 0.1) < 1e-12
        # Midpoint corner (j=2, k=2) is a 50/50 blend → between the two extremes.
        x_mid = nodes[:, 2, 2, 0]
        assert x_pinned[1] < x_mid[1] < x_uniform[1]

    def test_edge_spec_stamped_onto_block(self):
        from shore.primitive import BoxMesh

        spec = Spacing1D(law="tanh", beta=4.0)
        box = BoxMesh.build(
            vmin=(0, 0, 0),
            vmax=(1, 1, 1),
            ni=5,
            nj=5,
            nk=5,
            edge_pins={"j_lo_k_lo": spec},
        )
        assert box.blocks[0].edge("j_lo_k_lo").spacing is spec

    def test_pinned_first_step(self):
        """Spacing1D(law='pinned', first_step=...) on an edge produces a
        block whose first interval along that edge equals the requested
        first_step."""
        from shore.primitive import BoxMesh

        box = BoxMesh.build(
            vmin=(0, 0, 0),
            vmax=(1, 1, 1),
            ni=10,
            nj=5,
            nk=5,
            edge_pins={
                "j_lo_k_lo": Spacing1D(
                    law="pinned",
                    first_step=0.05,
                    last_step=0.05,
                    interior="uniform",
                )
            },
        )
        x_edge = box.blocks[0].nodes[:, 0, 0, 0]
        assert abs((x_edge[1] - x_edge[0]) - 0.05) < 1e-12
        assert abs((x_edge[-1] - x_edge[-2]) - 0.05) < 1e-12


# ---------------------------------------------------------------------------
# Per-edge i/j on AnnularCylinderMesh
# ---------------------------------------------------------------------------


class TestAnnulusPerEdge:
    def test_periodic_full_annulus_rejects_j_pins(self):
        from shore.primitive import AnnularCylinderMesh

        with pytest.raises(ParameterError, match="j-edges"):
            AnnularCylinderMesh.build(
                r_in=0.5,
                r_out=2.0,
                z0=0.0,
                z1=4.0,
                ni=10,
                nj=64,
                nk=20,
                edge_pins={"i_lo_k_lo": Spacing1D(law="tanh", beta=4.0)},
            )

    def test_partial_sector_accepts_j_pin(self):
        from shore.primitive import AnnularCylinderMesh

        ann = AnnularCylinderMesh.build(
            r_in=0.5,
            r_out=2.0,
            z0=0.0,
            z1=4.0,
            ni=10,
            nj=20,
            nk=10,
            theta_start=0.0,
            theta_end=np.pi,  # half annulus → non-periodic j
            edge_pins={"i_lo_k_lo": Spacing1D(law="tanh", beta=4.0)},
        )
        assert ann.blocks[0].periodic_j is False
        assert ann.blocks[0].edge("i_lo_k_lo").spacing.law == "tanh"

    def test_radial_edge_pin_blends_through_interior(self):
        from shore.primitive import AnnularCylinderMesh

        ann = AnnularCylinderMesh.build(
            r_in=0.5,
            r_out=2.0,
            z0=0.0,
            z1=4.0,
            ni=11,
            nj=20,
            nk=5,
            theta_start=0.0,
            theta_end=np.pi,
            edge_pins={"j_lo_k_lo": Spacing1D(law="tanh", beta=4.0)},
        )
        nodes = ann.blocks[0].nodes
        # Pinned edge: i-axis at j=0, k=0.  Compute radius along i.
        r_pinned = np.linalg.norm(nodes[:, 0, 0, :2], axis=-1)
        # Endpoints match r_in and r_out in both cases.
        assert abs(r_pinned[0] - 0.5) < 1e-9
        assert abs(r_pinned[-1] - 2.0) < 1e-9
        # Uniform-i would give r[1] = 0.5 + (2.0 - 0.5)/10 = 0.65; tanh
        # clusters near r_in so r[1] should be much closer to 0.5.
        assert r_pinned[1] < 0.6


# ---------------------------------------------------------------------------
# FlatCapsCylinder rejects edge_pins
# ---------------------------------------------------------------------------


class TestFlatCapsRejectsEdgePins:
    def test_any_pin_is_rejected(self):
        from shore.primitive import FlatCapsCylinder

        with pytest.raises(ParameterError, match="not supported"):
            FlatCapsCylinder.build(
                r_out=2.0,
                z0=0.0,
                z1=4.0,
                ni=8,
                nc=10,
                nk=8,
                edge_pins={"i_lo_j_lo": Spacing1D(law="tanh", beta=4.0)},
            )

    def test_unset_edge_pins_still_works(self):
        from shore.primitive import FlatCapsCylinder

        fc = FlatCapsCylinder.build(
            r_out=2.0,
            z0=0.0,
            z1=4.0,
            ni=8,
            nc=10,
            nk=8,
            edge_pins=None,
        )
        assert len(fc.blocks) == 5


# ---------------------------------------------------------------------------
# OGridTopology relies on the existing Mesh.march resolver for per-edge k
# ---------------------------------------------------------------------------


class TestOGridPerEdgeKThroughMarch:
    def test_per_edge_k_override_changes_layer_thickness(self, sphere_stl):
        from shore.mesh.mesh import Mesh
        from shore.mesh.topology import OGridTopology

        m = Mesh.from_stl(
            sphere_stl,
            topology=OGridTopology(),
            ni=10,
            nj=20,
            nk=5,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        override = Spacing1D(law="geometric", ds=0.08, growth=1.1)
        for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
            m.blocks[0].edge(name).spacing = override
        m.march(ds=0.02, growth=1.1, smooth_strength=0.2, smooth_iters=2)
        blk = m.blocks[0]
        # First layer at the south pole row — should reflect override (~0.08),
        # not the call-site ds (0.02).  Laplacian smoothing perturbs the
        # raw step somewhat; we just need the result to be unambiguously
        # closer to the override than to the call-site value.
        first_layer_thickness = float(
            np.linalg.norm(blk.nodes[5, 10, 1, :] - blk.nodes[5, 10, 0, :])
        )
        assert 0.04 < first_layer_thickness < 0.10
        # And clearly NOT the call-site value (0.02).
        assert first_layer_thickness > 2.0 * 0.02
