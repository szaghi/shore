"""Tests for shore.surface.cgrid — C-grid surface construction."""

from __future__ import annotations

import warnings

import numpy as np
import pytest
import trimesh

from shore.exceptions import MeshInputError, ParameterError
from shore.mesh import CGridTopology, FaceBC
from shore.surface.cgrid import (
    build_c_curve,
    build_cgrid_k0,
    detect_te_edges,
    round_body_te_in_surface,
    slice_by_span,
)

# ---------------------------------------------------------------------------
# detect_te_edges — sharp TE
# ---------------------------------------------------------------------------


class TestDetectTESharp:
    def test_naca0012_sharp_te_single_strand(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = detect_te_edges(mesh, dihedral_deg=30.0, span_axis="z")
        assert len(r.polylines) == 1
        assert not r.is_blunt
        assert r.te_thickness == 0.0

    def test_naca0012_te_at_chord_extremum(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = detect_te_edges(mesh, dihedral_deg=30.0)
        # NACA 0012 TE is at x = 1.0 (chord = 1.0).
        p = r.polyline
        assert np.allclose(p[:, 0], 1.0, atol=1e-6)

    def test_te_polyline_sorted_by_span(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = detect_te_edges(mesh, dihedral_deg=30.0)
        z = r.polyline[:, 2]
        assert np.all(np.diff(z) >= 0)

    def test_te_spans_full_body(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = detect_te_edges(mesh, dihedral_deg=30.0)
        # Body spans 0..0.5 in z (default fixture).
        z = r.polyline[:, 2]
        assert z.min() <= 0.0 + 1e-6
        assert z.max() >= 0.5 - 1e-6


# ---------------------------------------------------------------------------
# detect_te_edges — blunt TE
# ---------------------------------------------------------------------------


class TestDetectTEBlunt:
    def test_naca0012_blunt_te_two_strands(self, naca0012_blunt_stl):
        mesh = trimesh.load(str(naca0012_blunt_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = detect_te_edges(mesh, dihedral_deg=20.0, span_axis="z")
        assert len(r.polylines) == 2
        assert r.is_blunt

    def test_blunt_thickness_matches_input(self, naca0012_blunt_stl):
        # Fixture is built with blunt_te_thickness=0.01 (1% chord).
        mesh = trimesh.load(str(naca0012_blunt_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = detect_te_edges(mesh, dihedral_deg=20.0)
        assert abs(r.te_thickness - 0.01) < 1e-6

    def test_blunt_strands_at_upper_and_lower(self, naca0012_blunt_stl):
        mesh = trimesh.load(str(naca0012_blunt_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = detect_te_edges(mesh, dihedral_deg=20.0)
        assert len(r.polylines) == 2
        y_means = sorted(p[:, 1].mean() for p in r.polylines)
        assert y_means[0] < 0.0 < y_means[1]
        # Symmetric NACA 0012: upper and lower strands equidistant from y=0.
        assert abs(y_means[0] + y_means[1]) < 1e-6


# ---------------------------------------------------------------------------
# detect_te_edges — error paths
# ---------------------------------------------------------------------------


class TestDetectTEErrors:
    def test_invalid_span_axis_rejected(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with pytest.raises(ParameterError, match="span_axis"):
            detect_te_edges(mesh, span_axis="w")  # type: ignore[arg-type]

    def test_invalid_dihedral_rejected(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with pytest.raises(ParameterError, match="dihedral_deg"):
            detect_te_edges(mesh, dihedral_deg=-5.0)
        with pytest.raises(ParameterError, match="dihedral_deg"):
            detect_te_edges(mesh, dihedral_deg=200.0)

    def test_smooth_body_rejected(self, sphere_stl):
        # An icosphere has no sharp dihedrals, so no TE candidates.
        mesh = trimesh.load(str(sphere_stl), force="mesh")
        with pytest.raises(MeshInputError, match="dihedral"):
            detect_te_edges(mesh, dihedral_deg=60.0)

    def test_high_dihedral_yields_empty(self, naca0012_stl):
        # The fixture's max dihedral peaks below 170°; asking for
        # candidates above that finds none -- error.
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with pytest.raises(MeshInputError, match="dihedral"):
            detect_te_edges(mesh, dihedral_deg=170.0)


# ---------------------------------------------------------------------------
# slice_by_span
# ---------------------------------------------------------------------------


class TestSliceBySpan:
    def _slice_naca(self, naca0012_stl, *, n_stations=5, ni_body=80):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            te = detect_te_edges(mesh, dihedral_deg=30.0)
        return slice_by_span(
            mesh,
            n_stations=n_stations,
            ni_body=ni_body,
            span_axis="z",
            te_result=te,
        )

    def test_returns_one_per_station(self, naca0012_stl):
        slices = self._slice_naca(naca0012_stl, n_stations=7)
        assert len(slices) == 7

    def test_polylines_have_common_ni(self, naca0012_stl):
        slices = self._slice_naca(naca0012_stl, ni_body=64)
        for s in slices:
            assert s.polyline_2d.shape == (64, 2)

    def test_starts_at_te(self, naca0012_stl):
        # NACA 0012 sharp TE at chord = 1.0, thickness = 0.0.
        slices = self._slice_naca(naca0012_stl, ni_body=80)
        for s in slices:
            assert np.isclose(s.polyline_2d[0, 0], 1.0, atol=1e-3)
            assert np.isclose(s.polyline_2d[0, 1], 0.0, atol=1e-3)

    def test_walks_suction_first(self, naca0012_stl):
        # After leaving the upper-TE, the walk should climb onto the
        # suction (positive-thickness) surface.
        slices = self._slice_naca(naca0012_stl, ni_body=20)
        for s in slices:
            assert s.polyline_2d[1, 1] > 0  # i=1 is on suction surface

    def test_le_at_polyline_midpoint(self, naca0012_stl):
        # NACA 0012: LE at chord = 0.0.  After arclength resampling,
        # the LE sits near the polyline midpoint (because the upper
        # and lower surfaces have the same arclength).
        slices = self._slice_naca(naca0012_stl, ni_body=80)
        mid_chord = slices[2].polyline_2d[40, 0]  # i = ni/2
        assert mid_chord < 0.05  # LE region

    def test_arclength_resampling_uniform(self, naca0012_stl):
        # Resampling at uniform arclength fractions should produce
        # roughly-uniform segment lengths.
        slices = self._slice_naca(naca0012_stl, ni_body=100)
        for s in slices:
            seg_lens = np.linalg.norm(np.diff(s.polyline_2d, axis=0), axis=1)
            # Coefficient of variation of segment lengths should be
            # small (< 0.5).  Won't be zero because the polyline is
            # sampled from triangle edges, not exact arclength.
            assert seg_lens.std() / seg_lens.mean() < 0.5

    def test_invalid_span_axis_rejected(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with pytest.raises(ParameterError, match="span_axis"):
            slice_by_span(
                mesh,
                n_stations=5,
                ni_body=80,
                span_axis="w",  # type: ignore[arg-type]
            )

    def test_too_few_stations_rejected(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with pytest.raises(ParameterError, match="n_stations"):
            slice_by_span(mesh, n_stations=1, ni_body=80, span_axis="z")

    def test_too_few_ni_rejected(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with pytest.raises(ParameterError, match="ni_body"):
            slice_by_span(mesh, n_stations=5, ni_body=4, span_axis="z")


# ---------------------------------------------------------------------------
# build_c_curve
# ---------------------------------------------------------------------------


class TestBuildCCurve:
    def _make_slice(self, naca0012_stl, *, ni_body=80):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            te = detect_te_edges(mesh, dihedral_deg=30.0)
        slices = slice_by_span(mesh, n_stations=5, ni_body=ni_body, span_axis="z", te_result=te)
        return slices[2]  # mid-span

    def test_total_length_matches_layout(self, naca0012_stl):
        s = self._make_slice(naca0012_stl, ni_body=80)
        c = build_c_curve(
            s,
            ni_wake=20,
            wake_length=20.0,
            freestream_2d=np.array([1.0, 0.0]),
        )
        # Layout: 20 + 80 + (20 - 1) - 1 = upper_wake + body[1:] + lower_wake[1:]
        # = 20 + 79 + 19 = 118
        assert c.polyline_2d.shape == (118, 2)
        assert c.ni_wake_upper == 20
        assert c.ni_wake_lower == 20

    def test_te_at_expected_indices(self, naca0012_stl):
        s = self._make_slice(naca0012_stl, ni_body=80)
        c = build_c_curve(
            s,
            ni_wake=20,
            wake_length=20.0,
            freestream_2d=np.array([1.0, 0.0]),
        )
        # Upper TE sits at i = ni_wake_upper - 1 = 19 (the wake terminus,
        # which is the body's i=0).  Lower TE sits at i = 19 + 80 - 1 = 98.
        upper_te = c.polyline_2d[19]
        lower_te = c.polyline_2d[19 + 80 - 1]
        # NACA 0012 sharp TE at chord = 1.0, thickness = 0.0.
        assert np.isclose(upper_te[0], 1.0, atol=1e-3)
        assert np.isclose(upper_te[1], 0.0, atol=1e-3)
        assert np.isclose(lower_te[0], 1.0, atol=1e-3)
        assert np.isclose(lower_te[1], 0.0, atol=1e-3)

    def test_wake_tips_at_freestream_distance(self, naca0012_stl):
        s = self._make_slice(naca0012_stl, ni_body=80)
        c = build_c_curve(
            s,
            ni_wake=20,
            wake_length=20.0,
            freestream_2d=np.array([1.0, 0.0]),
        )
        # i = 0 is the upper-wake far-field tip.
        upper_tip = c.polyline_2d[0]
        # i = -1 is the lower-wake far-field tip.
        lower_tip = c.polyline_2d[-1]
        # Both tips are 20 units downstream of the TE (at chord=1.0).
        assert np.isclose(upper_tip[0], 21.0, atol=1e-6)
        assert np.isclose(lower_tip[0], 21.0, atol=1e-6)
        # Sharp TE: both tips also at thickness = 0.
        assert np.isclose(upper_tip[1], 0.0, atol=1e-3)
        assert np.isclose(lower_tip[1], 0.0, atol=1e-3)

    def test_upper_wake_walks_toward_te(self, naca0012_stl):
        s = self._make_slice(naca0012_stl, ni_body=80)
        c = build_c_curve(
            s,
            ni_wake=20,
            wake_length=20.0,
            freestream_2d=np.array([1.0, 0.0]),
        )
        # Upper wake (indices 0..19): chord coordinate decreases.
        upper = c.polyline_2d[:20, 0]
        assert np.all(np.diff(upper) <= 0)

    def test_lower_wake_walks_away_from_te(self, naca0012_stl):
        s = self._make_slice(naca0012_stl, ni_body=80)
        c = build_c_curve(
            s,
            ni_wake=20,
            wake_length=20.0,
            freestream_2d=np.array([1.0, 0.0]),
        )
        # Lower wake (last 19 indices): chord coordinate increases.
        lower = c.polyline_2d[-19:, 0]
        assert np.all(np.diff(lower) >= 0)

    def test_wake_growth_geometric(self, naca0012_stl):
        s = self._make_slice(naca0012_stl, ni_body=80)
        c = build_c_curve(
            s,
            ni_wake=20,
            wake_length=20.0,
            freestream_2d=np.array([1.0, 0.0]),
            wake_growth=1.2,
        )
        # Upper wake cell sizes (from far-field walking toward TE)
        # should monotonically *decrease* as we approach TE.
        upper = c.polyline_2d[:20]
        cell_sizes = np.linalg.norm(np.diff(upper, axis=0), axis=1)
        assert np.all(np.diff(cell_sizes) <= 1e-9)

    def test_wake_growth_uniform(self, naca0012_stl):
        s = self._make_slice(naca0012_stl, ni_body=80)
        c = build_c_curve(
            s,
            ni_wake=20,
            wake_length=20.0,
            freestream_2d=np.array([1.0, 0.0]),
            wake_growth=1.0,
        )
        # Uniform wake: all cell sizes equal.
        upper = c.polyline_2d[:20]
        cell_sizes = np.linalg.norm(np.diff(upper, axis=0), axis=1)
        assert cell_sizes.std() / cell_sizes.mean() < 1e-9

    def test_invalid_ni_wake_rejected(self, naca0012_stl):
        s = self._make_slice(naca0012_stl)
        with pytest.raises(ParameterError, match="ni_wake"):
            build_c_curve(
                s,
                ni_wake=1,
                wake_length=20.0,
                freestream_2d=np.array([1.0, 0.0]),
            )

    def test_invalid_wake_length_rejected(self, naca0012_stl):
        s = self._make_slice(naca0012_stl)
        with pytest.raises(ParameterError, match="wake_length"):
            build_c_curve(
                s,
                ni_wake=20,
                wake_length=-1.0,
                freestream_2d=np.array([1.0, 0.0]),
            )

    def test_invalid_freestream_rejected(self, naca0012_stl):
        s = self._make_slice(naca0012_stl)
        with pytest.raises(ParameterError, match="freestream_2d"):
            build_c_curve(
                s,
                ni_wake=20,
                wake_length=20.0,
                freestream_2d=np.array([0.0, 0.0]),
            )


# ---------------------------------------------------------------------------
# build_cgrid_k0 (orchestrator)
# ---------------------------------------------------------------------------


class TestBuildCGridK0:
    def _build(self, naca0012_stl, **kwargs):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        defaults = {
            "ni_body": 80,
            "ni_wake": 20,
            "n_stations": 5,
            "wake_length": 20.0,
            "span_axis": "z",
        }
        defaults.update(kwargs)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return build_cgrid_k0(mesh, **defaults)

    def test_shape(self, naca0012_stl):
        r = self._build(naca0012_stl)
        # ni_total = 2 * (ni_wake - 1) + ni_body = 38 + 80 = 118
        assert r.nodes.shape == (118, 5, 3)
        assert r.ni_body == 80
        assert r.ni_wake == 20

    def test_span_axis_coord_constant_per_station(self, naca0012_stl):
        r = self._build(naca0012_stl)
        # All nodes at station k share span coordinate.
        for k in range(r.nodes.shape[1]):
            z = r.nodes[:, k, 2]
            assert z.std() < 1e-10

    def test_span_stations_monotonic(self, naca0012_stl):
        r = self._build(naca0012_stl)
        # Station span coords increase along k.
        z = r.nodes[0, :, 2]
        assert np.all(np.diff(z) > 0)

    def test_te_at_chord_extremum(self, naca0012_stl):
        r = self._build(naca0012_stl)
        # Upper-TE node sits at the boundary between wake and body:
        # i = ni_wake_upper - 1 = 19 (= ni_wake - 1).
        # Same for lower-TE at i = ni_wake - 1 + ni_body - 1 = 98.
        upper_te = r.nodes[r.ni_wake - 1, :, :]
        lower_te = r.nodes[r.ni_wake - 1 + r.ni_body - 1, :, :]
        # NACA 0012 sharp TE at x = 1.0, y = 0.
        assert np.allclose(upper_te[:, 0], 1.0, atol=1e-3)
        assert np.allclose(upper_te[:, 1], 0.0, atol=1e-3)
        assert np.allclose(lower_te[:, 0], 1.0, atol=1e-3)
        assert np.allclose(lower_te[:, 1], 0.0, atol=1e-3)

    def test_wake_tips_match_wake_length(self, naca0012_stl):
        r = self._build(naca0012_stl, wake_length=15.0)
        # Wake tips: i = 0 (upper) and i = -1 (lower).  Both at
        # chord = 1.0 (TE) + 15.0 = 16.0.
        upper_tip = r.nodes[0, :, 0]
        lower_tip = r.nodes[-1, :, 0]
        assert np.allclose(upper_tip, 16.0, atol=1e-6)
        assert np.allclose(lower_tip, 16.0, atol=1e-6)

    def test_returns_te_result(self, naca0012_stl):
        r = self._build(naca0012_stl)
        # Sanity: TE detection ran and produced one strand for sharp NACA.
        assert not r.te_result.is_blunt
        assert len(r.te_result.polylines) == 1

    def test_default_freestream_along_chord(self, naca0012_stl):
        # No freestream supplied: should default to +x (chord axis).
        r = self._build(naca0012_stl)
        # Wake tips lie at +x relative to TE.
        upper_te_x = r.nodes[r.ni_wake - 1, 0, 0]
        upper_tip_x = r.nodes[0, 0, 0]
        assert upper_tip_x > upper_te_x

    def test_invalid_span_axis_rejected(self, naca0012_stl):
        with pytest.raises(ParameterError, match="span_axis"):
            self._build(naca0012_stl, span_axis="w")

    def test_freestream_along_span_rejected(self, naca0012_stl):
        with pytest.raises(ParameterError, match="parallel to the span"):
            # Freestream purely along z (= span axis) leaves no in-plane
            # component for the wake.
            self._build(
                naca0012_stl,
                freestream_world=np.array([0.0, 0.0, 1.0]),
            )

    def test_zero_freestream_rejected(self, naca0012_stl):
        with pytest.raises(ParameterError, match="zero magnitude"):
            self._build(
                naca0012_stl,
                freestream_world=np.array([0.0, 0.0, 0.0]),
            )

    def test_per_block_layers_shapes(self, naca0012_stl):
        r = self._build(naca0012_stl)
        upper, body, lower = r.per_block_layers()
        # ni_wake=20, ni_body=80, n_stations=5
        assert upper.shape == (20, 5, 3)
        assert body.shape == (80, 5, 3)
        assert lower.shape == (20, 5, 3)
        assert r.ni_total == 118

    def test_per_block_layers_share_te_columns(self, naca0012_stl):
        r = self._build(naca0012_stl)
        upper, body, lower = r.per_block_layers()
        # Upper-TE column: upper[-1] == body[0] (numpy view aliasing).
        np.testing.assert_array_equal(upper[-1], body[0])
        # Lower-TE column: body[-1] == lower[0].
        np.testing.assert_array_equal(body[-1], lower[0])
        # And they share memory.
        assert np.shares_memory(upper, body)
        assert np.shares_memory(body, lower)

    def test_per_block_layers_aliasing_propagates(self, naca0012_stl):
        # Writing to a TE node through one view must update the other.
        r = self._build(naca0012_stl)
        upper, body, lower = r.per_block_layers()
        upper[-1, 0, :] = [42.0, 42.0, 42.0]
        np.testing.assert_array_equal(body[0, 0, :], [42.0, 42.0, 42.0])
        body[-1, 0, :] = [99.0, 99.0, 99.0]
        np.testing.assert_array_equal(lower[0, 0, :], [99.0, 99.0, 99.0])


# ---------------------------------------------------------------------------
# CGridTopology
# ---------------------------------------------------------------------------


class TestCGridTopologySingleBlock:
    """Legacy 1-block C-grid topology (opt-in via ``blocks=1``)."""

    def _build(self, naca0012_stl, *, nk=10, **kwargs):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        defaults = {
            "ni_body": 80,
            "ni_wake": 20,
            "n_stations": 5,
            "wake_length": 20.0,
            "span_axis": "z",
            "blocks": 1,
        }
        defaults.update(kwargs)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return CGridTopology().build(surface=None, nk=nk, mesh=mesh, **defaults)

    def test_returns_single_block(self, naca0012_stl):
        blocks = self._build(naca0012_stl)
        assert len(blocks) == 1
        assert blocks[0].label == "cgrid"

    def test_block_shape(self, naca0012_stl):
        blocks = self._build(naca0012_stl, nk=12)
        b = blocks[0]
        # ni_total = 2 * (ni_wake - 1) + ni_body = 38 + 80 = 118
        assert b.nodes.shape == (118, 5, 12, 3)

    def test_block_not_periodic(self, naca0012_stl):
        b = self._build(naca0012_stl)[0]
        assert not b.periodic_i
        assert not b.periodic_j

    def test_default_face_bcs(self, naca0012_stl):
        b = self._build(naca0012_stl)[0]
        # k_lo = WALL; everything else = FREE; nothing SHARED.
        assert b.face("k_lo").bc == FaceBC.WALL
        assert b.face("k_hi").bc == FaceBC.FREE
        assert b.face("i_lo").bc == FaceBC.FREE
        assert b.face("i_hi").bc == FaceBC.FREE
        assert b.face("j_lo").bc == FaceBC.FREE
        assert b.face("j_hi").bc == FaceBC.FREE

    def test_k0_matches_surface(self, naca0012_stl):
        # k=0 layer should match what build_cgrid_k0 produces directly.
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            k0 = build_cgrid_k0(
                mesh,
                ni_body=80,
                ni_wake=20,
                n_stations=5,
                wake_length=20.0,
                span_axis="z",
            )
        b = self._build(naca0012_stl)[0]
        np.testing.assert_array_equal(b.nodes[:, :, 0, :], k0.nodes)

    def test_pre_built_surface_accepted(self, naca0012_stl):
        # When surface is supplied, the topology skips the build_cgrid_k0
        # call and uses it directly.
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            k0 = build_cgrid_k0(
                mesh,
                ni_body=80,
                ni_wake=20,
                n_stations=5,
                wake_length=20.0,
                span_axis="z",
            )
        blocks = CGridTopology().build(surface=k0.nodes, nk=8, ni_body=80, ni_wake=20, blocks=1)
        b = blocks[0]
        assert b.nodes.shape == (118, 5, 8, 3)
        np.testing.assert_array_equal(b.nodes[:, :, 0, :], k0.nodes)

    def test_missing_mesh_rejected(self):
        # No surface, no mesh — must error.
        with pytest.raises(ValueError, match="surface"):
            CGridTopology().build(surface=None, nk=10, ni_body=80)

    def test_missing_required_param_rejected(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with pytest.raises(ValueError, match="ni_wake"):
            CGridTopology().build(
                surface=None,
                nk=10,
                mesh=mesh,
                ni_body=80,
                # ni_wake missing
                n_stations=5,
                wake_length=20.0,
            )


# ---------------------------------------------------------------------------
# End-to-end march
# ---------------------------------------------------------------------------


class TestCGridMarchSingleBlock:
    """Legacy 1-block C-grid marcher (opt-in via ``blocks=1``)."""

    def test_march_naca0012(self, naca0012_stl):
        from shore.mesh import Body, Mesh

        body = Body.from_stl(naca0012_stl, skip_hygiene=True)
        topo = CGridTopology()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks = topo.build(
                surface=None,
                nk=10,
                mesh=body.mesh,
                ni_body=80,
                ni_wake=20,
                n_stations=5,
                wake_length=15.0,
                blocks=1,
            )
            m = Mesh(body=body, blocks=blocks, topology=topo)
            m.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)
        b = m.blocks[0]
        assert b.nodes.shape == (118, 5, 10, 3)
        # Outer layer should sit further from the body than k=0.
        # Mid-span LE region: chord ~ 0.0, body at y ~ 0.06.
        # i = 19 + 40 = 59 is roughly LE.
        k0_y = b.nodes[59, 2, 0, 1]
        kN_y = b.nodes[59, 2, -1, 1]
        # The marched nodes should have moved outward (|y| larger).
        assert abs(kN_y) > abs(k0_y) + 1e-4

    def test_march_no_smoothing(self, naca0012_stl):
        from shore.mesh import Body, Mesh

        body = Body.from_stl(naca0012_stl, skip_hygiene=True)
        topo = CGridTopology()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks = topo.build(
                surface=None,
                nk=8,
                mesh=body.mesh,
                ni_body=60,
                ni_wake=15,
                n_stations=4,
                wake_length=10.0,
                blocks=1,
            )
            m = Mesh(body=body, blocks=blocks, topology=topo)
            m.march(ds=1e-3, growth=1.2, smooth_strength=0.0, smooth_iters=0)
        b = m.blocks[0]
        # 2 * (ni_wake - 1) + ni_body = 28 + 60 = 88
        assert b.nodes.shape == (88, 4, 8, 3)

    def test_march_raises_on_negative_jacobian(self, naca0012_stl):
        from shore.exceptions import MarchingError
        from shore.mesh import Body, Mesh

        body = Body.from_stl(naca0012_stl, skip_hygiene=True)
        topo = CGridTopology()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks = topo.build(
                surface=None,
                nk=10,
                mesh=body.mesh,
                ni_body=80,
                ni_wake=20,
                n_stations=5,
                wake_length=20.0,
                blocks=1,
            )
            m = Mesh(body=body, blocks=blocks, topology=topo)
            # Crank unlimited (legacy) smoothing so the LE inverts on first
            # march: this test exercises the Jacobian guard, not the smoother.
            with pytest.raises(MarchingError, match="Negative Jacobian"):
                m.march(
                    ds=1e-3,
                    growth=1.15,
                    smooth_strength=0.5,
                    smooth_iters=4,
                    smooth_mode="legacy",
                )


# ---------------------------------------------------------------------------
# 3-block sharp-TE C-grid (default topology in v0.2)
# ---------------------------------------------------------------------------


class TestThreeBlockCGrid:
    """Default 3-block C-grid: upper / body / lower with SHARED seams."""

    def _build(self, naca0012_stl, *, nk=10, **kwargs):
        from shore.mesh import Body

        body = Body.from_stl(naca0012_stl, skip_hygiene=True)
        topo = CGridTopology()
        defaults = {
            "ni_body": 80,
            "ni_wake": 20,
            "n_stations": 5,
            "wake_length": 15.0,
            "blocks": 3,
        }
        defaults.update(kwargs)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks = topo.build(surface=None, nk=nk, mesh=body.mesh, **defaults)
        return body, topo, blocks

    # ---- block list shape ---------------------------------------------------

    def test_returns_three_blocks_in_order(self, naca0012_stl):
        _, _, blocks = self._build(naca0012_stl)
        assert len(blocks) == 3
        assert [b.label for b in blocks] == ["upper", "body", "lower"]

    def test_block_shapes(self, naca0012_stl):
        _, _, blocks = self._build(naca0012_stl, nk=12)
        u, b, lo = blocks
        # ni_wake=20, ni_body=80, n_stations=5, nk=12
        assert u.nodes.shape == (20, 5, 12, 3)
        assert b.nodes.shape == (80, 5, 12, 3)
        assert lo.nodes.shape == (20, 5, 12, 3)

    # ---- TE column coincidence (numpy view aliasing) ------------------------

    def test_upper_te_column_coincident(self, naca0012_stl):
        _, _, blocks = self._build(naca0012_stl)
        u, b, _ = blocks
        # upper.i_hi (row -1) should equal body.i_lo (row 0) bit-exactly.
        np.testing.assert_array_equal(u.nodes[-1, :, :, :], b.nodes[0, :, :, :])

    def test_lower_te_column_coincident(self, naca0012_stl):
        _, _, blocks = self._build(naca0012_stl)
        _, b, lo = blocks
        np.testing.assert_array_equal(b.nodes[-1, :, :, :], lo.nodes[0, :, :, :])

    # ---- SHARED seam wiring ------------------------------------------------

    def test_upper_te_seam_shared(self, naca0012_stl):
        _, _, blocks = self._build(naca0012_stl)
        u, b, _ = blocks
        assert u.face("i_hi").bc == FaceBC.SHARED
        assert b.face("i_lo").bc == FaceBC.SHARED
        # Partner pointers cross-reference each other.
        assert u.face("i_hi").partner is b.face("i_lo")
        assert b.face("i_lo").partner is u.face("i_hi")

    def test_lower_te_seam_shared(self, naca0012_stl):
        _, _, blocks = self._build(naca0012_stl)
        _, b, lo = blocks
        assert b.face("i_hi").bc == FaceBC.SHARED
        assert lo.face("i_lo").bc == FaceBC.SHARED
        assert b.face("i_hi").partner is lo.face("i_lo")

    def test_wake_cut_seam_shared(self, naca0012_stl):
        _, _, blocks = self._build(naca0012_stl)
        u, _, lo = blocks
        assert u.face("k_lo").bc == FaceBC.SHARED
        assert lo.face("k_lo").bc == FaceBC.SHARED
        assert u.face("k_lo").partner is lo.face("k_lo")

    def test_wake_cut_axis_map_has_i_reversal(self, naca0012_stl):
        # The wake cut on upper.k_lo is the same physical line as
        # lower.k_lo, but upper's i runs far-tip→TE while lower's
        # i runs TE→far-tip.  AxisMap must record the i-flip.
        _, _, blocks = self._build(naca0012_stl)
        u, _, _ = blocks
        amap = u.face("k_lo").axis_map
        assert amap is not None
        # k_lo's along-axes are (i, j); we expect ((i, True), (j, False)).
        assert amap.along == (("i", True), ("j", False))

    def test_te_seam_axis_maps_are_identity(self, naca0012_stl):
        # Same-axis i↔i seams: along-seam (j, k) preserved, no flips.
        _, _, blocks = self._build(naca0012_stl)
        u, b, lo = blocks
        for face_owner, face_name in (
            (u, "i_hi"),
            (b, "i_lo"),
            (b, "i_hi"),
            (lo, "i_lo"),
        ):
            amap = face_owner.face(face_name).axis_map
            assert amap is not None
            assert amap.along == (("j", False), ("k", False))

    # ---- non-shared faces ---------------------------------------------------

    def test_outer_far_tips_are_free(self, naca0012_stl):
        _, _, blocks = self._build(naca0012_stl)
        u, _, lo = blocks
        # upper.i_lo = upper-wake far-field tip.
        assert u.face("i_lo").bc == FaceBC.FREE
        # lower.i_hi = lower-wake far-field tip.
        assert lo.face("i_hi").bc == FaceBC.FREE

    def test_body_k_lo_is_wall(self, naca0012_stl):
        _, _, blocks = self._build(naca0012_stl)
        _, b, _ = blocks
        # The body block's k_lo is the airfoil wall (not part of any seam).
        assert b.face("k_lo").bc == FaceBC.WALL

    # ---- end-to-end march ---------------------------------------------------

    def test_march_keeps_all_blocks_jacobian_positive(self, naca0012_stl):
        from shore.mesh import Mesh
        from shore.volume.jacobian import _jacobian_field_open

        body, topo, blocks = self._build(naca0012_stl, nk=10)
        m = Mesh(body=body, blocks=blocks, topology=topo)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)
        # Match the marcher's own tolerance: -1e-6 absorbs FP noise
        # in the per-cell scalar triple product without masking real
        # cell inversions (which are -1e-5 or larger).
        for blk in m.blocks:
            jac = _jacobian_field_open(blk.nodes)
            assert jac.min() > -1e-6, f"{blk.label}: jac.min() = {jac.min():.3e}"

    def test_march_preserves_te_column_coincidence(self, naca0012_stl):
        # After marching, the SHARED TE columns should still be
        # bit-coincident across all k layers — the marcher's
        # restitch step averages the upper-TE and lower-TE columns
        # back into agreement after every step.
        from shore.mesh import Mesh

        body, topo, blocks = self._build(naca0012_stl, nk=8)
        m = Mesh(body=body, blocks=blocks, topology=topo)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)
        u, b, lo = m.blocks
        np.testing.assert_array_equal(u.nodes[-1, :, :, :], b.nodes[0, :, :, :])
        np.testing.assert_array_equal(b.nodes[-1, :, :, :], lo.nodes[0, :, :, :])

    # ---- error paths --------------------------------------------------------

    def test_blunt_te_rejected(self, naca0012_blunt_stl):
        from shore.mesh import Body

        body = Body.from_stl(naca0012_blunt_stl, skip_hygiene=True)
        topo = CGridTopology()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ValueError, match="3-block layout requires a sharp TE"):
                topo.build(
                    surface=None,
                    nk=10,
                    mesh=body.mesh,
                    ni_body=80,
                    ni_wake=20,
                    n_stations=5,
                    wake_length=15.0,
                    blocks=3,
                )

    def test_invalid_blocks_value_rejected(self, naca0012_stl):
        from shore.mesh import Body

        body = Body.from_stl(naca0012_stl, skip_hygiene=True)
        topo = CGridTopology()
        with pytest.raises(ValueError, match="blocks"):
            topo.build(
                surface=None,
                nk=10,
                mesh=body.mesh,
                ni_body=80,
                ni_wake=20,
                n_stations=5,
                wake_length=15.0,
                blocks=2,
            )


# ---------------------------------------------------------------------------
# 4-block blunt-TE C-grid (with wake-base block W)
# ---------------------------------------------------------------------------


class TestFourBlockCGrid:
    """Blunt-TE 4-block C-grid: upper / body / lower / wake_base."""

    def _build(self, naca0012_blunt_stl, *, nk=10, nk_w=3, **kwargs):
        from shore.mesh import Body

        body = Body.from_stl(naca0012_blunt_stl, skip_hygiene=True)
        topo = CGridTopology()
        defaults = {
            "ni_body": 80,
            "ni_wake": 20,
            "n_stations": 5,
            "wake_length": 15.0,
            "blocks": 4,
            "nk_w": nk_w,
            "te_dihedral_deg": 20.0,  # blunt TE has shallower dihedral
        }
        defaults.update(kwargs)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks = topo.build(surface=None, nk=nk, mesh=body.mesh, **defaults)
        return body, topo, blocks

    # ---- block list shape ---------------------------------------------------

    def test_returns_four_blocks_in_order(self, naca0012_blunt_stl):
        _, _, blocks = self._build(naca0012_blunt_stl)
        assert len(blocks) == 4
        assert [b.label for b in blocks] == ["upper", "body", "lower", "wake_base"]

    def test_wake_base_shape(self, naca0012_blunt_stl):
        _, _, blocks = self._build(naca0012_blunt_stl, nk=12, nk_w=5)
        wb = blocks[3]
        # ni_wake=20, n_stations=5, nk_w=5
        assert wb.nodes.shape == (20, 5, 5, 3)

    def test_wake_base_nk_w_minimum_enforced(self, naca0012_blunt_stl):
        # nk_w >= 2 (need at least 2 layers to bridge the gap).
        from shore.mesh import Body

        body = Body.from_stl(naca0012_blunt_stl, skip_hygiene=True)
        topo = CGridTopology()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ValueError, match="nk_w"):
                topo.build(
                    surface=None,
                    nk=10,
                    mesh=body.mesh,
                    ni_body=80,
                    ni_wake=20,
                    n_stations=5,
                    wake_length=15.0,
                    blocks=4,
                    nk_w=1,
                    te_dihedral_deg=20.0,
                )

    # ---- W's geometry: bridges the TE thickness gap -------------------------

    def test_wake_base_k_lo_matches_lower_k0_reversed(self, naca0012_blunt_stl):
        # W's k=0 layer = lower's k=0 strip with i reversed (W's i and L's
        # i run opposite directions along the wake).
        _, _, blocks = self._build(naca0012_blunt_stl)
        _, _, lower, wb = blocks
        np.testing.assert_array_equal(wb.nodes[:, :, 0, :], lower.nodes[::-1, :, 0, :])

    def test_wake_base_k_hi_matches_upper_k0(self, naca0012_blunt_stl):
        # W's k=last layer = upper's k=0 strip directly (same i direction).
        _, _, blocks = self._build(naca0012_blunt_stl)
        upper, _, _, wb = blocks
        np.testing.assert_array_equal(wb.nodes[:, :, -1, :], upper.nodes[:, :, 0, :])

    def test_wake_base_thickness_matches_te_thickness(self, naca0012_blunt_stl):
        # The fixture is 1% chord blunt; W's k-extent (= TE thickness)
        # should be 0.01 along y (the thickness axis).
        _, _, blocks = self._build(naca0012_blunt_stl)
        wb = blocks[3]
        y_extent = wb.nodes[..., 1].max() - wb.nodes[..., 1].min()
        assert abs(y_extent - 0.01) < 1e-9

    # ---- SHARED seam wiring ------------------------------------------------

    def test_wake_base_kh_seam_to_upper_is_identity(self, naca0012_blunt_stl):
        _, _, blocks = self._build(naca0012_blunt_stl)
        upper, _, _, wb = blocks
        assert wb.face("k_hi").bc == FaceBC.SHARED
        assert upper.face("k_lo").bc == FaceBC.SHARED
        assert wb.face("k_hi").partner is upper.face("k_lo")
        # Same i direction → identity AxisMap.
        assert wb.face("k_hi").axis_map.along == (("i", False), ("j", False))

    def test_wake_base_kl_seam_to_lower_has_i_flip(self, naca0012_blunt_stl):
        _, _, blocks = self._build(naca0012_blunt_stl)
        _, _, lower, wb = blocks
        assert wb.face("k_lo").bc == FaceBC.SHARED
        assert lower.face("k_lo").bc == FaceBC.SHARED
        assert wb.face("k_lo").partner is lower.face("k_lo")
        # i flips because lower's i runs TE → far-tip while W's runs far-tip → TE.
        assert wb.face("k_lo").axis_map.along == (("i", True), ("j", False))

    def test_wake_base_te_face_is_wall(self, naca0012_blunt_stl):
        # W.i_hi sits at the TE base — the bluff-base wall of the airfoil.
        _, _, blocks = self._build(naca0012_blunt_stl)
        wb = blocks[3]
        assert wb.face("i_hi").bc == FaceBC.WALL

    def test_wake_base_far_tip_face_is_free(self, naca0012_blunt_stl):
        _, _, blocks = self._build(naca0012_blunt_stl)
        wb = blocks[3]
        assert wb.face("i_lo").bc == FaceBC.FREE

    # ---- end-to-end march --------------------------------------------------

    def test_march_keeps_all_four_blocks_jacobian_positive(self, naca0012_blunt_stl):
        from shore.mesh import Mesh
        from shore.volume.jacobian import _jacobian_field_open

        body, topo, blocks = self._build(naca0012_blunt_stl, nk=10)
        m = Mesh(body=body, blocks=blocks, topology=topo)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)
        for blk in m.blocks:
            jac = _jacobian_field_open(blk.nodes)
            assert jac.min() > -1e-6, f"{blk.label}: jac.min() = {jac.min():.3e}"

    def test_wake_base_unchanged_by_march(self, naca0012_blunt_stl):
        # W is a static block — geometry is fixed at build time.
        from shore.mesh import Mesh

        body, topo, blocks = self._build(naca0012_blunt_stl, nk=10)
        wb_before = blocks[3].nodes.copy()
        m = Mesh(body=body, blocks=blocks, topology=topo)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.march(ds=5e-4, growth=1.15, smooth_strength=0.1, smooth_iters=1)
        np.testing.assert_array_equal(m.blocks[3].nodes, wb_before)

    # ---- error paths --------------------------------------------------------

    def test_sharp_te_rejected_for_4_blocks(self, naca0012_stl):
        # A sharp TE has zero gap to bridge — W would be degenerate.
        from shore.mesh import Body

        body = Body.from_stl(naca0012_stl, skip_hygiene=True)
        topo = CGridTopology()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ValueError, match="4-block layout requires a blunt TE"):
                topo.build(
                    surface=None,
                    nk=10,
                    mesh=body.mesh,
                    ni_body=80,
                    ni_wake=20,
                    n_stations=5,
                    wake_length=15.0,
                    blocks=4,
                )


# ---------------------------------------------------------------------------
# Body-TE rounding (full-span)
# ---------------------------------------------------------------------------


class TestRoundBodyTEInSurface:
    def _build(self, naca0012_stl):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return build_cgrid_k0(
                mesh,
                ni_body=80,
                ni_wake=20,
                n_stations=5,
                wake_length=20.0,
                span_axis="z",
            )

    def test_rounds_in_place(self, naca0012_stl):
        r = self._build(naca0012_stl)
        # Verify input is sharp.
        i_te_upper = r.ni_wake - 1
        i_te_lower = i_te_upper + r.ni_body - 1
        for j in range(r.n_stations):
            np.testing.assert_allclose(
                r.nodes[i_te_upper, j, :], r.nodes[i_te_lower, j, :], atol=1e-9
            )
        # Round, with a known thickness.
        round_body_te_in_surface(
            r.nodes,
            ni_wake=r.ni_wake,
            ni_body=r.ni_body,
            te_thickness=0.005,
            span_axis="z",
        )
        # After rounding, the i=0/i=ni-1 columns are separated by 0.005 in y.
        for j in range(r.n_stations):
            sep = float(np.linalg.norm(r.nodes[i_te_upper, j, :] - r.nodes[i_te_lower, j, :]))
            np.testing.assert_allclose(sep, 0.005, atol=1e-12)
            # Antisymmetric split around the original sharp-TE location.
            midpoint = 0.5 * (r.nodes[i_te_upper, j, :] + r.nodes[i_te_lower, j, :])
            # The midpoint y should still be 0 (NACA 0012 sharp TE at y=0).
            assert abs(midpoint[1]) < 1e-9

    def test_wake_te_columns_follow_body(self, naca0012_stl):
        """Wake blocks' k=0 strips share memory with body's TE columns; rounding propagates."""
        r = self._build(naca0012_stl)
        upper, body, lower = r.per_block_layers()  # views
        # All three views share memory with r.nodes.
        round_body_te_in_surface(
            r.nodes,
            ni_wake=r.ni_wake,
            ni_body=r.ni_body,
            te_thickness=0.005,
            span_axis="z",
        )
        # Upper wake's i=ni_wake-1 column == body's i=0 column (upper-TE node).
        np.testing.assert_array_equal(upper[-1, :, :], body[0, :, :])
        # Lower wake's i=0 column == body's i=ni-1 column (lower-TE node).
        np.testing.assert_array_equal(lower[0, :, :], body[-1, :, :])
        # And both moved by ±te/2 in y.
        for j in range(r.n_stations):
            assert abs(upper[-1, j, 1] - 0.0025) < 1e-12  # +te/2
            assert abs(lower[0, j, 1] + 0.0025) < 1e-12  # -te/2

    def test_wake_strips_at_full_te_offset(self, naca0012_stl):
        """Whole upper / lower wake k=0 strips sit at y = ±te/2, not just at the TE.

        Regression for the singular-cells-after-TE bug: if only the
        TE-side endpoints get shifted, the wake strips converge to y=0
        downstream of the TE and produce a degenerate near-zero-thickness
        wake_base block.
        """
        r = self._build(naca0012_stl)
        upper, _body, lower = r.per_block_layers()
        round_body_te_in_surface(
            r.nodes,
            ni_wake=r.ni_wake,
            ni_body=r.ni_body,
            te_thickness=0.005,
            span_axis="z",
        )
        half = 0.0025
        # Every node of the upper wake at every station should have y = +half.
        for j in range(r.n_stations):
            np.testing.assert_allclose(upper[:, j, 1], +half, atol=1e-12)
            np.testing.assert_allclose(lower[:, j, 1], -half, atol=1e-12)

    def test_invalid_te_thickness_rejected(self, naca0012_stl):
        r = self._build(naca0012_stl)
        with pytest.raises(ParameterError, match="te_thickness"):
            round_body_te_in_surface(
                r.nodes,
                ni_wake=r.ni_wake,
                ni_body=r.ni_body,
                te_thickness=-0.01,
                span_axis="z",
            )

    def test_invalid_span_axis_rejected(self, naca0012_stl):
        r = self._build(naca0012_stl)
        with pytest.raises(ParameterError, match="span_axis"):
            round_body_te_in_surface(
                r.nodes,
                ni_wake=r.ni_wake,
                ni_body=r.ni_body,
                te_thickness=0.001,
                span_axis="bogus",
            )

    def test_already_blunt_body_rejected(self, naca0012_stl):
        """Calling on an already-rounded surface raises (avoids double-rounding)."""
        r = self._build(naca0012_stl)
        round_body_te_in_surface(
            r.nodes,
            ni_wake=r.ni_wake,
            ni_body=r.ni_body,
            te_thickness=0.001,
            span_axis="z",
        )
        with pytest.raises(ParameterError, match="non-sharp"):
            round_body_te_in_surface(
                r.nodes,
                ni_wake=r.ni_wake,
                ni_body=r.ni_body,
                te_thickness=0.001,
                span_axis="z",
            )


# ---------------------------------------------------------------------------
# Conformal caps (CGridTopology cap_tips path)
# ---------------------------------------------------------------------------


class TestCGridConformalCaps:
    """When cap_tips is set on CGridTopology.build, the C-grid is rounded
    over its full span and the cap blocks are emitted with airfoil-wall
    nodes bit-coincident with the body wrap.
    """

    def _build(self, naca0012_stl, *, cap_tips="both", **overrides):
        mesh = trimesh.load(str(naca0012_stl), force="mesh")
        topo = CGridTopology()
        kwargs: dict = {
            "surface": None,
            "nk": 8,
            "ni_body": 60,
            "ni_wake": 15,
            "n_stations": 5,
            "wake_length": 5.0,
            "span_axis": "z",
            "spacing_k": None,
            "mesh": mesh,
            "blocks": 3,  # auto-promote to 4 when caps are on
            "cap_tips": cap_tips,
            "cap_te_thickness_frac": 0.005,
            "cap_le_frac": 0.05,
            "cap_te_frac": 0.95,
            "cap_n_cap_j": 8,
        }
        kwargs.update(overrides)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return topo.build(**kwargs)

    def test_returns_4_cgrid_plus_12_cap_blocks_for_both(self, naca0012_stl):
        """4 C-grid + 6 cap_jlo + 6 cap_jhi = 16 blocks total."""
        blocks = self._build(naca0012_stl, cap_tips="both")
        assert len(blocks) == 16
        labels = [b.label for b in blocks]
        assert labels[:4] == ["upper", "body", "lower", "wake_base"]
        # 6 jlo cap blocks: le, upper, lower, te_upper, te_lower, wake_base.
        jlo_labels = sorted(b.label for b in blocks if b.label.startswith("cap_jlo_"))
        assert jlo_labels == [
            "cap_jlo_le",
            "cap_jlo_lower",
            "cap_jlo_te_lower",
            "cap_jlo_te_upper",
            "cap_jlo_upper",
            "cap_jlo_wake_base",
        ]

    def test_returns_4_plus_6_blocks_for_jlo_only(self, naca0012_stl):
        blocks = self._build(naca0012_stl, cap_tips="jlo")
        assert len(blocks) == 10
        labels = [b.label for b in blocks]
        assert labels[:4] == ["upper", "body", "lower", "wake_base"]
        for lbl in labels[4:]:
            assert lbl.startswith("cap_jlo_")

    def test_body_te_is_rounded(self, naca0012_stl):
        """After build, the body's TE is split apart by the requested thickness."""
        # Use an absolute thickness so we know exactly what to compare to.
        blocks = self._build(
            naca0012_stl,
            cap_tips="both",
            cap_te_thickness=0.005,
        )
        body = blocks[1]
        # body[0] = upper-TE column (i=0); body[-1] = lower-TE column (i=ni-1).
        for j in range(body.nodes.shape[1]):
            sep = float(np.linalg.norm(body.nodes[0, j, 0, :] - body.nodes[-1, j, 0, :]))
            np.testing.assert_allclose(sep, 0.005, atol=1e-12)

    def test_cap_airfoil_wall_matches_body_bit_exact(self, naca0012_stl):
        """Cap airfoil-wall nodes coincide bit-exactly with body wall nodes at the tip."""
        blocks = self._build(naca0012_stl, cap_tips="jlo")
        body = blocks[1]
        # cap_jlo_le: the LE-wrap.  Its (j_hi, k_hi) airfoil-wall edge
        # should equal body.nodes[le_up : le_lo + 1, 0, 0, :].
        labels = {b.label: b for b in blocks}
        cap_le = labels["cap_jlo_le"]
        cap_upper = labels["cap_jlo_upper"]
        cap_lower = labels["cap_jlo_lower"]
        cap_te_upper = labels["cap_jlo_te_upper"]
        cap_te_lower = labels["cap_jlo_te_lower"]

        # The tip is jlo (j=0).
        body_wall_at_tip = body.nodes[:, 0, 0, :]  # (ni_body, 3)

        # cap_upper: i sweeps body i = le_up → te_up (reversed); the wall
        # edge equals body_wall_at_tip[le_up, ..., te_up].  We verify that
        # cap_upper[0, -1, -1, :] sits on the body wall and cap_upper[-1, -1, -1, :]
        # also does — but to do that we need to know the LE/TE cut indices
        # which the topology computed internally.  We invert the geometry:
        # for each i, cap_upper[i, -1, -1, :] should coincide with some
        # body wall point.  Easier: assert that every cap airfoil-wall node
        # is bit-exactly equal to *some* body wall node by membership.
        body_set = {tuple(p) for p in body_wall_at_tip}
        for cap_blk in (cap_le, cap_upper, cap_lower, cap_te_upper, cap_te_lower):
            wall_edge = cap_blk.nodes[:, -1, -1, :]
            for p in wall_edge:
                assert tuple(p) in body_set, (
                    f"Cap airfoil-wall node {p} on {cap_blk.label} is not "
                    f"bit-coincident with any body wall node at the tip."
                )

    def test_blocks_1_plus_caps_rejected(self, naca0012_stl):
        with pytest.raises(ValueError, match="blocks=1"):
            self._build(naca0012_stl, cap_tips="both", blocks=1)

    def test_invalid_cap_tips_rejected(self, naca0012_stl):
        with pytest.raises(ValueError, match="cap_tips"):
            self._build(naca0012_stl, cap_tips="bogus")

    # Note: a blunt-body conformal-cap test would go here
    # (`naca0012_blunt_stl` fixture + ``blocks=4`` + ``cap_tips="both"``,
    # asserting that the cap uses the body's existing TE thickness
    # without re-rounding).  The blunt-body C-grid build pipeline isn't
    # exercised in the test suite today (the `build_cgrid_k0` orchestrator
    # is only tested against the sharp NACA fixture), so this test would
    # depend on an untested upstream path; deferred until the blunt
    # `build_cgrid_k0` path is itself covered.
