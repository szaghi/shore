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


# ---------------------------------------------------------------------------
# CGridTopology
# ---------------------------------------------------------------------------


class TestCGridTopology:
    def _build(self, naca0012_stl, *, nk=10, **kwargs):
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
        blocks = CGridTopology().build(surface=k0.nodes, nk=8)
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


class TestCGridMarch:
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
            )
            m = Mesh(body=body, blocks=blocks, topology=topo)
            # Crank smoothing so the LE inverts on first march.
            with pytest.raises(MarchingError, match="Negative Jacobian"):
                m.march(
                    ds=1e-3,
                    growth=1.15,
                    smooth_strength=0.5,
                    smooth_iters=4,
                )
