"""Tests for shore.shapes — STL primitives + section extrusion."""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError
from shore.profiles import naca
from shore.shapes import box, capsule, cylinder, extrude_section, sphere, torus

# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------


class TestSphere:
    def test_default(self):
        m = sphere()
        assert m.is_watertight
        assert abs(m.bounds[1, 0] - 1.0) < 1e-3
        assert abs(m.bounds[0, 0] + 1.0) < 1e-3

    def test_custom_radius(self):
        m = sphere(radius=2.5)
        assert abs(m.bounds[1, 0] - 2.5) < 1e-3

    def test_subdivisions_changes_face_count(self):
        m0 = sphere(subdivisions=0)
        m3 = sphere(subdivisions=3)
        assert len(m3.faces) > len(m0.faces) * 50  # 4^3 = 64x more

    def test_invalid_radius_rejected(self):
        with pytest.raises(ParameterError, match="radius"):
            sphere(radius=-1.0)

    def test_invalid_subdivisions_rejected(self):
        with pytest.raises(ParameterError, match="subdivisions"):
            sphere(subdivisions=-1)


class TestCylinder:
    def test_default(self):
        m = cylinder()
        assert m.is_watertight
        # Default height=1, axis=z
        assert abs(m.bounds[1, 2] - 0.5) < 1e-6
        assert abs(m.bounds[0, 2] + 0.5) < 1e-6

    def test_custom(self):
        m = cylinder(radius=2.0, height=3.0)
        assert abs(m.bounds[1, 0] - 2.0) < 1e-3  # radius
        assert abs(m.extents[2] - 3.0) < 1e-6  # height

    def test_invalid_rejected(self):
        with pytest.raises(ParameterError, match="radius"):
            cylinder(radius=0)
        with pytest.raises(ParameterError, match="height"):
            cylinder(height=-1)
        with pytest.raises(ParameterError, match="sections"):
            cylinder(sections=2)


class TestTorus:
    def test_default(self):
        m = torus()
        assert m.is_watertight
        # Default r_major=1, r_minor=0.25 -> outer radius 1.25
        assert abs(m.bounds[1, 0] - 1.25) < 1e-2

    def test_genus_1(self):
        # A ring torus has Euler number 0 (genus 1).
        m = torus(r_major=2.0, r_minor=0.5)
        assert m.euler_number == 0

    def test_minor_geq_major_rejected(self):
        with pytest.raises(ParameterError, match="self-non-intersecting"):
            torus(r_major=1.0, r_minor=1.0)

    def test_invalid_rejected(self):
        with pytest.raises(ParameterError, match="r_major"):
            torus(r_major=-1.0)
        with pytest.raises(ParameterError, match="sections"):
            torus(major_sections=2)


class TestBox:
    def test_default_unit_cube(self):
        m = box()
        assert m.is_watertight
        assert np.allclose(m.extents, [1.0, 1.0, 1.0])

    def test_extents(self):
        m = box(extents=(2.0, 3.0, 5.0))
        assert np.allclose(m.extents, [2.0, 3.0, 5.0])

    def test_invalid_rejected(self):
        with pytest.raises(ParameterError, match="extents"):
            box(extents=(1.0, -1.0, 1.0))
        with pytest.raises(ParameterError, match="extents"):
            box(extents=(1.0, 1.0))  # type: ignore[arg-type]


class TestCapsule:
    def test_default_watertight(self):
        m = capsule()
        assert m.is_watertight

    def test_total_length_includes_caps(self):
        # Length along z = height + 2 * radius (both caps).
        m = capsule(radius=0.5, height=1.0)
        # trimesh.creation.capsule places the cylinder's bottom at z=0
        # with hemispherical caps, so total z-extent ~= height + 2*radius.
        assert abs(m.extents[2] - (1.0 + 2 * 0.5)) < 1e-6


# ---------------------------------------------------------------------------
# extrude_section
# ---------------------------------------------------------------------------


class TestExtrudeSection:
    def test_naca_open_ends(self):
        sec = naca("0012", n_pts=40)
        wing = extrude_section(sec, span=1.0, n_stations=5)
        # No end caps -> not watertight by design.
        assert not wing.is_watertight
        # Span direction (z) extent matches.
        assert abs(wing.extents[2] - 1.0) < 1e-6

    def test_naca_with_end_caps_watertight(self):
        sec = naca("0012", n_pts=40)
        wing = extrude_section(sec, span=0.5, n_stations=8, end_caps=True)
        assert wing.is_watertight

    def test_taper(self):
        sec = naca("0012", n_pts=40)
        # taper=0.5: tip chord half of root chord.
        wing = extrude_section(sec, span=1.0, n_stations=5, taper=0.5)
        # Root chord at z=0 should span x ~= [0, 1].
        # Tip chord at z=1 should span x ~= [0, 0.5].
        z0_pts = wing.vertices[wing.vertices[:, 2] < 1e-6]
        z1_pts = wing.vertices[wing.vertices[:, 2] > 1.0 - 1e-6]
        root_chord = z0_pts[:, 0].max() - z0_pts[:, 0].min()
        tip_chord = z1_pts[:, 0].max() - z1_pts[:, 0].min()
        assert abs(root_chord - 1.0) < 1e-3
        assert abs(tip_chord - 0.5) < 1e-3

    def test_sweep(self):
        sec = naca("0012", n_pts=40)
        wing = extrude_section(sec, span=1.0, n_stations=5, sweep_deg=45.0)
        # Sweep 45 deg -> tip LE shifted by tan(45)*1 = 1 along +x.
        z0_pts = wing.vertices[wing.vertices[:, 2] < 1e-6]
        z1_pts = wing.vertices[wing.vertices[:, 2] > 1.0 - 1e-6]
        root_le = z0_pts[:, 0].min()
        tip_le = z1_pts[:, 0].min()
        assert abs((tip_le - root_le) - 1.0) < 1e-3

    def test_twist(self):
        sec = naca("2412", n_pts=40)
        # Twist 10 deg at the tip, rotation about the LE.  The TE sits
        # one chord aft of the LE, so its vertical displacement under
        # twist is sin(10°) * 1 ≈ 0.1736.
        wing = extrude_section(sec, span=1.0, n_stations=5, twist_deg=10.0)
        z0_pts = wing.vertices[wing.vertices[:, 2] < 1e-6]
        z1_pts = wing.vertices[wing.vertices[:, 2] > 1.0 - 1e-6]
        root_y_at_te = z0_pts[np.argmax(z0_pts[:, 0])][1]
        tip_y_at_te = z1_pts[np.argmax(z1_pts[:, 0])][1]
        assert abs((tip_y_at_te - root_y_at_te) - np.sin(np.radians(10.0))) < 1e-3

    def test_invalid_section_shape_rejected(self):
        with pytest.raises(ParameterError, match="shape"):
            extrude_section(np.zeros((10,)), span=1.0, n_stations=5)

    def test_invalid_n_stations_rejected(self):
        sec = naca("0012", n_pts=40)
        with pytest.raises(ParameterError, match="n_stations"):
            extrude_section(sec, span=1.0, n_stations=1)

    def test_invalid_span_rejected(self):
        sec = naca("0012", n_pts=40)
        with pytest.raises(ParameterError, match="span"):
            extrude_section(sec, span=-1.0, n_stations=5)

    def test_invalid_taper_rejected(self):
        sec = naca("0012", n_pts=40)
        with pytest.raises(ParameterError, match="taper"):
            extrude_section(sec, span=1.0, n_stations=5, taper=0.0)

    def test_invalid_span_axis_rejected(self):
        sec = naca("0012", n_pts=40)
        with pytest.raises(ParameterError, match="span_axis"):
            extrude_section(sec, span=1.0, n_stations=5, span_axis="w")

    def test_returns_trimesh(self):
        sec = naca("0012", n_pts=40)
        wing = extrude_section(sec, span=1.0, n_stations=5)
        assert isinstance(wing, trimesh.Trimesh)

    def test_span_axis_x(self):
        sec = naca("0012", n_pts=40)
        wing = extrude_section(sec, span=1.0, n_stations=5, span_axis="x")
        assert abs(wing.extents[0] - 1.0) < 1e-6
