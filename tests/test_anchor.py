"""Tests for shore.surface.anchor — Phase 1."""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import MeshInputError, ParameterError
from shore.surface.anchor import resolve_anchor


@pytest.fixture(scope="module")
def unit_sphere() -> trimesh.Trimesh:
    return trimesh.creation.icosphere(subdivisions=3, radius=1.0)


@pytest.fixture(scope="module")
def ellipsoid() -> trimesh.Trimesh:
    """Axis-aligned ellipsoid (a=2, b=1, c=1) — star-shaped, non-spherical."""
    sphere = trimesh.creation.icosphere(subdivisions=3, radius=1.0)
    verts = sphere.vertices.copy()
    verts[:, 0] *= 2.0
    return trimesh.Trimesh(vertices=verts, faces=sphere.faces, process=False)


# ---------------------------------------------------------------------------
# Strategy: centroid
# ---------------------------------------------------------------------------

class TestCentroid:
    def test_centroid_inside_sphere(self, unit_sphere):
        anchor = resolve_anchor(unit_sphere, strategy="centroid")
        assert anchor.shape == (3,)
        assert np.linalg.norm(anchor) < 0.5

    def test_centroid_is_near_origin_for_sphere(self, unit_sphere):
        anchor = resolve_anchor(unit_sphere, strategy="centroid")
        assert anchor == pytest.approx([0.0, 0.0, 0.0], abs=0.1)


# ---------------------------------------------------------------------------
# Strategy: bbox
# ---------------------------------------------------------------------------

class TestBbox:
    def test_bbox_inside_sphere(self, unit_sphere):
        anchor = resolve_anchor(unit_sphere, strategy="bbox")
        assert anchor.shape == (3,)
        inside = unit_sphere.contains(anchor.reshape(1, 3))
        assert bool(inside[0])

    def test_bbox_near_centre(self, unit_sphere):
        anchor = resolve_anchor(unit_sphere, strategy="bbox")
        assert anchor == pytest.approx([0.0, 0.0, 0.0], abs=0.1)

    def test_bbox_ellipsoid(self, ellipsoid):
        anchor = resolve_anchor(ellipsoid, strategy="bbox")
        inside = ellipsoid.contains(anchor.reshape(1, 3))
        assert bool(inside[0])


# ---------------------------------------------------------------------------
# Strategy: incenter
# ---------------------------------------------------------------------------

class TestIncenter:
    def test_incenter_inside_sphere(self, unit_sphere):
        anchor = resolve_anchor(unit_sphere, strategy="incenter")
        inside = unit_sphere.contains(anchor.reshape(1, 3))
        assert bool(inside[0])

    def test_incenter_near_centre(self, unit_sphere):
        anchor = resolve_anchor(unit_sphere, strategy="incenter")
        assert anchor == pytest.approx([0.0, 0.0, 0.0], abs=0.3)


# ---------------------------------------------------------------------------
# Strategy: auto
# ---------------------------------------------------------------------------

class TestAuto:
    def test_auto_returns_inside_point(self, unit_sphere):
        anchor = resolve_anchor(unit_sphere, strategy="auto")
        inside = unit_sphere.contains(anchor.reshape(1, 3))
        assert bool(inside[0])

    def test_auto_agrees_with_bbox_on_sphere(self, unit_sphere):
        auto = resolve_anchor(unit_sphere, strategy="auto")
        assert auto == pytest.approx([0.0, 0.0, 0.0], abs=0.3)

    def test_auto_ellipsoid_valid(self, ellipsoid):
        anchor = resolve_anchor(ellipsoid, strategy="auto")
        inside = ellipsoid.contains(anchor.reshape(1, 3))
        assert bool(inside[0])


# ---------------------------------------------------------------------------
# Strategy: user
# ---------------------------------------------------------------------------

class TestUser:
    def test_user_valid_point(self, unit_sphere):
        pt = np.array([0.0, 0.0, 0.0])
        anchor = resolve_anchor(unit_sphere, strategy="user", user_point=pt)
        assert anchor == pytest.approx(pt)

    def test_user_outside_raises(self, unit_sphere):
        pt = np.array([2.0, 0.0, 0.0])
        with pytest.raises(MeshInputError, match="outside"):
            resolve_anchor(unit_sphere, strategy="user", user_point=pt)

    def test_user_without_point_raises_parameter_error(self, unit_sphere):
        with pytest.raises(ParameterError, match="user_point"):
            resolve_anchor(unit_sphere, strategy="user")


# ---------------------------------------------------------------------------
# Unknown strategy
# ---------------------------------------------------------------------------

class TestUnknown:
    def test_unknown_strategy_raises(self, unit_sphere):
        with pytest.raises(ParameterError, match="Unknown anchor strategy"):
            resolve_anchor(unit_sphere, strategy="magic")


# ---------------------------------------------------------------------------
# Multi-strategy agreement on sphere
# ---------------------------------------------------------------------------

class TestStrategyAgreement:
    def test_all_strategies_agree_on_sphere(self, unit_sphere):
        """For a symmetric unit sphere all strategies should give near-origin."""
        for strat in ("centroid", "bbox", "incenter", "auto"):
            anchor = resolve_anchor(unit_sphere, strategy=strat)
            assert anchor == pytest.approx([0.0, 0.0, 0.0], abs=0.3), (
                f"strategy={strat!r} gave anchor far from origin: {anchor}"
            )
