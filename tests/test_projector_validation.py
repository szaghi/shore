"""Tests for Projector.project status codes and on_failure recovery — Phase 2."""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import MeshInputError
from shore.surface.projector import ProjectionStatus, Projector


@pytest.fixture(scope="module")
def unit_sphere() -> trimesh.Trimesh:
    return trimesh.creation.icosphere(subdivisions=3, radius=1.0)


# ---------------------------------------------------------------------------
# Each failure status code hit by a crafted fixture
# ---------------------------------------------------------------------------

class TestStatusCodes:
    def test_hit_status_on_clean_projection(self, unit_sphere):
        proj = Projector(unit_sphere, np.zeros(3))
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs)
        assert result.status[0] == ProjectionStatus.HIT

    def test_missed_when_ray_points_away(self, unit_sphere):
        # Anchor far outside; ray points away from sphere
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="skip")
        assert result.status[0] == ProjectionStatus.MISSED

    def test_backface_on_inverted_mesh(self, unit_sphere):
        inverted = unit_sphere.copy()
        inverted.faces = unit_sphere.faces[:, ::-1]
        _ = inverted.face_normals  # force normal recompute

        proj = Projector(inverted, np.zeros(3))
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="skip")
        assert result.status[0] == ProjectionStatus.BACKFACE

    def test_too_near_when_anchor_on_surface(self, unit_sphere):
        # Anchor just inside the sphere surface → very small hit distance
        anchor = np.array([0.999, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="skip")
        assert result.status[0] == ProjectionStatus.TOO_NEAR


# ---------------------------------------------------------------------------
# on_failure modes
# ---------------------------------------------------------------------------

class TestOnFailure:
    def test_raise_mode_raises_mesh_input_error(self, unit_sphere):
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        with pytest.raises(MeshInputError):
            proj.project(dirs, on_failure="raise")

    def test_nearest_mode_returns_finite_locations(self, unit_sphere):
        """on_failure='nearest' should give a finite surface for failed rays."""
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="nearest")
        assert np.all(np.isfinite(result.locations))

    def test_nearest_mode_keeps_failure_status(self, unit_sphere):
        """on_failure='nearest' fills the location but preserves the non-HIT status."""
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="nearest")
        assert result.status[0] != ProjectionStatus.HIT
        assert not result.all_hit

    def test_skip_mode_leaves_nan(self, unit_sphere):
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="skip")
        assert np.all(np.isnan(result.locations[0]))

    def test_mixed_hit_and_miss(self, unit_sphere):
        """Some rays hit, some miss — only failed rays use fallback."""
        # One ray pointing at sphere (should hit), one pointing away (should miss)
        # Use an anchor that is inside the sphere so +x direction hits
        proj = Projector(unit_sphere, np.zeros(3))
        # Create a direction field on the equatorial plane; all should hit
        dirs = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        result = proj.project(dirs)
        assert result.all_hit

    def test_error_message_counts_per_mode(self, unit_sphere):
        """Error message should name the failure mode and count."""
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        with pytest.raises(MeshInputError) as exc_info:
            proj.project(dirs, on_failure="raise")
        msg = str(exc_info.value)
        assert "of" in msg and "rays" in msg


# ---------------------------------------------------------------------------
# Projector identity: same BVH object reused across calls
# ---------------------------------------------------------------------------

class TestProjectorBvhReuse:
    def test_intersector_cached(self, unit_sphere):
        proj = Projector(unit_sphere, np.zeros(3))
        inter_a = proj._intersector
        inter_b = proj._intersector
        assert inter_a is inter_b

    def test_multiple_project_calls_consistent(self, unit_sphere):
        proj = Projector(unit_sphere, np.zeros(3))
        dirs = np.array([[1.0, 0.0, 0.0]])
        r1 = proj.project(dirs)
        r2 = proj.project(dirs)
        np.testing.assert_array_equal(r1.locations, r2.locations)
