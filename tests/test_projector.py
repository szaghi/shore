"""Tests for shore.surface.projector — Phase 0 + Phase 2."""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import MeshInputError
from shore.surface.projector import ProjectionResult, ProjectionStatus, Projector

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def unit_sphere() -> trimesh.Trimesh:
    return trimesh.creation.icosphere(subdivisions=3, radius=1.0)


@pytest.fixture(scope="module")
def sphere_proj(unit_sphere) -> Projector:
    """Projector with anchor at origin (inside unit sphere)."""
    return Projector(unit_sphere, np.zeros(3))


@pytest.fixture(scope="module")
def inverted_sphere(unit_sphere) -> trimesh.Trimesh:
    """Unit sphere with all face windings reversed (normals pointing inward)."""
    inv = unit_sphere.copy()
    inv.faces = unit_sphere.faces[:, ::-1]
    _ = inv.face_normals  # force normal recompute
    return inv


# ---------------------------------------------------------------------------
# Basic hit behaviour
# ---------------------------------------------------------------------------


class TestProjectorHits:
    def test_single_axis_hit(self, sphere_proj):
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = sphere_proj.project(dirs)
        assert result.all_hit
        assert result.locations[0] == pytest.approx([1.0, 0.0, 0.0], abs=0.05)

    def test_all_latlon_rays_hit(self, sphere_proj):
        ni, nj = 8, 12
        lats = np.linspace(0.2, np.pi - 0.2, ni)
        lons = np.linspace(0, 2 * np.pi, nj, endpoint=False)
        la, lo = np.meshgrid(lats, lons, indexing="ij")
        dirs = np.stack(
            [
                np.sin(la) * np.cos(lo),
                np.sin(la) * np.sin(lo),
                np.cos(la),
            ],
            axis=-1,
        ).reshape(-1, 3)
        result = sphere_proj.project(dirs)
        assert result.all_hit

    def test_result_shapes(self, sphere_proj):
        n = 20
        dirs = np.random.default_rng(42).standard_normal((n, 3))
        dirs /= np.linalg.norm(dirs, axis=-1, keepdims=True)
        result = sphere_proj.project(dirs)
        assert result.locations.shape == (n, 3)
        assert result.status.shape == (n,)
        assert result.distances.shape == (n,)

    def test_hit_mask_all_true(self, sphere_proj):
        dirs = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        result = sphere_proj.project(dirs)
        assert result.hit_mask.all()

    def test_hit_distances_near_one(self, sphere_proj):
        dirs = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = sphere_proj.project(dirs)
        assert np.all(result.distances == pytest.approx(1.0, abs=0.05))

    def test_hit_locations_on_unit_sphere(self, sphere_proj):
        rng = np.random.default_rng(7)
        dirs = rng.standard_normal((16, 3))
        dirs /= np.linalg.norm(dirs, axis=-1, keepdims=True)
        result = sphere_proj.project(dirs)
        radii = np.linalg.norm(result.locations, axis=-1)
        assert radii == pytest.approx(np.ones(16), abs=0.05)


# ---------------------------------------------------------------------------
# Missed rays
# ---------------------------------------------------------------------------


class TestProjectorMissed:
    def test_raises_by_default(self, unit_sphere):
        # Anchor outside: rays aimed away from sphere all miss
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        with pytest.raises(MeshInputError, match="rays failed"):
            proj.project(dirs)

    def test_nearest_fallback_returns_finite(self, unit_sphere):
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="nearest")
        assert np.all(np.isfinite(result.locations))

    def test_skip_leaves_nan(self, unit_sphere):
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="skip")
        assert result.status[0] == ProjectionStatus.MISSED
        assert np.all(np.isnan(result.locations[0]))

    def test_error_message_contains_count(self, unit_sphere):
        anchor = np.array([5.0, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        dirs = np.array([[1.0, 0.0, 0.0]])
        with pytest.raises(MeshInputError, match=r"\d+ of \d+ rays"):
            proj.project(dirs)


# ---------------------------------------------------------------------------
# Backface detection
# ---------------------------------------------------------------------------


class TestProjectorBackface:
    def test_inverted_mesh_produces_backface(self, inverted_sphere):
        # Anchor at origin, inverted mesh → normals point inward.
        # For inside-out casting (outward ray, inward normal), dot <= 0 → BACKFACE.
        proj = Projector(inverted_sphere, np.zeros(3))
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="skip")
        assert result.status[0] == ProjectionStatus.BACKFACE

    def test_normal_mesh_no_backface(self, sphere_proj):
        dirs = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = sphere_proj.project(dirs)
        for s in result.status:
            assert s == ProjectionStatus.HIT


# ---------------------------------------------------------------------------
# TOO_NEAR detection
# ---------------------------------------------------------------------------


class TestProjectorTooNear:
    def test_anchor_on_surface_triggers_too_near(self, unit_sphere):
        # Anchor just inside the surface; hit distance < 0.05 * bbox_diag
        anchor = np.array([0.999, 0.0, 0.0])
        proj = Projector(unit_sphere, anchor)
        # bbox_diag ≈ 2 * sqrt(3) ≈ 3.46, threshold = 0.05 * 3.46 ≈ 0.17
        # hit distance ≈ 1 - 0.999 = 0.001 < 0.17 → TOO_NEAR
        dirs = np.array([[1.0, 0.0, 0.0]])
        result = proj.project(dirs, on_failure="skip")
        assert result.status[0] == ProjectionStatus.TOO_NEAR


# ---------------------------------------------------------------------------
# ProjectionResult helpers
# ---------------------------------------------------------------------------


class TestProjectionResult:
    def test_all_hit_true_when_all_hit(self):
        n = 5
        locs = np.ones((n, 3))
        status = np.full(n, ProjectionStatus.HIT, dtype=object)
        dist = np.ones(n)
        r = ProjectionResult(locs, status, dist)
        assert r.all_hit

    def test_all_hit_false_when_any_miss(self):
        n = 5
        locs = np.ones((n, 3))
        status = np.full(n, ProjectionStatus.HIT, dtype=object)
        status[2] = ProjectionStatus.MISSED
        dist = np.ones(n)
        r = ProjectionResult(locs, status, dist)
        assert not r.all_hit

    def test_hit_mask_shape(self):
        n = 4
        locs = np.ones((n, 3))
        status = np.full(n, ProjectionStatus.HIT, dtype=object)
        dist = np.ones(n)
        r = ProjectionResult(locs, status, dist)
        assert r.hit_mask.shape == (n,)


# ---------------------------------------------------------------------------
# Phase 0 equivalence — ported functions produce same output
# ---------------------------------------------------------------------------


class TestEquivalence:
    def test_project_sphere_to_surface_shape(self, sphere_stl):
        import trimesh as _tm

        from shore.surface.remesh import project_sphere_to_surface

        mesh = _tm.load(str(sphere_stl), force="mesh")
        result = project_sphere_to_surface(mesh, ni=8, nj=12)
        assert result.shape == (8, 12, 3)

    def test_project_sphere_to_surface_on_unit_sphere(self, sphere_stl):
        import trimesh as _tm

        from shore.surface.remesh import project_sphere_to_surface

        mesh = _tm.load(str(sphere_stl), force="mesh")
        result = project_sphere_to_surface(mesh, ni=8, nj=12)
        radii = np.linalg.norm(result.reshape(-1, 3), axis=-1)
        assert radii.min() == pytest.approx(1.0, abs=0.05)
        assert radii.max() == pytest.approx(1.0, abs=0.05)

    def test_project_caps_to_surface_shape(self, sphere_stl):
        import trimesh as _tm

        from shore.surface.remesh import project_caps_to_surface

        mesh = _tm.load(str(sphere_stl), force="mesh")
        surface, nc_cap = project_caps_to_surface(mesh, ni=8, nj=12)
        assert surface.ndim == 3
        assert surface.shape[1] == 12
        assert nc_cap >= 2

    def test_project_caps_on_unit_sphere(self, sphere_stl):
        import trimesh as _tm

        from shore.surface.remesh import project_caps_to_surface

        mesh = _tm.load(str(sphere_stl), force="mesh")
        surface, _ = project_caps_to_surface(mesh, ni=8, nj=12)
        radii = np.linalg.norm(surface.reshape(-1, 3), axis=-1)
        assert radii.min() == pytest.approx(1.0, abs=0.05)

    def test_external_projector_gives_same_result(self, sphere_stl):
        """Passing a pre-built Projector should give identical output."""
        import trimesh as _tm

        from shore.surface.projector import Projector
        from shore.surface.remesh import project_sphere_to_surface

        mesh = _tm.load(str(sphere_stl), force="mesh")
        proj = Projector(mesh, mesh.centroid)
        result_with = project_sphere_to_surface(mesh, ni=8, nj=12, projector=proj)
        result_without = project_sphere_to_surface(mesh, ni=8, nj=12)
        np.testing.assert_allclose(result_with, result_without)
