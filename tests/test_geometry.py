"""Tests for shore.geometry — Geometry protocol and the STL backend (ADR 0002)."""

from __future__ import annotations

import sys

import numpy as np
import pytest
import trimesh

from shore.exceptions import MeshInputError, ParameterError
from shore.geometry import (
    BODY,
    EntityKind,
    EntityRef,
    Geometry,
    Projection,
    StlGeometry,
    load_geometry,
)
from shore.mesh.body import Body
from shore.surface.projector import Projector
from shore.volume.hyperbolic import _stl_normals

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def unit_sphere() -> trimesh.Trimesh:
    return trimesh.creation.icosphere(subdivisions=3, radius=1.0)


@pytest.fixture(scope="module")
def geom(unit_sphere) -> StlGeometry:
    return StlGeometry(unit_sphere)


def _radial_dirs(n: int, seed: int = 0) -> np.ndarray:
    v = np.random.default_rng(seed).normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


# ---------------------------------------------------------------------------
# Protocol conformance
# ---------------------------------------------------------------------------


class TestProtocol:
    def test_stl_geometry_is_a_geometry(self, geom):
        assert isinstance(geom, Geometry)

    def test_unrelated_object_is_not_a_geometry(self):
        assert not isinstance(object(), Geometry)

    def test_entities_body_only(self, geom):
        assert geom.entities(EntityKind.BODY) == (BODY,)
        for kind in (EntityKind.SURFACE, EntityKind.CURVE, EntityKind.POINT):
            assert geom.entities(kind) == ()

    def test_bounds_shape_and_copy(self, geom, unit_sphere):
        b = geom.bounds
        assert b.shape == (2, 3)
        np.testing.assert_array_equal(b, unit_sphere.bounds)
        b[0, 0] = 99.0  # caller mutation must not leak into the geometry
        assert geom.bounds[0, 0] != 99.0


# ---------------------------------------------------------------------------
# Ray casting: bit-identical to the legacy Projector
# ---------------------------------------------------------------------------


class TestRayCast:
    def test_matches_projector_bit_for_bit(self, geom, unit_sphere):
        dirs = _radial_dirs(200)
        legacy = Projector(unit_sphere, np.zeros(3)).project(dirs)
        new = geom.ray_cast(np.zeros(3), dirs)
        np.testing.assert_array_equal(new.locations, legacy.locations)
        np.testing.assert_array_equal(new.distances, legacy.distances)
        assert new.all_hit

    def test_projector_cached_per_anchor(self, geom):
        a = geom.projector(np.zeros(3))
        assert geom.projector(np.zeros(3)) is a
        assert geom.projector(np.array([0.1, 0.0, 0.0])) is not a

    def test_on_failure_forwarded(self, geom):
        # anchor outside the sphere, pointing away: every ray misses
        res = geom.ray_cast(np.array([5.0, 0.0, 0.0]), np.array([[1.0, 0.0, 0.0]]), "skip")
        assert not res.all_hit
        assert np.all(np.isnan(res.locations))
        with pytest.raises(MeshInputError):
            geom.ray_cast(np.array([5.0, 0.0, 0.0]), np.array([[1.0, 0.0, 0.0]]))


# ---------------------------------------------------------------------------
# Closest-point projection
# ---------------------------------------------------------------------------


class TestProject:
    def test_radial_points_land_on_surface(self, geom):
        dirs = _radial_dirs(100, seed=1)
        proj = geom.project(2.0 * dirs)
        assert isinstance(proj, Projection)
        r = np.linalg.norm(proj.points, axis=-1)
        # subdivision-3 icosphere: facet sag below 1 %
        assert np.all(r > 0.99)
        assert np.all(r <= 1.0 + 1e-12)
        # closest point is not radially below the query on a faceted sphere
        np.testing.assert_allclose(
            proj.distances, np.linalg.norm(2.0 * dirs - proj.points, axis=-1), atol=1e-12
        )
        assert np.all(proj.distances >= 2.0 - r - 1e-12)
        assert proj.params is None

    def test_normals_unit_and_outward(self, geom):
        dirs = _radial_dirs(100, seed=2)
        proj = geom.project(1.5 * dirs)
        np.testing.assert_allclose(np.linalg.norm(proj.normals, axis=-1), 1.0, atol=1e-12)
        assert np.all(np.einsum("ij,ij->i", proj.normals, dirs) > 0.99)

    def test_matches_trimesh_closest_point(self, geom, unit_sphere):
        pts = 0.5 * _radial_dirs(50, seed=3)  # inside the body too
        ref, ref_d, _ = trimesh.proximity.closest_point(unit_sphere, pts)
        proj = geom.project(pts)
        np.testing.assert_array_equal(proj.points, ref)
        np.testing.assert_array_equal(proj.distances, ref_d)

    def test_normals_match_stl_normals_without_creases(self, geom, unit_sphere):
        # no edge of the icosphere is sharper than 30 deg: the crease-aware
        # normals (ADR 0014) equal the legacy vertex normals here
        pts = 1.2 * _radial_dirs(50, seed=4)
        proj = geom.project(pts)
        np.testing.assert_allclose(
            proj.normals, _stl_normals(unit_sphere, pts, interpolate=True), atol=1e-14
        )

    def test_leading_shape_preserved(self, geom):
        pts = 2.0 * _radial_dirs(24, seed=5).reshape(4, 6, 3)
        proj = geom.project(pts)
        assert proj.points.shape == (4, 6, 3)
        assert proj.distances.shape == (4, 6)
        assert proj.normals.shape == (4, 6, 3)

    def test_exact_flag_is_accepted(self, geom):
        pts = 2.0 * _radial_dirs(10, seed=6)
        np.testing.assert_array_equal(
            geom.project(pts, exact=True).points, geom.project(pts).points
        )

    def test_body_entity_accepted(self, geom):
        pts = 2.0 * _radial_dirs(5, seed=7)
        np.testing.assert_array_equal(geom.project(pts, BODY).points, geom.project(pts).points)

    def test_other_entity_rejected(self, geom):
        with pytest.raises(ParameterError):
            geom.project(np.zeros((1, 3)), EntityRef(EntityKind.SURFACE, 0))

    def test_evaluate_rejected(self, geom):
        with pytest.raises(ParameterError):
            geom.evaluate(BODY, np.zeros((1, 2)))


# ---------------------------------------------------------------------------
# Body integration and loading
# ---------------------------------------------------------------------------


class TestBodyIntegration:
    def test_body_geometry_is_stl_geometry(self, unit_sphere):
        body = Body(mesh=unit_sphere)
        assert isinstance(body.geometry, StlGeometry)
        assert body.geometry is body.geometry
        assert body.geometry.mesh is unit_sphere

    def test_body_projector_shared_with_geometry(self, unit_sphere):
        body = Body(mesh=unit_sphere)
        assert body.projector is body.geometry.projector(body.anchor)
        assert body.projector is body.projector


class TestLoadGeometry:
    def test_stl(self, sphere_stl):
        g = load_geometry(sphere_stl)
        assert isinstance(g, StlGeometry)
        assert isinstance(g, Geometry)

    def test_stl_suffix_case_insensitive(self, sphere_stl, tmp_path):
        upper = tmp_path / "SPHERE.STL"
        upper.write_bytes(sphere_stl.read_bytes())
        assert isinstance(load_geometry(upper), StlGeometry)

    @pytest.mark.parametrize("suffix", [".step", ".stp", ".iges", ".igs", ".brep"])
    def test_cad_without_occt_backend(self, tmp_path, suffix, monkeypatch):
        # simulate a missing [cad] extra: importing the backend fails
        monkeypatch.setitem(sys.modules, "shore.geometry.occt", None)
        with pytest.raises(MeshInputError, match=r"shore-mesh\[cad\]"):
            load_geometry(tmp_path / f"wing{suffix}")

    def test_unknown_suffix(self, tmp_path):
        with pytest.raises(MeshInputError, match="unknown geometry format"):
            load_geometry(tmp_path / "wing.obj")


class TestCornerNormals:
    """Crease-aware STL normals (ADR 0014) beyond the frozen spec."""

    def test_chunking_does_not_change_the_result(self, monkeypatch):
        import shore.geometry.stl as stl

        mesh = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=24)
        whole = stl._corner_normals(mesh, 30.0)
        monkeypatch.setattr(stl, "_CHUNK", 7)
        np.testing.assert_array_equal(stl._corner_normals(mesh, 30.0), whole)

    def test_planar_boundary_found_and_line_skipped(self):
        from shore.geometry.stl import _planar_boundaries

        cyl = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=24)
        cyl.update_faces(cyl.triangles_center[:, 1] > 0.0)
        cyl.remove_unreferenced_vertices()
        [(_, normal)] = _planar_boundaries(cyl)
        np.testing.assert_allclose(np.abs(normal), [0.0, 1.0, 0.0], atol=1e-12)
        assert _planar_boundaries(trimesh.creation.icosphere(2)) == []

    def test_legacy_marching_normals_untouched(self):
        # the legacy shore mesh path keeps trimesh's vertex normals (ADR 0014)
        mesh = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
        points = np.array([[0.5, 0.3, 0.5]])
        legacy = _stl_normals(mesh, points)
        assert abs(legacy[0, 2] - 1.0) > 1e-3  # the corner-averaged tilt is still there
        np.testing.assert_allclose(StlGeometry(mesh).project(points).normals, [[0, 0, 1]])


class TestBodyNormalsForBlending:
    """The marcher's normal blending keeps trimesh's vertex normals on STL (ADR 0015)."""

    def test_stl_body_normals_are_the_vertex_normal_blend(self):
        mesh = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
        points = trimesh.sample.sample_surface(mesh, 500, seed=1)[0].reshape(25, 20, 3)
        normals = Body(mesh=mesh).normals(points)
        assert normals.shape == (25, 20, 3)
        np.testing.assert_array_equal(
            normals.reshape(-1, 3), _stl_normals(mesh, points.reshape(-1, 3))
        )

    def test_blended_cubed_sphere_march_across_a_sharp_edge(self, tmp_path):
        # with the crease-aware normals the first layer folded at the cylinder's rim
        import warnings

        from shore.mesh.mesh import Mesh
        from shore.mesh.topology import CubedSphereTopology

        trimesh.creation.cylinder(radius=0.5, height=1.6, sections=96).export(
            str(tmp_path / "cyl.stl")
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mesh = Mesh.from_stl(tmp_path / "cyl.stl", CubedSphereTopology(), ni=33, nj=64, nk=4)
            mesh.march(ds=0.01, growth=1.15, blend_normals_k=5)
        assert len(mesh.blocks) == 6


class TestProjectCache:
    """Closest points remembered per query point (ADR 0032): same answers bit for bit."""

    def test_batches_and_repeats_are_bit_identical(self, unit_sphere):
        rng = np.random.default_rng(3)
        pts = rng.normal(size=(200, 3)) * 1.5
        pts[50:60] = pts[0:10]  # repeats inside one batch
        cold = StlGeometry(unit_sphere).project(pts)
        warm = StlGeometry(unit_sphere)
        order = rng.permutation(len(pts))
        warm.project(pts[order[:120]])  # fills part of the cache
        again = warm.project(pts)
        for field in ("points", "distances", "normals"):
            assert np.array_equal(getattr(cold, field), getattr(again, field))
        warm._closest.clear()
        assert np.array_equal(warm.project(pts).points, cold.points)

    def test_cache_is_bounded(self, unit_sphere, monkeypatch):
        g = StlGeometry(unit_sphere)
        monkeypatch.setattr(StlGeometry, "_CACHE_SIZE", 50)
        g.project(np.random.default_rng(0).normal(size=(40, 3)))
        g.project(np.random.default_rng(1).normal(size=(40, 3)))
        assert len(g._closest) <= 50
