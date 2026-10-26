"""Tests for shore.geometry.occt — OCCT backend of the Geometry protocol (ADR 0002).

Skipped when build123d (``shore-mesh[cad]``) is not installed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

bd = pytest.importorskip("build123d")

from shore.exceptions import MeshInputError, ParameterError  # noqa: E402
from shore.geometry import (  # noqa: E402
    BODY,
    EntityKind,
    EntityRef,
    Geometry,
    load_geometry,
)
from shore.geometry.occt import OcctGeometry  # noqa: E402

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _write_iges(shape: object, path: Path) -> None:
    from OCP.IGESControl import IGESControl_Writer

    writer = IGESControl_Writer()
    writer.AddShape(shape.wrapped)
    writer.ComputeModel()
    assert writer.Write(str(path))


@pytest.fixture(scope="module")
def cad_dir(tmp_path_factory) -> Path:
    d = tmp_path_factory.mktemp("cad")
    bd.export_step(bd.Box(1, 1, 1), str(d / "box.step"))
    bd.export_brep(bd.Box(1, 1, 1), str(d / "box.brep"))
    _write_iges(bd.Box(1, 1, 1), d / "box.igs")
    bd.export_step(bd.Sphere(1), str(d / "sphere.step"))
    return d


@pytest.fixture(scope="module")
def box() -> OcctGeometry:
    return OcctGeometry(bd.Box(1, 1, 1).wrapped)


@pytest.fixture(scope="module")
def sphere() -> OcctGeometry:
    return OcctGeometry(bd.Sphere(1).wrapped)


@pytest.fixture(scope="module")
def far_sphere() -> OcctGeometry:
    return OcctGeometry((bd.Pos(1.0e3, 0, 0) * bd.Sphere(1)).wrapped)


@pytest.fixture(scope="module")
def cylinder() -> OcctGeometry:
    return OcctGeometry(bd.Cylinder(1, 2).wrapped)


def _dirs(n: int, seed: int = 0) -> np.ndarray:
    v = np.random.default_rng(seed).normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _face_with_normal(g: OcctGeometry, axis: np.ndarray) -> EntityRef:
    """The planar box face whose outward normal is *axis*."""
    for ref in g.entities(EntityKind.SURFACE):
        centre = g.project(np.zeros((1, 3)) + 2.0 * axis, ref, exact=True)
        if np.allclose(centre.normals[0], axis):
            return ref
    raise AssertionError(f"no face with normal {axis}")


# ---------------------------------------------------------------------------
# Loading and topology
# ---------------------------------------------------------------------------


class TestLoad:
    @pytest.mark.parametrize("name", ["box.step", "box.brep", "box.igs"])
    def test_formats(self, cad_dir, name):
        g = load_geometry(cad_dir / name)
        assert isinstance(g, OcctGeometry)
        assert isinstance(g, Geometry)
        assert len(g.entities(EntityKind.SURFACE)) == 6
        assert len(g.entities(EntityKind.CURVE)) == 12
        assert len(g.entities(EntityKind.POINT)) == 8
        np.testing.assert_allclose(g.bounds, [[-0.5] * 3, [0.5] * 3], atol=1e-7)

    def test_body_entity(self, box):
        assert box.entities(EntityKind.BODY) == (BODY,)

    def test_missing_file(self, tmp_path):
        with pytest.raises(MeshInputError, match="not found"):
            OcctGeometry.from_file(tmp_path / "nope.step")

    def test_corrupt_step(self, tmp_path):
        bad = tmp_path / "bad.step"
        bad.write_text("not a step file\n")
        with pytest.raises(MeshInputError):
            load_geometry(bad)

    @pytest.mark.parametrize(
        "part",
        [lambda: bd.Sphere(1), lambda: bd.fillet(bd.Box(1, 1, 1).edges(), 0.1)],
        ids=["sphere", "filleted-box"],
    )
    def test_merged_proxy_watertight(self, part):
        # pole / corner-patch fans collapse on merging; degenerate faces must go
        merged = OcctGeometry(part().wrapped).merged_proxy
        assert merged.is_watertight
        assert merged.volume > 0.0

    def test_deflection_default_relative(self, sphere):
        assert sphere.deflection == pytest.approx(1e-3 * sphere.diag)


# ---------------------------------------------------------------------------
# Body projection
# ---------------------------------------------------------------------------


class TestBodyProjection:
    def test_exact_on_sphere(self, sphere):
        proj = sphere.project(2.0 * _dirs(200), exact=True)
        np.testing.assert_allclose(np.linalg.norm(proj.points, axis=-1), 1.0, atol=1e-12)
        np.testing.assert_allclose(proj.distances, 1.0, atol=1e-12)
        np.testing.assert_allclose(np.einsum("ij,ij->i", proj.normals, proj.points), 1.0, atol=1e-9)

    def test_exact_far_from_origin(self, far_sphere):
        c = np.array([1.0e3, 0.0, 0.0])
        proj = far_sphere.project(c + 2.0 * _dirs(100, seed=1), exact=True)
        np.testing.assert_allclose(np.linalg.norm(proj.points - c, axis=-1), 1.0, atol=1e-9)

    def test_tier1_near_surface_within_deflection(self, sphere):
        on = sphere.project(2.0 * _dirs(500, seed=2), exact=True).points
        q = on * (1.0 + 0.5 * sphere.deflection)  # just off the surface
        tier1 = sphere.project(q)
        exact = sphere.project(q, exact=True)
        err = np.linalg.norm(tier1.points - exact.points, axis=-1)
        assert err.max() <= 2.0 * sphere.deflection

    def test_normals_outward_on_reversed_faces(self, box):
        # a box has faces of both orientations; every normal must point outward
        q = 2.0 * _dirs(300, seed=3)
        for exact in (False, True):
            proj = box.project(q, exact=exact)
            assert np.all(np.einsum("ij,ij->i", proj.normals, proj.points) > 0.0)

    def test_box_exact_and_axis_normals(self, box):
        q = np.array([[2.0, 0.1, -0.2], [0.1, -3.0, 0.2], [0.0, 0.2, 5.0]])
        proj = box.project(q, exact=True)
        np.testing.assert_allclose(
            proj.points, [[0.5, 0.1, -0.2], [0.1, -0.5, 0.2], [0.0, 0.2, 0.5]]
        )
        np.testing.assert_allclose(proj.normals, [[1, 0, 0], [0, -1, 0], [0, 0, 1]], atol=1e-12)

    def test_entity_index_names_the_face(self, box):
        q = np.array([[2.0, 0.1, -0.2]])
        body = box.project(q)
        face = EntityRef(EntityKind.SURFACE, int(body.entity_index[0]))
        np.testing.assert_allclose(box.project(q, face).points, body.points)
        assert body.params is None

    def test_leading_shape_preserved(self, sphere):
        proj = sphere.project(2.0 * _dirs(12, seed=4).reshape(3, 4, 3), exact=True)
        assert proj.points.shape == (3, 4, 3)
        assert proj.distances.shape == (3, 4)
        assert proj.entity_index.shape == (3, 4)


# ---------------------------------------------------------------------------
# Surface and curve entities
# ---------------------------------------------------------------------------


class TestSurfaceEntity:
    def test_params_evaluate_round_trip(self, sphere):
        face = EntityRef(EntityKind.SURFACE, 0)
        proj = sphere.project(2.0 * _dirs(50, seed=5), face, exact=True)
        assert proj.params.shape == (50, 2)
        np.testing.assert_allclose(sphere.evaluate(face, proj.params), proj.points, atol=1e-12)

    def test_trimmed_face_falls_back_to_boundary(self, box):
        # foot of (2, 2, 0) on the +x plane is (0.5, 2, 0): outside the trimmed
        # face, so the exact answer is the nearest boundary point (0.5, 0.5, 0)
        face = _face_with_normal(box, np.array([1.0, 0.0, 0.0]))
        proj = box.project(np.array([[2.0, 2.0, 0.0]]), face, exact=True)
        np.testing.assert_allclose(proj.points[0], [0.5, 0.5, 0.0], atol=1e-9)
        np.testing.assert_allclose(box.evaluate(face, proj.params)[0], proj.points[0], atol=1e-9)

    def test_tier1_params_near_exact(self, sphere):
        face = EntityRef(EntityKind.SURFACE, 0)
        on = sphere.project(2.0 * _dirs(100, seed=6), face, exact=True).points
        tier1 = sphere.project(on, face)
        back = sphere.evaluate(face, tier1.params)
        # interpolated (u, v) land within a few deflections, poles included
        assert np.linalg.norm(back - on, axis=-1).max() <= 4.0 * sphere.deflection

    def test_tier1_params_at_pole(self, sphere):
        # u is singular at the poles; its arbitrary pole value must not leak in
        face = EntityRef(EntityKind.SURFACE, 0)
        theta = np.linspace(0.0, 2.0 * np.pi, 64, endpoint=False)
        ring = np.stack([0.03 * np.cos(theta), 0.03 * np.sin(theta), np.ones(64)], axis=1)
        on = sphere.project(ring, face, exact=True).points
        back = sphere.evaluate(face, sphere.project(on, face).params)
        assert np.linalg.norm(back - on, axis=-1).max() <= 4.0 * sphere.deflection


class TestCurveEntity:
    def test_straight_edge_exact(self, box):
        q = np.array([[3.0, 0.1, 2.0]])
        for ref in box.entities(EntityKind.CURVE):
            exact = box.project(q, ref, exact=True)
            tier1 = box.project(q, ref)
            np.testing.assert_allclose(tier1.points, exact.points, atol=1e-12)
            np.testing.assert_allclose(box.evaluate(ref, exact.params), exact.points, atol=1e-12)
            assert np.all(np.isnan(exact.normals))

    def test_circular_edge_exact(self, cylinder):
        q = 3.0 * _dirs(40, seed=7)
        rims = 0
        for ref in cylinder.entities(EntityKind.CURVE):
            exact = cylinder.project(q, ref, exact=True)
            r = np.linalg.norm(exact.points[:, :2], axis=-1)
            if np.allclose(np.abs(exact.points[:, 2]), 1.0):  # the two rim circles
                rims += 1
                np.testing.assert_allclose(r, 1.0, atol=1e-10)
                # tier 1 is specified for queries near the entity
                near = exact.points * (1.0 + 0.5 * cylinder.deflection)
                tier1 = cylinder.project(near, ref)
                err = np.linalg.norm(tier1.points - exact.points, axis=-1)
                assert err.max() <= 2.0 * cylinder.deflection
        assert rims == 2

    def test_vertex(self, box):
        ref = EntityRef(EntityKind.POINT, 0)
        proj = box.project(np.zeros((2, 3)), ref)
        np.testing.assert_allclose(np.abs(proj.points), 0.5)
        np.testing.assert_allclose(box.evaluate(ref, np.zeros(2)), proj.points)


# ---------------------------------------------------------------------------
# Ray casting and errors
# ---------------------------------------------------------------------------


class TestRayCastAndErrors:
    def test_ray_cast_sphere(self, sphere):
        res = sphere.ray_cast(np.zeros(3), _dirs(300, seed=8))
        assert res.all_hit
        r = np.linalg.norm(res.locations, axis=-1)
        assert np.all(r <= 1.0 + 1e-12)
        # BRepMesh bounds edge sag by the deflection; interiors sag up to ~2x
        assert np.all(r >= 1.0 - 2.0 * sphere.deflection)

    def test_entity_out_of_range(self, box):
        with pytest.raises(ParameterError, match="out of range"):
            box.project(np.zeros((1, 3)), EntityRef(EntityKind.SURFACE, 6))

    def test_evaluate_body_rejected(self, box):
        with pytest.raises(ParameterError):
            box.evaluate(BODY, np.zeros((1, 2)))
