"""CAD input through the meshing pipeline (ADR 0002, roadmap phase 1 acceptance).

A STEP sphere and a STEP NACA 0012 wing must reproduce the STL k=0 surfaces
to within the STL facet error. Skipped when build123d (``shore-mesh[cad]``)
is not installed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import trimesh
from typer.testing import CliRunner

bd = pytest.importorskip("build123d")

from shore.cli import app  # noqa: E402
from shore.geometry import AnchoredCaster, StlGeometry  # noqa: E402
from shore.geometry.occt import OcctGeometry  # noqa: E402
from shore.io.geo import read_geo  # noqa: E402
from shore.mesh import Body, CubedSphereTopology, Mesh, OGridTopology  # noqa: E402
from shore.surface.projector import Projector  # noqa: E402
from tests.conftest import _make_naca0012_stl, _naca0012_section  # noqa: E402

runner = CliRunner()

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def cad(tmp_path_factory) -> Path:
    d = tmp_path_factory.mktemp("cad_pipeline")
    bd.export_step(bd.Sphere(1), str(d / "sphere.step"))
    trimesh.creation.icosphere(subdivisions=3).export(str(d / "sphere.stl"))
    # the same polygonal NACA 0012 wing as the STL fixture, as a STEP solid
    section = _naca0012_section(80)
    face = bd.make_face(bd.Polyline(*[tuple(p) for p in section], close=True))
    bd.export_step(bd.extrude(face, amount=0.5), str(d / "wing.step"))
    _make_naca0012_stl(d / "wing.stl")
    return d


def _k0(mesh: Mesh) -> np.ndarray:
    return np.concatenate([b.nodes[:, :, 0].reshape(-1, 3) for b in mesh.blocks])


def _stl_sag(path: Path) -> float:
    """Largest inward deviation of a unit-sphere STL from the true sphere."""
    m = trimesh.load(path)
    return float(1.0 - np.linalg.norm(m.triangles_center, axis=1).min())


# ---------------------------------------------------------------------------
# Body
# ---------------------------------------------------------------------------


class TestBody:
    def test_from_file_step(self, cad):
        body = Body.from_file(cad / "sphere.step")
        assert body.is_cad
        assert isinstance(body.geometry, OcctGeometry)
        assert isinstance(body.projector, AnchoredCaster)
        assert body.projector.snap
        assert body.mesh.is_watertight  # merged tessellation for mesh consumers
        assert np.linalg.norm(body.anchor) < 0.1

    def test_from_file_stl_unchanged(self, cad):
        body = Body.from_file(cad / "sphere.stl")
        assert not body.is_cad
        assert isinstance(body.geometry, StlGeometry)
        assert isinstance(body.projector, Projector)

    def test_normals_exact_geometry(self, cad):
        body = Body.from_file(cad / "sphere.step")
        d = np.random.default_rng(0).normal(size=(50, 3))
        d /= np.linalg.norm(d, axis=1, keepdims=True)
        # tier-1 normals: interpolated from exact nodal normals
        assert np.all(np.einsum("ij,ij->i", body.normals(1.2 * d), d) > 0.999)


class TestAnchoredCaster:
    def test_no_snap_matches_projector(self):
        mesh = trimesh.creation.icosphere(subdivisions=3)
        d = np.random.default_rng(1).normal(size=(100, 3))
        d /= np.linalg.norm(d, axis=1, keepdims=True)
        ref = Projector(mesh, np.zeros(3)).project(d)
        got = AnchoredCaster(StlGeometry(mesh), np.zeros(3)).project(d)
        np.testing.assert_array_equal(got.locations, ref.locations)
        np.testing.assert_array_equal(got.distances, ref.distances)

    def test_snap_puts_hits_on_exact_surface(self, cad):
        g = OcctGeometry.from_file(cad / "sphere.step")
        d = np.random.default_rng(2).normal(size=(100, 3))
        d /= np.linalg.norm(d, axis=1, keepdims=True)
        raw = AnchoredCaster(g, np.zeros(3)).project(d)
        snapped = AnchoredCaster(g, np.zeros(3), snap=True).project(d)
        assert np.abs(np.linalg.norm(raw.locations, axis=1) - 1.0).max() > 1e-4
        np.testing.assert_allclose(np.linalg.norm(snapped.locations, axis=1), 1.0, atol=1e-12)
        np.testing.assert_allclose(snapped.distances, 1.0, atol=1e-12)


# ---------------------------------------------------------------------------
# Acceptance: STEP reproduces STL k=0 within the facet error
# ---------------------------------------------------------------------------


class TestSphereAcceptance:
    @pytest.mark.parametrize("topology", [CubedSphereTopology, OGridTopology])
    def test_k0_exact_on_step(self, cad, topology):
        mesh = Mesh.from_file(cad / "sphere.step", topology(), ni=17, nj=32, nk=6)
        np.testing.assert_allclose(np.linalg.norm(_k0(mesh), axis=1), 1.0, atol=1e-12)

    def test_k0_step_vs_stl_within_facet_sag(self, cad):
        step = Mesh.from_file(cad / "sphere.step", CubedSphereTopology(), ni=17, nj=32, nk=6)
        stl = Mesh.from_file(cad / "sphere.stl", CubedSphereTopology(), ni=17, nj=32, nk=6)
        diff = np.abs(_k0(step) - _k0(stl)).max()
        assert diff <= _stl_sag(cad / "sphere.stl")

    def test_march_with_exact_normals(self, cad):
        mesh = Mesh.from_file(cad / "sphere.step", CubedSphereTopology(), ni=17, nj=32, nk=6)
        mesh.march(ds=1e-2, growth=1.1, blend_normals_k=3)
        for block in mesh.blocks:
            r = np.linalg.norm(block.nodes, axis=-1)
            assert np.all(np.diff(r, axis=2) > 0.0)  # every layer moves outward


class TestCli:
    def test_mesh_step_sphere(self, cad, tmp_path):
        out = tmp_path / "s"
        args = ["--ni", "10", "--nj", "12", "--nk", "6", "--ds", "0.05", "--growth", "1.1"]
        r = runner.invoke(app, ["mesh", str(cad / "sphere.step"), *args, "-o", str(out)])
        assert r.exit_code == 0, r.output
        for label in ("sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"):
            nodes = np.asarray(read_geo(tmp_path / f"s_{label}.geo"))
            np.testing.assert_allclose(np.linalg.norm(nodes[:, :, 0], axis=-1), 1.0, atol=1e-9)

    def test_project_step_sphere(self, cad, tmp_path):
        out = tmp_path / "k0.npz"
        r = runner.invoke(
            app, ["project", str(cad / "sphere.step"), "--ni", "10", "--nj", "12", "-o", str(out)]
        )
        assert r.exit_code == 0, r.output
        k0 = np.load(out)["k0"]
        np.testing.assert_allclose(np.linalg.norm(k0, axis=-1), 1.0, atol=1e-12)

    def test_cgrid_step_wing_matches_stl(self, cad, tmp_path):
        args = [
            "--topology", "cgrid", "--ni-body", "60", "--ni-wake", "15",
            "--n-stations", "4", "--wake-length", "10.0", "--nk", "8",
            "--ds", "1e-3", "--smooth", "0.0", "--smooth-iters", "0",
        ]  # fmt: skip
        for ext in ("stl", "step"):
            r = runner.invoke(
                app, ["mesh", str(cad / f"wing.{ext}"), *args, "-o", str(tmp_path / ext)]
            )
            assert r.exit_code == 0, r.output
        for label in ("upper", "body", "lower"):
            a = np.asarray(read_geo(tmp_path / f"stl_{label}.geo"))
            b = np.asarray(read_geo(tmp_path / f"step_{label}.geo"))
            assert a.shape == b.shape
            # same polygonal wing in both formats: agreement to .geo print precision
            np.testing.assert_allclose(b, a, atol=1e-8)


class TestSurfaceGridSeam:
    def test_annulus_face_across_the_step_seam_is_rejected(self, tmp_path):
        from shore.exceptions import ParameterError
        from shore.geometry import BODY
        from shore.project import Project

        bd.export_step(bd.Cylinder(radius=0.5, height=4.0), str(tmp_path / "cyl.step"))
        p = Project()
        p.load_geometry(tmp_path / "cyl.step")
        p.build_primitive(
            "annulus", r_in=0.5, r_out=1.5, z0=-1.0, z1=1.0, ni=5, nj=13, nk=9,
            theta_start=-0.6, theta_end=0.6,
        )  # fmt: skip
        k = int(p.body.geometry.project(np.array([0.5, 0.0, 0.0]), BODY).entity_index)
        with pytest.raises(ParameterError, match="seam"):
            p.surface_grid(p.mesh.blocks[0].label, "i_lo", surface=k)
        assert len(p.history) == 2  # rejected before any snapshot or history entry
