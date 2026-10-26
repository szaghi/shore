"""FROZEN acceptance spec: hyperbolic volume marching (ADR 0005, roadmap phase 4).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py`` (the metric
definitions the thresholds below are written in). The ONE permitted edit is
deleting the ``pytestmark`` line below once the implementation passes; strict
xfail makes a passing-but-still-marked suite fail. SHA-256 of this file and
of ``quality.py`` are recorded in ADR 0005.

Verdict: ``pytest tests/spec`` must report 13 passed, 0 skipped, 0 xfailed.

Interface under test
--------------------
- ``Mesh.march(..., march_scheme="hyperbolic")``; the default stays
  ``"algebraic"`` (existing grids unchanged); an unknown value raises
  ``ParameterError``.
- ``Project.march(..., march_scheme=...)`` journals it; CLI
  ``shore mesh --march-scheme hyperbolic``.

Thresholds and their reference (algebraic, limited smoothing, measured
2026-10-26 with ``shore.quality`` on the same cases):

| Case | Metric | Algebraic | Required (hyperbolic) |
|---|---|---|---|
| H1 coarse sphere, ds 1e-3 | wall orthogonality mean / max (deg) | 12.15 / 18.4 | <= 3 / <= 10 |
| H1 | first cell / ds | 0.751 .. 1.000 | 0.95 .. 1.05 |
| H1 | outer-layer body distance / sum(steps) | 0.719 .. 1.000 | 0.90 .. 1.10 |
| H2 60x60x10 regression grid | wall orthogonality mean / max | 2.01 / 5.9 | <= 2.0 / <= 5.9 |
| H2 | min scaled Jacobian | 0.754 | >= 0.70 |
| H2 | outer-layer distance / sum(steps) | 0.932 .. 1.000 | 0.95 .. 1.05 |
| H3 plus prism (4 concave edges) | completes nk=15 | fails k=13 (legacy k=1) | yes, 0 folded |
| H3 | outer-layer distance / sum(steps) | n/a | 0.50 .. 1.50 |
| H4 3-block C-grid NACA 0012 | wall orthogonality mean / max | 0.48 / 18.2 | <= 1.0 / <= 20 |

Out of scope of this spec: blunt-TE C-grids, CH/OH channel blocks (they do
not march), stretched / non-star-shaped CAD bodies, grids finer than the
listed cases, and run time. Hyperbolic-specific parameters (dissipation,
volume smoothing) may be chosen freely but must not be passed per case: every
case runs with the scheme's defaults.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh
from typer.testing import CliRunner

from shore.cli import app
from shore.exceptions import ParameterError
from shore.mesh import Body, CGridTopology, CubedSphereTopology, Mesh, OGridTopology
from shore.mesh.spacing import Spacing1D
from shore.project import Project
from shore.quality import k_orthogonality_deg, k_spacing, summary

pytestmark = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

HYP = {"march_scheme": "hyperbolic"}


# ---------------------------------------------------------------------------
# Spec helpers (frozen with the spec)
# ---------------------------------------------------------------------------


def plus_prism(arm: float = 1.0, half_width: float = 0.3, half_height: float = 0.3):
    """Plus-shaped prism: two bars crossing in xy, extruded in z.

    Star-shaped w.r.t. the origin; four concave 90-degree edges along z.
    """
    a, w = arm, half_width
    poly = np.array(
        [[a, -w], [a, w], [w, w], [w, a], [-w, a], [-w, w],
         [-a, w], [-a, -w], [-w, -w], [-w, -a], [w, -a], [w, -w]],
        dtype=float,
    )  # fmt: skip
    n = len(poly)
    bot = np.c_[poly, np.full(n, -half_height)]
    top = np.c_[poly, np.full(n, half_height)]
    verts = np.vstack([bot, top, [[0.0, 0.0, -half_height], [0.0, 0.0, half_height]]])
    cb, ct = 2 * n, 2 * n + 1
    faces = []
    for i in range(n):
        j = (i + 1) % n
        faces += [[i, j, n + j], [i, n + j, n + i], [cb, j, i], [ct, n + i, n + j]]
    return trimesh.Trimesh(verts, np.array(faces), process=True)


def wall_orthogonality(mesh: Mesh) -> np.ndarray:
    return np.concatenate([k_orthogonality_deg(b.nodes, 0).ravel() for b in mesh.blocks])


def first_cell_ratio(mesh: Mesh, ds: float) -> np.ndarray:
    return np.concatenate([k_spacing(b.nodes)[:, :, 0].ravel() for b in mesh.blocks]) / ds


def outer_distance_ratio(mesh: Mesh, total: float) -> np.ndarray:
    outer = np.concatenate([b.nodes[:, :, -1].reshape(-1, 3) for b in mesh.blocks])
    return mesh.body.geometry.project(outer).distances / total


def assert_seams_coincide(mesh: Mesh) -> None:
    """Every SHARED / DIRICHLET face pair holds the same nodes, to round-off.

    Each node of a seam face must have a node of the partner face within
    ``1e-12 * diag`` (``diag`` = the mesh bounding-box diagonal). Not
    bitwise: the existing algebraic C-grid wake cut already differs at the
    1e-17 level, so bitwise would demand more than the reference meets.
    """
    from scipy.spatial import cKDTree

    from shore.mesh.face import FaceBC, _face_slice

    allnodes = np.concatenate([b.nodes.reshape(-1, 3) for b in mesh.blocks])
    tol = 1e-12 * float(np.linalg.norm(allnodes.max(axis=0) - allnodes.min(axis=0)))
    owner = {id(f): b for b in mesh.blocks for f in b.faces()}
    for b in mesh.blocks:
        for f in b.faces():
            if f.bc in (FaceBC.SHARED, FaceBC.DIRICHLET) and id(f.partner) in owner:
                mine = _face_slice(b, f.name).reshape(-1, 3)
                theirs = _face_slice(owner[id(f.partner)], f.partner.name).reshape(-1, 3)
                gap = cKDTree(theirs).query(mine)[0].max()
                assert gap <= tol, f"{b.label}.{f.name}: seam gap {gap:.3e} > {tol:.3e}"


@pytest.fixture(scope="module")
def plus_stl(tmp_path_factory):
    path = tmp_path_factory.mktemp("plus") / "plus.stl"
    plus_prism().export(str(path))
    return path


# ---------------------------------------------------------------------------
# H1: coarse sphere, small first step (the smoothing-defect regime)
# ---------------------------------------------------------------------------


class TestH1CoarseSphere:
    DS, GROWTH, NK = 1e-3, 1.15, 8

    @pytest.fixture(scope="class")
    def mesh(self, sphere_stl):
        m = Mesh.from_file(sphere_stl, CubedSphereTopology(), ni=10, nj=12, nk=self.NK)
        m.march(ds=self.DS, growth=self.GROWTH, **HYP)
        return m

    def test_valid(self, mesh):
        assert summary(mesh.blocks)["folded_cells"] == 0

    def test_wall_orthogonality(self, mesh):
        ortho = wall_orthogonality(mesh)
        assert ortho.mean() <= 3.0
        assert ortho.max() <= 10.0

    def test_first_cell_and_thickness(self, mesh):
        first = first_cell_ratio(mesh, self.DS)
        assert first.min() >= 0.95 and first.max() <= 1.05
        total = (self.DS * self.GROWTH ** np.arange(self.NK - 1)).sum()
        outer = outer_distance_ratio(mesh, total)
        assert outer.min() >= 0.90 and outer.max() <= 1.10

    def test_seams_and_determinism(self, mesh, sphere_stl):
        assert_seams_coincide(mesh)
        again = Mesh.from_file(sphere_stl, CubedSphereTopology(), ni=10, nj=12, nk=self.NK)
        again.march(ds=self.DS, growth=self.GROWTH, **HYP)
        for a, b in zip(mesh.blocks, again.blocks, strict=True):
            np.testing.assert_array_equal(a.nodes, b.nodes)


# ---------------------------------------------------------------------------
# H2: the 60x60x10 regression grid (must not be worse than algebraic)
# ---------------------------------------------------------------------------


class TestH2RegressionGrid:
    @pytest.fixture(scope="class")
    def mesh(self, sphere_stl):
        m = Mesh.from_file(
            sphere_stl,
            CubedSphereTopology(),
            ni=60,
            nj=60,
            nk=10,
            theta_cap=float(np.pi / 6),
            spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
        )
        m.march(ds=0.02, growth=1.1, **HYP)
        return m

    def test_orthogonality_and_validity(self, mesh):
        ortho = wall_orthogonality(mesh)
        assert ortho.mean() <= 2.0
        assert ortho.max() <= 5.9
        assert summary(mesh.blocks)["scaled_jacobian_min"] >= 0.70

    def test_thickness(self, mesh):
        outer = outer_distance_ratio(mesh, (0.02 * 1.1 ** np.arange(9)).sum())
        assert outer.min() >= 0.95 and outer.max() <= 1.05
        assert_seams_coincide(mesh)


# ---------------------------------------------------------------------------
# H3: concave robustness (plus prism; algebraic fails here)
# ---------------------------------------------------------------------------


class TestH3ConcavePlusPrism:
    DS, GROWTH, NK = 0.01, 1.2, 15

    @pytest.mark.parametrize("topology", [CubedSphereTopology, OGridTopology])
    def test_completes_without_folds(self, plus_stl, topology):
        mesh = Mesh.from_file(plus_stl, topology(), ni=24, nj=32, nk=self.NK)
        mesh.march(ds=self.DS, growth=self.GROWTH, **HYP)
        assert summary(mesh.blocks)["folded_cells"] == 0
        total = (self.DS * self.GROWTH ** np.arange(self.NK - 1)).sum()
        outer = outer_distance_ratio(mesh, total)
        assert outer.min() >= 0.50 and outer.max() <= 1.50

    def test_algebraic_still_fails_here(self, plus_stl):
        # keeps the case discriminating: if algebraic ever passes, revisit H3
        from shore.exceptions import MarchingError

        mesh = Mesh.from_file(plus_stl, CubedSphereTopology(), ni=24, nj=32, nk=self.NK)
        with pytest.raises(MarchingError):
            mesh.march(ds=self.DS, growth=self.GROWTH)
        hyp = Mesh.from_file(plus_stl, CubedSphereTopology(), ni=24, nj=32, nk=self.NK)
        hyp.march(ds=self.DS, growth=self.GROWTH, **HYP)


# ---------------------------------------------------------------------------
# H4: 3-block sharp-TE C-grid
# ---------------------------------------------------------------------------


class TestH4CGrid:
    def test_wing(self, naca0012_stl):
        body = Body.from_stl(naca0012_stl, skip_hygiene=True)
        topo = CGridTopology()
        blocks = topo.build(
            surface=None, nk=12, mesh=body.mesh, ni_body=80, ni_wake=20,
            n_stations=5, wake_length=10.0, blocks=3,
        )  # fmt: skip
        mesh = Mesh(body=body, blocks=blocks, topology=topo)
        mesh.march(ds=1e-3, growth=1.2, **HYP)
        assert summary(mesh.blocks)["folded_cells"] == 0
        wing = [b for b in mesh.blocks if b.label == "body"]
        ortho = np.concatenate([k_orthogonality_deg(b.nodes, 0).ravel() for b in wing])
        assert ortho.mean() <= 1.0
        assert ortho.max() <= 20.0
        assert_seams_coincide(mesh)


# ---------------------------------------------------------------------------
# H5: interface
# ---------------------------------------------------------------------------


class TestH5Interface:
    def test_default_and_validation(self, sphere_stl):
        import inspect

        assert inspect.signature(Mesh.march).parameters["march_scheme"].default == "algebraic"
        mesh = Mesh.from_file(sphere_stl, CubedSphereTopology(), ni=10, nj=12, nk=4)
        with pytest.raises(ParameterError):
            mesh.march(ds=0.05, march_scheme="elliptic")

    def test_project_journal(self, sphere_stl):
        p = Project()
        p.load_geometry(sphere_stl)
        p.build_topology("cubed_sphere", nk=6, ds=0.05, growth=1.1, ni=10, nj=12)
        p.march(march_scheme="hyperbolic")
        assert "march_scheme='hyperbolic'" in p.journal()

    def test_cli_flag(self, sphere_stl, tmp_path):
        r = CliRunner().invoke(
            app,
            ["mesh", str(sphere_stl), "--ni", "10", "--nj", "12", "--nk", "6", "--ds", "0.05",
             "--march-scheme", "hyperbolic", "-o", str(tmp_path / "s")],
        )  # fmt: skip
        assert r.exit_code == 0, r.output
