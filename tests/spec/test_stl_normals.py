"""FROZEN acceptance spec: crease-aware STL surface normals (ADR 0014).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_stl_normals.py``
reports 17 passed, 0 skipped, 0 xfailed. At freeze: 2 passed, 15 xfailed (N2
icosphere and N6 are regression guards that hold before and after).

Interface under test
--------------------
``shore.geometry.StlGeometry(mesh, crease_deg=30.0)`` (and
``StlGeometry.from_file(path, skip_hygiene=False, crease_deg=30.0)``):
``project`` returns, at the closest point, the barycentric blend of the
closest triangle's three **corner normals**:

- a corner normal is the corner-angle-weighted mean of the normals of the
  faces around its vertex whose normal is within *crease_deg* of the
  triangle's own (trimesh's vertex-normal weighting, Thuerrner & Wuethrich
  1998, restricted to one side of every sharper edge);
- on an open boundary whose vertices lie in a plane (the cut of a half
  model), the corner normals of its vertices are projected into that plane;
  other open boundaries are left alone;
- ``crease_deg = 180`` on a closed mesh is the previous normal (the blend of
  trimesh's vertex normals);
- ``ParameterError`` for *crease_deg* outside (0, 180].

Points and distances are unchanged; ``volume.hyperbolic._stl_normals`` (the
legacy ``shore mesh`` marching) is not touched.

Cases (normals compared as vectors, to 1e-12; ``old`` is the previous
formula, defined here):

- N1 box: every sampled surface point gets its face normal (old: up to
  54.5 deg off).
- N2 curved: a 96-section cylinder's side normals are the exact radial
  normals of the points and its caps' are +-z (old: up to 45.6 deg off); an
  icosphere (no edge sharper than 30 deg) keeps the old normals.
- N3 regression: ``crease_deg=180`` on closed meshes (icosphere, box,
  cylinder) is the old normal.
- N4 half models: the half post and the half box (cut at y = 0) give the
  whole meshes' normals at points of the half off its crease edges (old:
  up to 1.9 deg off), and every cut vertex's normals lie in the plane.
- N5 downstream: ``march_surface(plane=)`` on the STL half model equals the
  y >= 0 half of the periodic march of the mirrored curve on the whole
  meshes within 1e-6 x the first step (ADR 0012's E3 could not require
  this; prototype 2e-13 box, 7e-8 post, where trimesh splits every side
  quad along the same diagonal so the whole cylinder is not exactly
  mirror-symmetric; old 0.47, 3.9e-3).
- N6 an open boundary that is not planar keeps the crease normals: an
  icosphere opened by a curved cut gives the old normals.
- N7 interface: ``crease_deg`` checked; ``from_file`` passes it on.

Out of scope of this spec: points on a crease edge (the closest triangle
there is either side's, and so is the normal), feature detection beyond the
dihedral angle, the legacy ``_stl_normals`` path, CAD tessellation proxies
(they inherit the change; checked by the CAD job), run time.
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh
import trimesh.proximity
import trimesh.triangles

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

TOL = 1e-12


def stl(mesh, **kwargs):
    from shore.geometry import StlGeometry

    return StlGeometry(mesh, **kwargs)


def old_normals(mesh: trimesh.Trimesh, points: np.ndarray) -> np.ndarray:
    """The previous normal: blend of trimesh's vertex normals over the closest triangle."""
    closest, _, tri = trimesh.proximity.closest_point(mesh, points)
    bary = trimesh.triangles.points_to_barycentric(mesh.triangles[tri], closest)
    n = (bary[..., None] * mesh.vertex_normals[mesh.faces[tri]]).sum(axis=1)
    return n / np.linalg.norm(n, axis=-1, keepdims=True)


def sample(mesh: trimesh.Trimesh, count: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    points, faces = trimesh.sample.sample_surface(mesh, count, seed=seed)
    return np.asarray(points), np.asarray(faces)


def box() -> trimesh.Trimesh:
    return trimesh.creation.box(extents=(4.0, 4.0, 1.0))


def cylinder() -> trimesh.Trimesh:
    return trimesh.creation.cylinder(radius=0.5, height=1.6, sections=96)


def half(mesh: trimesh.Trimesh, keep) -> trimesh.Trimesh:
    out = mesh.copy()
    out.update_faces(keep(out))
    out.remove_unreferenced_vertices()
    return out


def half_post() -> trimesh.Trimesh:
    return half(cylinder(), lambda m: m.triangles_center[:, 1] > 0.0)


def half_box() -> trimesh.Trimesh:
    m = trimesh.creation.box(extents=(4.0, 2.0, 1.0))
    m.apply_translation((0.0, 1.0, 0.0))
    return half(m, lambda m: ~np.all(np.isclose(m.vertices[m.faces][..., 1], 0.0), axis=1))


# ---------------------------------------------------------------------------
# N1-N3: closed meshes
# ---------------------------------------------------------------------------


class TestClosedMeshes:
    def test_n1_box_faces_are_flat(self):
        mesh = box()
        points, faces = sample(mesh, 2000, seed=1)
        normals = stl(mesh).project(points).normals
        assert np.abs(normals - mesh.face_normals[faces]).max() <= TOL

    def test_n2_cylinder_side_radial_caps_axial(self):
        mesh = cylinder()
        points, _ = sample(mesh, 4000, seed=2)
        normals = stl(mesh).project(points).normals
        side = np.abs(np.abs(points[:, 2]) - 0.8) > 1e-9
        radial = points[side] * [1.0, 1.0, 0.0]
        radial /= np.linalg.norm(radial, axis=-1, keepdims=True)
        assert np.abs(normals[side] - radial).max() <= TOL
        axial = np.c_[0 * points[~side, :2], np.sign(points[~side, 2])]
        assert np.abs(normals[~side] - axial).max() <= TOL

    def test_n2_icosphere_unchanged(self):
        mesh = trimesh.creation.icosphere(subdivisions=3)
        points, _ = sample(mesh, 2000, seed=3)
        assert np.abs(stl(mesh).project(points).normals - old_normals(mesh, points)).max() <= TOL

    @pytest.mark.parametrize("make", [lambda: trimesh.creation.icosphere(3), box, cylinder])
    def test_n3_crease_180_is_the_old_normal(self, make):
        mesh = make()
        points, _ = sample(mesh, 2000, seed=4)
        normals = stl(mesh, crease_deg=180.0).project(points).normals
        assert np.abs(normals - old_normals(mesh, points)).max() <= TOL


# ---------------------------------------------------------------------------
# N4-N5: half models
# ---------------------------------------------------------------------------


class TestHalfModels:
    @pytest.mark.parametrize(
        ("make_half", "make_whole"),
        [(half_post, cylinder), (half_box, box)],
    )
    def test_n4_half_gives_the_whole_meshes_normals(self, make_half, make_whole):
        part, whole = make_half(), make_whole()
        points, _ = sample(part, 4000, seed=5)
        # stay off the crease edges, where the closest triangle is a tie
        edges = part.face_adjacency_edges[part.face_adjacency_angles > np.radians(30.0)]
        if len(edges):
            a, b = part.vertices[edges[:, 0]], part.vertices[edges[:, 1]]
            ab = b - a
            t = np.clip(np.einsum("pej,ej->pe", points[:, None] - a, ab) / (ab * ab).sum(-1), 0, 1)
            d = np.linalg.norm(points[:, None] - (a + t[..., None] * ab), axis=-1).min(axis=1)
            points = points[d > 1e-6]
        got, want = stl(part).project(points).normals, stl(whole).project(points).normals
        assert np.abs(got - want).max() <= TOL

    def test_n4_cut_vertex_normals_in_the_plane(self):
        part = half_post()
        cut = part.vertices[np.abs(part.vertices[:, 1]) < 1e-12]
        assert len(cut) > 0
        assert np.abs(stl(part).project(cut).normals[:, 1]).max() <= TOL

    @pytest.mark.parametrize("side", ["box", "post"])
    def test_n5_half_model_march_is_mirror_symmetric(self, side):
        from shore.surface.intersect import intersection_curve
        from shore.surface.march import march_surface

        hbox, hpost = half_box(), half_post()
        hbox.apply_translation((0.0, 0.0, -0.5))
        hpost.apply_translation((0.0, 0.0, 0.2))
        wbox = box().apply_translation((0.0, 0.0, -0.5))
        wpost = cylinder().apply_translation((0.0, 0.0, 0.2))
        a, b = stl(hbox), stl(hpost)
        plane = ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))
        curve = intersection_curve(a, b, np.array([0.45, 0.2, 0.02]), 33, plane=plane)
        steps = 0.01 * 1.12 ** np.arange(15 if side == "box" else 12)
        if side == "box":
            body, whole, start = a, stl(wbox), curve
        else:
            body, whole, start = b, stl(wpost), curve[::-1]
        grid = march_surface(body, start, steps, plane=plane)
        closed = np.vstack([start, start[-2:0:-1] * [1.0, -1.0, 1.0]])
        mirror = march_surface(whole, closed, steps, periodic=True)[: len(start)]
        assert np.abs(grid - mirror).max() <= 1e-6 * steps[0]


# ---------------------------------------------------------------------------
# N6: non-planar open boundaries
# ---------------------------------------------------------------------------


class TestOpenBoundary:
    def test_n6_curved_cut_keeps_crease_normals(self):
        sphere = trimesh.creation.icosphere(subdivisions=3)
        centres = sphere.triangles_center
        part = half(sphere, lambda m: centres[:, 2] > 0.4 * centres[:, 0] ** 2 - 0.2)
        points, _ = sample(part, 2000, seed=6)
        assert np.abs(stl(part).project(points).normals - old_normals(part, points)).max() <= TOL


# ---------------------------------------------------------------------------
# N7: interface
# ---------------------------------------------------------------------------


class TestInterface:
    @pytest.mark.parametrize("crease", [0.0, -1.0, 181.0, float("nan")])
    def test_n7_bad_crease_angle(self, crease):
        with pytest.raises(ParameterError):
            stl(box(), crease_deg=crease)

    def test_n7_from_file_passes_crease_angle(self, tmp_path):
        from shore.geometry import StlGeometry

        mesh = box()
        mesh.export(str(tmp_path / "box.stl"))
        geometry = StlGeometry.from_file(tmp_path / "box.stl", crease_deg=180.0)
        points, _ = sample(geometry.mesh, 500, seed=7)
        normals = geometry.project(points).normals
        assert np.abs(normals - old_normals(geometry.mesh, points)).max() <= TOL
