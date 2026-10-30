"""FROZEN acceptance spec: feature curves and points of STL geometry (ADR 0049).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_stl_curves.py``
reports 7 passed, 0 skipped, 0 xfailed. At freeze: 2 passed (C5: guards, the
refusals hold before the implementation), 5 xfailed.

Interface under test
--------------------
``StlGeometry`` (``crease_deg``, default 30):
- feature edges: the mesh edges whose two facets' normals differ by more
  than ``crease_deg``, and the open-boundary edges;
- ``entities(CURVE)``: the feature edges chained into polylines, split at
  every vertex where other than two feature edges meet; ``entities(POINT)``:
  those vertices; ``features()`` gives ``(curves, points)`` as float64
  arrays, curves in the order of their lowest vertex index;
- ``project(points, EntityRef(CURVE, k))``: the closest points of the
  polyline, NaN normals; ``EntityRef(POINT, k)``: the point;
- ``fingerprint``: a curve by its length-weighted centroid, length and box,
  a point by its position, so ``match_keys`` re-finds them in an edited STL;
- ``display()``: the curves and points drawn with the triangles;
- a CURVE or POINT index out of range: ``ParameterError``.

Cases
-----
- C1 a smooth sphere: no curves, no points (as ADR 0038 draws it).
- C2 a box: 12 two-point curves between its 8 corners; projection; display.
- C3 a blunt-TE wing (tip cap, no root cap): its two TE corner lines, the
  root and tip outlines, the TE bases; the 4 TE corners as points.
- C4 a mesh edge associated with the TE line (``associate``) and snapped:
  on the line; the journal replayed on a longer wing: on its TE line.
- C5 (x 2) out-of-range entities refused (``project``, ``associate``).
- C6 the box's points: a vertex associated with a corner and snapped.

Out of scope of this spec: plane sections and smooth-feature curves (the LE
of a smooth wing has none), curves' effect on premesh quality (none measured
on ADR 0046's wing), curve parameters (``evaluate`` stays refused).
"""

from __future__ import annotations

import numpy as np
import pytest
import trimesh

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")


def wing_mesh(span=2.0):
    from shore.profiles import naca
    from shore.shapes import extrude_section

    w = extrude_section(naca("0012", 60, closed_te=False), span=span, n_stations=21,
                        span_axis="z", end_caps=True)  # fmt: skip
    root = np.all(np.abs(w.vertices[w.faces][:, :, 2]) < 1e-12, axis=1)
    w.update_faces(~root)
    w.remove_unreferenced_vertices()
    return w


def te_line(curves, sign):
    """The TE corner line on the side *sign* (+1 upper, -1 lower)."""
    out = []
    for k, c in enumerate(curves):
        if np.allclose(c[:, 0], 1.0) and np.all(np.sign(c[:, 1]) == sign) and np.ptp(c[:, 2]) > 1:
            out.append(k)
    assert len(out) == 1
    return out[0]


def on_polyline(curve, points):
    a, b = curve[:-1], curve[1:]
    ab = b - a
    t = np.clip(np.einsum("nsj,sj->ns", points[:, None] - a, ab) / (ab * ab).sum(1), 0, 1)
    return np.linalg.norm(a + t[..., None] * ab - points[:, None], axis=-1).min(axis=1)


class TestFeatures:
    @PENDING
    def test_c1_sphere(self):
        from shore.geometry.base import EntityKind
        from shore.geometry.stl import StlGeometry

        g = StlGeometry(trimesh.creation.icosphere(subdivisions=3))
        curves, points = g.features()
        assert curves == [] and points.shape == (0, 3)
        assert g.entities(EntityKind.CURVE) == () and g.entities(EntityKind.POINT) == ()
        assert g.display().curves == () and len(g.display().vertices) == 0

    @PENDING
    def test_c2_box(self):
        from shore.geometry.base import EntityKind, EntityRef
        from shore.geometry.stl import StlGeometry

        g = StlGeometry(trimesh.creation.box(extents=(2.0, 2.0, 2.0)))
        curves, points = g.features()
        assert len(curves) == 12 and all(c.shape == (2, 3) for c in curves)
        assert points.shape == (8, 3) and np.allclose(np.abs(points), 1.0)
        corners = {tuple(p) for p in points}
        assert all(tuple(c[0]) in corners and tuple(c[-1]) in corners for c in curves)
        k = next(k for k, c in enumerate(curves)
                 if np.allclose(c[:, 0], 1.0) and np.allclose(c[:, 1], 1.0))  # fmt: skip
        pr = g.project(
            np.array([[3.0, 2.0, 0.25], [1.0, 1.0, 5.0]]), EntityRef(EntityKind.CURVE, k)
        )
        np.testing.assert_allclose(pr.points, [[1.0, 1.0, 0.25], [1.0, 1.0, 1.0]], atol=1e-15)
        np.testing.assert_allclose(pr.distances, [np.sqrt(5.0), 4.0])
        assert np.isnan(pr.normals).all()
        d = g.display()
        assert len(d.curves) == 12 and np.array_equal(d.vertices, points)
        keys = {g.fingerprint(r) for r in g.entities(EntityKind.CURVE)}
        assert len(keys) == 12  # distinct fingerprints

    @PENDING
    def test_c3_wing(self):
        from shore.geometry.stl import StlGeometry

        g = StlGeometry(wing_mesh())
        curves, points = g.features()
        assert len(curves) == 6 and points.shape == (4, 3)
        for sign in (1, -1):
            c = curves[te_line(curves, sign)]
            assert np.ptp(c[:, 2]) == pytest.approx(2.0)
        outlines = [c for c in curves if len(c) > 50]
        assert sorted(float(c[0, 2]) for c in outlines) == [0.0, 2.0]  # root and tip
        assert all(np.ptp(c[:, 2]) == 0.0 and np.ptp(c[:, 0]) == pytest.approx(1.0)
                   for c in outlines)  # fmt: skip
        np.testing.assert_allclose(np.sort(points[:, 2]), [0, 0, 2, 2])
        np.testing.assert_allclose(points[:, 0], 1.0)

    @PENDING
    def test_c4_associate_and_replay(self, tmp_path):
        from shore.project import Project
        from shore.project.project import run_journal

        wing_mesh().export(tmp_path / "wing.stl")
        p = Project()
        p.load_geometry(tmp_path / "wing.stl", name="wing", skip_hygiene=True)
        curves = p.geometries["wing"].geometry.features()[0]
        upper = te_line(curves, 1)
        p.build_primitive("box", vmin=(1.0, 0.0, 0.1), vmax=(1.5, 0.5, 1.9), ni=3, nj=3, nk=9)
        p.associate("box", "i_lo_j_lo", "wing", "curve", [upper])
        p.snap()
        edge = p.mesh.blocks[0].nodes[0, 0, :]
        assert on_polyline(curves[upper], edge).max() < 1e-12
        journal = p.journal().replace(repr(str(tmp_path / "wing.stl")), repr("wing.stl"))
        (tmp_path / "in").mkdir()
        wing_mesh(span=2.4).export(tmp_path / "in" / "wing.stl")
        (tmp_path / "journal.py").write_text(journal)
        q = run_journal(tmp_path / "journal.py", input_dir=tmp_path / "in")
        new = q.geometries["wing"].geometry.features()[0]
        te = new[te_line(new, 1)]
        assert np.ptp(te[:, 2]) == pytest.approx(2.4)
        assert on_polyline(te, q.mesh.blocks[0].nodes[0, 0, :]).max() < 1e-12

    @pytest.mark.parametrize("where", ["project", "associate"])
    def test_c5_out_of_range(self, tmp_path, where):
        from shore.geometry.base import EntityKind, EntityRef
        from shore.geometry.stl import StlGeometry
        from shore.project import Project

        if where == "project":
            g = StlGeometry(trimesh.creation.box())
            with pytest.raises(ParameterError):
                g.project(np.zeros((1, 3)), EntityRef(EntityKind.CURVE, 12))
            with pytest.raises(ParameterError):
                g.project(np.zeros((1, 3)), EntityRef(EntityKind.POINT, 8))
            return
        trimesh.creation.box().export(tmp_path / "box.stl")
        p = Project()
        p.load_geometry(tmp_path / "box.stl", name="cube")
        p.build_primitive("box", vmin=(0, 0, 0), vmax=(1, 1, 1), ni=3, nj=3, nk=3)
        journal = p.journal()
        with pytest.raises(ParameterError):
            p.associate("box", "i_lo_j_lo", "cube", "curve", [12])
        assert p.journal() == journal

    @PENDING
    def test_c6_point(self, tmp_path):
        from shore.project import Project

        trimesh.creation.box(extents=(2.0, 2.0, 2.0)).export(tmp_path / "box.stl")
        p = Project()
        p.load_geometry(tmp_path / "box.stl", name="cube")
        points = p.geometries["cube"].geometry.features()[1]
        k = int(np.argmin(np.linalg.norm(points - [1.0, 1.0, 1.0], axis=1)))
        p.build_primitive("box", vmin=(0.5, 0.5, 0.5), vmax=(1.5, 1.5, 1.5), ni=3, nj=3, nk=3)
        p.associate("box", "v0", "cube", "point", [k])
        p.snap()
        np.testing.assert_array_equal(p.mesh.blocks[0].nodes[0, 0, 0], [1.0, 1.0, 1.0])
