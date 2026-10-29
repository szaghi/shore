"""FROZEN acceptance spec: GUI slice 2, geometry display and slice planes (ADR 0038).

Frozen 2026-10-29, before any implementation. Implementing work may not edit
this file, ``tests_gui/conftest.py``, the pytest configuration in
``pyproject.toml`` or ``src/shore/quality.py``. The ONE permitted edit is
deleting an ``@PENDING`` line once that case passes (strict xfail turns a
passing pending case into a failure). ``tests/spec/SHA256SUMS`` holds the
hashes, taken with ``@PENDING`` lines removed; ``scripts/check_spec.py``
checks them.

Verdict for acceptance (GUI CI job, which installs ``.[dev,gui,cad]``):
``pytest tests_gui/spec/test_gui_inspect.py`` reports 12 passed, 0 skipped,
0 xfailed. At freeze: 0 passed, 12 xfailed.

Interface under test
--------------------
- ``Geometry.display() -> shore.geometry.base.Display``: triangles (float64
  ``points``, ``triangles``), each owned by entity ``EntityRef(surface_kind,
  triangle_entity[t])`` (the body for an STL, a face for CAD), ``curves``
  (one polyline per CURVE) and ``vertices`` (one per POINT). STL and OCCT
  backends; core, no GUI dependency.
- ``PygfxViewport``:
  - ``set_geometries([(name, Display)])``, ``geometry_names``;
  - ``set_slices([(SliceKey(block, axis, index), nodes)])``, ``slice_keys``;
  - ``show(kind, visible)`` for ``"block"``, ``"geometry"``, ``"slice"``;
  - ``uploaded_positions(key)``, ``camera_position()``, ``pixel_size(xyz)``.
  - ``pick`` returns ``GeometryPick(geometry, entity, xyz)`` on a geometry,
    and ``PickResult(SliceKey, "k=2", cell, xyz)`` on a slice. Curves and
    vertices are drawn over the faces they bound, so they stay pickable.
- ``MainWindow``:
  - ``geometry_tree`` (one row per geometry: name, entity summary);
  - ``set_visible(kind, visible)``;
  - ``add_slice(block_key, axis, index) -> SliceKey`` (``KeyError`` for an
    unknown block, ``IndexError`` out of range), ``remove_slice``,
    ``slices``.
  - Quality colours slices too.
  - View state (visibility, slices) is not a command: the journal and the
    project are unchanged. Slices of blocks that disappear are dropped on
    refresh.

Cases
-----
- V1 every project geometry is shown and listed; the uploaded positions
  are its display points.
- V2 STL picks (blocks hidden): the body, at the picked point.
- V3 CAD face picks: the face that owns the triangle, at the picked point.
  V2, V3 and V5 pick only points that no triangle hides from the camera (a
  CPU ray cast).
- V4 (x 2) CAD edge and vertex picks.
- V5 (x 2) slice picks: slice key, cell, exact centre; invalid slices
  refused.
- V6 quality colours on a slice.
- V7 (x 2) view state is not journaled; slices follow their blocks.
- V8 (x 2) ``display()`` of the STL and OCCT backends.

Out of scope of this spec: the look (colours, line widths), on-screen
mouse interaction, association editing, the performance of large CAD
models.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-29; implementation pending")

CUBED = {"topology": "cubed_sphere", "nk": "5", "ds": "0.05", "growth": "1.1",
         "params": "{'ni': 9, 'nj': 16}"}  # fmt: skip
FLAT_CAPS = {"kind": "flat_caps",
             "params": "{'r_out': 3.0, 'z0': -1.0, 'z1': 4.0, 'ni': 5, 'nc': 6, 'nk': 5, "
                       "'inner_frac': 0.4}"}  # fmt: skip


@pytest.fixture(scope="module")
def step_part(tmp_path_factory) -> Path:
    import build123d as bd

    path = tmp_path_factory.mktemp("cad") / "part.step"
    bd.export_step(bd.Box(2.0, 1.0, 1.0) - bd.Cylinder(0.3, 2.0), str(path))
    return path


def window(qapp):
    from shore.gui.window import MainWindow

    return MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)


def sphere_session(w, sphere_stl):
    w.run_form("load_geometry", {"path": str(sphere_stl)})
    w.run_form("build_topology", CUBED)
    w.wait_idle()
    w.run_form("march", {})
    w.wait_idle()
    assert w.last_error is None


def only(w, kind: str) -> None:
    for k in ("block", "geometry", "slice"):
        w.set_visible(k, k == kind)
    w.viewport.show_all()
    w.viewport.draw_now()


def unoccluded(eye, xyz, points, triangles) -> bool:
    """No triangle hides *xyz* from *eye* (a CPU ray cast, independent of the GPU)."""
    import trimesh

    mesh = trimesh.Trimesh(points, triangles, process=False)
    direction = xyz - eye
    dist = float(np.linalg.norm(direction))
    hits, _, _ = mesh.ray.intersects_location([eye], [direction / dist])
    return not len(hits) or float(np.linalg.norm(hits - eye, axis=1).min()) >= dist * (1 - 1e-6)


def facing_triangles(vp, display, step: int = 1):
    """(triangle, centroid) of triangles facing the camera, inside the image."""
    eye = vp.camera_position()
    width, height = (int(v) for v in vp.renderer.logical_size)
    for t in range(0, len(display.triangles), step):
        corners = display.points[display.triangles[t]]
        c = corners.mean(axis=0)
        n = np.cross(corners[1] - corners[0], corners[2] - corners[0])
        view = eye - c
        if n @ view <= 0.3 * np.linalg.norm(n) * np.linalg.norm(view):
            continue
        x, y = vp.project(c)
        inside = 10 <= x < width - 10 and 10 <= y < height - 10
        if inside and unoccluded(eye, c, display.points, display.triangles):
            yield t, c, x, y


# ---------------------------------------------------------------------------
# V1-V4: geometry
# ---------------------------------------------------------------------------


class TestGeometry:
    def test_v1_shown_and_listed(self, qapp, sphere_stl, step_part):
        w = window(qapp)
        w.run_form("load_geometry", {"path": str(sphere_stl)})
        w.run_form("load_geometry", {"path": str(step_part), "name": "part"})
        assert w.last_error is None
        names = list(w.session.project.geometries)
        assert w.viewport.geometry_names == names == ["body", "part"]
        rows = [w.geometry_tree.topLevelItem(i) for i in range(w.geometry_tree.topLevelItemCount())]
        assert [r.text(0) for r in rows] == names
        assert rows[0].text(1) == "1280 triangles"
        assert rows[1].text(1) == "7 faces, 15 edges, 10 vertices"
        for name in names:
            d = w.session.project.geometries[name].geometry.display()
            extent = float(np.ptp(d.points, axis=0).max())
            err = np.abs(w.viewport.uploaded_positions(name) - d.points).max()
            assert err <= 1e-6 * extent

    def test_v2_stl_picks(self, qapp, sphere_stl):
        from shore.geometry.base import BODY
        from shore.gui.viewport import GeometryPick

        w = window(qapp)
        sphere_session(w, sphere_stl)
        only(w, "geometry")
        vp = w.viewport
        d = w.session.project.geometries["body"].geometry.display()
        checked = 0
        for _, c, x, y in facing_triangles(vp, d, step=29):
            r = vp.pick(x, y)
            assert isinstance(r, GeometryPick) and (r.geometry, r.entity) == ("body", BODY)
            assert np.linalg.norm(r.xyz - c) <= 2.0 * vp.pixel_size(c)
            checked += 1
        assert checked >= 5

    def test_v3_cad_face_picks(self, qapp, step_part):
        from shore.geometry.base import EntityKind, EntityRef
        from shore.gui.viewport import GeometryPick

        w = window(qapp)
        w.run_form("load_geometry", {"path": str(step_part), "name": "part"})
        only(w, "geometry")
        vp = w.viewport
        d = w.session.project.geometries["part"].geometry.display()
        faces = set()
        for t, c, x, y in facing_triangles(vp, d):
            r = vp.pick(x, y)
            if isinstance(r, GeometryPick) and r.entity.kind is not EntityKind.SURFACE:
                continue  # an edge or a vertex drawn over this triangle
            expected = EntityRef(EntityKind.SURFACE, int(d.triangle_entity[t]))
            assert isinstance(r, GeometryPick) and r.entity == expected, (t, r)
            assert np.linalg.norm(r.xyz - c) <= 2.0 * vp.pixel_size(c)
            faces.add(expected.index)
        assert len(faces) >= 3

    @pytest.mark.parametrize("kind", ["CURVE", "POINT"])
    def test_v4_cad_edge_and_vertex_picks(self, qapp, step_part, kind):
        from shore.geometry.base import EntityKind
        from shore.gui.viewport import GeometryPick

        w = window(qapp)
        w.run_form("load_geometry", {"path": str(step_part), "name": "part"})
        only(w, "geometry")
        vp = w.viewport
        d = w.session.project.geometries["part"].geometry.display()
        width, height = (int(v) for v in vp.renderer.logical_size)
        if kind == "CURVE":
            # the polyline's middle (straight edges have two points: their ends are vertices)
            targets = [
                (c, 0.5 * (curve[(len(curve) - 1) // 2] + curve[len(curve) // 2]))
                for c, curve in enumerate(d.curves)
            ]
        else:
            targets = list(enumerate(d.vertices))
        hits = 0
        for index, xyz in targets:
            x, y = vp.project(xyz)
            if not (10 <= x < width - 10 and 10 <= y < height - 10):
                continue
            r = vp.pick(x, y)
            if isinstance(r, GeometryPick) and r.entity.kind is EntityKind[kind]:
                hits += r.entity.index == index
        assert hits >= 3


# ---------------------------------------------------------------------------
# V5-V7: slices and view state
# ---------------------------------------------------------------------------


class TestSlices:
    def test_v5_slice_picks(self, qapp, sphere_stl):
        from shore.gui.viewport import PickResult, SliceKey, slice_grid

        w = window(qapp)
        sphere_session(w, sphere_stl)
        key = w.session.blocks()[0][0]
        sk = w.add_slice(key, "k", 2)
        assert sk == SliceKey(key, "k", 2) and w.viewport.slice_keys == [sk]
        only(w, "slice")
        vp = w.viewport
        grid = slice_grid(w.session.block(key).nodes, "k", 2)
        from shore.gui.viewport import quad_indices

        na, nb, _ = grid.shape
        quads = quad_indices(na, nb)
        tris = np.concatenate([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
        flat = grid.reshape(-1, 3)
        eye = vp.camera_position()
        checked = 0
        for p in range(0, grid.shape[0] - 1, 2):
            for q in range(0, grid.shape[1] - 1):
                centre = grid[p : p + 2, q : q + 2].reshape(-1, 3).mean(axis=0)
                if not unoccluded(eye, centre, flat, tris):
                    continue  # the curved slice hides part of itself
                x, y = vp.project(centre)
                r = vp.pick(x, y)
                assert isinstance(r, PickResult) and r.key == sk and r.face == "k=2"
                assert r.cell == (p, q, 2) and np.array_equal(r.xyz, centre)
                checked += 1
        assert checked >= 6

    def test_v5_invalid_slices(self, qapp, sphere_stl):
        from shore.gui.session import BlockKey

        w = window(qapp)
        sphere_session(w, sphere_stl)
        key = w.session.blocks()[0][0]
        with pytest.raises(KeyError):
            w.add_slice(BlockKey("main", "nope"), "k", 1)
        with pytest.raises(IndexError):
            w.add_slice(key, "k", 99)
        assert w.slices == []

    def test_v6_quality_on_slices(self, qapp, sphere_stl):
        from shore.gui.viewport import colormap, slice_values
        from shore.quality import scaled_jacobian

        w = window(qapp)
        sphere_session(w, sphere_stl)
        key, block = w.session.blocks()[2]
        sk = w.add_slice(key, "i", 3)
        w.quality_action.setChecked(True)
        w.refresh_quality()
        expected = colormap(slice_values(scaled_jacobian(block.nodes), "i", 3), (-1.0, 1.0))
        assert np.array_equal(w.viewport.face_colours(sk), expected)

    def test_v7_view_state_not_journaled(self, qapp, sphere_stl):
        w = window(qapp)
        sphere_session(w, sphere_stl)
        lines = w.session.journal_lines()
        history = list(w.session.project.history)
        key = w.session.blocks()[0][0]
        sk = w.add_slice(key, "j", 1)
        w.set_visible("block", False)
        w.set_visible("geometry", False)
        w.remove_slice(sk)
        w.add_slice(key, "k", 3)
        assert w.session.journal_lines() == lines
        assert list(w.session.project.history) == history
        assert not w.viewport.visible("block") and not w.viewport.visible("geometry")

    def test_v7_slices_follow_blocks(self, qapp, sphere_stl):
        w = window(qapp)
        sphere_session(w, sphere_stl)
        w.run_form("use_component", {"name": "background"})
        w.run_form("build_primitive", FLAT_CAPS)
        w.wait_idle()
        main = w.session.blocks()[0][0]
        background = w.session.blocks()[-1][0]
        kept = w.add_slice(main, "k", 2)
        w.add_slice(background, "k", 1)
        w.undo()  # the background's primitive goes
        assert w.slices == [kept]
        assert w.viewport.slice_keys == [kept]


# ---------------------------------------------------------------------------
# V8: the display of the geometry backends
# ---------------------------------------------------------------------------


class TestDisplay:
    def test_v8_stl_display(self, sphere_stl):
        from shore.geometry import load_geometry
        from shore.geometry.base import EntityKind

        g = load_geometry(sphere_stl)
        d = g.display()
        assert d.surface_kind is EntityKind.BODY
        assert d.points.dtype == np.float64 and np.array_equal(d.points, g.mesh.vertices)
        assert np.array_equal(d.triangles, g.mesh.faces)
        assert np.array_equal(d.triangle_entity, np.zeros(len(g.mesh.faces)))
        assert d.curves == () and len(d.vertices) == 0

    def test_v8_cad_display(self, step_part):
        from shore.geometry import load_geometry
        from shore.geometry.base import EntityKind

        g = load_geometry(step_part)
        d = g.display()
        assert d.surface_kind is EntityKind.SURFACE
        assert d.triangles.max() < len(d.points)
        assert set(np.unique(d.triangle_entity)) == set(range(len(g.entities(EntityKind.SURFACE))))
        assert len(d.curves) == len(g.entities(EntityKind.CURVE))
        assert len(d.vertices) == len(g.entities(EntityKind.POINT))
        for curve in d.curves:  # every curve point lies on the part's tessellation
            dist = np.min(np.linalg.norm(d.points[None, :, :] - curve[:, None, :], axis=-1), axis=1)
            assert dist.max() <= 0.05
