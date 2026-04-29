"""Tests for VTK export (shore.io.vtk)."""

from __future__ import annotations

import numpy as np
import pytest

pv = pytest.importorskip("pyvista", reason="pyvista not installed")

from shore.io.vtk import (  # noqa: E402
    geo_to_pyvista,
    surface_to_pyvista,
    write_surface_vts,
    write_vts,
)
from tests.conftest import sphere_layer  # noqa: E402


def _make_volume(ni: int = 6, nj: int = 8, nk: int = 4) -> np.ndarray:
    """Synthetic O-grid volume: concentric sphere shells at geometric radii."""
    radii = 1.0 * 1.1 ** np.arange(nk)
    layers = [sphere_layer(ni, nj, radius=float(r)) for r in radii]
    return np.stack(layers, axis=2)


class TestGeotoPyvista:
    def test_dimensions(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid)
        ni, nj, nk, _ = grid.shape
        assert sg.dimensions == (ni, nj + 1, nk)  # j is closed

    def test_point_count(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid)
        ni, nj, nk, _ = grid.shape
        assert sg.n_points == ni * (nj + 1) * nk

    def test_jacobian_cell_data_present(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid, with_jacobian=True)
        assert "Jacobian" in sg.cell_data

    def test_jacobian_cell_count(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid, with_jacobian=True)
        ni, nj, nk, _ = grid.shape
        expected = (ni - 1) * nj * (nk - 1)
        assert sg.cell_data["Jacobian"].shape[0] == expected

    def test_jacobian_all_positive(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid, with_jacobian=True)
        assert sg.cell_data["Jacobian"].min() > 0.0

    def test_no_jacobian(self):
        grid = _make_volume(6, 8, 4)
        sg = geo_to_pyvista(grid, with_jacobian=False)
        assert "Jacobian" not in sg.cell_data


class TestSurfaceToPyvista:
    def test_dimensions(self):
        k0 = sphere_layer(8, 10, radius=1.0)
        sg = surface_to_pyvista(k0)
        ni, nj, _ = k0.shape
        assert sg.dimensions == (ni, nj + 1, 1)

    def test_point_count(self):
        k0 = sphere_layer(8, 10, radius=1.0)
        sg = surface_to_pyvista(k0)
        ni, nj, _ = k0.shape
        assert sg.n_points == ni * (nj + 1) * 1


class TestWriteVts:
    def test_creates_file(self, tmp_path):
        grid = _make_volume(6, 8, 4)
        out = tmp_path / "test.vts"
        write_vts(out, grid)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_appends_extension(self, tmp_path):
        grid = _make_volume(6, 8, 4)
        out = tmp_path / "test"
        write_vts(out, grid)
        assert (tmp_path / "test.vts").exists()

    def test_readable_by_pyvista(self, tmp_path):
        grid = _make_volume(6, 8, 4)
        out = tmp_path / "test.vts"
        write_vts(out, grid)
        loaded = pv.read(str(out))
        assert loaded.n_points > 0
        assert "Jacobian" in loaded.cell_data


class TestWriteSurfaceVts:
    def test_creates_file(self, tmp_path):
        k0 = sphere_layer(8, 10, radius=1.0)
        out = tmp_path / "surf.vts"
        write_surface_vts(out, k0)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_readable_by_pyvista(self, tmp_path):
        k0 = sphere_layer(8, 10, radius=1.0)
        out = tmp_path / "surf.vts"
        write_surface_vts(out, k0)
        loaded = pv.read(str(out))
        assert loaded.n_points > 0


# ---------------------------------------------------------------------------
# Marker boxes → VTK (write_boxes_vtm + write_caps_vtm boxes= layer)
# ---------------------------------------------------------------------------


from shore.io.boxes import Box  # noqa: E402
from shore.io.vtk import write_boxes_vtm, write_caps_vtm  # noqa: E402


def _aabb_box(label: str = "sub0") -> Box:
    return Box.from_aabb((-0.5, -0.5, -0.5), (0.5, 0.5, 0.5), label)


class TestWriteBoxesVtm:
    def test_empty_list_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="must not be empty"):
            write_boxes_vtm([], tmp_path / "boxes.vtm")

    def test_extension_added(self, tmp_path):
        out = write_boxes_vtm([_aabb_box()], tmp_path / "no_ext")
        assert out.suffix == ".vtm"

    def test_one_block_per_box(self, tmp_path):
        boxes = [_aabb_box("a"), _aabb_box("b")]
        out = write_boxes_vtm(boxes, tmp_path / "boxes.vtm")
        mb = pv.read(str(out))
        assert sorted(mb.keys()) == ["box_1", "box_2"]

    def test_each_block_is_one_hexahedron_cell(self, tmp_path):
        out = write_boxes_vtm([_aabb_box()], tmp_path / "boxes.vtm")
        mb = pv.read(str(out))
        grid = mb["box_1"]
        assert grid.n_cells == 1
        assert grid.n_points == 8
        # VTK_HEXAHEDRON has cell type id 12.
        assert int(grid.celltypes[0]) == int(pv.CellType.HEXAHEDRON)

    def test_cell_data_carries_associated_block_and_index(self, tmp_path):
        boxes = [_aabb_box("wall.sub0"), _aabb_box("wall.cap_north")]
        out = write_boxes_vtm(boxes, tmp_path / "boxes.vtm")
        mb = pv.read(str(out))
        b1 = mb["box_1"]
        b2 = mb["box_2"]
        assert b1.cell_data["associated_block"][0] == "wall.sub0"
        assert b1.cell_data["box_index"][0] == 1
        assert b2.cell_data["associated_block"][0] == "wall.cap_north"
        assert b2.cell_data["box_index"][0] == 2

    def test_vertex_remap_preserves_world_coordinates(self, tmp_path):
        """SHORE → VTK reordering must not move any vertex.  The set of
        8 written points must equal the input vertex set as an
        unordered point cloud."""
        box = _aabb_box()
        out = write_boxes_vtm([box], tmp_path / "boxes.vtm")
        mb = pv.read(str(out))
        grid = mb["box_1"]
        written = np.asarray(grid.points)
        # Sort both for set comparison (reorder is allowed; coordinate
        # values must match exactly).
        original_sorted = box.vertices[np.lexsort(box.vertices.T)]
        written_sorted = written[np.lexsort(written.T)]
        np.testing.assert_array_equal(written_sorted, original_sorted)

    def test_skewed_hexahedron_round_trips(self, tmp_path):
        """A non-axis-aligned (skewed) box must also round-trip
        without mangling — vertex reorder is independent of orientation."""
        # 30°-rotated cube around z.
        c, s = np.cos(np.pi / 6), np.sin(np.pi / 6)
        rot = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        verts = (Box.from_aabb((-1, -1, -1), (1, 1, 1), "x").vertices) @ rot.T
        box = Box(vertices=verts, associated_block="x")
        out = write_boxes_vtm([box], tmp_path / "boxes.vtm")
        mb = pv.read(str(out))
        grid = mb["box_1"]
        assert grid.n_cells == 1
        # Volume of a unit-2 cube = 8; rotated cube has the same volume.
        # PyVista's volume property gives us a sanity check that the
        # cell connectivity is right (a bow-tie connectivity would
        # produce a different volume).
        assert grid.volume == pytest.approx(8.0, rel=1e-9)


class TestWriteCapsVtmBoxesLayer:
    """write_caps_vtm with boxes=... layers boxes onto the structured
    blocks so ParaView shows everything in one MultiBlock."""

    def _make_caps_volume(self):
        return _make_volume(6, 8, 4)

    def test_no_boxes_keeps_block_count(self, tmp_path):
        block = self._make_caps_volume()
        out = tmp_path / "caps.vtm"
        write_caps_vtm(out, [block], labels=["only"])
        mb = pv.read(str(out))
        assert sorted(mb.keys()) == ["only"]

    def test_boxes_appended_with_box_keys(self, tmp_path):
        block = self._make_caps_volume()
        boxes = [_aabb_box("only"), _aabb_box("only")]
        out = tmp_path / "caps.vtm"
        write_caps_vtm(out, [block], labels=["only"], boxes=boxes)
        mb = pv.read(str(out))
        assert sorted(mb.keys()) == ["box_1", "box_2", "only"]

    def test_box_cells_inherit_metadata(self, tmp_path):
        block = self._make_caps_volume()
        out = tmp_path / "caps.vtm"
        write_caps_vtm(
            out,
            [block],
            labels=["only"],
            boxes=[_aabb_box("wall.cap_north")],
        )
        mb = pv.read(str(out))
        assert mb["box_1"].cell_data["associated_block"][0] == "wall.cap_north"


class TestCLIBoxesExport:
    def test_basic_invocation(self, tmp_path):
        # Build input JSON.
        import json

        from typer.testing import CliRunner

        from shore.cli import app
        from shore.io.boxes import box_to_json

        boxes_json = tmp_path / "boxes.json"
        boxes_json.write_text(json.dumps([box_to_json(_aabb_box("a"))]))

        runner = CliRunner()
        out = tmp_path / "boxes.vtm"
        r = runner.invoke(
            app,
            [
                "boxes-export",
                str(boxes_json),
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        mb = pv.read(str(out))
        assert "box_1" in mb.keys()

    def test_missing_file_exits_1(self, tmp_path):
        from typer.testing import CliRunner

        from shore.cli import app

        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "boxes-export",
                str(tmp_path / "ghost.json"),
                "-o",
                str(tmp_path / "out.vtm"),
            ],
        )
        assert r.exit_code == 1
        assert "not found" in r.output.lower()

    def test_empty_list_exits_1(self, tmp_path):
        from typer.testing import CliRunner

        from shore.cli import app

        boxes_json = tmp_path / "empty.json"
        boxes_json.write_text("[]")
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "boxes-export",
                str(boxes_json),
                "-o",
                str(tmp_path / "out.vtm"),
            ],
        )
        assert r.exit_code == 1
        assert "zero boxes" in r.output

    def test_non_list_json_exits_1(self, tmp_path):
        from typer.testing import CliRunner

        from shore.cli import app

        boxes_json = tmp_path / "bad.json"
        boxes_json.write_text('{"not": "a list"}')
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "boxes-export",
                str(boxes_json),
                "-o",
                str(tmp_path / "out.vtm"),
            ],
        )
        assert r.exit_code == 1
        assert "list" in r.output.lower()

    def test_invalid_json_exits_1(self, tmp_path):
        from typer.testing import CliRunner

        from shore.cli import app

        boxes_json = tmp_path / "broken.json"
        boxes_json.write_text("this is not json {{{")
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "boxes-export",
                str(boxes_json),
                "-o",
                str(tmp_path / "out.vtm"),
            ],
        )
        assert r.exit_code == 1
        assert "not valid JSON" in r.output
