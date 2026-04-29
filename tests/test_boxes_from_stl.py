"""Tests for shore.boxes.boxes_from_stl (layer 2 — STL → marker boxes).

The default method ``"voxel-fill"`` voxelises the body's AABB and
emits each inside-body voxel as a :class:`Box`, guaranteeing
body-volume coverage.  Single-box methods (``aabb``, ``obb``) are
best-effort and emit a UserWarning when their measured outward gap
exceeds the budget.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
import trimesh
from typer.testing import CliRunner

from shore.boxes import boxes_from_stl
from shore.cli import app
from shore.exceptions import MeshInputError, ParameterError

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def ellipsoid_stl(tmp_path_factory) -> Path:
    """4x-elongated ellipsoid (sphere stretched along x)."""
    tmp = tmp_path_factory.mktemp("ellipsoid")
    mesh = trimesh.creation.icosphere(subdivisions=3, radius=1.0)
    mesh.vertices[:, 0] *= 4.0
    path = tmp / "ellipsoid.stl"
    mesh.export(str(path))
    return path


@pytest.fixture(scope="module")
def cube_stl(tmp_path_factory) -> Path:
    """Watertight 1x1x1 cube — exact known volume for coverage tests."""
    tmp = tmp_path_factory.mktemp("cube")
    mesh = trimesh.creation.box(extents=(1.0, 1.0, 1.0))
    path = tmp / "cube.stl"
    mesh.export(str(path))
    return path


@pytest.fixture(scope="module")
def open_strip_stl(tmp_path_factory) -> Path:
    """Non-watertight STL: a cube with one face removed (5 quads,
    each split into 2 triangles → 10 triangles, 8 vertices).  Used
    to test the voxel-fill watertightness check."""
    tmp = tmp_path_factory.mktemp("open")
    cube = trimesh.creation.box(extents=(1.0, 1.0, 1.0))
    # Drop the +z face triangles so the cube is open on top.
    # +z face has both triangle normals pointing in +z.
    keep = cube.face_normals[:, 2] < 0.99
    cube = cube.submesh([np.where(keep)[0]], append=True)
    path = tmp / "open.stl"
    cube.export(str(path))
    return path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _all_inside_union(points: np.ndarray, boxes) -> bool:
    """True iff every point lies inside at least one box (parallelepiped coords)."""
    bbox_diag = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0)))
    tol = max(1e-9, bbox_diag * 1e-6)

    covered = np.zeros(len(points), dtype=bool)
    for box in boxes:
        v0 = box.vertices[0]
        edge_i = box.vertices[1] - v0
        edge_j = box.vertices[2] - v0
        edge_k = box.vertices[4] - v0
        basis = np.column_stack([edge_i, edge_j, edge_k])
        try:
            local = np.linalg.solve(basis, (points - v0).T).T
        except np.linalg.LinAlgError:
            continue
        inside = (local >= -tol).all(axis=1) & (local <= 1.0 + tol).all(axis=1)
        covered |= inside

    return bool(covered.all())


def _max_box_vertex_outward_gap(boxes, mesh: trimesh.Trimesh) -> float:
    """Max outward gap from any box vertex to the STL surface.

    Used by tests to verify voxel-fill respects the gap budget at
    every box's corners.
    """
    all_verts = np.vstack([b.vertices for b in boxes])
    prox = trimesh.proximity.ProximityQuery(mesh)
    signed = prox.signed_distance(all_verts)
    return float(np.maximum(0.0, -signed).max())


# ---------------------------------------------------------------------------
# Method = voxel-fill (default; guarantees body-volume coverage)
# ---------------------------------------------------------------------------


class TestVoxelFill:
    def test_default_method_is_voxel_fill(self, sphere_stl):
        """Implicit method= → voxel-fill."""
        boxes = boxes_from_stl(sphere_stl, "sub0", max_outward_gap=0.20)
        # Voxel-fill on a sphere produces many boxes, not 1.
        assert len(boxes) > 10

    def test_returns_multiple_boxes(self, sphere_stl):
        boxes = boxes_from_stl(sphere_stl, "sub0", max_outward_gap=0.20, method="voxel-fill")
        assert len(boxes) > 1

    def test_box_associated_block_propagated(self, sphere_stl):
        boxes = boxes_from_stl(
            sphere_stl,
            "wall.cap_north",
            max_outward_gap=0.20,
            method="voxel-fill",
        )
        assert all(b.associated_block == "wall.cap_north" for b in boxes)

    def test_covers_body_interior(self, cube_stl):
        """Random points sampled deep inside a unit cube must all be
        covered by the box union."""
        boxes = boxes_from_stl(cube_stl, "x", max_outward_gap=0.20, method="voxel-fill")
        # Sample 1000 points in the [-0.4, 0.4]^3 deep interior.
        rng = np.random.default_rng(0)
        pts = rng.uniform(-0.4, 0.4, size=(1000, 3))
        assert _all_inside_union(pts, boxes), "voxel-fill did not cover the body interior"

    def test_covers_sphere_interior(self, sphere_stl):
        """Random points sampled in the unit sphere of radius 0.4
        (well inside the unit-radius STL sphere) must all be covered."""
        boxes = boxes_from_stl(sphere_stl, "x", max_outward_gap=0.20, method="voxel-fill")
        rng = np.random.default_rng(0)
        # Uniform points inside r=0.4 sphere.
        dirs = rng.standard_normal((1000, 3))
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        radii = 0.4 * rng.random(1000) ** (1.0 / 3.0)
        pts = dirs * radii[:, None]
        assert _all_inside_union(pts, boxes)

    def test_box_vertices_within_gap_of_body(self, sphere_stl):
        """Every box vertex's outward distance to the body surface
        must not exceed max_outward_gap (with FP tolerance)."""
        gap = 0.20
        boxes = boxes_from_stl(sphere_stl, "x", max_outward_gap=gap, method="voxel-fill")
        mesh = trimesh.load(sphere_stl, force="mesh")
        max_gap = _max_box_vertex_outward_gap(boxes, mesh)
        assert max_gap <= gap + 1e-9, f"voxel-fill produced a corner gap {max_gap} > budget {gap}"

    def test_covers_body_surface(self, sphere_stl):
        """Every STL surface vertex must lie inside the union of boxes.

        This guards the bug where a body-shell of width ~side/2 was
        left uncovered because only voxels with centre strictly inside
        the body were emitted.
        """
        boxes = boxes_from_stl(sphere_stl, "x", max_outward_gap=0.20, method="voxel-fill")
        mesh = trimesh.load(sphere_stl, force="mesh")
        verts = np.asarray(mesh.vertices, dtype=np.float64)
        assert _all_inside_union(verts, boxes), (
            "voxel-fill left part of the body surface outside the box union"
        )

    def test_greedy_merge_reduces_box_count(self, sphere_stl):
        """The greedy AABB merge must collapse the cubic seed grid
        into far fewer boxes.  At gap=0.30 the seed grid has thousands
        of voxels (~3888 on a unit-radius icosphere); after merging
        the count must be at least 10x smaller.
        """
        boxes = boxes_from_stl(sphere_stl, "x", max_outward_gap=0.30, method="voxel-fill")
        # Seed grid: side = gap / sqrt(3) ≈ 0.173, AABB extent ≈ 2.0,
        # so ~12^3 ≈ 1700 candidate voxels and ~3800 inside-the-cover.
        # Merge target: ≤ 200 (>=10x reduction).  Empirically: ~50.
        assert len(boxes) <= 200, f"expected <=200 merged boxes, got {len(boxes)}"

    def test_merged_boxes_keep_corner_budget(self, sphere_stl):
        """Greedy merge must never emit a merged AABB whose corners
        violate the outward-gap budget."""
        gap = 0.30
        boxes = boxes_from_stl(sphere_stl, "x", max_outward_gap=gap, method="voxel-fill")
        mesh = trimesh.load(sphere_stl, force="mesh")
        max_gap = _max_box_vertex_outward_gap(boxes, mesh)
        assert max_gap <= gap + 1e-9, f"merged boxes exceeded corner gap: {max_gap} > {gap}"

    def test_merged_boxes_still_cover_body_surface(self, sphere_stl):
        """Merging must preserve full body surface coverage —
        a merged AABB is a strict superset of its parts, so this is
        a sanity check that the merge doesn't drop voxels."""
        boxes = boxes_from_stl(sphere_stl, "x", max_outward_gap=0.30, method="voxel-fill")
        mesh = trimesh.load(sphere_stl, force="mesh")
        verts = np.asarray(mesh.vertices, dtype=np.float64)
        assert _all_inside_union(verts, boxes)

    def test_tighter_gap_produces_more_boxes(self, sphere_stl):
        """Halving the gap budget should grow the box count.
        After merging, the count tracks body surface area rather
        than volume, so the scaling is ~4x (2D), not ~8x (3D)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            n_loose = len(
                boxes_from_stl(sphere_stl, "x", max_outward_gap=0.40, method="voxel-fill")
            )
            n_tight = len(
                boxes_from_stl(sphere_stl, "x", max_outward_gap=0.20, method="voxel-fill")
            )
        assert n_tight > 2 * n_loose, (
            f"halving gap should grow count: loose={n_loose}, tight={n_tight}"
        )

    def test_non_watertight_rejected(self, open_strip_stl):
        """voxel-fill needs a watertight STL for inside/outside
        classification.  An open mesh must be rejected."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(MeshInputError, match="watertight"):
                boxes_from_stl(
                    open_strip_stl,
                    "x",
                    max_outward_gap=0.10,
                    method="voxel-fill",
                )

    def test_axis_aligned(self, sphere_stl):
        """Voxel-fill emits axis-aligned boxes only."""
        boxes = boxes_from_stl(sphere_stl, "x", max_outward_gap=0.20, method="voxel-fill")
        for box in boxes:
            v = box.vertices
            edge_i = v[1] - v[0]
            edge_j = v[2] - v[0]
            edge_k = v[4] - v[0]
            # i-edge is along world x: y, z components zero.
            assert np.isclose(edge_i[1], 0.0) and np.isclose(edge_i[2], 0.0)
            assert np.isclose(edge_j[0], 0.0) and np.isclose(edge_j[2], 0.0)
            assert np.isclose(edge_k[0], 0.0) and np.isclose(edge_k[1], 0.0)


# ---------------------------------------------------------------------------
# Method = aabb (single-box, best-effort)
# ---------------------------------------------------------------------------


class TestAABB:
    def test_returns_one_box(self, sphere_stl):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            boxes = boxes_from_stl(sphere_stl, "sub0", max_outward_gap=1.0, method="aabb")
        assert len(boxes) == 1

    def test_world_axis_aligned(self, sphere_stl):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            boxes = boxes_from_stl(sphere_stl, "sub0", max_outward_gap=1.0, method="aabb")
        v = boxes[0].vertices
        edge_i = v[1] - v[0]
        assert np.isclose(edge_i[1], 0.0) and np.isclose(edge_i[2], 0.0)

    def test_unit_sphere_tight_extents(self, sphere_stl):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            boxes = boxes_from_stl(sphere_stl, "sub0", max_outward_gap=1.0, method="aabb")
        extents = boxes[0].vertices.max(axis=0) - boxes[0].vertices.min(axis=0)
        assert np.allclose(extents, 2.0, atol=1e-2)

    def test_warns_on_tight_budget(self, sphere_stl):
        """A unit sphere's AABB has corner gap r*(sqrt(3)-1) ≈ 0.73;
        a tight budget triggers the best-effort warning."""
        with pytest.warns(UserWarning, match="max outward gap.*max_outward_gap"):
            boxes_from_stl(sphere_stl, "sub0", max_outward_gap=0.05, method="aabb")

    def test_no_warning_on_loose_budget(self, sphere_stl):
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            boxes_from_stl(sphere_stl, "sub0", max_outward_gap=2.0, method="aabb")


# ---------------------------------------------------------------------------
# Method = obb (single-box, best-effort)
# ---------------------------------------------------------------------------


class TestOBB:
    def test_returns_one_box(self, sphere_stl):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            boxes = boxes_from_stl(sphere_stl, "sub0", max_outward_gap=1.0, method="obb")
        assert len(boxes) == 1

    def test_obb_aligns_with_principal_axis_on_ellipsoid(self, ellipsoid_stl):
        """OBB on a 4x-elongated ellipsoid must be longest along ~x."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            boxes = boxes_from_stl(ellipsoid_stl, "sub0", max_outward_gap=2.0, method="obb")
        extents = boxes[0].vertices.max(axis=0) - boxes[0].vertices.min(axis=0)
        sorted_ext = np.sort(extents)[::-1]
        assert sorted_ext[0] > 7.5
        assert sorted_ext[1] < 2.5
        assert sorted_ext[2] < 2.5

    def test_warns_on_tight_budget(self, sphere_stl):
        with pytest.warns(UserWarning, match="max outward gap.*max_outward_gap"):
            boxes_from_stl(sphere_stl, "sub0", max_outward_gap=0.05, method="obb")


# ---------------------------------------------------------------------------
# Parameter validation
# ---------------------------------------------------------------------------


class TestValidation:
    def test_zero_gap_rejected(self, sphere_stl):
        with pytest.raises(ParameterError, match="max_outward_gap must be > 0"):
            boxes_from_stl(sphere_stl, "x", max_outward_gap=0.0)

    def test_negative_gap_rejected(self, sphere_stl):
        with pytest.raises(ParameterError, match="max_outward_gap must be > 0"):
            boxes_from_stl(sphere_stl, "x", max_outward_gap=-0.05)

    def test_unknown_method_rejected(self, sphere_stl):
        with pytest.raises(ParameterError, match="method must be"):
            boxes_from_stl(
                sphere_stl,
                "x",
                max_outward_gap=0.05,
                method="cylinder",  # type: ignore[arg-type]
            )

    def test_bad_face_samples_rejected(self, sphere_stl):
        with pytest.raises(ParameterError, match="face_samples_per_axis"):
            boxes_from_stl(
                sphere_stl,
                "sub0",
                max_outward_gap=0.10,
                method="obb",
                face_samples_per_axis=0,
            )

    def test_missing_stl_rejected(self, tmp_path):
        with pytest.raises(MeshInputError, match="not found"):
            boxes_from_stl(
                tmp_path / "ghost.stl",
                "x",
                max_outward_gap=0.05,
            )


# ---------------------------------------------------------------------------
# CLI: shore boxes-from-stl
# ---------------------------------------------------------------------------


class TestCLI:
    def test_basic_invocation(self, sphere_stl, tmp_path):
        runner = CliRunner()
        out = tmp_path / "boxes.json"
        r = runner.invoke(
            app,
            [
                "boxes-from-stl",
                str(sphere_stl),
                "--block",
                "sub0",
                "--max-outward-gap",
                "0.20",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        import json

        payload = json.loads(out.read_text())
        assert isinstance(payload, list)
        assert len(payload) >= 1
        assert payload[0]["associated_block"] == "sub0"
        assert len(payload[0]["vertices"]) == 8

    def test_default_method_is_voxel_fill(self, sphere_stl, tmp_path):
        """Without --method, voxel-fill is the default → many boxes
        on a sphere with a moderate budget."""
        runner = CliRunner()
        out = tmp_path / "boxes.json"
        r = runner.invoke(
            app,
            [
                "boxes-from-stl",
                str(sphere_stl),
                "--block",
                "sub0",
                "--max-outward-gap",
                "0.20",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        assert "voxel-fill" in r.output
        import json

        payload = json.loads(out.read_text())
        assert len(payload) > 10

    def test_method_aabb(self, sphere_stl, tmp_path):
        runner = CliRunner()
        out = tmp_path / "boxes.json"
        r = runner.invoke(
            app,
            [
                "boxes-from-stl",
                str(sphere_stl),
                "--block",
                "sub0",
                "--max-outward-gap",
                "1.0",
                "--method",
                "aabb",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        import json

        assert len(json.loads(out.read_text())) == 1

    def test_method_obb(self, sphere_stl, tmp_path):
        runner = CliRunner()
        out = tmp_path / "boxes.json"
        r = runner.invoke(
            app,
            [
                "boxes-from-stl",
                str(sphere_stl),
                "--block",
                "sub0",
                "--max-outward-gap",
                "1.0",
                "--method",
                "obb",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output

    def test_unknown_method_exits_1(self, sphere_stl, tmp_path):
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "boxes-from-stl",
                str(sphere_stl),
                "--block",
                "sub0",
                "--max-outward-gap",
                "0.10",
                "--method",
                "kdtree-obb",
                "-o",
                str(tmp_path / "out.json"),
            ],
        )
        assert r.exit_code == 1
        assert "method" in r.output.lower()

    def test_missing_stl_exits_1(self, tmp_path):
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "boxes-from-stl",
                str(tmp_path / "ghost.stl"),
                "--block",
                "sub0",
                "--max-outward-gap",
                "0.10",
                "-o",
                str(tmp_path / "out.json"),
            ],
        )
        assert r.exit_code == 1
        assert "not found" in r.output.lower()

    def test_output_pipes_into_cc_par_boxes(self, sphere_stl, tmp_path):
        """End-to-end: shore boxes-from-stl → shore cc-par --boxes."""
        runner = CliRunner()
        import json

        adj = tmp_path / "adj.json"
        adj.write_text(
            json.dumps(
                {
                    "version": 2,
                    "blocks": [
                        {
                            "label": "sub0",
                            "shape": [4, 4, 4],
                            "faces": {
                                f: {
                                    "bc": "WALL" if f == "k_lo" else "FREE",
                                    "partner": None,
                                    "axis_map": None,
                                }
                                for f in (
                                    "i_lo",
                                    "i_hi",
                                    "j_lo",
                                    "j_hi",
                                    "k_lo",
                                    "k_hi",
                                )
                            },
                        }
                    ],
                }
            )
        )
        boxes_json = tmp_path / "boxes.json"
        r1 = runner.invoke(
            app,
            [
                "boxes-from-stl",
                str(sphere_stl),
                "--block",
                "sub0",
                "--max-outward-gap",
                "0.30",
                "-o",
                str(boxes_json),
            ],
        )
        assert r1.exit_code == 0, r1.output

        ccpar = tmp_path / "out.cc.par"
        r2 = runner.invoke(
            app,
            [
                "cc-par",
                str(adj),
                str(ccpar),
                "--grd",
                "x.grd",
                "--boxes",
                str(boxes_json),
            ],
        )
        assert r2.exit_code == 0, r2.output
        text = ccpar.read_text()
        n_boxes = json.loads(boxes_json.read_text())
        assert f"{len(n_boxes)} ! boxes number" in text
