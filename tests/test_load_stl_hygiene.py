"""Tests for load_stl mesh hygiene — Phase 1."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
import trimesh

from shore.exceptions import MeshInputError
from shore.surface.io import load_stl

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _export_stl(mesh: trimesh.Trimesh, path: Path) -> Path:
    mesh.export(str(path))
    return path


# ---------------------------------------------------------------------------
# Basic load
# ---------------------------------------------------------------------------


class TestLoadStlBasic:
    def test_loads_clean_stl(self, sphere_stl):
        mesh = load_stl(sphere_stl)
        assert isinstance(mesh, trimesh.Trimesh)
        assert len(mesh.faces) > 0

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(MeshInputError, match="not found"):
            load_stl(tmp_path / "nonexistent.stl")


# ---------------------------------------------------------------------------
# Hygiene: degenerate face removal
# ---------------------------------------------------------------------------


class TestDegenerateFaceRemoval:
    def test_degenerate_faces_removed(self, tmp_path):
        """A mesh with a duplicated zero-area face should have it removed."""
        sphere = trimesh.creation.icosphere(subdivisions=2, radius=1.0)
        # Add a degenerate triangle (three identical vertices → zero area)
        v0 = sphere.vertices[0]
        new_verts = np.vstack([sphere.vertices, [v0, v0, v0]])
        degenerate_face = np.array(
            [[len(sphere.vertices), len(sphere.vertices) + 1, len(sphere.vertices) + 2]]
        )
        new_faces = np.vstack([sphere.faces, degenerate_face])
        dirty = trimesh.Trimesh(vertices=new_verts, faces=new_faces, process=False)
        path = _export_stl(dirty, tmp_path / "dirty.stl")

        n_dirty = len(dirty.faces)
        mesh = load_stl(path)
        # Degenerate face should have been removed
        assert len(mesh.faces) < n_dirty

    def test_skip_hygiene_preserves_faces(self, tmp_path):
        """skip_hygiene=True leaves the mesh unmodified."""
        sphere = trimesh.creation.icosphere(subdivisions=2, radius=1.0)
        path = _export_stl(sphere, tmp_path / "sphere.stl")
        mesh = load_stl(path, skip_hygiene=True)
        # Should load without any modification
        assert isinstance(mesh, trimesh.Trimesh)
        assert len(mesh.faces) > 0


# ---------------------------------------------------------------------------
# Hygiene: duplicate vertex merging
# ---------------------------------------------------------------------------


class TestVertexMerge:
    def test_duplicate_vertices_merged(self, tmp_path):
        """A mesh exported with per-face vertices (duplicated) should have vertices merged."""
        sphere = trimesh.creation.icosphere(subdivisions=2, radius=1.0)
        # Explicitly split every face (creates 3x faces x 3 duplicated vertices)
        verts_split = sphere.vertices[sphere.faces.reshape(-1)]
        faces_split = np.arange(len(verts_split)).reshape(-1, 3)
        split_mesh = trimesh.Trimesh(vertices=verts_split, faces=faces_split, process=False)
        n_vert_before = len(split_mesh.vertices)
        path = _export_stl(split_mesh, tmp_path / "split.stl")

        mesh = load_stl(path)
        # Merged vertex count should be substantially less
        assert len(mesh.vertices) < n_vert_before


# ---------------------------------------------------------------------------
# Non-watertight warning
# ---------------------------------------------------------------------------


class TestNonWatertightWarning:
    def test_open_mesh_emits_warning(self, tmp_path):
        """A mesh with a hole should trigger a UserWarning."""
        sphere = trimesh.creation.icosphere(subdivisions=2, radius=1.0)
        # Remove one face to create a hole
        open_mesh = trimesh.Trimesh(
            vertices=sphere.vertices,
            faces=sphere.faces[1:],  # drop face 0
            process=False,
        )
        assert not open_mesh.is_watertight
        path = _export_stl(open_mesh, tmp_path / "open.stl")

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            load_stl(path)

        msgs = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
        assert any("watertight" in m.lower() for m in msgs), (
            f"Expected watertight warning; got: {msgs}"
        )

    def test_watertight_mesh_no_warning(self, sphere_stl):
        """A clean watertight mesh should produce no watertight warning."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            load_stl(sphere_stl)
        msgs = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
        assert not any("watertight" in m.lower() for m in msgs)
