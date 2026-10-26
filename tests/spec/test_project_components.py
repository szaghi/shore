"""FROZEN acceptance spec: multi-component, portable projects (ADR 0016).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job):
``pytest tests/spec/test_project_components.py`` reports 14 passed, 0 skipped,
0 xfailed. At freeze: 1 passed (P3, a regression guard: format-1 projects load
today), 13 xfailed.

Interface under test (``shore.project.Project``)
------------------------------------------------
- ``load_geometry(path, anchor_strategy="auto", skip_hygiene=False, *,
  name="body")``: a named geometry; reloading a name replaces it and
  discards the components built on it. ``geometries``: mapping name -> Body.
- ``use_component(name)``: makes *name* the current component (created
  empty if new); journaled and undoable. Every mesh-building or
  mesh-changing command acts on the current component; ``mesh`` / ``build``
  are the current component's. ``components``: mapping name -> Mesh or
  None, in creation order; ``current_component``. The default component is
  ``"main"``.
- ``build_topology(..., geometry="body")``, ``surface_grid(...,
  geometry="body")``, ``build_collar(other, ..., geometry="body")``: the
  geometry by name. ``build_collar``'s *other* is a geometry name when one
  is loaded under it, otherwise a file path (format-1 journals).
- ``write_xall(directory, *, grd="mesh.grd", np=1)``: for every component
  with a mesh, in creation order, into ``directory/<component>/``:
  ``write_geo(<component>, per_block=True)``,
  ``write_adjacency(<component>.adjacency.json)``,
  ``write_grd(<component>.grd)`` from those ``.geo`` files in block order,
  ``write_cc_par(<component>.adjacency.json, cc.par,
  grd_filename=<component>.grd)``; then, in *directory*, *grd* =
  ``merge_grd`` of the component ``.grd`` files, ``cc.par`` =
  ``merge_cc_par`` of the component ``cc.par`` files with
  ``grd_filename=grd``, and ``proc.input`` = ``write_proc_input`` of the
  merged ``.grd`` with *np* ranks (one load balance over the whole run).
- ``save(directory)`` copies every input file the journal reads
  (geometries, ``build_collar`` files, ``.geo`` / sidecar inputs) into
  ``directory/inputs/`` and writes the journal against those copies;
  ``project.toml`` (format 2) records each copy's SHA-256. ``load`` replays
  with those inputs from wherever the directory now is, fails with
  ``MeshInputError`` when a copy no longer matches its hash, and still
  loads format-1 projects.
- ``ParameterError`` for an unknown geometry name, an empty component or
  geometry name, and ``write_xall`` with no mesh in any component.

Cases (small meshes: an icosphere STL O-grid, a box primitive, an STL box
and post collar):

- P1 a saved project moved elsewhere, with the original STL deleted, loads
  bit-identically.
- P2 a changed input copy is detected on load.
- P3 a format-1 project (journal with an absolute STL path) still loads.
- P4 two named geometries; ``build_collar`` by name equals ``build_collar``
  by file path, bit for bit.
- P5 a body-fitted O-grid and a collar coexist as two components; the first
  is unchanged by building the second.
- P6 ``write_xall``: the merged ``.grd``, ``cc.par`` and ``proc.input`` are
  byte-identical to the accepted tools applied as described above.
- P7 undo and redo across components.
- P8 save / load of a two-geometry, two-component project.
- P9 reloading a geometry name discards only the components built on it.
- P10 errors (x 4) and ``geometries`` / ``components`` reporting.

Out of scope of this spec: per-component ``cc.par`` options in
``write_xall`` (block metadata, marker boxes, natural BCs), relocating
*output* paths of journaled write commands (replay writes them where the
journal says), CLI flags, the GUI.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest
import trimesh

from shore.exceptions import MeshInputError, ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

OGRID = {"nk": 5, "ds": 0.05, "growth": 1.1, "ni": 9, "nj": 16}
COLLAR = {
    "seed": (0.55, 0.0, 0.02), "n_loop": 16, "ds_surface": 0.02, "growth_surface": 1.15,
    "m_a": 5, "m_b": 5, "ds_wall": 1e-3, "growth_wall": 1.25, "nk": 8,
}  # fmt: skip


@pytest.fixture
def inputs(tmp_path) -> dict[str, Path]:
    """A sphere (O-grid body), and a box with a post (collar), as STL files."""
    d = tmp_path / "in"
    d.mkdir()
    trimesh.creation.icosphere(subdivisions=2, radius=1.0).export(str(d / "sphere.stl"))
    box = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
    box.apply_translation((0.0, 0.0, -0.5))
    box.export(str(d / "box.stl"))
    post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=48)
    post.apply_translation((0.0, 0.0, 0.2))
    post.export(str(d / "post.stl"))
    return {"sphere": d / "sphere.stl", "box": d / "box.stl", "post": d / "post.stl"}


def project():
    from shore.project import Project

    return Project()


def nodes(mesh) -> list[np.ndarray]:
    return [b.nodes.copy() for b in mesh.blocks]


def same(a: list[np.ndarray], b: list[np.ndarray]) -> bool:
    return len(a) == len(b) and all(np.array_equal(x, y) for x, y in zip(a, b, strict=True))


def two_components(inputs):
    """'main': O-grid on the sphere; 'collar': box x post collar, by names."""
    p = project()
    p.load_geometry(inputs["sphere"])
    p.build_topology("ogrid", **OGRID)
    p.march()
    p.load_geometry(inputs["box"], name="box")
    p.load_geometry(inputs["post"], name="post")
    p.use_component("collar")
    p.build_collar("post", geometry="box", **COLLAR)
    return p


# ---------------------------------------------------------------------------
# P1-P3: portable save / load
# ---------------------------------------------------------------------------


class TestPortable:
    @PENDING
    def test_p1_moved_project_loads_without_the_original_inputs(self, inputs, tmp_path):
        from shore.project import Project

        p = project()
        p.load_geometry(inputs["sphere"])
        p.build_topology("ogrid", **OGRID)
        p.march()
        built = nodes(p.mesh)
        p.save(tmp_path / "saved")
        shutil.copytree(tmp_path / "saved", tmp_path / "elsewhere" / "moved")
        shutil.rmtree(tmp_path / "saved")
        inputs["sphere"].unlink()
        q = Project.load(tmp_path / "elsewhere" / "moved")
        assert same(nodes(q.mesh), built)

    @PENDING
    def test_p2_changed_input_is_detected(self, inputs, tmp_path):
        from shore.project import Project

        p = project()
        p.load_geometry(inputs["sphere"])
        p.build_topology("ogrid", **OGRID)
        p.save(tmp_path / "saved")
        copies = list((tmp_path / "saved" / "inputs").iterdir())
        assert len(copies) == 1
        trimesh.creation.icosphere(subdivisions=2, radius=1.1).export(str(copies[0]))
        with pytest.raises(MeshInputError):
            Project.load(tmp_path / "saved")

    def test_p3_format_1_project_still_loads(self, inputs, tmp_path):
        from shore.project import Project, run_journal

        d = tmp_path / "format1"
        (d / "grids").mkdir(parents=True)
        journal = (
            "# SHORE journal: replay with `shore run <this file>`.\n"
            "from shore.project import Project\n\n"
            "p = Project()\n"
            f"p.load_geometry({str(inputs['sphere'])!r})\n"
            "p.build_topology('ogrid', nk=5, ds=0.05, growth=1.1, ni=9, nj=16)\n"
            "p.march()\n"
        )
        (d / "journal.py").write_text(journal)
        reference = run_journal(d / "journal.py")
        labels = [b.label for b in reference.mesh.blocks]
        for i, block in enumerate(reference.mesh.blocks):
            np.savez(d / "grids" / f"{i:04d}.npz", nodes=block.nodes)
        (d / "project.toml").write_text(
            'format = 1\nshore_version = "0"\ncommands = 3\nlabels = ['
            + ", ".join(f'"{lbl}"' for lbl in labels)
            + "]\n"
        )
        q = Project.load(d)
        assert same(nodes(q.mesh), nodes(reference.mesh))


# ---------------------------------------------------------------------------
# P4-P5: named geometries, components
# ---------------------------------------------------------------------------


class TestComponents:
    @PENDING
    def test_p4_collar_by_name_equals_collar_by_path(self, inputs):
        by_path = project()
        by_path.load_geometry(inputs["box"])
        by_path.build_collar(str(inputs["post"]), **COLLAR)
        by_name = project()
        by_name.load_geometry(inputs["box"], name="box")
        by_name.load_geometry(inputs["post"], name="post")
        by_name.build_collar("post", geometry="box", **COLLAR)
        assert sorted(by_name.geometries) == ["box", "post"]
        assert same(nodes(by_name.mesh), nodes(by_path.mesh))

    @PENDING
    def test_p5_components_coexist(self, inputs):
        p = project()
        p.load_geometry(inputs["sphere"])
        p.build_topology("ogrid", **OGRID)
        p.march()
        main = nodes(p.mesh)
        p.load_geometry(inputs["box"], name="box")
        p.load_geometry(inputs["post"], name="post")
        p.use_component("collar")
        assert p.mesh is None
        p.build_collar("post", geometry="box", **COLLAR)
        assert list(p.components) == ["main", "collar"]
        assert p.current_component == "collar"
        assert [b.label for b in p.components["collar"].blocks] == ["collar"]
        assert same(nodes(p.components["main"]), main)


# ---------------------------------------------------------------------------
# P6: Xall output of the whole run
# ---------------------------------------------------------------------------


class TestXall:
    @PENDING
    def test_p6_write_xall_equals_the_accepted_tools(self, inputs, tmp_path):
        from shore.io.cc_par import merge_cc_par
        from shore.io.grd import merge_grd
        from shore.io.proc_input import write_proc_input

        p = two_components(inputs)
        p.write_xall(tmp_path / "run", grd="mesh.grd", np=3)

        ref = tmp_path / "ref"
        grds, ccs = [], []
        for name in ("main", "collar"):
            c = ref / name
            p.use_component(name)
            p.write_geo(str(c / name), per_block=True)
            p.write_adjacency(str(c / f"{name}.adjacency.json"))
            geos = [str(c / f"{name}_{b.label}.geo") for b in p.mesh.blocks]
            p.write_grd(str(c / f"{name}.grd"), geos)
            p.write_cc_par(str(c / f"{name}.adjacency.json"), str(c / "cc.par"),
                           grd_filename=f"{name}.grd")  # fmt: skip
            grds.append(c / f"{name}.grd")
            ccs.append(c / "cc.par")
        merge_grd(grds, ref / "mesh.grd")
        merge_cc_par(ccs, ref / "cc.par", grd_filename="mesh.grd")
        write_proc_input(ref / "mesh.grd", ref / "proc.input", np=3)
        for name in ("mesh.grd", "cc.par", "proc.input"):
            assert (tmp_path / "run" / name).read_bytes() == (ref / name).read_bytes(), name
        for name in ("main", "collar"):
            assert (tmp_path / "run" / name / "cc.par").read_bytes() == (
                ref / name / "cc.par"
            ).read_bytes()


# ---------------------------------------------------------------------------
# P7-P9: history
# ---------------------------------------------------------------------------


class TestHistory:
    @PENDING
    def test_p7_undo_redo_across_components(self, inputs):
        p = two_components(inputs)
        collar = nodes(p.mesh)
        p.undo()  # the collar
        assert p.current_component == "collar" and p.mesh is None
        p.undo()  # use_component
        assert p.current_component == "main"
        assert list(p.components) == ["main"]
        p.redo()
        p.redo()
        assert p.current_component == "collar"
        assert same(nodes(p.mesh), collar)

    @PENDING
    def test_p8_save_load_two_components(self, inputs, tmp_path):
        from shore.project import Project

        p = two_components(inputs)
        p.save(tmp_path / "saved")
        shutil.move(str(tmp_path / "saved"), str(tmp_path / "moved"))
        for path in inputs.values():
            path.unlink()
        q = Project.load(tmp_path / "moved")
        assert list(q.components) == ["main", "collar"]
        assert q.current_component == "collar"
        for name in ("main", "collar"):
            assert same(nodes(q.components[name]), nodes(p.components[name]))

    @PENDING
    def test_p9_reloading_a_geometry_discards_only_its_components(self, inputs):
        p = two_components(inputs)
        main = nodes(p.components["main"])
        p.load_geometry(inputs["post"], name="post")
        assert p.components["collar"] is None
        assert same(nodes(p.components["main"]), main)


# ---------------------------------------------------------------------------
# P10: errors and reporting
# ---------------------------------------------------------------------------


class TestErrors:
    @PENDING
    def test_p10_unknown_geometry(self, inputs):
        p = project()
        p.load_geometry(inputs["box"], name="box")
        with pytest.raises(ParameterError):
            p.build_topology("ogrid", geometry="nope", **OGRID)

    @PENDING
    @pytest.mark.parametrize("call", ["component", "geometry"])
    def test_p10_empty_names(self, inputs, call):
        p = project()
        with pytest.raises(ParameterError):
            if call == "component":
                p.use_component("")
            else:
                p.load_geometry(inputs["box"], name="")

    @PENDING
    def test_p10_write_xall_without_a_mesh(self, tmp_path):
        p = project()
        p.use_component("empty")
        with pytest.raises(ParameterError):
            p.write_xall(tmp_path / "run")

    @PENDING
    def test_p10_reporting(self, inputs):
        p = project()
        assert dict(p.geometries) == {} and list(p.components) == ["main"]
        assert p.current_component == "main" and p.components["main"] is None
        p.load_geometry(inputs["sphere"])
        assert list(p.geometries) == ["body"]
