"""FROZEN acceptance spec: per-component Xall options in ``write_xall`` (ADR 0017).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_xall_options.py``
reports 9 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 9 xfailed.

Interface under test
--------------------
``Project.write_xall(directory, *, grd="mesh.grd", np=1, components=None,
output_basename="cc.", multigrid_levels=(4, 1), boundary_layer_thickness=-0.1,
numerical_beach=-1.0)`` (ADR 0016 plus):

- the run-wide ``cc.par`` header options (``output_basename``,
  ``multigrid_levels``, ``boundary_layer_thickness``, ``numerical_beach``)
  are passed to every component's ``write_cc_par``: ``merge_cc_par`` takes
  the header from the first component and requires the others to agree;
- ``components`` maps a component name to its options:

  - ``free_face_bc``, ``meta`` (block label -> ``level`` / ``group`` /
    ``priority`` / ``comment`` / ``free_face_bc``) and ``boxes`` (marker-box
    JSON entries): passed to that component's ``write_cc_par``, as
    ``write_cc_par``'s own arguments of the same names;
  - ``body``: the Xall body number of every block of the component in
    ``proc.input``;
  - ``proc_meta`` (block label -> ``proc`` / ``group`` / ``body``):
    per-block ``proc.input`` overrides, ``body`` winning over the
    component's;

  ``proc.input`` is ``write_proc_input`` of the merged ``.grd`` with these
  assignments at the merged 1-based block indices (components in creation
  order, blocks in order);
- ``ParameterError`` for a component name that does not exist or has no
  mesh, an unknown option key, and a block label (in ``meta`` or
  ``proc_meta``) that the component does not have; checked before any file
  is written.

Cases (the two components of ADR 0016's spec: an O-grid ``main`` on an
icosphere STL, and a box x post ``collar``):

- X1 per-component cc.par options: merged and per-component ``cc.par``
  byte-identical to the accepted writers with the same options, merged.
- X2 run-wide header options: likewise, with every header option changed.
- X3 proc.input: component bodies and per-block overrides, byte-identical
  to ``write_proc_input`` with the assignments at the merged indices.
- X4 the journal line of a ``write_xall`` with options replays to the same
  files.
- X5 errors (x 5).

Out of scope of this spec: options for ``.geo`` / ``.grd`` writing,
re-balancing with per-component rank sets, CLI flags.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import trimesh

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

OGRID = {"nk": 5, "ds": 0.05, "growth": 1.1, "ni": 9, "nj": 16}
COLLAR = {
    "seed": (0.55, 0.0, 0.02), "n_loop": 16, "ds_surface": 0.02, "growth_surface": 1.15,
    "m_a": 5, "m_b": 5, "ds_wall": 1e-3, "growth_wall": 1.25, "nk": 8,
}  # fmt: skip
BOX = [[x, y, z] for z in (-2.0, 2.0) for y in (-2.0, 2.0) for x in (-2.0, 2.0)]


@pytest.fixture
def project(tmp_path):
    """'main': O-grid on an icosphere; 'collar': box x post collar."""
    from shore.project import Project

    d = tmp_path / "in"
    d.mkdir()
    trimesh.creation.icosphere(subdivisions=2, radius=1.0).export(str(d / "sphere.stl"))
    box = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
    box.apply_translation((0.0, 0.0, -0.5))
    box.export(str(d / "box.stl"))
    post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=48)
    post.apply_translation((0.0, 0.0, 0.2))
    post.export(str(d / "post.stl"))
    p = Project()
    p.load_geometry(d / "sphere.stl")
    p.build_topology("ogrid", **OGRID)
    p.march()
    p.load_geometry(d / "box.stl", name="box")
    p.load_geometry(d / "post.stl", name="post")
    p.use_component("collar")
    p.build_collar("post", geometry="box", **COLLAR)
    return p


def reference(p, out: Path, *, header, cc_options, assignments, np_ranks) -> None:
    """The accepted per-component writers with the options, then the merges."""
    from shore.io.cc_par import merge_cc_par
    from shore.io.grd import merge_grd
    from shore.io.proc_input import write_proc_input

    grds, ccs = [], []
    for name in ("main", "collar"):
        c = out / name
        p.use_component(name)
        p.write_geo(str(c / name), per_block=True)
        p.write_adjacency(str(c / f"{name}.adjacency.json"))
        geos = [str(c / f"{name}_{b.label}.geo") for b in p.mesh.blocks]
        p.write_grd(str(c / f"{name}.grd"), geos)
        options = {"grd_filename": f"{name}.grd", **header, **cc_options.get(name, {})}
        p.write_cc_par(str(c / f"{name}.adjacency.json"), str(c / "cc.par"), **options)
        grds.append(c / f"{name}.grd")
        ccs.append(c / "cc.par")
    merge_grd(grds, out / "mesh.grd")
    merge_cc_par(ccs, out / "cc.par", grd_filename="mesh.grd")
    write_proc_input(out / "mesh.grd", out / "proc.input", np=np_ranks, **assignments)


def assert_same(run: Path, ref: Path, names=("mesh.grd", "cc.par", "proc.input")) -> None:
    for name in names:
        assert (run / name).read_bytes() == (ref / name).read_bytes(), name
    for c in ("main", "collar"):
        assert (run / c / "cc.par").read_bytes() == (ref / c / "cc.par").read_bytes(), c


# ---------------------------------------------------------------------------
# X1-X3: options reach the right writer
# ---------------------------------------------------------------------------


class TestOptions:
    @PENDING
    def test_x1_per_component_cc_par_options(self, project, tmp_path):
        main_label = project.components["main"].blocks[0].label
        options = {
            "main": {
                "free_face_bc": 9,
                "meta": {main_label: {"level": 2, "group": 1, "priority": 3, "comment": "hull"}},
                "boxes": [{"associated_block": main_label, "vertices": BOX}],
            },
            "collar": {"meta": {"collar": {"free_face_bc": 3}}},
        }
        project.write_xall(tmp_path / "run", np=2, components=options)
        cc_options = {
            "main": {**options["main"], "boxes": tuple(options["main"]["boxes"])},
            "collar": options["collar"],
        }
        reference(project, tmp_path / "ref", header={}, cc_options=cc_options,
                  assignments={}, np_ranks=2)  # fmt: skip
        assert_same(tmp_path / "run", tmp_path / "ref")

    @PENDING
    def test_x2_run_wide_header_options(self, project, tmp_path):
        header = {
            "output_basename": "run.", "multigrid_levels": (3, 1),
            "boundary_layer_thickness": -0.2, "numerical_beach": 0.5,
        }  # fmt: skip
        project.write_xall(tmp_path / "run", **header)
        reference(project, tmp_path / "ref", header=header, cc_options={},
                  assignments={}, np_ranks=1)  # fmt: skip
        assert_same(tmp_path / "run", tmp_path / "ref")

    @PENDING
    def test_x3_proc_input_bodies_and_overrides(self, project, tmp_path):
        main_blocks = len(project.components["main"].blocks)
        main_label = project.components["main"].blocks[0].label
        options = {
            "main": {"body": 1, "proc_meta": {main_label: {"group": 5}}},
            "collar": {"body": 2, "proc_meta": {"collar": {"proc": 1, "body": 7}}},
        }
        project.write_xall(tmp_path / "run", np=2, components=options)
        collar_index = main_blocks + 1
        body = dict.fromkeys(range(1, main_blocks + 1), 1) | {collar_index: 7}
        assignments = {
            "proc_assignment": {collar_index: 1},
            "group_assignment": {1: 5},
            "body_assignment": body,
        }
        reference(project, tmp_path / "ref", header={}, cc_options={},
                  assignments=assignments, np_ranks=2)  # fmt: skip
        assert_same(tmp_path / "run", tmp_path / "ref")


# ---------------------------------------------------------------------------
# X4: journal
# ---------------------------------------------------------------------------


class TestJournal:
    @PENDING
    def test_x4_journal_line_replays(self, project, tmp_path):
        from shore.project import run_journal

        options = {"collar": {"body": 2, "meta": {"collar": {"level": 2}}}}
        project.write_xall(tmp_path / "run", np=2, components=options, multigrid_levels=(3, 1))
        line = project.journal().splitlines()[-1]
        assert line.startswith("p.write_xall(")
        replay = project.journal().replace(str(tmp_path / "run"), str(tmp_path / "again"))
        (tmp_path / "j.py").write_text(replay)
        run_journal(tmp_path / "j.py")
        for name in ("mesh.grd", "cc.par", "proc.input"):
            again, run = tmp_path / "again" / name, tmp_path / "run" / name
            assert again.read_bytes() == run.read_bytes(), name


# ---------------------------------------------------------------------------
# X5: errors
# ---------------------------------------------------------------------------


class TestErrors:
    @PENDING
    @pytest.mark.parametrize(
        "components",
        [
            {"nope": {"body": 1}},
            {"collar": {"colour": 3}},
            {"collar": {"meta": {"not_a_block": {"level": 2}}}},
            {"collar": {"proc_meta": {"not_a_block": {"proc": 0}}}},
            {"empty": {"body": 1}},
        ],
        ids=["unknown-component", "unknown-key", "meta-label", "proc-meta-label", "no-mesh"],
    )
    def test_x5_errors(self, project, tmp_path, components):
        project.use_component("empty")
        project.use_component("collar")
        with pytest.raises(ParameterError):
            project.write_xall(tmp_path / "run", components=components)
        assert not (tmp_path / "run").exists()  # checked before anything is written
