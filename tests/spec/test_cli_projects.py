"""FROZEN acceptance spec: CLI for collars and multi-component projects (ADR 0018).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_cli_projects.py``
reports 13 passed, 0 skipped, 0 xfailed. At freeze: 0 passed, 13 xfailed.

Interface under test (``shore`` CLI)
------------------------------------
- ``shore collar BODY_A BODY_B --seed X Y Z --n-loop N --ds-surface D
  --growth-surface G --m-a M --m-b M --ds-wall D --growth-wall G --nk N
  [--label L] [--plane PX PY PZ NX NY NZ] [--name-a A] [--name-b B]
  [-o STEM] [--adjacency-json PATH] [--project DIR [--component C]]``: the
  collar (ADR 0010; with ``--plane`` the half-collar, ADR 0013) of the two
  geometry files, as ``Project.build_collar(name_b, geometry=name_a, ...)``
  after ``load_geometry`` of each file under its name (default: the file
  stem) and ``use_component(C)`` (default ``"collar"``). ``-o`` writes the
  ``.geo`` file(s) as ``write_geo``.
- ``--project DIR [--component C]`` on ``shore mesh``, ``shore collar`` and
  ``shore primitive *``: the project in *DIR* is loaded (or created when
  *DIR* does not exist), the command builds component *C* (``shore mesh``
  and ``shore primitive``: default ``"main"``), and the project is saved
  back to *DIR*. ``shore mesh`` loads its input under ``--geometry-name``
  (default: the file stem). A geometry name already loaded from a file with
  the same SHA-256 is reused, not reloaded, so the components built on it
  survive.
- ``shore xall PROJECT_DIR --out RUN [--options xall.toml]``: loads the
  project and runs ``write_xall(RUN, ...)``. The options file has a
  ``[run]`` table (``grd``, ``np``, ``output_basename``,
  ``multigrid_levels``, ``boundary_layer_thickness``, ``numerical_beach``)
  and ``[components.<name>]`` tables (``free_face_bc``, ``body``,
  ``[components.<name>.meta.<label>]``,
  ``[components.<name>.proc_meta.<label>]``, and ``boxes``, a path to a
  boxes JSON file as ``shore boxes-from-stl`` writes, relative to the
  options file).
- Exit code 1 with a message for: an unknown key in the options file, a
  component in it that does not exist, a directory that is not a project,
  ``--component`` without ``--project``, and ``shore collar`` with neither
  ``-o`` nor ``--project``.

Cases (small: an icosphere STL O-grid, an STL box with a post, its half
model):

- L1 ``shore collar --project``: the saved collar equals the Python API's,
  bit for bit.
- L2 half collar with ``--plane``: equals the Python API's; SYMMETRY j-faces.
- L3 ``shore collar -o``: the ``.geo`` file equals ``write_geo`` of the
  Python API's collar, byte for byte.
- L4 ``shore mesh --project``: equals the replay of the same run's
  ``--journal``.
- L5 mesh, collar and primitive into one project: three components, the
  earlier ones unchanged.
- L6 two collars on the same files into one project: both survive
  (geometries reused).
- L7 ``shore xall --options``: byte-identical to ``Project.load(...)
  .write_xall`` with the options the file describes.
- L8 errors (x 6).

Out of scope of this spec: a CLI for every Project method, progress output,
CAD inputs (the same code path; CAD job), run time.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import trimesh
from typer.testing import CliRunner

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

runner = CliRunner()
COLLAR = {
    "seed": (0.55, 0.0, 0.02), "n_loop": 16, "ds_surface": 0.02, "growth_surface": 1.15,
    "m_a": 5, "m_b": 5, "ds_wall": 1e-3, "growth_wall": 1.25, "nk": 8,
}  # fmt: skip
HALF = {**COLLAR, "seed": (0.45, 0.2, 0.02), "plane": ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))}
MESH = ["--topology", "ogrid", "--ni", "9", "--nj", "16", "--nk", "5", "--ds", "0.05",
        "--growth", "1.1"]  # fmt: skip


def collar_args(options: dict) -> list[str]:
    args = ["--seed", *map(str, options["seed"])]
    keys = ("n_loop", "ds_surface", "growth_surface", "m_a", "m_b", "ds_wall", "growth_wall", "nk")
    for key in keys:
        args += ["--" + key.replace("_", "-"), str(options[key])]
    if "plane" in options:
        args += ["--plane", *map(str, [*options["plane"][0], *options["plane"][1]])]
    return args


def shore(*args: object):
    from shore.cli import app

    return runner.invoke(app, [str(a) for a in args])


@pytest.fixture
def files(tmp_path) -> dict[str, Path]:
    d = tmp_path / "in"
    d.mkdir()
    trimesh.creation.icosphere(subdivisions=2, radius=1.0).export(str(d / "sphere.stl"))
    box = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
    box.apply_translation((0.0, 0.0, -0.5))
    box.export(str(d / "box.stl"))
    post = trimesh.creation.cylinder(radius=0.5, height=1.6, sections=48)
    post.apply_translation((0.0, 0.0, 0.2))
    post.export(str(d / "post.stl"))
    hbox = trimesh.creation.box(extents=(4.0, 2.0, 1.0))
    hbox.apply_translation((0.0, 1.0, -0.5))
    hbox.update_faces(~np.all(np.isclose(hbox.vertices[hbox.faces][..., 1], 0.0), axis=1))
    hpost = post.copy()
    hpost.update_faces(hpost.triangles_center[:, 1] > 0.0)
    for mesh in (hbox, hpost):
        mesh.remove_unreferenced_vertices()
    hbox.export(str(d / "hbox.stl"))
    hpost.export(str(d / "hpost.stl"))
    return {k: d / f"{k}.stl" for k in ("sphere", "box", "post", "hbox", "hpost")}


def api_collar(a: Path, b: Path, options: dict):
    from shore.project import Project

    p = Project()
    p.load_geometry(a, name=a.stem)
    p.load_geometry(b, name=b.stem)
    p.use_component("collar")
    p.build_collar(b.stem, geometry=a.stem, **options)
    return p


def nodes(mesh) -> list[np.ndarray]:
    return [blk.nodes for blk in mesh.blocks]


def same(a, b) -> bool:
    a, b = nodes(a), nodes(b)
    return len(a) == len(b) and all(np.array_equal(x, y) for x, y in zip(a, b, strict=True))


def load(d: Path):
    from shore.project import Project

    return Project.load(d)


# ---------------------------------------------------------------------------
# L1-L3: shore collar
# ---------------------------------------------------------------------------


class TestCollar:
    def test_l1_collar_into_a_project(self, files, tmp_path):
        r = shore("collar", files["box"], files["post"], *collar_args(COLLAR),
                  "--project", tmp_path / "proj")  # fmt: skip
        assert r.exit_code == 0, r.output
        q = load(tmp_path / "proj")
        assert list(q.components) == ["main", "collar"]
        assert same(q.components["collar"], api_collar(files["box"], files["post"], COLLAR).mesh)

    def test_l2_half_collar(self, files, tmp_path):
        from shore.mesh.face import FaceBC

        r = shore("collar", files["hbox"], files["hpost"], *collar_args(HALF),
                  "--project", tmp_path / "proj")  # fmt: skip
        assert r.exit_code == 0, r.output
        collar = load(tmp_path / "proj").components["collar"]
        assert same(collar, api_collar(files["hbox"], files["hpost"], HALF).mesh)
        assert collar.blocks[0].face("j_lo").bc is FaceBC.SYMMETRY

    def test_l3_collar_geo_output(self, files, tmp_path):
        r = shore("collar", files["box"], files["post"], *collar_args(COLLAR),
                  "-o", tmp_path / "out" / "c")  # fmt: skip
        assert r.exit_code == 0, r.output
        p = api_collar(files["box"], files["post"], COLLAR)
        p.write_geo(tmp_path / "ref" / "c")
        written = sorted(f.name for f in (tmp_path / "out").iterdir())
        assert written == sorted(f.name for f in (tmp_path / "ref").iterdir())
        for name in written:
            assert (tmp_path / "out" / name).read_bytes() == (tmp_path / "ref" / name).read_bytes()


# ---------------------------------------------------------------------------
# L4-L6: components from several commands
# ---------------------------------------------------------------------------


class TestProjects:
    def test_l4_mesh_into_a_project(self, files, tmp_path):
        from shore.project import run_journal

        r = shore("mesh", files["sphere"], *MESH, "-o", tmp_path / "plain" / "m",
                  "--journal", tmp_path / "plain.py")  # fmt: skip
        assert r.exit_code == 0, r.output
        r = shore("mesh", files["sphere"], *MESH, "-o", tmp_path / "out" / "m",
                  "--project", tmp_path / "proj")  # fmt: skip
        assert r.exit_code == 0, r.output
        q = load(tmp_path / "proj")
        assert list(q.components) == ["main"] and list(q.geometries) == ["sphere"]
        assert same(q.components["main"], run_journal(tmp_path / "plain.py").mesh)

    def test_l5_three_components(self, files, tmp_path):
        proj = tmp_path / "proj"
        r = shore("mesh", files["sphere"], *MESH, "-o", tmp_path / "o1" / "m", "--project", proj)
        assert r.exit_code == 0, r.output
        main = [n.copy() for n in nodes(load(proj).components["main"])]
        r = shore("collar", files["box"], files["post"], *collar_args(COLLAR), "--project", proj)
        assert r.exit_code == 0, r.output
        r = shore("primitive", "box", "--min", -3, -3, -3, "--max", 3, 3, 3, "--ni", 5, "--nj", 5,
                  "--nk", 5, "-o", tmp_path / "o2" / "bg.geo", "--project", proj,
                  "--component", "bg")  # fmt: skip
        assert r.exit_code == 0, r.output
        q = load(proj)
        assert list(q.components) == ["main", "collar", "bg"]
        kept = zip(nodes(q.components["main"]), main, strict=True)
        assert all(np.array_equal(a, b) for a, b in kept)
        assert same(q.components["collar"], api_collar(files["box"], files["post"], COLLAR).mesh)

    def test_l6_geometries_reused(self, files, tmp_path):
        proj = tmp_path / "proj"
        for component in ("c1", "c2"):
            r = shore("collar", files["box"], files["post"], *collar_args(COLLAR),
                      "--project", proj, "--component", component)  # fmt: skip
            assert r.exit_code == 0, r.output
        q = load(proj)
        assert q.components["c1"] is not None and q.components["c2"] is not None
        assert sorted(q.geometries) == ["box", "post"]


# ---------------------------------------------------------------------------
# L7: shore xall
# ---------------------------------------------------------------------------

OPTIONS = """\
[run]
grd = "run.grd"
np = 2
multigrid_levels = [3, 1]
numerical_beach = 0.5

[components.main]
free_face_bc = 9
body = 1

[components.main.meta.ogrid]
level = 2
comment = "hull"

[components.collar]
body = 2
boxes = "boxes.json"

[components.collar.proc_meta.collar]
proc = 1
"""


class TestXall:
    def test_l7_xall_with_options(self, files, tmp_path):
        proj = tmp_path / "proj"
        assert shore("mesh", files["sphere"], *MESH, "-o", tmp_path / "o" / "m",
                     "--project", proj).exit_code == 0  # fmt: skip
        assert shore("collar", files["box"], files["post"], *collar_args(COLLAR),
                     "--project", proj).exit_code == 0  # fmt: skip
        box = [[x, y, z] for z in (-1.0, 1.0) for y in (-1.0, 1.0) for x in (0.0, 1.0)]
        (tmp_path / "opts").mkdir()
        (tmp_path / "opts" / "boxes.json").write_text(
            json.dumps([{"associated_block": "collar", "vertices": box}])
        )
        (tmp_path / "opts" / "xall.toml").write_text(OPTIONS)
        r = shore("xall", proj, "--out", tmp_path / "run", "--options",
                  tmp_path / "opts" / "xall.toml")  # fmt: skip
        assert r.exit_code == 0, r.output
        load(proj).write_xall(
            tmp_path / "ref", grd="run.grd", np=2, multigrid_levels=(3, 1), numerical_beach=0.5,
            components={
                "main": {"free_face_bc": 9, "body": 1,
                         "meta": {"ogrid": {"level": 2, "comment": "hull"}}},
                "collar": {"body": 2, "boxes": [{"associated_block": "collar", "vertices": box}],
                           "proc_meta": {"collar": {"proc": 1}}},
            },
        )  # fmt: skip
        for name in ("run.grd", "cc.par", "proc.input"):
            assert (tmp_path / "run" / name).read_bytes() == (tmp_path / "ref" / name).read_bytes()


# ---------------------------------------------------------------------------
# L8: errors
# ---------------------------------------------------------------------------


class TestErrors:
    @pytest.mark.parametrize(
        "options",
        ["[run]\ncolour = 3\n", "[components.nope]\nbody = 1\n"],
        ids=["unknown-key", "unknown-component"],
    )
    def test_l8_bad_options_file(self, files, tmp_path, options):
        proj = tmp_path / "proj"
        assert shore("collar", files["box"], files["post"], *collar_args(COLLAR),
                     "--project", proj).exit_code == 0  # fmt: skip
        (tmp_path / "xall.toml").write_text(options)
        r = shore("xall", proj, "--out", tmp_path / "run", "--options", tmp_path / "xall.toml")
        assert r.exit_code == 1
        assert not (tmp_path / "run" / "proc.input").exists()

    def test_l8_not_a_project(self, tmp_path):
        (tmp_path / "empty").mkdir()
        assert shore("xall", tmp_path / "empty", "--out", tmp_path / "run").exit_code == 1

    def test_l8_component_without_project(self, files, tmp_path):
        r = shore("collar", files["box"], files["post"], *collar_args(COLLAR),
                  "-o", tmp_path / "c", "--component", "x")  # fmt: skip
        assert r.exit_code == 1

    def test_l8_collar_with_nothing_to_write(self, files):
        assert shore("collar", files["box"], files["post"], *collar_args(COLLAR)).exit_code == 1

    def test_l8_mesh_component_without_project(self, files, tmp_path):
        r = shore("mesh", files["sphere"], *MESH, "-o", tmp_path / "m", "--component", "x")
        assert r.exit_code == 1
