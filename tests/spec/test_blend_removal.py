"""FROZEN acceptance spec: removal of the marcher's normal blending (ADR 0034).

Frozen 2026-10-28, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, ``src/shore/quality.py`` or the data
files below. The ONE permitted edit is deleting an ``@PENDING`` line once
that case passes (strict xfail turns a passing pending case into a failure).
``tests/spec/SHA256SUMS`` holds the hashes, taken with ``@PENDING`` lines
removed; ``scripts/check_spec.py`` checks them.

It replaces ADR 0019's ``test_blend_deprecation.py`` (re-freeze approved
2026-10-28).

Verdict for acceptance (main CI job): ``pytest tests/spec/test_blend_removal.py``
reports 18 passed, 0 skipped, 0 xfailed. At freeze: 7 passed (R7, R8 x 5,
R9a: regression guards), 11 xfailed.

Interface under test
--------------------
- ``blend_normals_k > 0`` is refused everywhere with ``ParameterError``
  (Python) or exit code 1 (CLI), before any layer is marched, with a
  message naming ``blend_normals_k`` (``--blend-normals-k`` on the CLI) and
  the word ``removed``: ``Mesh.march``, ``Project.march``, a journal
  replaying it, ``shore mesh --blend-normals-k K`` and ``[march]
  blend_normals_k = K`` in a ``--config`` file. The mesh is unchanged and
  nothing is journaled.
- ``blend_normals_k = 0`` given explicitly is still accepted, with a notice
  naming it and ``removed`` (``DeprecationWarning`` in Python, a console line
  on the CLI); the result is the same as without it (grids, journal, files).
- Without the option, nothing warns, and every march is what it was before
  the removal (golden ``data/blend_removal_default.npz``, made at the
  removal's parent commit; atol 1e-10 for runner rounding).
- A project saved with blending (``data/blended_project.zip``, made before
  the removal: an O-grid sphere marched with ``blend_normals_k=3``) still
  loads without replay (ADR 0033): its journal keeps the line, its grids are
  the saved ones. Its replayed load is refused with ``ParameterError``.

Cases: the icosphere and NACA 0012 STLs of ``tests/conftest.py``.

- R1 (x 3: O-grid, cubed sphere, C-grid) ``Mesh.march`` refuses k > 0.
- R2 ``Project.march`` and a journal replay refuse k > 0.
- R3 the CLI flag and R4 the config key exit 1.
- R5 (x 2) k = 0 in ``Mesh.march`` / ``Project.march``: notice, same result.
- R6 (x 2) k = 0 as CLI flag / config key: notice, same files.
- R7 no warning without the option.
- R8 (x 5) unblended marches equal the golden: O-grid and cubed sphere,
  algebraic and hyperbolic schemes; C-grid algebraic.
- R9a the blended project loads without replay; R9b its replay is refused.

Out of scope of this spec: ``Body.normals`` (kept: ADR 0005's spec measures
k-line angles against it), OH and CH (they share the cubed-sphere and
C-grid marchers), the message's exact wording.
"""

from __future__ import annotations

import filecmp
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-28; implementation pending")

DATA = Path(__file__).parent / "data"
runner = CliRunner()
OGRID = {"nk": 5, "ds": 0.05, "growth": 1.1, "ni": 9, "nj": 16}
CGRID = {"nk": 6, "ds": 1e-3, "ni_body": 40, "ni_wake": 10, "n_stations": 4, "wake_length": 5.0}
MESH_ARGS = ["--topology", "ogrid", "--ni", "9", "--nj", "16", "--nk", "5", "--ds", "0.05",
             "--growth", "1.1"]  # fmt: skip
REMOVED = r"(?s)blend_normals_k.*removed|removed.*blend_normals_k"


def project(stl: Path, topology: str):
    from shore.project import Project

    p = Project()
    p.load_geometry(stl, skip_hygiene=topology == "cgrid")
    p.build_topology(topology, **(CGRID if topology == "cgrid" else OGRID))
    return p


def stl_for(topology, sphere_stl, naca0012_stl) -> Path:
    return naca0012_stl if topology == "cgrid" else sphere_stl


def notices(record) -> list[str]:
    return [str(w.message) for w in record if issubclass(w.category, DeprecationWarning)]


def removal_notice(messages: list[str]) -> bool:
    return any("blend_normals_k" in m and "removed" in m for m in messages)


def cli_says_removed(output: str) -> bool:
    return "--blend-normals-k" in output and "removed" in output.lower()


def same_tree(a: Path, b: Path) -> bool:
    cmp = filecmp.dircmp(a, b)
    pending = [cmp]
    while pending:
        c = pending.pop()
        if c.left_only or c.right_only:
            return False
        _, mismatch, errors = filecmp.cmpfiles(c.left, c.right, c.common_files, shallow=False)
        if mismatch or errors:
            return False
        pending.extend(c.subdirs.values())
    return True


@pytest.fixture
def blended(tmp_path) -> Path:
    with zipfile.ZipFile(DATA / "blended_project.zip") as z:
        z.extractall(tmp_path / "blended")
    return tmp_path / "blended"


# ---------------------------------------------------------------------------
# R1-R4: k > 0 is refused
# ---------------------------------------------------------------------------


class TestRefused:
    @pytest.mark.parametrize("topology", ["ogrid", "cubed_sphere", "cgrid"])
    def test_r1_mesh_march_refuses(self, sphere_stl, naca0012_stl, topology):
        mesh = project(stl_for(topology, sphere_stl, naca0012_stl), topology).mesh
        before = [b.nodes.copy() for b in mesh.blocks]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ParameterError, match=REMOVED):
                mesh.march(ds=0.05, growth=1.1, blend_normals_k=3)
        assert all(np.array_equal(a, b.nodes) for a, b in zip(before, mesh.blocks, strict=True))

    def test_r2_project_and_journal_refuse(self, sphere_stl, tmp_path):
        from shore.project import run_journal

        p = project(sphere_stl, "ogrid")
        before = [b.nodes.copy() for b in p.mesh.blocks]
        journal = p.journal()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ParameterError, match=REMOVED):
                p.march(blend_normals_k=3)
        assert p.journal() == journal
        assert all(np.array_equal(a, b.nodes) for a, b in zip(before, p.mesh.blocks, strict=True))
        p.save_journal(tmp_path / "j.py")
        source = (tmp_path / "j.py").read_text() + "p.march(blend_normals_k=3)\n"
        (tmp_path / "j.py").write_text(source)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ParameterError, match=REMOVED):
                run_journal(tmp_path / "j.py")

    def test_r3_cli_flag_exits(self, sphere_stl, tmp_path):
        from shore.cli import app

        r = runner.invoke(app, ["mesh", str(sphere_stl), *MESH_ARGS, "-o", str(tmp_path / "m"),
                                "--blend-normals-k", "3"])  # fmt: skip
        assert r.exit_code == 1, r.output
        assert cli_says_removed(r.output), r.output

    def test_r4_config_key_exits(self, sphere_stl, tmp_path):
        from shore.cli import app

        config = f"version = 1\ninput = {str(sphere_stl)!r}\n\n[march]\nblend_normals_k = 3\n"
        (tmp_path / "c.toml").write_text(config)
        r = runner.invoke(app, ["mesh", str(sphere_stl), *MESH_ARGS, "-o", str(tmp_path / "m"),
                                "--config", str(tmp_path / "c.toml")])  # fmt: skip
        assert r.exit_code == 1, r.output
        assert cli_says_removed(r.output), r.output


# ---------------------------------------------------------------------------
# R5-R6: k = 0 is accepted with a notice
# ---------------------------------------------------------------------------


class TestZero:
    def test_r5_mesh_march_zero(self, sphere_stl):
        plain = project(sphere_stl, "ogrid").mesh
        plain.march(ds=0.05, growth=1.1)
        mesh = project(sphere_stl, "ogrid").mesh
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            mesh.march(ds=0.05, growth=1.1, blend_normals_k=0)
        assert removal_notice(notices(record)), notices(record)
        for a, b in zip(plain.blocks, mesh.blocks, strict=True):
            assert np.array_equal(a.nodes, b.nodes)

    def test_r5_project_march_zero(self, sphere_stl):
        plain = project(sphere_stl, "ogrid")
        plain.march()
        p = project(sphere_stl, "ogrid")
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            p.march(blend_normals_k=0)
        assert removal_notice(notices(record)), notices(record)
        assert p.journal() == plain.journal()
        for a, b in zip(plain.mesh.blocks, p.mesh.blocks, strict=True):
            assert np.array_equal(a.nodes, b.nodes)

    def test_r6_cli_flag_zero(self, sphere_stl, tmp_path):
        from shore.cli import app

        base = ["mesh", str(sphere_stl), *MESH_ARGS]
        r = runner.invoke(app, [*base, "-o", str(tmp_path / "a" / "m")])
        assert r.exit_code == 0, r.output
        r = runner.invoke(app, [*base, "-o", str(tmp_path / "b" / "m"), "--blend-normals-k", "0"])
        assert r.exit_code == 0, r.output
        assert cli_says_removed(r.output), r.output
        assert same_tree(tmp_path / "a", tmp_path / "b")

    def test_r6_config_key_zero(self, sphere_stl, tmp_path):
        from shore.cli import app

        config = f"version = 1\ninput = {str(sphere_stl)!r}\n\n[march]\nblend_normals_k = 0\n"
        (tmp_path / "c.toml").write_text(config)
        base = ["mesh", str(sphere_stl), *MESH_ARGS]
        r = runner.invoke(app, [*base, "-o", str(tmp_path / "a" / "m")])
        assert r.exit_code == 0, r.output
        r = runner.invoke(
            app, [*base, "-o", str(tmp_path / "b" / "m"), "--config", str(tmp_path / "c.toml")]
        )
        assert r.exit_code == 0, r.output
        assert cli_says_removed(r.output), r.output
        assert same_tree(tmp_path / "a", tmp_path / "b")


# ---------------------------------------------------------------------------
# R7-R8: without the option nothing changes
# ---------------------------------------------------------------------------

GOLDEN = [
    ("ogrid_algebraic", "ogrid", "algebraic"),
    ("ogrid_hyperbolic", "ogrid", "hyperbolic"),
    ("cubed_sphere_algebraic", "cubed_sphere", "algebraic"),
    ("cubed_sphere_hyperbolic", "cubed_sphere", "hyperbolic"),
    ("cgrid_algebraic", "cgrid", "algebraic"),
]


class TestUnchanged:
    def test_r7_no_warning_without_option(self, sphere_stl):  # already holds: a guard
        p = project(sphere_stl, "ogrid")
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            p.march()
        assert not [m for m in notices(record) if "blend" in m]

    @pytest.mark.parametrize(("name", "topology", "scheme"), GOLDEN)
    def test_r8_marches_equal_golden(self, sphere_stl, naca0012_stl, name, topology, scheme):
        p = project(stl_for(topology, sphere_stl, naca0012_stl), topology)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.march(march_scheme=scheme)
        with np.load(DATA / "blend_removal_default.npz") as gold:
            keys = sorted(k for k in gold.files if k.startswith(f"{name}/"))
            assert len(keys) == len(p.mesh.blocks)
            for key, block in zip(keys, p.mesh.blocks, strict=True):
                np.testing.assert_allclose(block.nodes, gold[key], rtol=0, atol=1e-10)


# ---------------------------------------------------------------------------
# R9: a project saved with blending
# ---------------------------------------------------------------------------


class TestBlendedProject:
    def test_r9a_loads_without_replay(self, blended):  # already holds: a guard
        from shore.project import Project

        p = Project.load(blended, replay=False)
        assert "blend_normals_k=3" in p.journal()
        with np.load(blended / "grids" / "main" / "0000.npz") as g:
            nodes = g["nodes"]
        assert np.array_equal(p.mesh.blocks[0].nodes, nodes)
        with np.load(DATA / "blend_removal_default.npz") as gold:
            assert not np.allclose(nodes, gold["ogrid_algebraic/00"])  # really blended

    def test_r9b_replay_refused(self, blended):
        from shore.project import Project

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ParameterError, match=REMOVED):
                Project.load(blended)
