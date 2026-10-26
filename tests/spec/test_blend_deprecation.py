"""FROZEN acceptance spec: deprecation of the marcher's normal blending (ADR 0019).

Frozen 2026-10-26, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job):
``pytest tests/spec/test_blend_deprecation.py`` reports 8 passed, 0 skipped,
0 xfailed. At freeze: 1 passed (B5, a regression guard: no warning without
blending), 7 xfailed.

Interface under test
--------------------
- ``Mesh.march(..., blend_normals_k=k)`` with ``k > 0`` emits a
  ``DeprecationWarning`` whose message names ``blend_normals_k`` and
  ``ADR 0019``, and otherwise marches as before (ADR 0015's vertex normals).
  So do ``Project.march`` and a journal replaying it.
- ``shore mesh --blend-normals-k K`` (K > 0), and a ``--config`` file with
  ``[march] blend_normals_k = K``, print a deprecation notice naming
  ``--blend-normals-k`` on the console and still mesh (exit code 0).
  (``DeprecationWarning`` from library code is hidden by Python's default
  filters, so the CLI says it itself.)
- On C-grid topologies (``cgrid``, and ``ch``, which marches its wrap the
  same way) ``blend_normals_k > 0`` raises ``ParameterError`` before any
  layer is marched: the blend gives the wake-cut nodes, which are off the
  body, the normals of the nearest trailing-edge point, and every normal
  source tried folds there (ADR 0019). The project's mesh is unchanged.
- ``blend_normals_k = 0`` (the default) warns about nothing.

Cases: an icosphere STL (cubed sphere, O-grid) and the NACA 0012 STL of
``tests/conftest.py`` (C-grid).

- B1 ``Mesh.march`` warns; B2 ``Project.march`` and a journal replay warn.
- B3 the CLI flag and B4 the config key print a notice; the run succeeds.
- B5 no warning at 0 (regression guard).
- B6 C-grid with blending raises (``Project.march``), mesh unchanged; B7
  the CLI exits 1 with the reason.

Out of scope of this spec: removing the option (a later release), the
direction quality of the blend (ADR 0019 measured it), CH beyond sharing
the C-grid code path (unit test).
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-26; implementation pending")

runner = CliRunner()
OGRID = {"nk": 5, "ds": 0.05, "growth": 1.1, "ni": 9, "nj": 16}
CGRID = {"nk": 6, "ds": 1e-3, "ni_body": 40, "ni_wake": 10, "n_stations": 4, "wake_length": 5.0}


def sphere_project(sphere_stl, topology="ogrid"):
    from shore.project import Project

    p = Project()
    p.load_geometry(sphere_stl)
    if topology == "ogrid":
        p.build_topology("ogrid", **OGRID)
    else:
        p.build_topology("cubed_sphere", nk=5, ds=0.05, growth=1.1, ni=9, nj=16)
    return p


def deprecations(record) -> list[str]:
    return [str(w.message) for w in record if issubclass(w.category, DeprecationWarning)]


# ---------------------------------------------------------------------------
# B1-B2: the warning
# ---------------------------------------------------------------------------


class TestWarning:
    @PENDING
    @pytest.mark.parametrize("topology", ["ogrid", "cubed_sphere"])
    def test_b1_mesh_march_warns(self, sphere_stl, topology):
        mesh = sphere_project(sphere_stl, topology).mesh
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            mesh.march(ds=0.05, growth=1.1, blend_normals_k=3)
        messages = deprecations(record)
        assert any("blend_normals_k" in m and "ADR 0019" in m for m in messages), messages

    @PENDING
    def test_b2_project_march_and_replay_warn(self, sphere_stl, tmp_path):
        from shore.project import run_journal

        p = sphere_project(sphere_stl)
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            p.march(blend_normals_k=3)
        assert any("ADR 0019" in m for m in deprecations(record))
        p.save_journal(tmp_path / "j.py")
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            q = run_journal(tmp_path / "j.py")
        assert any("ADR 0019" in m for m in deprecations(record))
        for a, b in zip(q.mesh.blocks, p.mesh.blocks, strict=True):
            assert np.array_equal(a.nodes, b.nodes)


# ---------------------------------------------------------------------------
# B3-B4: the CLI says it
# ---------------------------------------------------------------------------

MESH_ARGS = ["--topology", "ogrid", "--ni", "9", "--nj", "16", "--nk", "5", "--ds", "0.05",
             "--growth", "1.1"]  # fmt: skip


class TestCli:
    @PENDING
    def test_b3_cli_flag_notice(self, sphere_stl, tmp_path):
        from shore.cli import app

        r = runner.invoke(app, ["mesh", str(sphere_stl), *MESH_ARGS, "-o", str(tmp_path / "m"),
                                "--blend-normals-k", "3"])  # fmt: skip
        assert r.exit_code == 0, r.output
        assert "--blend-normals-k" in r.output and "deprecated" in r.output.lower()

    @PENDING
    def test_b4_config_key_notice(self, sphere_stl, tmp_path):
        from shore.cli import app

        config = f"version = 1\ninput = {str(sphere_stl)!r}\n\n[march]\nblend_normals_k = 3\n"
        (tmp_path / "c.toml").write_text(config)
        r = runner.invoke(app, ["mesh", str(sphere_stl), *MESH_ARGS, "-o", str(tmp_path / "m"),
                                "--config", str(tmp_path / "c.toml")])  # fmt: skip
        assert r.exit_code == 0, r.output
        assert "--blend-normals-k" in r.output and "deprecated" in r.output.lower()


# ---------------------------------------------------------------------------
# B5: nothing at 0
# ---------------------------------------------------------------------------


class TestDefault:
    def test_b5_no_warning_without_blending(self, sphere_stl):
        p = sphere_project(sphere_stl)
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            p.march()
        assert not [m for m in deprecations(record) if "blend" in m]


# ---------------------------------------------------------------------------
# B6-B7: C-grid refuses
# ---------------------------------------------------------------------------


class TestCgrid:
    @PENDING
    def test_b6_cgrid_blend_raises(self, naca0012_stl):
        from shore.project import Project

        p = Project()
        p.load_geometry(naca0012_stl, skip_hygiene=True)
        p.build_topology("cgrid", **CGRID)
        before = [b.nodes.copy() for b in p.mesh.blocks]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ParameterError, match="blend_normals_k"):
                p.march(blend_normals_k=3)
        assert all(np.array_equal(a, b.nodes) for a, b in zip(before, p.mesh.blocks, strict=True))

    @PENDING
    def test_b7_cgrid_cli_exits(self, naca0012_stl, tmp_path):
        from shore.cli import app

        args = ["mesh", str(naca0012_stl), "--topology", "cgrid", "--ni-body", "40",
                "--ni-wake", "10", "--n-stations", "4", "--wake-length", "5.0", "--nk", "6",
                "--ds", "1e-3", "-o", str(tmp_path / "w"), "--blend-normals-k", "3"]  # fmt: skip
        r = runner.invoke(app, args)
        assert r.exit_code == 1
        assert "blend_normals_k" in r.output or "--blend-normals-k" in r.output
