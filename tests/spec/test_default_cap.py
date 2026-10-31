"""FROZEN acceptance spec: a 30 deg default cap for the cubed sphere (ADR 0054).

Frozen 2026-10-31, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_default_cap.py``
reports 7 passed, 0 skipped, 0 xfailed. At freeze: 2 passed, 5 xfailed.

Interface under test
--------------------
The cubed sphere's default ``theta_cap`` is pi / 6 (30 deg) in
``Project.build_topology``, ``shore mesh`` and ``Mesh.from_stl``; the
O-grid's stays 0.05 pi. ``Project.build_topology("cubed_sphere", ...)``
without ``theta_cap`` journals the default explicitly
(``theta_cap=0.5235987755982988``), so a journal never depends on a default.
Projects are saved as format 4; a project saved before (format 3 or older)
replays its cubed spheres built without ``theta_cap`` with the old default,
0.05 pi, and journals that value when saved again.

Cases
-----
- C1 the Python default: journal line explicit, the march fold-free and
  without the cap/equator C1 fallback warning (sphere, ni 33, nj 64, nk 25).
- C2 the CLI default: ``shore mesh`` exits 0, journals the cap explicitly,
  and writes C1's grid (to the .geo files' 11 digits).
- C3 ``Mesh.from_stl`` with ``CubedSphereTopology()``: the same grid.
- C4 the O-grid default is unchanged: equal to ``theta_cap=0.05 pi``
  (a guard; holds at freeze).
- C5 a format-3 project (``data/cap9_project.zip``, made at the base
  commit with the old default) loads, its grids checked bit for bit
  (a guard; holds at freeze).
- C6 saved again, that project is format 4 and journals
  ``theta_cap=0.15707963267948966``; it still loads bit for bit.
- C7 a new project saves ``format = 4``.

Out of scope of this spec: bare journal scripts (``shore run``) written
before ADR 0054 without ``theta_cap`` replay with the new default (they
carry no format); OH pins its cap at 45 deg as before.
"""

from __future__ import annotations

import re
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-31; implementation pending")

DATA = Path(__file__).parent / "data"
SIZE = {"nk": 25, "ds": 0.01, "growth": 1.15, "ni": 33, "nj": 64}
PI6 = repr(float(np.pi / 6))
OLD = repr(float(0.05 * np.pi))


@pytest.fixture(scope="module")
def ball(tmp_path_factory):
    from shore.shapes import sphere

    path = tmp_path_factory.mktemp("cap") / "ball.stl"
    sphere(radius=1.0, subdivisions=4).export(path)
    return path


def default_project(stl, topology="cubed_sphere", **extra):
    from shore.project import Project

    p = Project()
    p.load_geometry(stl)
    p.build_topology(topology, **{**SIZE, **extra})
    return p


def nodes(blocks):
    return {b.label: b.nodes.copy() for b in blocks}


def same(a, b):
    assert a.keys() == b.keys()
    for k in a:
        assert np.array_equal(a[k], b[k]), k


@pytest.fixture(scope="module")
def reference(ball):
    """The cubed sphere marched at an explicit 30 deg cap."""
    p = default_project(ball, theta_cap=float(np.pi / 6))
    p.march()
    return nodes(p.mesh.blocks)


class TestCubedSphereDefault:
    @PENDING
    def test_c1_python_default(self, ball, reference):
        from shore.quality import scaled_jacobian

        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            p = default_project(ball)
            p.march()
        assert f"theta_cap={PI6}" in p.journal()
        assert not [w for w in record if "incompatible" in str(w.message)]
        assert min(float(scaled_jacobian(b.nodes).min()) for b in p.mesh.blocks) > 0.0
        same(nodes(p.mesh.blocks), reference)

    @PENDING
    def test_c2_cli_default(self, ball, reference, tmp_path):
        from shore.cli import app
        from shore.io.geo import read_geo

        args = ["--ni", "33", "--nj", "64", "--nk", "25", "--ds", "0.01", "--growth", "1.15"]
        out = ["-o", str(tmp_path / "s"), "--journal", str(tmp_path / "j.py")]
        r = CliRunner().invoke(app, ["mesh", str(ball), *args, *out])
        assert r.exit_code == 0, r.output
        assert f"theta_cap={PI6}" in (tmp_path / "j.py").read_text()
        for label, grid in reference.items():  # .geo holds 11 significant digits
            g = read_geo(tmp_path / f"s_{label}.geo")
            np.testing.assert_allclose(g, grid, rtol=0, atol=1e-10)

    @PENDING
    def test_c3_mesh_from_stl(self, ball, reference):
        from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

        m = Mesh.from_stl(ball, topology=CubedSphereTopology(), ni=33, nj=64, nk=25,
                          spacing_k=Spacing1D(law="geometric", ds=0.01, growth=1.15))  # fmt: skip
        m.march(ds=0.01, growth=1.15)
        same(nodes(m.blocks), reference)

    def test_c4_ogrid_unchanged(self, ball):
        a = default_project(ball, "ogrid", nk=5, nj=17)
        b = default_project(ball, "ogrid", nk=5, nj=17, theta_cap=float(0.05 * np.pi))
        a.march()
        b.march()
        same(nodes(a.mesh.blocks), nodes(b.mesh.blocks))


class TestFormat:
    @pytest.fixture
    def legacy(self, tmp_path):
        with zipfile.ZipFile(DATA / "cap9_project.zip") as z:
            z.extractall(tmp_path)
        return tmp_path / "cap9"

    def test_c5_format_3_loads(self, legacy):
        from shore.project import Project

        assert re.search(r"^format = 3$", (legacy / "project.toml").read_text(), re.M)
        p = Project.load(legacy)  # replays, checks the saved grids bit for bit
        assert len(p.mesh.blocks) == 6

    @PENDING
    def test_c6_resaved_explicit(self, legacy, tmp_path):
        from shore.project import Project

        p = Project.load(legacy)
        p.save(tmp_path / "again")
        text = (tmp_path / "again" / "journal.py").read_text()
        assert f"theta_cap={OLD}" in text
        assert re.search(r"^format = 4$", (tmp_path / "again" / "project.toml").read_text(), re.M)
        q = Project.load(tmp_path / "again")
        same(nodes(q.mesh.blocks), nodes(p.mesh.blocks))

    @PENDING
    def test_c7_new_format(self, ball, tmp_path):
        p = default_project(ball, nk=5, ds=0.05, ni=9, nj=16)
        p.save(tmp_path / "new")
        assert re.search(r"^format = 4$", (tmp_path / "new" / "project.toml").read_text(), re.M)
