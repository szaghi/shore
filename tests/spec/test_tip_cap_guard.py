"""FROZEN acceptance spec: no silent folds in the wing-tip caps (ADR 0058).

Frozen 2026-10-31, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_tip_cap_guard.py``
reports 6 passed, 0 skipped, 0 xfailed. At freeze: 1 passed, 5 xfailed.

Interface under test
--------------------
The C-grid's wing-tip caps (``cap_tips`` on ``CGridTopology.build``, the
``--cap-tips`` flag of ``shore mesh``; ``CGridTipCap``, ``Project.build_tip_cap``
and ``shore cap-tip``) are built, not marched, so no marcher guard sees them.
A cap build either returns blocks whose every cell is valid (minimum corner
scaled Jacobian > 0, ``shore.quality``; at least one cell in every block) or
raises ``ParameterError`` naming a cap block and the fold, or its missing
cells; a refused build leaves the project as it was, and the CLI exits 1
without writing ``.geo`` files. The property holds for any construction: the
caps of the base commit always fold (ADR 0057), so they are refused; a valid
construction passes the same cases.

Cases
-----
- P1 (x 2) ``build_topology("cgrid", cap_tips="both")`` on a sharp and a
  blunt NACA 0012 wing.
- P2 ``Project.build_tip_cap`` from a C-grid body block ``.geo``.
- P3 ``shore cap-tip`` on the same file.
- P4 ``shore mesh --topology cgrid --cap-tips both``.
- P5 a C-grid without caps is unchanged, bit for bit
  (``data/tip_cap_guard_clean.npz``, made at the base commit).

Out of scope of this spec: what a valid cap looks like (ADR 0059).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.exceptions import ParameterError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-31; implementation pending")

DATA = Path(__file__).parent / "data"
CGRID = {"nk": 12, "ds": 2e-3, "ni_body": 40, "ni_wake": 10, "n_stations": 4, "wake_length": 5.0}
CLI_CGRID = [
    "--topology", "cgrid", "--nk", "12", "--ds", "2e-3", "--ni-body", "40", "--ni-wake", "10",
    "--n-stations", "4", "--wake-length", "5",
]  # fmt: skip


@pytest.fixture(scope="module")
def wings(tmp_path_factory):
    from shore.profiles import naca
    from shore.shapes import extrude_section

    d = tmp_path_factory.mktemp("tip_cap_guard")
    for name, closed in (("sharp", True), ("blunt", False)):
        wing = extrude_section(naca("0012", 60, closed_te=closed), span=3.0, n_stations=5)
        wing.export(d / f"{name}.stl")
    return d


@pytest.fixture(scope="module")
def body_geo(wings):
    """The body block of a sharp C-grid, written as ``shore mesh`` writes it."""
    from shore.project import Project

    p = Project()
    p.load_geometry(wings / "sharp.stl", skip_hygiene=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.build_topology("cgrid", **CGRID)
        p.march()
    p.write_geo(wings / "w", per_block=True)
    (path,) = wings.glob("w*body*.geo")
    return path


def valid_or_refused(call, project):
    """The cap blocks are valid, or the build was refused and changed nothing."""
    from shore.quality import scaled_jacobian

    before = project.mesh
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            call()
    except ParameterError as exc:
        msg = str(exc)
        assert "cap" in msg
        assert any(w in msg.lower() for w in ("fold", "no cells", "scaled jacobian"))
        assert project.mesh is before
        return
    caps = [b for b in project.mesh.blocks if b.label.startswith("cap")]
    assert caps
    for b in caps:
        assert min(b.nodes.shape[:3]) >= 2, b.label
        assert float(scaled_jacobian(b.nodes).min()) > 0.0, b.label


def geo_files_valid(directory):
    from shore.io.geo import read_geo
    from shore.quality import scaled_jacobian

    files = sorted(Path(directory).glob("*.geo"))
    assert files
    for f in files:
        nodes = read_geo(f)
        assert min(nodes.shape[:3]) >= 2, f.name
        assert float(scaled_jacobian(nodes).min()) > 0.0, f.name


class TestCaps:
    @pytest.mark.parametrize("te", ["sharp", "blunt"])
    def test_p1_conformal_caps(self, wings, te):
        from shore.project import Project

        p = Project()
        p.load_geometry(wings / f"{te}.stl", skip_hygiene=True)
        extra = {"blocks": 4} if te == "blunt" else {}
        valid_or_refused(lambda: p.build_topology("cgrid", cap_tips="both", **CGRID, **extra), p)

    def test_p2_tip_cap_api(self, body_geo):
        from shore.project import Project

        p = Project()
        valid_or_refused(lambda: p.build_tip_cap(body_geo, tip="jlo"), p)

    def test_p3_tip_cap_cli(self, body_geo, tmp_path):
        from shore.cli import app

        out = tmp_path / "cap"
        out.mkdir()
        args = ["cap-tip", str(body_geo), "--tip", "jhi", "-o", str(out / "c")]
        r = CliRunner().invoke(app, args)
        if r.exit_code == 0:
            geo_files_valid(out)
        else:
            assert r.exit_code == 1
            assert not list(out.glob("*.geo"))

    def test_p4_mesh_cli(self, wings, tmp_path):
        from shore.cli import app

        out = tmp_path / "mesh"
        out.mkdir()
        args = ["mesh", str(wings / "sharp.stl"), *CLI_CGRID, "--cap-tips", "both"]
        r = CliRunner().invoke(app, [*args, "-o", str(out / "w")])
        if r.exit_code == 0:
            geo_files_valid(out)
        else:
            assert r.exit_code == 1
            assert not list(out.glob("*.geo"))


class TestUnchanged:
    def test_p5_no_caps_bit_for_bit(self, wings):
        from shore.project import Project

        p = Project()
        p.load_geometry(wings / "sharp.stl", skip_hygiene=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.build_topology("cgrid", **CGRID)
            p.march()
        gold = np.load(DATA / "tip_cap_guard_clean.npz")
        assert sorted(gold.files) == sorted(b.label for b in p.mesh.blocks)
        for b in p.mesh.blocks:
            assert np.array_equal(b.nodes, gold[b.label]), b.label
