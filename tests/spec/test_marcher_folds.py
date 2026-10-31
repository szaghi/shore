"""FROZEN acceptance spec: no silent folds in the O-grid and C-grid marches (ADR 0056).

Frozen 2026-10-31, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_marcher_folds.py``
reports 7 passed, 0 skipped, 0 xfailed. At freeze: 4 passed, 3 xfailed.

Interface under test
--------------------
ADR 0052's guard, for the other marched topologies: the O-grid march
(``Mesh._march_ogrid``, cells across the j wrap included) and the C-grid
marches (``_march_cgrid_three``, which CH's wrap uses, and the 1-block
``_march_cgrid``) refuse a layer with a folded cell (minimum corner scaled
Jacobian <= 0, ``shore.quality``) with ``MarchingError`` instead of writing
it. Their old guard tests one corner of each cell against an absolute
tolerance, which let the folds below through. Marches without folds are
unchanged, bit for bit.

Cases
-----
- G1 an O-grid on a box (2 x 1 x 0.6; ni 25, nj 33, nk 20, ds 0.01): 4
  folded cells before (min -0.004), now refused.
- G2 a sharp-TE C-grid on a NACA 0012 wing (nk 20, ds 2e-3, growth 1.2):
  4 folded cells before (min -0.98), now refused.
- G3 the same from ``shore mesh --topology ogrid``: exit 1, no ``.geo``.
- G4 (x 3) clean cases march, every cell of the marched blocks positive:
  O-grid on a sphere, C-grid with growth 1.1, CH (its wrap).
- G5 a clean O-grid and a clean C-grid equal their grids marched before
  ADR 0056 (``data/marcher_folds_clean.npz``, made at the base commit).

Out of scope of this spec: CH's channel and outflow blocks (filled, not
marched; their handedness is ADR 0057's), mending folds instead of refusing.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.exceptions import MarchingError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-31; implementation pending")

DATA = Path(__file__).parent / "data"
OGRID_BOX = {"nk": 20, "ds": 0.01, "growth": 1.15, "ni": 25, "nj": 33}
CGRID = {"ni_body": 60, "ni_wake": 15, "n_stations": 4, "wake_length": 5.0, "ds": 2e-3}
CH = {"chord_min": -2.0, "chord_max": 15.0, "normal_min": -2.5, "normal_max": 2.5, "nk_channel": 4}
MARCHED = {"upper", "body", "lower", "ogrid", "main"}


@pytest.fixture(scope="module")
def bodies(tmp_path_factory):
    from shore.profiles import naca
    from shore.shapes import box, extrude_section, sphere

    d = tmp_path_factory.mktemp("marcher_folds")
    box(extents=(2.0, 1.0, 0.6)).export(d / "box.stl")
    sphere(radius=1.0, subdivisions=3).export(d / "sphere.stl")
    wing = extrude_section(naca("0012", 60, closed_te=True), span=3.0, n_stations=5)
    wing.export(d / "wing.stl")
    return d


def project(stl, topology, **params):
    from shore.project import Project

    p = Project()
    p.load_geometry(stl, skip_hygiene=topology in ("cgrid", "ch"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.build_topology(topology, **params)
    return p


def march(p):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.march()


def refused(p):
    before = [b.nodes.copy() for b in p.mesh.blocks]
    with pytest.raises(MarchingError, match=r"(?i)fold"):
        march(p)
    for b, old in zip(p.mesh.blocks, before, strict=True):
        assert np.array_equal(b.nodes, old)


def min_marched(p):
    from shore.quality import scaled_jacobian

    blocks = [b for b in p.mesh.blocks if not b.label.startswith("chan_")]
    return min(float(scaled_jacobian(b.nodes).min()) for b in blocks)


class TestRefused:
    @PENDING
    def test_g1_ogrid_box(self, bodies):
        refused(project(bodies / "box.stl", "ogrid", **OGRID_BOX))

    @PENDING
    def test_g2_cgrid_wing(self, bodies):
        refused(project(bodies / "wing.stl", "cgrid", nk=20, growth=1.2, **CGRID))

    @PENDING
    def test_g3_cli(self, bodies, tmp_path):
        from shore.cli import app

        r = CliRunner().invoke(app, [
            "mesh", str(bodies / "box.stl"), "--topology", "ogrid", "--ni", "25", "--nj", "33",
            "--nk", "20", "--ds", "0.01", "--growth", "1.15", "-o", str(tmp_path / "b"),
        ])  # fmt: skip
        assert r.exit_code == 1
        assert "fold" in r.output.lower()
        assert not list(tmp_path.glob("*.geo"))


class TestAccepted:
    @pytest.mark.parametrize("name", ["ogrid", "cgrid", "ch"])
    def test_g4_clean_cases(self, bodies, name):
        if name == "ogrid":
            p = project(bodies / "sphere.stl", "ogrid", nk=12, ds=0.02, growth=1.15, ni=17, nj=33)
        elif name == "cgrid":
            p = project(bodies / "wing.stl", "cgrid", nk=20, growth=1.1, **CGRID)
        else:
            p = project(bodies / "wing.stl", "ch", nk=12, growth=1.1, options=CH, **CGRID)
        march(p)
        assert min_marched(p) > 0.0

    def test_g5_unchanged(self, bodies):
        og = project(bodies / "sphere.stl", "ogrid", nk=12, ds=0.02, growth=1.15, ni=17, nj=33)
        cg = project(bodies / "wing.stl", "cgrid", nk=20, growth=1.1, **CGRID)
        march(og)
        march(cg)
        with np.load(DATA / "marcher_folds_clean.npz") as gold:
            got = {f"ogrid/{b.label}": b.nodes for b in og.mesh.blocks}
            got |= {f"cgrid/{b.label}": b.nodes for b in cg.mesh.blocks}
            assert set(gold.files) == set(got)
            for key, nodes in got.items():
                assert np.array_equal(nodes, gold[key]), key
