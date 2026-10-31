"""FROZEN acceptance spec: no silent folds in the cubed-sphere march (ADR 0052).

Frozen 2026-10-30, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_cap_folds.py``
reports 10 passed, 0 skipped, 0 xfailed. At freeze: 6 passed, 4 xfailed.

Interface under test
--------------------
``Project.march()`` (and ``shore mesh``) on a cubed sphere: a march that
would leave a folded cell (minimum corner scaled Jacobian <= 0, the frozen
``shore.quality.scaled_jacobian``), anywhere in any block, the cap-corner
cells at the 3-block junctions included, raises ``MarchingError`` instead of
returning the grid. The message says the cell is folded and names the
remedies (``theta_cap``, ``hyperbolic``). The failed march is not recorded
and leaves the mesh as it was (``March`` restores on failure); the CLI exits
1 and writes no ``.geo``. Marches without folds are unchanged, bit for bit.

Cases
-----
- F1 the default cap (theta_cap 0.05 pi) on a sphere, algebraic scheme: 72
  folded cap-corner cells before ADR 0052 (min -0.83), now refused.
- F2 the same from ``shore mesh``: exit 1, no ``.geo`` written.
- F3 a capsule with a 45 deg cap, hyperbolic scheme: 8 folds (min -0.006).
- F4 a sphere at 1e-3 first cells: only 4 folded cells (min -0.085).
- F5 (x 4) clean cases still march, with every cell positive.
- F6 a clean case equals the grid marched before ADR 0052
  (``data/cap_folds_clean.npz``, made at the base commit), bit for bit.
- F7 the OH hybrid (cubed-sphere wrap, cap pinned at 45 deg) on a sphere
  still marches, every cell positive.

Out of scope of this spec: the default ``theta_cap`` (unchanged here; ADR
0052 proposes 30 deg separately, a re-freeze of ADR 0034's R8 golden), the
other topologies' guards (single-corner, absolute tolerance), smoothing
that keeps cells valid instead of refusing.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.exceptions import MarchingError

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-30; implementation pending")

DATA = Path(__file__).parent / "data"
FOLD = r"(?i)fold"


@pytest.fixture(scope="module")
def bodies(tmp_path_factory):
    from shore.shapes import capsule, sphere

    d = tmp_path_factory.mktemp("cap_folds")
    sphere(radius=1.0, subdivisions=4).export(d / "sphere.stl")
    capsule(radius=0.5, height=2.0).export(d / "capsule.stl")
    return d


def project(stl, theta_deg=None, ni=33, nj=64, nk=25, ds=0.01):
    from shore.project import Project

    p = Project()
    p.load_geometry(stl)
    params = {"nk": nk, "ds": ds, "growth": 1.15, "ni": ni, "nj": nj}
    if theta_deg is not None:
        params["theta_cap"] = float(np.deg2rad(theta_deg))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the cap/equator C1 fallback note
        p.build_topology("cubed_sphere", **params)
    return p


def march(p, scheme="algebraic"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.march(march_scheme=scheme)


def min_scaled(p):
    from shore.quality import scaled_jacobian

    return min(float(scaled_jacobian(b.nodes).min()) for b in p.mesh.blocks)


def refused(p, scheme):
    before = [b.nodes.copy() for b in p.mesh.blocks]
    n = len(p.history)
    with pytest.raises(MarchingError, match=FOLD) as err:
        march(p, scheme)
    assert len(p.history) == n  # not recorded
    for b, old in zip(p.mesh.blocks, before, strict=True):
        assert np.array_equal(b.nodes, old)  # restored
    return str(err.value)


class TestRefused:
    def test_f1_default_cap_algebraic(self, bodies):
        msg = refused(project(bodies / "sphere.stl"), "algebraic")
        assert "theta_cap" in msg and "hyperbolic" in msg

    def test_f2_cli(self, bodies, tmp_path):
        from shore.cli import app

        r = CliRunner().invoke(
            app,
            ["mesh", str(bodies / "sphere.stl"), "--ni", "33", "--nj", "64", "--nk", "25",
             "--ds", "0.01", "--growth", "1.15", "-o", str(tmp_path / "s")],
        )  # fmt: skip
        assert r.exit_code == 1
        assert "fold" in r.output.lower()
        assert not list(tmp_path.glob("*.geo"))

    def test_f3_capsule_45_hyperbolic(self, bodies):
        refused(project(bodies / "capsule.stl", 45), "hyperbolic")

    def test_f4_few_folds(self, bodies):
        refused(project(bodies / "sphere.stl", None, ni=41, nj=64, nk=31, ds=1e-3), "algebraic")


class TestAccepted:
    @pytest.mark.parametrize(
        ("body", "theta", "scheme"),
        [("sphere", 30, "algebraic"), ("sphere", 35, "hyperbolic"),
         ("capsule", 30, "algebraic"), ("capsule", 35, "hyperbolic")],
    )  # fmt: skip
    def test_f5_clean_cases_march(self, bodies, body, theta, scheme):
        p = project(bodies / f"{body}.stl", theta)
        march(p, scheme)
        assert min_scaled(p) > 0.0

    def test_f6_unchanged(self, bodies):
        p = project(bodies / "sphere.stl", 35, ni=25, nj=48, nk=21, ds=0.02)
        march(p)
        with np.load(DATA / "cap_folds_clean.npz") as gold:
            assert set(gold.files) == {b.label for b in p.mesh.blocks}
            for b in p.mesh.blocks:
                assert np.array_equal(b.nodes, gold[b.label])

    def test_f7_oh_hybrid(self, bodies):
        from shore.project import Project

        p = Project()
        p.load_geometry(bodies / "sphere.stl")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            box = {"box_vmin": (-4.0,) * 3, "box_vmax": (4.0,) * 3, "nk_channel": 5}
            p.build_topology("oh", nk=9, ds=0.02, growth=1.15, ni=17, nj=32, options=box)
            p.march()
        assert min_scaled(p) > 0.0
