"""FROZEN acceptance spec: the smoothed collar volume as the default (ADR 0026).

Frozen 2026-10-28, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml``, or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_collar_default.py``
reports 3 passed, 0 skipped, 0 xfailed. At freeze: 1 passed (a project saved
before the switch still loads: a regression guard), 2 xfailed.

Interface under test
--------------------
``CollarParams.volume`` defaults to ``"smoothed"`` (ADR 0025's construction,
smoothing only near loop corners at the collar's scale); the default is
pinned by ADR 0025's V3 and ADR 0024's S6 (re-frozen with this ADR). Saved
projects become format 3. A ``BuildCollar`` without a volume uses the
default, except while a project saved before format 3 replays: its collars
were built with ``"stations"``, so they replay with it and the saved grids
still match bit for bit.

Cases
-----
A plate (an STL box) pierced by an STL box prism 1 x 0.1 (a loop with four
90 deg corners), n_loop 48, surface steps ``0.01 * 1.2**j`` (8), wall-normal
``1e-4 * 1.3**k`` (12).

- D1 ``Project.save`` writes ``format = 4`` (re-frozen 2026-10-31: format 4
  since ADR 0054; format 3 was this spec's).
- D2 a project saved before format 3 (a collar journalled without a volume,
  ``format = 2``) loads, and its collar is the ``"stations"`` block.
- D3 a new project's collar without a volume is the ``"smoothed"`` block,
  and saves and loads (bit-identical grids).

Out of scope of this spec: the construction itself (ADR 0025), its quality
(ADR 0010, 0013, 0020-0025 specs, run with the new default).
"""

from __future__ import annotations

import re

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-10-28; implementation pending")

STEPS_S = 0.01 * 1.2 ** np.arange(8)
STEPS_K = 1e-4 * 1.3 ** np.arange(12)
SEED = (0.0, 0.06, 0.0)


def write_bodies(d):
    import trimesh

    plate = trimesh.creation.box(extents=(4.0, 4.0, 1.0))
    plate.apply_translation((0.0, 0.0, -0.5))
    plate.export(str(d / "plate.stl"))
    trimesh.creation.box(extents=(1.0, 0.1, 1.0)).export(str(d / "prism.stl"))


def kernel_block(d, volume: str):
    from shore.collar import CollarParams, build_collar
    from shore.geometry import load_geometry

    return build_collar(
        load_geometry(d / "plate.stl"), load_geometry(d / "prism.stl"), np.asarray(SEED),
        n_loop=48, steps_a=STEPS_S, steps_b=STEPS_S, steps_k=STEPS_K,
        params=CollarParams(volume=volume),
    ).nodes  # fmt: skip


def project(d, **kw):
    from shore.project import Project

    p = Project()
    p.load_geometry(d / "plate.stl", name="plate")
    p.build_collar(
        d / "prism.stl", geometry="plate", seed=SEED, n_loop=48, ds_surface=0.01,
        growth_surface=1.2, m_a=8, m_b=8, ds_wall=1e-4, growth_wall=1.3, nk=12, **kw,
    )  # fmt: skip
    return p


class TestDefault:
    @PENDING
    def test_d1_format_3(self, tmp_path):
        write_bodies(tmp_path)
        project(tmp_path).save(tmp_path / "p")
        assert re.search(r"^format = 4$", (tmp_path / "p" / "project.toml").read_text(), re.M)

    def test_d2_old_projects_load(self, tmp_path):  # holds before the switch: a guard
        from shore.project import Project

        write_bodies(tmp_path)
        project(tmp_path, volume="stations").save(tmp_path / "p")
        journal = tmp_path / "p" / "journal.py"
        journal.write_text(journal.read_text().replace(", volume='stations'", ""))
        meta = tmp_path / "p" / "project.toml"
        meta.write_text(re.sub(r"^format = \d+$", "format = 2", meta.read_text(), flags=re.M))
        assert "volume=" not in journal.read_text()
        loaded = Project.load(tmp_path / "p")
        assert np.array_equal(loaded.mesh.blocks[0].nodes, kernel_block(tmp_path, "stations"))

    def test_d3_new_projects_smoothed(self, tmp_path):
        from shore.project import Project

        write_bodies(tmp_path)
        p = project(tmp_path)
        smoothed = kernel_block(tmp_path, "smoothed")
        assert np.array_equal(p.mesh.blocks[0].nodes, smoothed)
        assert not np.array_equal(smoothed, kernel_block(tmp_path, "stations"))
        p.save(tmp_path / "p")
        assert np.array_equal(Project.load(tmp_path / "p").mesh.blocks[0].nodes, smoothed)
