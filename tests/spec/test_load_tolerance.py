"""FROZEN acceptance spec: projects load to a tolerance, not bit for bit (ADR 0061).

Frozen 2026-11-01, before any implementation. Implementing work may not edit
this file, ``tests/spec/__init__.py``, ``tests/conftest.py``, the pytest
configuration in ``pyproject.toml`` or ``src/shore/quality.py``. The ONE
permitted edit is deleting an ``@PENDING`` line once that case passes (strict
xfail turns a passing pending case into a failure). ``tests/spec/SHA256SUMS``
holds the hashes, taken with ``@PENDING`` lines removed;
``scripts/check_spec.py`` checks them.

Verdict for acceptance (main CI job): ``pytest tests/spec/test_load_tolerance.py``
reports 4 passed, 0 skipped, 0 xfailed. At freeze: 2 passed, 2 xfailed.

Interface under test
--------------------
``Project.load`` replays the journal and checks each block against its saved
grid. The check allows each node to differ from the saved one by up to 1e-10
of the block's size (its largest bounding-box extent): the last bits of a
replay depend on the machine (the BLAS kernel picked per CPU), so a project
saved on one machine must load on another. A larger difference is refused
with ``MeshInputError``, and the message gives the largest deviation.

Cases
-----
- L1 a saved project whose grids are moved by 1e-13 of their size loads.
- L2 moved by 1e-6 of their size, it is refused (a guard: refused before).
- L3 the refusal names the block and the largest deviation.
- L4 an unchanged project still loads (a guard).

Out of scope of this spec: loading without replay (ADR 0033), which reads the
saved grids as they are.
"""

from __future__ import annotations

import re
import warnings

import numpy as np
import pytest

PENDING = pytest.mark.xfail(strict=True, reason="spec frozen 2026-11-01; implementation pending")


@pytest.fixture
def saved(tmp_path):
    from shore.project import Project
    from shore.shapes import sphere

    sphere(radius=1.0, subdivisions=3).export(tmp_path / "ball.stl")
    p = Project()
    p.load_geometry(tmp_path / "ball.stl")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.build_topology("cubed_sphere", nk=6, ds=0.05, growth=1.1, ni=9, nj=16)
        p.march()
    p.save(tmp_path / "proj")
    return tmp_path / "proj"


def move(project, fraction):
    """Shift every saved node by *fraction* of its block's size."""
    for grid in sorted((project / "grids").rglob("*.npz")):
        with np.load(grid) as g:
            arrays = dict(g)
        nodes = arrays["nodes"]
        size = float(np.ptp(nodes.reshape(-1, 3), axis=0).max())
        arrays["nodes"] = nodes + fraction * size
        np.savez(grid, **arrays)


def load(project):
    from shore.project import Project

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Project.load(project)


class TestLoad:
    def test_l1_tiny_differences_load(self, saved):
        move(saved, 1e-13)
        assert len(load(saved).mesh.blocks) == 6

    def test_l2_real_differences_refused(self, saved):
        from shore.exceptions import MeshInputError

        move(saved, 1e-6)
        with pytest.raises(MeshInputError, match="differs"):
            load(saved)

    def test_l3_message(self, saved):
        from shore.exceptions import MeshInputError

        move(saved, 1e-6)
        with pytest.raises(MeshInputError) as err:
            load(saved)
        msg = str(err.value)
        assert "sub0" in msg
        found = re.search(r"([0-9.]+e-0?[0-9]+)", msg)
        assert found and 1e-7 < float(found.group(1)) < 1e-5

    def test_l4_unchanged_loads(self, saved):
        assert len(load(saved).mesh.blocks) == 6
