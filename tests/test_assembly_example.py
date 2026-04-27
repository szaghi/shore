"""Smoke-test for the ``assembly_sphere_in_cylinder.py`` example.

Imports the orchestrator's helper functions directly (no subprocess) and
verifies that:
- the in-memory icosphere is written;
- the near-body cubed-sphere mesh marches cleanly to the expected outer
  radius (~0.7 with the script defaults);
- the background flat-caps cylinder spans the expected domain (r=10,
  z in [-10, 50]);
- the bundled ``assembly.vtm`` ends up on disk with the right block keys.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from shore.io.geo import read_geo

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(scope="module")
def assembly_module():
    """Import ``examples/assembly_sphere_in_cylinder.py`` as a module.

    The example file isn't a package; we load it from its path so the
    test doesn't depend on cwd or PYTHONPATH gymnastics.
    """
    path = EXAMPLES / "assembly_sphere_in_cylinder.py"
    spec = importlib.util.spec_from_file_location("assembly_example", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["assembly_example"] = module
    spec.loader.exec_module(module)
    return module


class TestAssemblyDefaults:
    def test_pipeline_runs_to_completion(self, assembly_module, tmp_path):
        m = assembly_module
        stl = m.make_sphere_stl(tmp_path)
        nb = m.build_near_body(stl, tmp_path / "fuselage")
        bg = m.build_background(tmp_path / "background")
        # 6 near-body files + 5 background files.
        assert len(nb) == 6
        assert len(bg) == 5
        for p in nb + bg:
            assert p.exists()

    def test_near_body_outer_shell_radius(self, assembly_module, tmp_path):
        m = assembly_module
        stl = m.make_sphere_stl(tmp_path)
        m.build_near_body(stl, tmp_path / "fuselage")
        # sub0 outer shell radius should be ~0.7 (= 0.5 + 0.2 thickness).
        nb_sub0 = read_geo(tmp_path / "fuselage_sub0.geo")
        r_outer = np.linalg.norm(nb_sub0[:, :, -1, :], axis=-1)
        assert 0.65 < r_outer.min() < 0.75
        assert 0.65 < r_outer.max() < 0.75

    def test_background_extent(self, assembly_module, tmp_path):
        m = assembly_module
        m.build_background(tmp_path / "background")
        # Outer arc on r=10.0 exactly; axial extent -10..50.
        bg_e = read_geo(tmp_path / "background_sub_e.geo")
        outer = bg_e[-1, :, :, :2].reshape(-1, 2)
        radii = np.linalg.norm(outer, axis=-1)
        assert np.allclose(radii, 10.0)
        z = bg_e[..., 2]
        assert np.isclose(z.min(), -10.0)
        assert np.isclose(z.max(), 50.0)

    def test_assembly_vtm_written(self, assembly_module, tmp_path):
        m = assembly_module
        stl = m.make_sphere_stl(tmp_path)
        nb = m.build_near_body(stl, tmp_path / "fuselage")
        bg = m.build_background(tmp_path / "background")
        try:
            vtm = m.write_assembly_vtm(nb, bg, tmp_path / "assembly.vtm")
        except ImportError:
            pytest.skip("pyvista not installed")
        assert vtm.exists()
        # PyVista's MultiBlock.save also creates a sister directory of
        # per-block .vts files.
        sister = tmp_path / "assembly"
        assert sister.is_dir()
        # 11 blocks → 11 .vts files inside the sister directory.
        vts_files = list(sister.glob("*.vts"))
        assert len(vts_files) == 11
