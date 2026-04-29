"""Smoke-test for the ``split_cubed_sphere.py`` example.

Imports the orchestrator's helpers directly (no subprocess) and verifies:
- the cubed-sphere mesh builds and marches cleanly;
- ``apply_splits`` produces 10 blocks (8 sub-block chunks + 2 caps);
- the adjacency JSON sidecar lands on disk with the expected schema
  (10 blocks, intact j-periodic equatorial ring across the cuts,
  new k-seams between paired chunks, caps unchanged).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(scope="module")
def split_module():
    """Import ``examples/04-split/split_cubed_sphere.py`` as a module."""
    path = EXAMPLES / "04-split" / "split_cubed_sphere.py"
    spec = importlib.util.spec_from_file_location("split_example", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["split_example"] = module
    spec.loader.exec_module(module)
    return module


class TestSplitExampleDefaults:
    def test_pipeline_runs_to_completion(self, split_module, tmp_path):
        m = split_module
        stl = m.make_sphere_stl(tmp_path)
        mesh = m.build_cubed_sphere(stl)
        new_blocks = m.apply_splits(mesh.blocks)
        assert len(mesh.blocks) == 6
        assert len(new_blocks) == 10
        labels = {b.label for b in new_blocks}
        for q in range(4):
            assert f"sub{q}__k0" in labels
            assert f"sub{q}__k1" in labels
        assert "cap_north" in labels
        assert "cap_south" in labels

    def test_adjacency_json_round_trip(self, split_module, tmp_path):
        from shore.io.adjacency import read_adjacency_json

        m = split_module
        stl = m.make_sphere_stl(tmp_path)
        mesh = m.build_cubed_sphere(stl)
        new_blocks = m.apply_splits(mesh.blocks)

        m.write_geo_files(new_blocks, tmp_path / "sphere")
        adj_path = m.write_adjacency(new_blocks, tmp_path / "sphere.adjacency.json")
        assert adj_path.exists()

        payload = read_adjacency_json(adj_path)
        from shore.io.adjacency import ADJACENCY_VERSION

        assert payload["version"] == ADJACENCY_VERSION
        assert len(payload["blocks"]) == 10

        # j-periodic equatorial ring: every sub_q__k0 must share j_lo
        # with sub_(q-1)__k0.j_hi.  This is the cross-block adjacency
        # the splitter is supposed to preserve.
        by_label = {b["label"]: b for b in payload["blocks"]}
        for q in range(4):
            for k_chunk in (0, 1):
                blk = by_label[f"sub{q}__k{k_chunk}"]
                qprev = (q - 1) % 4
                j_lo = blk["faces"]["j_lo"]
                assert j_lo["bc"] == "SHARED"
                assert j_lo["partner"] == {
                    "label": f"sub{qprev}__k{k_chunk}",
                    "face": "j_hi",
                }

        # Internal k-seam: sub_q__k0.k_hi SHARED with sub_q__k1.k_lo.
        for q in range(4):
            k_hi_face = by_label[f"sub{q}__k0"]["faces"]["k_hi"]
            k_lo_face = by_label[f"sub{q}__k1"]["faces"]["k_lo"]
            assert k_hi_face["bc"] == "SHARED"
            assert k_hi_face["partner"] == {
                "label": f"sub{q}__k1",
                "face": "k_lo",
            }
            assert k_lo_face["partner"] == {
                "label": f"sub{q}__k0",
                "face": "k_hi",
            }

        # Caps pass through unchanged: same shape as the input mesh.
        for cap in ("cap_north", "cap_south"):
            shape = by_label[cap]["shape"]
            # Cap blocks are (nc, nc, nk).  Doesn't matter what the
            # exact numbers are; the point is they aren't split.
            assert shape[2] == m._NK

    def test_assembly_vtm_written(self, split_module, tmp_path):
        m = split_module
        stl = m.make_sphere_stl(tmp_path)
        mesh = m.build_cubed_sphere(stl)
        new_blocks = m.apply_splits(mesh.blocks)
        geo_paths = m.write_geo_files(new_blocks, tmp_path / "sphere")
        try:
            vtm = m.write_assembly_vtm(geo_paths, tmp_path / "sphere.vtm")
        except ImportError:
            pytest.skip("pyvista not installed")
        if vtm is None:
            pytest.skip("VTK export skipped at runtime (pyvista import failed)")
        assert vtm.exists()
        # PyVista's MultiBlock.save creates a sister directory of per-block .vts.
        sister = tmp_path / "sphere"
        assert sister.is_dir()
        vts_files = list(sister.glob("*.vts"))
        assert len(vts_files) == 10


class TestSplitExampleConfigDriven:
    def test_loads_canonical_config(self, split_module):
        """The shipped TOML loads to the same 4 paired k-cuts the script
        uses by default."""
        directives = split_module._load_split_directives_from_toml(
            EXAMPLES / "04-split" / "config_split_cubed_sphere.toml"
        )
        assert len(directives) == 4
        labels = sorted(d.label for d in directives)
        assert labels == ["sub0", "sub1", "sub2", "sub3"]
        for d in directives:
            assert d.axis == "k"
            assert d.at == (6,)
