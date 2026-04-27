"""Tests for ``shore.split.split_blocks`` and the JSON adjacency sidecar."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.cli import app
from shore.exceptions import ParameterError
from shore.io.adjacency import (
    ADJACENCY_VERSION,
    read_adjacency_json,
    write_adjacency_json,
)
from shore.io.geo import read_geo
from shore.mesh.face import FaceBC
from shore.primitive import BoxMesh
from shore.split import BlockSplit, split_blocks

runner = CliRunner()


# ---------------------------------------------------------------------------
# BlockSplit dataclass validation
# ---------------------------------------------------------------------------


class TestBlockSplitValidation:
    def test_axis_must_be_i_j_k(self):
        with pytest.raises(ParameterError, match="axis must be"):
            BlockSplit(label="x", axis="t", at=[4])

    def test_at_must_be_non_empty(self):
        with pytest.raises(ParameterError, match="not a split"):
            BlockSplit(label="x", axis="i", at=[])

    def test_at_must_be_positive(self):
        with pytest.raises(ParameterError, match="must be >= 1"):
            BlockSplit(label="x", axis="i", at=[0])

    def test_dedupes_and_sorts(self):
        sp = BlockSplit(label="x", axis="i", at=[16, 8, 16, 4])
        assert sp.at == (4, 8, 16)


# ---------------------------------------------------------------------------
# split_blocks: single-block cases
# ---------------------------------------------------------------------------


def _box(ni=12, nj=8, nk=6) -> list:
    return BoxMesh.build(vmin=(0, 0, 0), vmax=(10, 5, 3), ni=ni, nj=nj, nk=nk).blocks


class TestSplitSingleBlock:
    def test_single_cut_produces_two_chunks(self):
        new = split_blocks(_box(), [BlockSplit("box", "i", [6])])
        assert len(new) == 2
        assert new[0].label == "box__i0"
        assert new[1].label == "box__i1"
        # Shape: 7 + 6 = 13; the shared vertex at i=6 is on both chunks.
        assert new[0].nodes.shape == (7, 8, 6, 3)
        assert new[1].nodes.shape == (6, 8, 6, 3)

    def test_two_cuts_produce_three_chunks(self):
        # nj=8 → vertices 0..7.  Cuts at [3, 5] give chunks
        #   j0: vertices 0..3 (4 nodes)
        #   j1: vertices 3..5 (3 nodes; the cut vertex is shared)
        #   j2: vertices 5..7 (3 nodes)
        new = split_blocks(_box(), [BlockSplit("box", "j", [3, 5])])
        assert [b.label for b in new] == ["box__j0", "box__j1", "box__j2"]
        assert new[0].nodes.shape == (12, 4, 6, 3)
        assert new[1].nodes.shape == (12, 3, 6, 3)
        assert new[2].nodes.shape == (12, 3, 6, 3)

    def test_node_coords_preserved_bit_exactly(self):
        original = _box()
        new = split_blocks(original, [BlockSplit("box", "k", [3])])
        # The two chunks together must reproduce the original block.
        rebuilt = np.concatenate([new[0].nodes, new[1].nodes[:, :, 1:, :]], axis=2)
        np.testing.assert_array_equal(rebuilt, original[0].nodes)

    def test_internal_seam_shared_with_partner(self):
        new = split_blocks(_box(), [BlockSplit("box", "i", [6])])
        a, b = new
        a_hi = a.face("i_hi")
        b_lo = b.face("i_lo")
        assert a_hi.bc == FaceBC.SHARED
        assert b_lo.bc == FaceBC.SHARED
        assert a_hi.partner is b_lo
        assert b_lo.partner is a_hi
        # Bit-exact node coincidence at the seam.
        np.testing.assert_array_equal(a.nodes[-1, :, :, :], b.nodes[0, :, :, :])

    def test_external_face_bcs_preserved_per_chunk(self):
        new = split_blocks(_box(), [BlockSplit("box", "i", [6])])
        a, b = new
        # Wall (k_lo) and free (k_hi, j_lo, j_hi) inherited on every chunk.
        for blk in (a, b):
            assert blk.face("k_lo").bc == FaceBC.WALL
            assert blk.face("j_lo").bc == FaceBC.FREE
            assert blk.face("j_hi").bc == FaceBC.FREE
            assert blk.face("k_hi").bc == FaceBC.FREE
        # i_lo only on the first chunk; i_hi only on the last (these are
        # the two "external" faces along the split axis).
        assert a.face("i_lo").bc == FaceBC.FREE
        assert b.face("i_hi").bc == FaceBC.FREE
        # Internal seam: a.i_hi and b.i_lo are SHARED.
        assert a.face("i_hi").bc == FaceBC.SHARED
        assert b.face("i_lo").bc == FaceBC.SHARED


class TestSplitErrorPaths:
    def test_unknown_label(self):
        with pytest.raises(ParameterError, match="no block named"):
            split_blocks(_box(), [BlockSplit("missing", "i", [4])])

    def test_duplicate_split_label_in_one_call(self):
        with pytest.raises(ParameterError, match="appears in multiple"):
            split_blocks(
                _box(),
                [
                    BlockSplit("box", "i", [4]),
                    BlockSplit("box", "j", [3]),
                ],
            )

    def test_cut_index_out_of_range(self):
        with pytest.raises(ParameterError, match="out of range"):
            split_blocks(_box(ni=12), [BlockSplit("box", "i", [12])])


# ---------------------------------------------------------------------------
# split_blocks: cross-block adjacency rewiring
# ---------------------------------------------------------------------------


class TestSplitAdjacencyAcrossBlocks:
    def test_unsplit_partner_rejects_perpendicular_split(self, sphere_stl):
        """Splitting sub0 along k while leaving sub3 unsplit should fail
        (the j-periodic seam between sub0 and sub3 crosses the k axis on
        sub0's side only)."""
        from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=10,
            nj=12,
            nk=4,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        with pytest.raises(ParameterError, match="crosses a split axis"):
            split_blocks(m.blocks, [BlockSplit("sub0", "k", [2])])

    def test_paired_perpendicular_split_works(self, sphere_stl):
        """Splitting all 4 sub-blocks along k at the same index produces
        a coherent 8-chunk equatorial ring with all SHARED j-seams
        rewired."""
        from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=10,
            nj=12,
            nk=4,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        new = split_blocks(
            m.blocks,
            [BlockSplit(f"sub{q}", "k", [2]) for q in range(4)],
        )
        # 8 sub-block chunks + 2 caps = 10 blocks.
        assert len(new) == 10
        labels = [b.label for b in new]
        for q in range(4):
            assert f"sub{q}__k0" in labels
            assert f"sub{q}__k1" in labels
        assert "cap_north" in labels
        assert "cap_south" in labels

        # Every sub-block chunk should have its j_lo SHARED with the
        # previous sub's j_hi at the same k chunk.
        chunks = {b.label: b for b in new}
        for q in range(4):
            for k in (0, 1):
                blk = chunks[f"sub{q}__k{k}"]
                qprev = (q - 1) % 4
                expected_partner = chunks[f"sub{qprev}__k{k}"]
                assert blk.face("j_lo").bc == FaceBC.SHARED
                assert blk.face("j_lo").partner is expected_partner.face("j_hi")


# ---------------------------------------------------------------------------
# Adjacency JSON sidecar
# ---------------------------------------------------------------------------


class TestAdjacencyJSON:
    def test_round_trip_box_split(self, tmp_path):
        new = split_blocks(_box(), [BlockSplit("box", "i", [6])])
        out = write_adjacency_json(new, tmp_path / "adj")
        assert out.suffix == ".json"
        assert out.exists()

        payload = read_adjacency_json(out)
        assert payload["version"] == ADJACENCY_VERSION
        assert len(payload["blocks"]) == 2
        labels = {b["label"] for b in payload["blocks"]}
        assert labels == {"box__i0", "box__i1"}

        # Verify the SHARED partner reference points at the right block / face.
        i0 = next(b for b in payload["blocks"] if b["label"] == "box__i0")
        assert i0["faces"]["i_hi"]["bc"] == "SHARED"
        assert i0["faces"]["i_hi"]["partner"] == {
            "label": "box__i1",
            "face": "i_lo",
        }

    def test_extension_added_when_missing(self, tmp_path):
        new = split_blocks(_box(), [BlockSplit("box", "i", [6])])
        out = write_adjacency_json(new, tmp_path / "no_ext")
        assert out.name == "no_ext.json"

    def test_read_rejects_wrong_version(self, tmp_path):
        from shore.exceptions import MeshInputError

        p = tmp_path / "bad.json"
        p.write_text(json.dumps({"version": 99, "blocks": []}))
        with pytest.raises(MeshInputError, match="unsupported adjacency version"):
            read_adjacency_json(p)


# ---------------------------------------------------------------------------
# CLI: shore split
# ---------------------------------------------------------------------------


class TestSplitCLI:
    def _make_input(self, tmp_path: Path) -> Path:
        """Write a single 12x8x6 box .geo to disk and return its path."""
        from shore.io.geo import write_geo

        box = BoxMesh.build(vmin=(0, 0, 0), vmax=(10, 5, 3), ni=12, nj=8, nk=6)
        p = tmp_path / "background_box.geo"
        write_geo(p, box.blocks[0].nodes)
        return p

    def _make_config(self, tmp_path: Path, label: str = "box") -> Path:
        cfg = tmp_path / "split.toml"
        cfg.write_text(
            textwrap.dedent(
                f"""
                version = 1
                kind    = "split"

                [[splits]]
                label = "{label}"
                axis  = "i"
                at    = [6]
                """
            )
        )
        return cfg

    def test_basic_invocation(self, tmp_path):
        geo = self._make_input(tmp_path)
        cfg = self._make_config(tmp_path)
        out_stem = tmp_path / "out"

        r = runner.invoke(
            app,
            ["split", str(geo), "-c", str(cfg), "-o", str(out_stem)],
        )
        assert r.exit_code == 0, r.output
        assert (tmp_path / "out_box__i0.geo").exists()
        assert (tmp_path / "out_box__i1.geo").exists()
        assert (tmp_path / "out.adjacency.json").exists()

        # Verify node counts on the chunks.
        a = read_geo(tmp_path / "out_box__i0.geo")
        b = read_geo(tmp_path / "out_box__i1.geo")
        assert a.shape == (7, 8, 6, 3)
        assert b.shape == (6, 8, 6, 3)

    def test_missing_input_exits_1(self, tmp_path):
        cfg = self._make_config(tmp_path)
        r = runner.invoke(
            app,
            ["split", str(tmp_path / "ghost.geo"), "-c", str(cfg)],
        )
        assert r.exit_code == 1
        assert "file not found" in r.output

    def test_unknown_label_in_config_exits_1(self, tmp_path):
        geo = self._make_input(tmp_path)
        cfg = self._make_config(tmp_path, label="not_in_input")
        r = runner.invoke(
            app,
            ["split", str(geo), "-c", str(cfg), "-o", str(tmp_path / "out")],
        )
        assert r.exit_code == 1
        assert "no block named" in r.output
