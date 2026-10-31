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

    def test_internal_seam_carries_identity_axis_map(self):
        """Splitting along i creates an identity (j↔j, k↔k, no flips)
        axis_map on each side of the new seam — the chunks share node
        memory, so there is no axis swap or reversal."""
        from shore.mesh.face import AxisMap

        new = split_blocks(_box(), [BlockSplit("box", "i", [6])])
        a, b = new
        assert a.face("i_hi").axis_map == AxisMap.identity("i_hi")
        assert b.face("i_lo").axis_map == AxisMap.identity("i_lo")

    def test_external_face_bcs_preserved_per_chunk(self):
        blocks = _box()
        blocks[0]._set_face_bc("k_lo", FaceBC.WALL)  # a box is all FREE (ADR 0055)
        new = split_blocks(blocks, [BlockSplit("box", "i", [6])])
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

    def test_multi_axis_same_block_one_call(self):
        """A block listed with two axes in one call is cut on both (a chunk
        grid), unioning per-axis cuts — the refactor supports this directly."""
        new = split_blocks(
            _box(ni=12, nj=8, nk=6),
            [BlockSplit("box", "i", [6]), BlockSplit("box", "j", [4])],
        )
        labels = sorted(b.label for b in new)
        assert labels == ["box__i0__j0", "box__i0__j1", "box__i1__j0", "box__i1__j1"]

    def test_cut_index_out_of_range(self):
        with pytest.raises(ParameterError, match="out of range"):
            split_blocks(_box(ni=12), [BlockSplit("box", "i", [12])])


# ---------------------------------------------------------------------------
# split_blocks: cross-block adjacency rewiring
# ---------------------------------------------------------------------------


class TestSplitAdjacencyAcrossBlocks:
    def test_lone_k_split_auto_propagates_across_seams(self, sphere_stl):
        """Splitting only sub0 on k now AUTO-PROPAGATES the k-cut to every
        block sharing a k-perpendicular seam with it (the j-periodic ring
        sub1/sub3 and the DIRICHLET caps), so the result is conformal — no
        manual paired directives, no rejection."""
        from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

        m = Mesh.from_stl(
            sphere_stl,
            topology=CubedSphereTopology(),
            ni=10,
            nj=12,
            nk=4,
            spacing_k=Spacing1D(law="geometric", ds=0.05, growth=1.1),
        )
        new = split_blocks(m.blocks, [BlockSplit("sub0", "k", [2])])
        labels = {b.label for b in new}
        # sub0's k-cut propagated to its j-ring partners and the caps:
        # every cubed-sphere block ends up k-split into __k0 / __k1.
        for base in ("sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"):
            assert f"{base}__k0" in labels, f"{base} not k-split"
            assert f"{base}__k1" in labels
        # and no connection face is left dangling.
        chunks = {b.label: b for b in new}
        for b in chunks.values():
            for f in b.faces():
                if f.bc.name in ("SHARED", "DIRICHLET"):
                    assert f.partner is not None

    def test_paired_perpendicular_split_works(self, sphere_stl):
        """Splitting all 6 cubed-sphere blocks (4 subs + 2 caps) along k at
        the same index produces a coherent ring with SHARED j-seams AND the
        DIRICHLET cap<->sub seams rewired onto the matching k-chunks.

        The caps MUST be co-split: the cap<->sub DIRICHLET seam is
        perpendicular to k, so a k-split of the sub without a matching k-split
        of the cap leaves that seam unmatched (split_blocks rejects it).
        """
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
            [
                BlockSplit(lbl, "k", [2])
                for lbl in ("sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south")
            ],
        )
        # 8 sub chunks + 4 cap chunks = 12 blocks.
        assert len(new) == 12
        labels = [b.label for b in new]
        for q in range(4):
            assert f"sub{q}__k0" in labels
            assert f"sub{q}__k1" in labels
        for cap in ("cap_north", "cap_south"):
            assert f"{cap}__k0" in labels
            assert f"{cap}__k1" in labels

        chunks = {b.label: b for b in new}
        # SHARED j-seams: each sub chunk's j_lo SHARED with the previous sub.
        for q in range(4):
            for k in (0, 1):
                blk = chunks[f"sub{q}__k{k}"]
                qprev = (q - 1) % 4
                expected_partner = chunks[f"sub{qprev}__k{k}"]
                assert blk.face("j_lo").bc == FaceBC.SHARED
                assert blk.face("j_lo").partner is expected_partner.face("j_hi")
        # DIRICHLET cap<->sub seams: each sub chunk keeps DIRICHLET i-faces
        # with a live partner on a matching cap chunk (regression: these used
        # to silently drop to partner=None through the split).
        for q in range(4):
            for k in (0, 1):
                blk = chunks[f"sub{q}__k{k}"]
                for fn in ("i_lo", "i_hi"):
                    assert blk.face(fn).bc == FaceBC.DIRICHLET
                    assert blk.face(fn).partner is not None


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

    def test_read_rejects_v1_with_regenerate_message(self, tmp_path):
        """v1 sidecars (no axis_map) are unambiguous data loss; reject with
        a clear regenerate-it message rather than silently mis-emitting
        cc.par on cross-axis seams."""
        from shore.exceptions import MeshInputError

        p = tmp_path / "v1.json"
        p.write_text(json.dumps({"version": 1, "blocks": []}))
        with pytest.raises(MeshInputError, match=r"version 1.*[Rr]egenerate"):
            read_adjacency_json(p)

    def test_axis_map_round_trip(self, tmp_path):
        """SHARED faces must serialise their axis_map and read it back
        bit-identical."""
        new = split_blocks(_box(), [BlockSplit("box", "i", [6])])
        out = write_adjacency_json(new, tmp_path / "adj")

        payload = read_adjacency_json(out)
        i0 = next(b for b in payload["blocks"] if b["label"] == "box__i0")
        # Internal seam (i_hi between consecutive chunks) is same-axis no-flip:
        # along (j, k) -> partner along (j, k), no flips.
        assert i0["faces"]["i_hi"]["axis_map"] == [["j", False], ["k", False]]
        # Non-connection faces have axis_map=None.
        assert i0["faces"]["k_lo"]["axis_map"] is None

    def test_read_rejects_malformed_axis_map(self, tmp_path):
        """A 2-element list with bad axis labels must fail at read time."""
        from shore.exceptions import MeshInputError

        bad = {
            "version": 2,
            "blocks": [
                {
                    "label": "x",
                    "shape": [3, 3, 3],
                    "faces": {
                        "i_lo": {"bc": "WALL", "partner": None, "axis_map": None},
                        "i_hi": {
                            "bc": "SHARED",
                            "partner": {"label": "y", "face": "i_lo"},
                            "axis_map": [["bad", False], ["k", False]],
                        },
                        "j_lo": {"bc": "FREE", "partner": None, "axis_map": None},
                        "j_hi": {"bc": "FREE", "partner": None, "axis_map": None},
                        "k_lo": {"bc": "WALL", "partner": None, "axis_map": None},
                        "k_hi": {"bc": "FREE", "partner": None, "axis_map": None},
                    },
                }
            ],
        }
        p = tmp_path / "bad.json"
        p.write_text(json.dumps(bad))
        with pytest.raises(MeshInputError, match="axis_map axis label"):
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

    def _make_two_block_geo(self, tmp_path: Path) -> tuple[Path, Path, Path]:
        """Two boxes sharing an i_hi/i_lo SHARED seam + their adjacency sidecar.

        Returns (geo_a, geo_b, adjacency_json).
        """
        from shore.io.adjacency import write_adjacency_json
        from shore.io.geo import write_geo
        from shore.mesh.block import HexBlock
        from shore.mesh.face import AxisMap

        a = BoxMesh.build(vmin=(0, 0, 0), vmax=(8, 8, 8), ni=9, nj=9, nk=9)
        b = BoxMesh.build(vmin=(8, 0, 0), vmax=(16, 8, 8), ni=9, nj=9, nk=9)
        ba = HexBlock(label="a", nodes=a.blocks[0].nodes)
        bb = HexBlock(label="b", nodes=b.blocks[0].nodes)
        ba._set_face_bc(
            "i_hi", FaceBC.SHARED, partner=bb.face("i_lo"), axis_map=AxisMap.identity("i_hi")
        )
        bb._set_face_bc(
            "i_lo", FaceBC.SHARED, partner=ba.face("i_hi"), axis_map=AxisMap.identity("i_lo")
        )
        ga = tmp_path / "bg_a.geo"
        gb = tmp_path / "bg_b.geo"
        write_geo(ga, ba.nodes)
        write_geo(gb, bb.nodes)
        adj = write_adjacency_json([ba, bb], tmp_path / "bg.adjacency.json")
        return ga, gb, adj

    def test_adjacency_preserves_seam_through_split(self, tmp_path):
        """--adjacency wires the input a<->b SHARED seam so it survives the
        split (without it, the .geo blocks are FREE and the seam is lost).

        Co-split both blocks on k (perpendicular to the a.i_hi<->b.i_lo seam)
        with matching indices -- split_blocks then re-pairs the matching
        k-chunks across the seam.
        """
        ga, gb, adj = self._make_two_block_geo(tmp_path)
        cfg = tmp_path / "split.toml"
        cfg.write_text(
            'version = 1\nkind = "split"\n\n'
            '[[splits]]\nlabel = "a"\naxis = "k"\nat = [4]\n'
            '[[splits]]\nlabel = "b"\naxis = "k"\nat = [4]\n'
        )
        out = tmp_path / "out"
        r = runner.invoke(
            app,
            ["split", str(ga), str(gb), "-c", str(cfg), "-o", str(out), "--adjacency", str(adj)],
        )
        assert r.exit_code == 0, r.output
        payload = read_adjacency_json(tmp_path / "out.adjacency.json")
        faces = {blk["label"]: blk["faces"] for blk in payload["blocks"]}
        # each b k-chunk's i_lo stays SHARED with the matching a k-chunk.
        b0 = faces["b__k0"]["i_lo"]
        assert b0["bc"] == "SHARED"
        assert b0["partner"]["label"] == "a__k0"

    def test_adjacency_omitted_loses_seam(self, tmp_path):
        """Control: without --adjacency the same co-split leaves b's i_lo FREE."""
        ga, gb, _ = self._make_two_block_geo(tmp_path)
        cfg = tmp_path / "split.toml"
        cfg.write_text(
            'version = 1\nkind = "split"\n\n'
            '[[splits]]\nlabel = "a"\naxis = "k"\nat = [4]\n'
            '[[splits]]\nlabel = "b"\naxis = "k"\nat = [4]\n'
        )
        out = tmp_path / "out"
        r = runner.invoke(app, ["split", str(ga), str(gb), "-c", str(cfg), "-o", str(out)])
        assert r.exit_code == 0, r.output
        payload = read_adjacency_json(tmp_path / "out.adjacency.json")
        faces = {blk["label"]: blk["faces"] for blk in payload["blocks"]}
        assert faces["b__k0"]["i_lo"]["bc"] == "FREE"

    def test_adjacency_suffix_label_match(self):
        """A sidecar label ('sub0') resolves to a full-stem block label
        ('sphere_o_sub0') via unique suffix match -- the combined multi-
        component split case where shore split keeps full stems."""
        import numpy as np

        from shore.io.adjacency import wire_blocks_from_adjacency
        from shore.mesh.block import HexBlock
        from shore.mesh.face import AxisMap, FaceBC

        a = HexBlock(label="sphere_o_sub0", nodes=np.zeros((5, 5, 5, 3)))
        b = HexBlock(label="sphere_o_sub1", nodes=np.zeros((5, 5, 5, 3)))

        def _shared(partner_label, partner_face, along):
            return {
                "bc": "SHARED",
                "partner": {"label": partner_label, "face": partner_face},
                "axis_map": [[along[0], False], [along[1], False]],
            }

        payload = {
            "version": ADJACENCY_VERSION,
            "blocks": [
                {  # sidecar uses topology labels, not full stems
                    "label": "sub0",
                    "shape": [5, 5, 5],
                    "faces": {"j_hi": _shared("sub1", "j_lo", ("i", "k"))},
                },
                {
                    "label": "sub1",
                    "shape": [5, 5, 5],
                    "faces": {"j_lo": _shared("sub0", "j_hi", ("i", "k"))},
                },
            ],
        }
        wire_blocks_from_adjacency([a, b], payload, source="t")
        assert a.face("j_hi").bc is FaceBC.SHARED
        assert a.face("j_hi").partner is b.face("j_lo")
        assert a.face("j_hi").axis_map == AxisMap.identity("j_hi")

    def test_adjacency_unmatched_label_errors(self):
        """A sidecar block matching no input block is a hard error (not a
        silent skip -- the regression that lost the cubed-sphere seams)."""
        import numpy as np

        from shore.exceptions import MeshInputError
        from shore.io.adjacency import wire_blocks_from_adjacency
        from shore.mesh.block import HexBlock

        blk = HexBlock(label="bg_000", nodes=np.zeros((5, 5, 5, 3)))
        payload = {
            "version": ADJACENCY_VERSION,
            "blocks": [{"label": "nonexistent", "shape": [5, 5, 5], "faces": {}}],
        }
        with pytest.raises(MeshInputError, match="matches no input block"):
            wire_blocks_from_adjacency([blk], payload, source="t")

    def test_multi_axis_split_one_block(self, tmp_path):
        """A block listed with two axes is cut on both (iterative rounds)."""
        from shore.io.geo import write_geo

        box = BoxMesh.build(vmin=(0, 0, 0), vmax=(8, 8, 8), ni=9, nj=9, nk=9)
        geo = tmp_path / "bg_box.geo"
        write_geo(geo, box.blocks[0].nodes)
        cfg = tmp_path / "split.toml"
        cfg.write_text(
            'version = 1\nkind = "split"\n\n'
            '[[splits]]\nlabel = "box"\naxis = "i"\nat = [4]\n'
            '[[splits]]\nlabel = "box"\naxis = "k"\nat = [4]\n'
        )
        out = tmp_path / "out"
        r = runner.invoke(app, ["split", str(geo), "-c", str(cfg), "-o", str(out)])
        assert r.exit_code == 0, r.output
        # i-cut then k-cut => 2x2 = 4 chunks, labelled __i{0,1}__k{0,1}.
        payload = read_adjacency_json(tmp_path / "out.adjacency.json")
        labels = sorted(blk["label"] for blk in payload["blocks"])
        assert labels == ["box__i0__k0", "box__i0__k1", "box__i1__k0", "box__i1__k1"]
