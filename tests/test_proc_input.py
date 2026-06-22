"""Tests for shore.io.proc_input — write_proc_input + merge_proc_input + CLI."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.cli import app
from shore.exceptions import ParameterError
from shore.io.geo import write_geo
from shore.io.grd import read_grd_metadata, write_grd
from shore.io.proc_input import (
    BalanceSummary,
    merge_proc_input,
    write_proc_input,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_grd(tmp: Path, shapes: list[tuple[int, int, int]], stem: str = "test") -> Path:
    """Build a .grd from a list of (ni_node, nj_node, nk_node) shapes."""
    geo_paths = []
    for i, (ni, nj, nk) in enumerate(shapes):
        path = tmp / f"{stem}_b{i}.geo"
        write_geo(path, np.zeros((ni, nj, nk, 3)))
        geo_paths.append(path)
    return write_grd(tmp / f"{stem}.grd", geo_paths)


def _expected_cells(shape_node: tuple[int, int, int]) -> int:
    """Cell count = product of (node - 1)."""
    ni, nj, nk = shape_node
    return (ni - 1) * (nj - 1) * (nk - 1)


# ---------------------------------------------------------------------------
# read_grd_metadata
# ---------------------------------------------------------------------------


class TestReadGrdMetadata:
    def test_returns_cell_counts_per_block(self, tmp_path):
        shapes = [(8, 6, 4), (12, 12, 8), (4, 4, 4)]
        grd = _make_grd(tmp_path, shapes)
        meta = read_grd_metadata(grd)
        assert meta == [(7, 5, 3), (11, 11, 7), (3, 3, 3)]

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="not found"):
            read_grd_metadata(tmp_path / "ghost.grd")

    def test_truncated_file_raises(self, tmp_path):
        bad = tmp_path / "bad.grd"
        bad.write_bytes(b"\x00\x00\x00")
        with pytest.raises(ValueError, match="too short"):
            read_grd_metadata(bad)


# ---------------------------------------------------------------------------
# write_proc_input — defaults / serial
# ---------------------------------------------------------------------------


class TestWriteSerial:
    def test_serial_default_assigns_all_to_rank_0(self, tmp_path):
        grd = _make_grd(tmp_path, [(8, 6, 4), (12, 12, 8), (4, 4, 4)])
        out_path, summary = write_proc_input(grd, tmp_path / "p.input", np=1)
        assert out_path.exists()
        assert summary.n_blocks == 3
        assert summary.np == 1
        assert summary.per_rank_blocks == [3]
        assert summary.imbalance == 0.0

    def test_serial_writes_valid_format(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4), (4, 4, 4)])
        write_proc_input(grd, tmp_path / "p.input", np=1)
        text = (tmp_path / "p.input").read_text().splitlines()
        assert len(text) == 7 + 2  # 7 header + 2 records
        # Line 4 carries the block count.
        assert text[3].split()[0] == "2"
        # Records: block_idx group body proc.
        for i, ln in enumerate(text[7:]):
            head, _, _comment = ln.partition("!")
            toks = head.split()
            assert int(toks[0]) == i + 1
            assert int(toks[1]) == 0  # group default
            assert int(toks[2]) == 0  # body default
            assert int(toks[3]) == 0  # proc — serial


# ---------------------------------------------------------------------------
# write_proc_input — greedy weight-balanced assignment
# ---------------------------------------------------------------------------


class TestGreedyAssignment:
    def test_uniform_blocks_round_robin_balance(self, tmp_path):
        # 4 identical blocks across 2 ranks: 2 blocks per rank, perfect balance.
        grd = _make_grd(tmp_path, [(5, 5, 5)] * 4)
        _, summary = write_proc_input(grd, tmp_path / "p.input", np=2)
        assert summary.per_rank_blocks == [2, 2]
        assert summary.per_rank_cells[0] == summary.per_rank_cells[1]
        assert summary.imbalance == 0.0

    def test_dominant_block_lands_alone(self, tmp_path):
        # One huge + two small.  Greedy puts the huge alone on rank 0.
        grd = _make_grd(tmp_path, [(20, 20, 20), (4, 4, 4), (4, 4, 4)])
        _, summary = write_proc_input(grd, tmp_path / "p.input", np=2)
        assert summary.per_rank_blocks == [1, 2]
        # The 19^3 = 6859 cells dwarf 2 * 27 = 54 cells.
        assert summary.per_rank_cells[0] == 19**3
        assert summary.per_rank_cells[1] == 2 * 3**3

    def test_three_blocks_two_ranks(self, tmp_path):
        # Reproduces the smoke-test case: cell weights 105, 847, 27.
        grd = _make_grd(tmp_path, [(8, 6, 4), (12, 12, 8), (4, 4, 4)])
        _, summary = write_proc_input(grd, tmp_path / "p.input", np=2)
        # Heaviest (block 2) goes to one rank; the other two to the other.
        # No move/swap can improve the 847-vs-132 span, so the local
        # search leaves the greedy seed untouched.
        sorted_loads = sorted(summary.per_rank_cells)
        assert sorted_loads == [27 + 105, 847]


# ---------------------------------------------------------------------------
# write_proc_input — weighting metric (interior cells, not nodes)
# ---------------------------------------------------------------------------


class TestCellWeighting:
    def test_per_rank_cells_are_interior_cells(self, tmp_path):
        # read_grd_metadata returns node-1 per axis, so the per-block
        # weight ni*nj*nk is the interior-cell count, not the node count.
        node_shapes = [(8, 6, 4), (12, 12, 8), (4, 4, 4)]
        grd = _make_grd(tmp_path, node_shapes)
        _, summary = write_proc_input(grd, tmp_path / "p.input", np=1)
        # Serial: all cells land on rank 0; total must equal sum of
        # interior cells (node-1 product), NOT node products.
        total_cells = sum(_expected_cells(s) for s in node_shapes)
        total_nodes = sum(ni * nj * nk for ni, nj, nk in node_shapes)
        assert summary.per_rank_cells[0] == total_cells
        assert summary.per_rank_cells[0] != total_nodes

    def test_thin_slabs_balanced_by_cells_not_nodes(self, tmp_path):
        # Thin slabs vs a cube: node-count and cell-count weightings give
        # a different node/cell ratio per block, so balancing the wrong
        # metric mis-distributes work.  We assert the CELL loads balance.
        #
        #   slab 160x81x4 (nodes) -> 159*80*3   = 38160 cells
        #   slab 160x81x4 (nodes) -> 159*80*3   = 38160 cells
        #   cube  40x40x40 (nodes) -> 39*39*39  = 59319 cells
        #   cube  40x40x40 (nodes) -> 39*39*39  = 59319 cells
        node_shapes = [(160, 81, 4), (160, 81, 4), (40, 40, 40), (40, 40, 40)]
        grd = _make_grd(tmp_path, node_shapes)
        _, summary = write_proc_input(grd, tmp_path / "p.input", np=2)

        cells = [_expected_cells(s) for s in node_shapes]
        # The summary cell loads must be the real interior-cell loads.
        assert sum(summary.per_rank_cells) == sum(cells)
        # Balanced on cells: one cube + one slab per rank -> 59319+38160.
        assert sorted(summary.per_rank_cells) == [59319 + 38160, 59319 + 38160]
        # Sanity: the per-block node/cell ratio is in the 1.07-1.36 band
        # the issue cites, i.e. node-weighting would NOT be equivalent.
        ratios = [(ni * nj * nk) / c for (ni, nj, nk), c in zip(node_shapes, cells, strict=True)]
        assert any(r > 1.3 for r in ratios)  # the thin slabs
        assert any(r < 1.1 for r in ratios)  # the cubes


# ---------------------------------------------------------------------------
# write_proc_input — local-search refinement
# ---------------------------------------------------------------------------


class TestLocalSearchRefinement:
    def _greedy_only_span(self, weights: list[int], np: int) -> int:
        """Peak-to-trough span of a pure greedy-LPT pass (no refinement)."""
        rank = [0] * np
        order = sorted(((weights[i], i) for i in range(len(weights))), key=lambda x: (-x[0], x[1]))
        for w, _i in order:
            r = min(range(np), key=lambda r: (rank[r], r))
            rank[r] += w
        return max(rank) - min(rank)

    def test_refinement_never_worse_than_greedy(self, tmp_path):
        # A spread of block sizes; refined span must be <= greedy-only span.
        node_shapes = [(9, 9, 9), (8, 8, 8), (7, 7, 7), (6, 6, 6), (5, 5, 5)]
        grd = _make_grd(tmp_path, node_shapes)
        _, summary = write_proc_input(grd, tmp_path / "p.input", np=2)
        cells = [_expected_cells(s) for s in node_shapes]
        refined_span = max(summary.per_rank_cells) - min(summary.per_rank_cells)
        greedy_span = self._greedy_only_span(cells, 2)
        assert refined_span <= greedy_span

    def test_refinement_strictly_beats_greedy_on_crafted_case(self):
        # weights [8,7,6,5,4] over 2 ranks: greedy LPT -> [17,13] (span 4);
        # a single swap reaches [15,15] (span 0).  Drive the private path
        # directly so the case is exact and independent of .grd plumbing.
        from shore.io.proc_input import _build_summary, _resolve_proc_assignment

        weights = [8, 7, 6, 5, 4]
        greedy_span = self._greedy_only_span(weights, 2)
        assert greedy_span == 4  # guard the premise

        procs = _resolve_proc_assignment(weights, 2, None, len(weights))
        summary = _build_summary(weights, procs, 2)
        refined_span = max(summary.per_rank_cells) - min(summary.per_rank_cells)
        assert refined_span < greedy_span
        assert refined_span == 0  # perfect 15/15 split

    def test_refinement_is_deterministic(self):
        from shore.io.proc_input import _resolve_proc_assignment

        weights = [8, 7, 6, 5, 4, 9, 3, 11]
        runs = [_resolve_proc_assignment(weights, 3, None, len(weights)) for _ in range(5)]
        assert all(r == runs[0] for r in runs)

    def test_pinned_blocks_never_moved_by_refinement(self):
        from shore.io.proc_input import _resolve_proc_assignment

        # Pin two blocks to rank 0; the refinement must leave them there
        # even though moving them off would balance better.
        weights = [100, 1, 1, 1]
        procs = _resolve_proc_assignment(weights, 2, {1: 0, 2: 0}, len(weights))
        assert procs[0] == 0  # pinned
        assert procs[1] == 0  # pinned

    def test_pinned_block_seeds_load_but_movables_refine_around_it(self):
        from shore.io.proc_input import _resolve_proc_assignment

        # Pin the single heavy block to rank 0; the light movable blocks
        # should pile onto rank 1 to balance — and the heavy one stays put.
        weights = [50, 10, 10, 10, 10, 10]
        procs = _resolve_proc_assignment(weights, 2, {1: 0}, len(weights))
        assert procs[0] == 0  # pinned heavy block unmoved
        loads = [0, 0]
        for w, p in zip(weights, procs, strict=True):
            loads[p] += w
        # rank 0 = 50 + maybe one 10; rank 1 = the rest.  Span small.
        assert max(loads) - min(loads) <= 10


# ---------------------------------------------------------------------------
# write_proc_input — explicit overrides
# ---------------------------------------------------------------------------


class TestExplicitOverrides:
    def test_explicit_proc_assignment(self, tmp_path):
        grd = _make_grd(tmp_path, [(5, 5, 5)] * 3)
        _, summary = write_proc_input(
            grd,
            tmp_path / "p.input",
            np=2,
            proc_assignment={1: 0, 2: 1, 3: 0},  # explicit triple
        )
        assert summary.per_rank_blocks == [2, 1]

    def test_partial_proc_assignment_greedy_fills_rest(self, tmp_path):
        grd = _make_grd(tmp_path, [(10, 10, 10), (5, 5, 5), (5, 5, 5)])
        # Pin block 1 to rank 1; greedy assigns the other two.
        _, summary = write_proc_input(
            grd,
            tmp_path / "p.input",
            np=2,
            proc_assignment={1: 1},
        )
        # Block 1 weight = 9^3 = 729 on rank 1; remaining 4^3 + 4^3 = 128 split.
        # Greedy puts both small blocks on rank 0 (lighter).
        assert summary.per_rank_blocks == [2, 1]
        assert summary.per_rank_cells[1] == 9**3

    def test_group_and_body_overrides(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)] * 2)
        write_proc_input(
            grd,
            tmp_path / "p.input",
            np=1,
            group_assignment={1: 5, 2: 5},
            body_assignment={1: 1, 2: 1},
        )
        text = (tmp_path / "p.input").read_text().splitlines()
        for ln in text[7:]:
            head, _, _ = ln.partition("!")
            toks = head.split()
            assert int(toks[1]) == 5  # group
            assert int(toks[2]) == 1  # body


# ---------------------------------------------------------------------------
# write_proc_input — error cases
# ---------------------------------------------------------------------------


class TestWriteErrors:
    def test_np_zero_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        with pytest.raises(ParameterError, match="np must be >= 1"):
            write_proc_input(grd, tmp_path / "p.input", np=0)

    def test_proc_out_of_range(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        with pytest.raises(ParameterError, match="outside range"):
            write_proc_input(
                grd,
                tmp_path / "p.input",
                np=2,
                proc_assignment={1: 5},
            )

    def test_unknown_block_index(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        with pytest.raises(ParameterError, match="references block"):
            write_proc_input(
                grd,
                tmp_path / "p.input",
                np=1,
                proc_assignment={5: 0},
            )


# ---------------------------------------------------------------------------
# merge_proc_input
# ---------------------------------------------------------------------------


class TestMergeProcInput:
    def test_concatenation_shifts_block_indices(self, tmp_path):
        # First component: 2 blocks on rank 0.
        grd_a = _make_grd(tmp_path, [(4, 4, 4)] * 2, stem="a")
        write_proc_input(grd_a, tmp_path / "a.proc.input", np=1)
        # Second component: 3 blocks on rank 1.
        grd_b = _make_grd(tmp_path, [(4, 4, 4)] * 3, stem="b")
        write_proc_input(
            grd_b,
            tmp_path / "b.proc.input",
            np=2,
            proc_assignment={1: 1, 2: 1, 3: 1},
        )

        out_path, summary = merge_proc_input(
            [tmp_path / "a.proc.input", tmp_path / "b.proc.input"],
            tmp_path / "merged.proc.input",
        )
        assert summary.n_blocks == 5

        text = out_path.read_text().splitlines()
        body = text[7:]
        assert len(body) == 5
        # Block indices: 1..5 (renumbered).
        assert [int(line.split()[0]) for line in body] == [1, 2, 3, 4, 5]
        # Procs preserved: 0,0,1,1,1.
        assert [int(line.split()[3]) for line in body] == [0, 0, 1, 1, 1]

    def test_one_input_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        write_proc_input(grd, tmp_path / "a.proc.input", np=1)
        with pytest.raises(ValueError, match="at least 2"):
            merge_proc_input([tmp_path / "a.proc.input"], tmp_path / "out.proc.input")

    def test_missing_input_rejected(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="not found"):
            merge_proc_input(
                [tmp_path / "ghost.proc.input", tmp_path / "ghost2.proc.input"],
                tmp_path / "out.proc.input",
            )

    def test_merge_with_grd_path_recomputes_real_cell_balance(self, tmp_path):
        # Build two components with very different cell counts so the
        # difference between block-count-only and cell-count summaries
        # is visible.
        grd_a = _make_grd(tmp_path, [(20, 20, 20)] * 2, stem="big")  # 19^3 cells each
        grd_b = _make_grd(tmp_path, [(4, 4, 4)] * 4, stem="small")  # 3^3 cells each
        write_proc_input(grd_a, tmp_path / "big.proc.input", np=2)
        write_proc_input(grd_b, tmp_path / "small.proc.input", np=2)

        from shore.io.grd import merge_grd

        merged_grd = merge_grd(
            [tmp_path / "big.grd", tmp_path / "small.grd"],
            tmp_path / "assembly.grd",
        )

        # Without grd_path: cells_known=False, no real imbalance number.
        _, sum_no_grd = merge_proc_input(
            [tmp_path / "big.proc.input", tmp_path / "small.proc.input"],
            tmp_path / "merged_a.proc.input",
        )
        assert sum_no_grd.cells_known is False
        assert sum_no_grd.imbalance == 0.0

        # With grd_path: real cell-count summary.
        _, sum_with_grd = merge_proc_input(
            [tmp_path / "big.proc.input", tmp_path / "small.proc.input"],
            tmp_path / "merged_b.proc.input",
            grd_path=merged_grd,
        )
        assert sum_with_grd.cells_known is True
        # Total cells = 2 * 19^3 + 4 * 3^3 = 13718 + 108 = 13826.
        assert sum(sum_with_grd.per_rank_cells) == 2 * 19**3 + 4 * 3**3
        # The cell-count summary differs from the block-count summary.
        # No-grd: per_rank_cells is all-zero.
        # With grd: per_rank_cells reflects real weights.
        assert sum(sum_no_grd.per_rank_cells) == 0
        assert max(sum_with_grd.per_rank_cells) > 0

    def test_merge_with_mismatched_grd_rejected(self, tmp_path):
        grd_a = _make_grd(tmp_path, [(4, 4, 4)] * 2, stem="a")
        grd_b = _make_grd(tmp_path, [(4, 4, 4)] * 2, stem="b")
        write_proc_input(grd_a, tmp_path / "a.proc.input", np=1)
        write_proc_input(grd_b, tmp_path / "b.proc.input", np=1)
        # Use grd_a as a stand-in for the "merged" .grd — wrong block
        # count (2 instead of 4) on purpose.
        with pytest.raises(ParameterError, match="does not match"):
            merge_proc_input(
                [tmp_path / "a.proc.input", tmp_path / "b.proc.input"],
                tmp_path / "merged.proc.input",
                grd_path=grd_a,
            )

    def test_comments_carry_input_stem(self, tmp_path):
        grd_a = _make_grd(tmp_path, [(4, 4, 4)], stem="wall")
        grd_b = _make_grd(tmp_path, [(4, 4, 4)], stem="bg")
        write_proc_input(grd_a, tmp_path / "wall.proc.input", np=1)
        write_proc_input(grd_b, tmp_path / "bg.proc.input", np=1)
        out, _ = merge_proc_input(
            [tmp_path / "wall.proc.input", tmp_path / "bg.proc.input"],
            tmp_path / "merged.proc.input",
        )
        text = out.read_text()
        assert "wall" in text and "bg" in text


# ---------------------------------------------------------------------------
# BalanceSummary.format
# ---------------------------------------------------------------------------


class TestBalanceSummary:
    def test_format_serial(self):
        s = BalanceSummary(
            n_blocks=3,
            np=1,
            per_rank_blocks=[3],
            per_rank_cells=[1000],
            imbalance=0.0,
        )
        out = s.format()
        assert "3 blocks across 1 ranks" in out
        assert "imbalance" not in out  # serial -> no imbalance line

    def test_format_parallel_includes_imbalance(self):
        s = BalanceSummary(
            n_blocks=4,
            np=2,
            per_rank_blocks=[2, 2],
            per_rank_cells=[1000, 800],
            imbalance=0.222,
        )
        out = s.format()
        assert "imbalance: 22.2%" in out

    def test_format_uses_M_suffix_for_million(self):
        s = BalanceSummary(
            n_blocks=1,
            np=1,
            per_rank_blocks=[1],
            per_rank_cells=[1_843_000],
            imbalance=0.0,
        )
        assert "1.84M cells" in s.format()

    def test_build_summary_cells_and_imbalance_hand_computed(self):
        # Hand-computed: blocks of cell-weights distributed so the loads
        # are deliberately uneven, then check per_rank_cells and the
        # (max-min)/mean imbalance match by hand.
        from shore.io.proc_input import _build_summary

        weights = [1000, 600, 400]
        procs = [0, 0, 1]  # rank0 = 1600, rank1 = 400
        s = _build_summary(weights, procs, 2)
        assert s.per_rank_cells == [1600, 400]
        assert s.per_rank_blocks == [2, 1]
        # mean = 2000/2 = 1000; imbalance = (1600-400)/1000 = 1.2.
        assert s.imbalance == pytest.approx(1.2)
        assert s.cells_known is True


# ---------------------------------------------------------------------------
# CLI: shore proc-input
# ---------------------------------------------------------------------------


class TestProcInputCLI:
    def test_basic_invocation(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)] * 2)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "proc-input",
                str(grd),
                str(tmp_path / "p.input"),
            ],
        )
        assert r.exit_code == 0, r.output
        assert (tmp_path / "p.input").exists()

    def test_np_2_prints_balance_summary(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)] * 4)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "proc-input",
                str(grd),
                str(tmp_path / "p.input"),
                "--np",
                "2",
            ],
        )
        assert r.exit_code == 0
        # Summary goes to stderr; CliRunner mixes stderr/stdout into r.output.
        # Whitespace can wrap so split-and-rejoin to be robust.
        out = " ".join(r.output.lower().split())
        assert "ranks" in out
        assert "imbalance" in out

    def test_np_zero_exits_1(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "proc-input",
                str(grd),
                str(tmp_path / "p.input"),
                "--np",
                "0",
            ],
        )
        assert r.exit_code == 1

    def test_missing_grd_exits_1(self, tmp_path):
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "proc-input",
                str(tmp_path / "ghost.grd"),
                str(tmp_path / "p.input"),
            ],
        )
        assert r.exit_code == 1
        assert "not found" in r.output.lower()

    def test_meta_inline_json(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)] * 2)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "proc-input",
                str(grd),
                str(tmp_path / "p.input"),
                "--np",
                "2",
                "--meta",
                json.dumps({"1": {"proc": 1, "body": 5}}),
            ],
        )
        assert r.exit_code == 0, r.output
        text = (tmp_path / "p.input").read_text().splitlines()
        # Block 1 should be on proc 1, body 5.
        toks = text[7].split("!")[0].split()
        assert int(toks[2]) == 5
        assert int(toks[3]) == 1

    def test_meta_invalid_block_key_exits_1(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "proc-input",
                str(grd),
                str(tmp_path / "p.input"),
                "--meta",
                json.dumps({"not_an_int": {"proc": 0}}),
            ],
        )
        assert r.exit_code == 1


# ---------------------------------------------------------------------------
# CLI: shore proc-input-merge
# ---------------------------------------------------------------------------


class TestProcInputMergeCLI:
    def test_basic_merge(self, tmp_path):
        grd_a = _make_grd(tmp_path, [(4, 4, 4)] * 2, stem="a")
        grd_b = _make_grd(tmp_path, [(4, 4, 4)] * 3, stem="b")
        write_proc_input(grd_a, tmp_path / "a.proc.input", np=1)
        write_proc_input(grd_b, tmp_path / "b.proc.input", np=1)

        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "proc-input-merge",
                str(tmp_path / "a.proc.input"),
                str(tmp_path / "b.proc.input"),
                "-o",
                str(tmp_path / "merged.proc.input"),
            ],
        )
        assert r.exit_code == 0, r.output

        text = (tmp_path / "merged.proc.input").read_text().splitlines()
        # 5 blocks total, indices 1..5.
        body = text[7:]
        assert [int(ln.split()[0]) for ln in body] == [1, 2, 3, 4, 5]

    def test_one_input_exits_1(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        write_proc_input(grd, tmp_path / "a.proc.input", np=1)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "proc-input-merge",
                str(tmp_path / "a.proc.input"),
                "-o",
                str(tmp_path / "out.proc.input"),
            ],
        )
        assert r.exit_code == 1
