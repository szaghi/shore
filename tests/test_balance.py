"""Tests for shore.balance — the multi-way equal-partition split planner.

These tests target the REWRITTEN planner (GitHub issue #6): equal multi-way
partition toward ``np * blocks_per_rank`` blocks with a multigrid cell floor,
seam-propagation-aware floor checking via ``shore.split.propagate_cuts``, and
best-effort (never-raises) behaviour on an unreachable tolerance.  The old
bisection / ``max_iterations`` API is gone.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from shore.balance import plan_balance
from shore.cli import app
from shore.exceptions import ParameterError
from shore.io.geo import write_geo
from shore.io.grd import write_grd

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_grd(tmp: Path, shapes: list[tuple[int, int, int]], stem: str = "test") -> Path:
    geo_paths = []
    for i, (ni, nj, nk) in enumerate(shapes):
        path = tmp / f"{stem}_b{i}.geo"
        write_geo(path, np.zeros((ni, nj, nk, 3)))
        geo_paths.append(path)
    return write_grd(tmp / f"{stem}.grd", geo_paths)


def _write_adjacency_stub(
    path: Path,
    n_blocks: int,
    *,
    labels: list[str] | None = None,
    seams: dict[str, dict] | None = None,
) -> Path:
    """Write a minimal adjacency JSON for testing without a full Mesh.

    `seams` maps block label -> {face_name: {"bc": ..., "partner": {...}, "axis_map": ...}}.
    """
    labels = labels or [f"block_{i + 1}" for i in range(n_blocks)]
    blocks = []
    for label in labels:
        face_dict = (seams or {}).get(label, {})
        # Fill in any unmentioned faces as WALL/FREE.
        for fn in ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi"):
            face_dict.setdefault(
                fn, {"bc": "WALL" if fn == "k_lo" else "FREE", "partner": None, "axis_map": None}
            )
        blocks.append({"label": label, "shape": [10, 10, 10], "faces": face_dict})
    payload = {"version": 2, "blocks": blocks}
    path.write_text(json.dumps(payload))
    return path


# ---------------------------------------------------------------------------
# Conformal Cartesian fixture — a row of equal boxes wired SHARED on i.
#
# The .grd holds CELL counts (vertices - 1); the HexBlocks built for the
# ``blocks=`` argument hold VERTEX counts, so the two must agree:
# grd cells = nodes vertices - 1.  These helpers keep them in lockstep so a
# plan can be FED BACK through ``split_blocks`` and floor-checked for real.
# ---------------------------------------------------------------------------


def _cartesian_row_adjacency(
    path: Path, labels: list[str], shape_vertices: tuple[int, int, int]
) -> Path:
    """Adjacency JSON for a 1-D row of equal boxes, SHARED across i-seams.

    Block n's ``i_hi`` is SHARED with block n+1's ``i_lo`` (identity axis_map
    over the in-plane (j, k) axes).  All other faces are WALL/FREE.  Equal
    extents on every axis make the tiling conformal by construction, so the
    planner's per-axis cut counts are seam-aligned.
    """
    ni, nj, nk = shape_vertices
    seams: dict[str, dict] = {}
    n = len(labels)
    for idx, label in enumerate(labels):
        faces: dict[str, dict] = {}
        if idx + 1 < n:
            faces["i_hi"] = {
                "bc": "SHARED",
                "partner": {"label": labels[idx + 1], "face": "i_lo"},
                "axis_map": [["j", False], ["k", False]],
            }
        if idx - 1 >= 0:
            faces["i_lo"] = {
                "bc": "SHARED",
                "partner": {"label": labels[idx - 1], "face": "i_hi"},
                "axis_map": [["j", False], ["k", False]],
            }
        seams[label] = faces
    p = _write_adjacency_stub(path, n, labels=labels, seams=seams)
    # Rewrite the (stub-fixed) shapes to the real vertex counts.
    payload = json.loads(p.read_text())
    for blk in payload["blocks"]:
        blk["shape"] = [ni, nj, nk]
    p.write_text(json.dumps(payload))
    return p


def _wired_blocks(adjacency_path: Path) -> list:
    """Build zero-node wired HexBlocks from an adjacency sidecar.

    Mirrors ``shore.cli._blocks_from_adjacency_shapes``: only extents and the
    seam graph matter for cut propagation / floor checking and for actually
    running ``split_blocks`` (node coordinates are irrelevant to shape/wiring
    assertions).
    """
    from shore.io.adjacency import read_adjacency_json, wire_blocks_from_adjacency
    from shore.mesh.block import HexBlock

    payload = read_adjacency_json(adjacency_path)
    blocks = []
    for b in payload["blocks"]:
        ni, nj, nk = b["shape"]
        blocks.append(HexBlock(label=b["label"], nodes=np.zeros((ni, nj, nk, 3))))
    wire_blocks_from_adjacency(blocks, payload, source=str(adjacency_path))
    return blocks


def _split_directives_from_plan(plan) -> list:
    """Convert a BalancePlan's grouped splits into BlockSplit directives."""
    from shore.split import BlockSplit

    grouped: dict[tuple[str, str], set[int]] = {}
    for label, axis, at in plan.splits:
        grouped.setdefault((label, axis), set()).add(at)
    return [
        BlockSplit(label=label, axis=axis, at=sorted(cuts))
        for (label, axis), cuts in grouped.items()
    ]


# ---------------------------------------------------------------------------
# plan_balance — basic behaviour
# ---------------------------------------------------------------------------


class TestPlanBasic:
    def test_serial_returns_trivial_plan(self, tmp_path):
        # np=1: imbalance is 0 by construction and nothing needs splitting.
        grd = _make_grd(tmp_path, [(20, 20, 20), (5, 5, 5), (5, 5, 5)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=1, tolerance=0.05)
        assert plan.splits == []
        assert plan.final_summary.imbalance == 0.0
        assert plan.initial_summary.imbalance == 0.0

    def test_already_balanced_needs_no_cuts(self, tmp_path):
        # Four identical small blocks across 4 ranks at blocks_per_rank=1:
        # target ~4 blocks, already 4, no cuts needed.
        grd = _make_grd(tmp_path, [(9, 9, 9)] * 4)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=4, blocks_per_rank=1.0, tolerance=0.05)
        assert plan.splits == []
        assert plan.converged

    def test_partitions_toward_np_times_blocks_per_rank(self, tmp_path):
        # One large block; the planner cuts it into ~np*blocks_per_rank chunks.
        # n_blocks must grow monotonically with blocks_per_rank.
        grd = _make_grd(tmp_path, [(65, 65, 65)])

        def n_blocks_for(bpr: float) -> int:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                plan = plan_balance(grd, np=4, blocks_per_rank=bpr, tolerance=0.01)
            return plan.final_summary.n_blocks

        n1 = n_blocks_for(1.0)
        n4 = n_blocks_for(4.0)
        n16 = n_blocks_for(16.0)
        assert n1 < n4 < n16, f"n_blocks should grow with blocks_per_rank: {n1}, {n4}, {n16}"
        # And blocks_per_rank=1 over a single block must at least cut it.
        assert n1 > 1

    def test_cut_count_is_pure_function_of_extent(self, tmp_path):
        # Two blocks with the SAME i-extent but different j/k get IDENTICAL
        # i-cuts (the conformal-by-construction property): _axis_chunk_count
        # depends only on the axis's own extent, not the 3-D volume.
        grd = _make_grd(tmp_path, [(65, 65, 9), (65, 9, 65)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=4, blocks_per_rank=4.0, tolerance=0.01)
        i_cuts = {}
        for label, axis, at in plan.splits:
            if axis == "i":
                i_cuts.setdefault(label, []).append(at)
        labels = sorted(i_cuts)
        assert len(labels) == 2, f"expected i-cuts on both blocks, got {i_cuts}"
        assert sorted(i_cuts[labels[0]]) == sorted(i_cuts[labels[1]]), (
            f"equal-i-extent blocks must get identical i-cuts: {i_cuts}"
        )


# ---------------------------------------------------------------------------
# plan_balance — multigrid floor (the load-bearing end-to-end check)
# ---------------------------------------------------------------------------


class TestPlanFloor:
    def test_floor_respected_through_real_split(self, tmp_path):
        # Build a conformal row of 2 equal boxes wired SHARED on i, plan with
        # a non-trivial multigrid level count, then ACTUALLY run split_blocks
        # with the planned directives and assert the real output honours the
        # floor, divisibility, has zero dangling SHARED/DIRICHLET partners,
        # and matches the plan's predicted block count EXACTLY.
        ngr = 3
        min_coarse = 4
        mg = 2 ** (ngr - 1)  # 4
        floor = min_coarse * mg  # 16 cells/axis
        # 65 vertices -> 64 cells/axis: divisible by 4, comfortably above the
        # 16-cell floor when cut.  ``_make_grd`` writes VERTEX-shaped geos and
        # the .grd reports cells = vertices - 1, so the grd's cell shape lines
        # up with the wired HexBlocks' (65-vertex) shape automatically.
        shape_v = (65, 65, 65)
        grd = _make_grd(tmp_path, [shape_v, shape_v])

        adj = _cartesian_row_adjacency(tmp_path / "row.json", ["A", "B"], shape_v)
        blocks = _wired_blocks(adj)

        with warnings.catch_warnings():
            warnings.simplefilter("error")  # no floor warning expected here
            plan = plan_balance(
                grd,
                np=8,
                blocks_per_rank=2.0,
                ngr=ngr,
                min_coarse_cells=min_coarse,
                tolerance=0.05,
                adjacency_path=adj,
                blocks=blocks,
            )

        assert plan.splits, "a 64^3 x2 mesh on np=8 should require cuts"

        from shore.split import split_blocks

        out = split_blocks(blocks, _split_directives_from_plan(plan))

        # 1. Predicted block count == actual split output count.
        assert len(out) == plan.final_summary.n_blocks, (
            f"plan predicted {plan.final_summary.n_blocks} blocks, split made {len(out)}"
        )

        # 2. Every axis of every output block honours the floor + divisibility.
        for b in out:
            ni, nj, nk, _ = b.nodes.shape
            cell_dims = (ni - 1, nj - 1, nk - 1)
            for ax, c in zip(("i", "j", "k"), cell_dims, strict=True):
                assert c >= floor, f"{b.label} {ax}={c} below floor {floor}"
                assert c % mg == 0, f"{b.label} {ax}={c} not divisible by {mg}"

        # 3. No connection face is left dangling after the real split.
        for b in out:
            for f in b.faces():
                if f.bc.name in ("SHARED", "DIRICHLET"):
                    assert f.partner is not None, f"{b.label}.{f.name} dangling {f.bc.name}"

    def test_input_below_floor_warns_and_not_converged(self, tmp_path):
        # An input block already below the multigrid floor cannot be made
        # conformal at this level count: the planner WARNS (best-effort) and
        # reports converged=False — it never raises.
        # ngr=4 -> mg=8, floor = 8*8 = 64 cells/axis.  A 30-cell block is
        # below floor on every axis.
        grd = _make_grd(tmp_path, [(31, 31, 31)])  # 30 cells/axis < 64
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            plan = plan_balance(grd, np=4, ngr=4, min_coarse_cells=8, tolerance=0.05)
        assert any("below the multigrid floor" in str(x.message) for x in w), (
            f"expected a floor warning, got {[str(x.message) for x in w]}"
        )
        assert not plan.converged


# ---------------------------------------------------------------------------
# plan_balance — parameter validation
# ---------------------------------------------------------------------------


class TestPlanErrors:
    def test_np_below_one_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(9, 9, 9)])
        with pytest.raises(ParameterError, match="np must be >= 1"):
            plan_balance(grd, np=0)

    def test_zero_tolerance_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(9, 9, 9)])
        with pytest.raises(ParameterError, match="tolerance must be > 0"):
            plan_balance(grd, np=2, tolerance=0.0)

    def test_ngr_below_one_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(9, 9, 9)])
        with pytest.raises(ParameterError, match="ngr must be >= 1"):
            plan_balance(grd, np=2, ngr=0)

    def test_min_coarse_cells_below_one_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(9, 9, 9)])
        with pytest.raises(ParameterError, match="min_coarse_cells must be >= 1"):
            plan_balance(grd, np=2, min_coarse_cells=0)

    def test_blocks_per_rank_non_positive_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(9, 9, 9)])
        with pytest.raises(ParameterError, match="blocks_per_rank must be > 0"):
            plan_balance(grd, np=2, blocks_per_rank=0.0)


# ---------------------------------------------------------------------------
# BalancePlan output — to_toml()
# ---------------------------------------------------------------------------


class TestPlanOutput:
    def test_to_toml_groups_same_label_axis(self, tmp_path):
        grd = _make_grd(tmp_path, [(65, 65, 65)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=8, blocks_per_rank=4.0, tolerance=0.01)
        toml = plan.to_toml()
        import re

        entries = re.findall(r'label = "([^"]+)"\naxis  = "([^"]+)"', toml)
        seen = set()
        for label, axis in entries:
            assert (label, axis) not in seen, (
                "[[splits]] entries should be grouped by (label, axis)"
            )
            seen.add((label, axis))

    def test_toml_is_valid_split_config(self, tmp_path):
        # End-to-end: planner output must parse through shore split's loader,
        # carry version/kind, and yield at least one directive.
        grd = _make_grd(tmp_path, [(65, 65, 65)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=8, blocks_per_rank=4.0, tolerance=0.01)
        toml_text = plan.to_toml()
        assert "version = 1" in toml_text
        assert 'kind    = "split"' in toml_text

        from shore.config import load_split_config

        toml_path = tmp_path / "splits.toml"
        toml_path.write_text(toml_text)
        cfg = load_split_config(toml_path)
        assert len(cfg.splits) >= 1
        # Each (label, axis) appears as exactly one grouped directive.
        keys = [(d.label, d.axis) for d in cfg.splits]
        assert len(keys) == len(set(keys)), f"directives not grouped by (label, axis): {keys}"

    def test_to_bsplit_par_header_and_lines(self, tmp_path):
        grd = _make_grd(tmp_path, [(65, 65, 65)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=8, blocks_per_rank=4.0, tolerance=0.01)
        par = plan.to_bsplit_par(ngr=3, grd_root="cc", icc_root="cc").splitlines()
        # 7-line header: sol, grd, icc, ngr, varn, varo, debug
        assert par[:7] == ["sol", "cc", "cc", "3", "0", "0", "0"]
        # one body line per cut: "<block_n> <1|2|3> <index>"
        body = par[7:]
        assert len(body) == len(plan.splits)
        for line in body:
            n, d, at = line.split()
            assert int(n) >= 1 and d in ("1", "2", "3") and int(at) >= 1

    def test_to_bsplit_par_matches_toml_cuts(self, tmp_path):
        # bsplit-par is just another rendering of the SAME plan: identical cut set.
        grd = _make_grd(tmp_path, [(65, 65, 65)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=8, blocks_per_rank=4.0, tolerance=0.01)
        axis_num = {"i": 1, "j": 2, "k": 3}
        toml_cuts = {(int(lbl.split("_")[1]), axis_num[ax], at) for lbl, ax, at in plan.splits}
        par_cuts = set()
        for line in plan.to_bsplit_par(ngr=1).splitlines()[7:]:
            n, d, at = line.split()
            par_cuts.add((int(n), int(d), int(at)))
        assert par_cuts == toml_cuts

    def test_to_bsplit_par_rejects_labeled_plan(self):
        from shore.balance import BalancePlan
        from shore.io.proc_input import BalanceSummary

        s = BalanceSummary(n_blocks=1, np=1, per_rank_blocks=[1], per_rank_cells=[1], imbalance=0.0)
        plan = BalancePlan(
            splits=[("sub0", "k", 4)],  # a labeled (non block_<n>) plan
            iterations=1,
            converged=True,
            initial_summary=s,
            final_summary=s,
            tolerance=0.05,
        )
        with pytest.raises(ParameterError, match="block_<n>"):
            plan.to_bsplit_par()


# ---------------------------------------------------------------------------
# CLI: shore balance
# ---------------------------------------------------------------------------


class TestBalanceCLI:
    def test_basic_invocation(self, tmp_path):
        grd = _make_grd(tmp_path, [(65, 65, 65)])
        runner = CliRunner()
        out = tmp_path / "splits.toml"
        r = runner.invoke(
            app,
            ["balance", str(grd), "-o", str(out), "--np", "8", "--blocks-per-rank", "4"],
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        text = out.read_text()
        assert "version = 1" in text
        assert 'kind    = "split"' in text

    def test_format_bsplit_par(self, tmp_path):
        grd = _make_grd(tmp_path, [(65, 65, 65)])
        runner = CliRunner()
        out = tmp_path / "bsplit.par"
        r = runner.invoke(
            app,
            [
                "balance",
                str(grd),
                "-o",
                str(out),
                "--np",
                "8",
                "--ngr",
                "3",
                "--blocks-per-rank",
                "4",
                "--format",
                "bsplit-par",
            ],
        )
        assert r.exit_code == 0, r.output
        lines = out.read_text().splitlines()
        assert lines[:7] == ["sol", "cc", "cc", "3", "0", "0", "0"]
        assert all(len(line.split()) == 3 for line in lines[7:])

    def test_format_bsplit_par_rejects_adjacency(self, tmp_path):
        shape_v = (65, 65, 65)
        grd = _make_grd(tmp_path, [shape_v, shape_v])
        adj = _cartesian_row_adjacency(tmp_path / "row.json", ["test_b0", "test_b1"], shape_v)
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "balance",
                str(grd),
                "-o",
                str(tmp_path / "x.par"),
                "--np",
                "8",
                "--format",
                "bsplit-par",
                "--adjacency",
                str(adj),
            ],
        )
        assert r.exit_code == 1
        assert "grd-only" in r.output

    def test_format_invalid_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(65, 65, 65)])
        runner = CliRunner()
        r = runner.invoke(
            app, ["balance", str(grd), "-o", str(tmp_path / "x"), "--np", "8", "--format", "bogus"]
        )
        assert r.exit_code == 1

    def test_all_flags_accepted_with_adjacency(self, tmp_path):
        shape_v = (65, 65, 65)
        grd = _make_grd(tmp_path, [shape_v, shape_v])
        adj = _cartesian_row_adjacency(tmp_path / "row.json", ["test_b0", "test_b1"], shape_v)
        runner = CliRunner()
        out = tmp_path / "splits.toml"
        r = runner.invoke(
            app,
            [
                "balance",
                str(grd),
                "-o",
                str(out),
                "--np",
                "8",
                "--ngr",
                "3",
                "--min-coarse-cells",
                "4",
                "--blocks-per-rank",
                "2",
                "--adjacency",
                str(adj),
            ],
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        assert "[[splits]]" in out.read_text()

    def test_missing_grd_exits_1(self, tmp_path):
        runner = CliRunner()
        r = runner.invoke(
            app,
            ["balance", str(tmp_path / "ghost.grd"), "-o", str(tmp_path / "out.toml"), "--np", "4"],
        )
        assert r.exit_code == 1
        assert "not found" in r.output.lower()

    def test_missing_adjacency_exits_1(self, tmp_path):
        grd = _make_grd(tmp_path, [(9, 9, 9), (9, 9, 9)])
        runner = CliRunner()
        r = runner.invoke(
            app,
            [
                "balance",
                str(grd),
                "-o",
                str(tmp_path / "out.toml"),
                "--np",
                "4",
                "--adjacency",
                str(tmp_path / "ghost.json"),
            ],
        )
        assert r.exit_code == 1
        assert "not found" in r.output.lower()
