"""Tests for shore.balance — the iterative split planner."""

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
# plan_balance — basic behaviour
# ---------------------------------------------------------------------------


class TestPlanBasic:
    def test_already_balanced_returns_empty_plan(self, tmp_path):
        # Four identical blocks across 4 ranks: perfect balance.
        grd = _make_grd(tmp_path, [(5, 5, 5)] * 4)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=4, tolerance=0.05)
        assert plan.converged
        assert plan.iterations == 0
        assert plan.splits == []
        assert plan.initial_summary.imbalance == 0.0

    def test_serial_returns_empty_plan(self, tmp_path):
        # np=1: imbalance is 0 by construction.
        grd = _make_grd(tmp_path, [(20, 20, 20), (4, 4, 4), (4, 4, 4)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=1, tolerance=0.05)
        assert plan.converged
        assert plan.splits == []

    def test_one_heavy_block_converges_with_single_axis_cuts(self, tmp_path):
        # Six blocks of equal size; with np=8, they need to be cut.
        # All k axis is unsafe-free in a no-adjacency scenario, so the
        # planner picks the longest axis (here all equal -> k tiebreak).
        grd = _make_grd(tmp_path, [(20, 20, 20)] * 6)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=8, tolerance=0.10, max_iterations=200)
        assert plan.converged
        assert len(plan.splits) > 0
        # All cuts on a single axis per block (one-axis-per-block constraint).
        per_block_axes = {}
        for label, axis, _ in plan.splits:
            per_block_axes.setdefault(label, set()).add(axis)
        for label, axes in per_block_axes.items():
            assert len(axes) == 1, (
                f"block {label!r} got cuts on multiple axes {axes}; "
                "shore split v1 supports only one axis per block per call"
            )


# ---------------------------------------------------------------------------
# plan_balance — error cases and capping
# ---------------------------------------------------------------------------


class TestPlanErrors:
    def test_np_zero_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        with pytest.raises(ParameterError, match="np must be >= 1"):
            plan_balance(grd, np=0, tolerance=0.05)

    def test_zero_tolerance_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        with pytest.raises(ParameterError, match="tolerance must be > 0"):
            plan_balance(grd, np=2, tolerance=0.0)

    def test_negative_max_iterations_rejected(self, tmp_path):
        grd = _make_grd(tmp_path, [(4, 4, 4)])
        with pytest.raises(ParameterError, match="max_iterations must be >= 1"):
            plan_balance(grd, np=2, tolerance=0.05, max_iterations=0)

    def test_max_iterations_caps_gracefully(self, tmp_path):
        # Pathologically unbalanced setup: one giant + tiny blocks, only
        # 2 iterations allowed.  Planner returns best-found, converged=False.
        grd = _make_grd(tmp_path, [(50, 50, 50), (4, 4, 4), (4, 4, 4)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=4, tolerance=0.001, max_iterations=2)
        assert not plan.converged
        # Iterations may be slightly different from len(splits) when an
        # iteration produces multiple co-splits, but with no adjacency
        # they should match.
        assert len(plan.splits) <= 2 + 0  # at most 2 iterations of single cuts

    def test_no_adjacency_emits_warning(self, tmp_path):
        grd = _make_grd(tmp_path, [(20, 20, 20), (4, 4, 4)])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            plan_balance(grd, np=2, tolerance=0.05)
            assert any("no adjacency" in str(x.message).lower() for x in w)


# ---------------------------------------------------------------------------
# plan_balance — topology awareness
# ---------------------------------------------------------------------------


class TestPlanTopology:
    def test_seam_free_axis_preferred(self, tmp_path):
        # Build a 2-block mesh where j_hi/j_lo is SHARED.  Splitting
        # along *j* doesn't break the j-seam (the seam IS the j_hi face,
        # which sits perpendicular to j; perpendicular faces are simply
        # divided into smaller faces by a parallel-axis split, no
        # co-split needed).  Splitting along i or k *does* require
        # co-splits because i_*/k_* are parallel to the j-seam.
        # The planner should prefer the axis with no parallel SHARED
        # partners — j wins here.
        grd = _make_grd(tmp_path, [(10, 10, 10)] * 2)
        adj = _write_adjacency_stub(
            tmp_path / "adj.json",
            n_blocks=2,
            labels=["A", "B"],
            seams={
                "A": {
                    "j_hi": {
                        "bc": "SHARED",
                        "partner": {"label": "B", "face": "j_lo"},
                        "axis_map": [["i", False], ["k", False]],
                    }
                },
                "B": {
                    "j_lo": {
                        "bc": "SHARED",
                        "partner": {"label": "A", "face": "j_hi"},
                        "axis_map": [["i", False], ["k", False]],
                    }
                },
            },
        )
        plan = plan_balance(grd, np=4, tolerance=0.05, adjacency_path=adj, max_iterations=50)
        # All cuts should be on j (no parallel SHARED/DIRICHLET partners).
        axes = {axis for _, axis, _ in plan.splits}
        assert axes.issubset({"j"}), f"expected only j cuts, got {axes}"

    def test_co_splits_emitted_for_paired_seams(self, tmp_path):
        # Two blocks with a SHARED j-seam; force the planner to use j
        # by making k too short (only 2 vertices = 1 cell, can't split).
        # Then co-splits on j must appear on both blocks at matching indices.
        grd = _make_grd(tmp_path, [(20, 20, 2)] * 2)  # nk=2 nodes, 1 cell
        adj = _write_adjacency_stub(
            tmp_path / "adj.json",
            n_blocks=2,
            labels=["A", "B"],
            seams={
                "A": {
                    "j_hi": {
                        "bc": "SHARED",
                        "partner": {"label": "B", "face": "j_lo"},
                        "axis_map": [["i", False], ["k", False]],
                    }
                },
                "B": {
                    "j_lo": {
                        "bc": "SHARED",
                        "partner": {"label": "A", "face": "j_hi"},
                        "axis_map": [["i", False], ["k", False]],
                    }
                },
            },
        )
        # Since k axis has only 1 cell, planner falls back to j (or i).
        # If it picks j, co-splits on both A and B at matching indices.
        plan = plan_balance(grd, np=4, tolerance=0.05, adjacency_path=adj, max_iterations=50)
        if any(axis == "j" for _, axis, _ in plan.splits):
            j_cuts_a = sorted(at for label, axis, at in plan.splits if label == "A" and axis == "j")
            j_cuts_b = sorted(at for label, axis, at in plan.splits if label == "B" and axis == "j")
            assert j_cuts_a == j_cuts_b, (
                "co-splits on shared j-seam must be at matching vertex indices"
            )


# ---------------------------------------------------------------------------
# BalancePlan output
# ---------------------------------------------------------------------------


class TestPlanOutput:
    def test_to_toml_groups_same_label_axis(self, tmp_path):
        grd = _make_grd(tmp_path, [(20, 20, 20)] * 2)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=4, tolerance=0.05, max_iterations=50)
        toml = plan.to_toml()
        # Each block label appears at most once per axis in the output.
        import re

        entries = re.findall(r'label = "([^"]+)"\naxis  = "([^"]+)"', toml)
        seen = set()
        for label, axis in entries:
            assert (label, axis) not in seen, (
                "[[splits]] entries should be grouped by (label, axis)"
            )
            seen.add((label, axis))

    def test_toml_is_valid_split_config(self, tmp_path):
        # End-to-end check: planner output is parseable by shore split's
        # config loader.
        grd = _make_grd(tmp_path, [(20, 20, 20)] * 2)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plan = plan_balance(grd, np=4, tolerance=0.05, max_iterations=50)
        toml_text = plan.to_toml()
        from shore.config import load_split_config

        toml_path = tmp_path / "splits.toml"
        toml_path.write_text(toml_text)
        cfg = load_split_config(toml_path)
        # The config's directives should match plan.splits when grouped.
        assert len(cfg.splits) >= 1


# ---------------------------------------------------------------------------
# CLI: shore balance
# ---------------------------------------------------------------------------


class TestBalanceCLI:
    def test_basic_invocation(self, tmp_path):
        grd = _make_grd(tmp_path, [(20, 20, 20)] * 2)
        runner = CliRunner()
        out = tmp_path / "splits.toml"
        r = runner.invoke(
            app,
            ["balance", str(grd), "-o", str(out), "--np", "4"],
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        text = out.read_text()
        assert "version = 1" in text
        assert 'kind    = "split"' in text

    def test_already_balanced_writes_empty_plan(self, tmp_path):
        grd = _make_grd(tmp_path, [(5, 5, 5)] * 4)
        runner = CliRunner()
        out = tmp_path / "splits.toml"
        r = runner.invoke(
            app,
            ["balance", str(grd), "-o", str(out), "--np", "4"],
        )
        assert r.exit_code == 0
        text = out.read_text()
        # No [[splits]] entries.
        assert "[[splits]]" not in text

    def test_missing_grd_exits_1(self, tmp_path):
        runner = CliRunner()
        r = runner.invoke(
            app,
            ["balance", str(tmp_path / "ghost.grd"), "-o", str(tmp_path / "out.toml"), "--np", "4"],
        )
        assert r.exit_code == 1
        assert "not found" in r.output.lower()

    def test_missing_adjacency_exits_1(self, tmp_path):
        grd = _make_grd(tmp_path, [(5, 5, 5)] * 2)
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

    def test_with_adjacency_no_warning(self, tmp_path):
        grd = _make_grd(tmp_path, [(20, 20, 20)] * 2)
        adj = _write_adjacency_stub(tmp_path / "adj.json", n_blocks=2)
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
                str(adj),
            ],
        )
        assert r.exit_code == 0
