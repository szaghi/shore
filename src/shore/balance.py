"""Iterative split planner for Xall MPI load balancing.

When :func:`shore.io.proc_input.write_proc_input`'s greedy assignment
exceeds the user's tolerance, the typical fix is to subdivide the
heaviest block(s) so the per-rank load can be distributed more
evenly.  This module plans those splits — without mutating any
files — and emits a ``split.toml`` that :func:`shore.split.split_blocks`
can then execute.

Algorithm
---------
1. Read ``.grd`` metadata to get the per-block ``(ni, nj, nk)`` cell
   counts.  Initial weight = ``ni * nj * nk``.
2. Greedy weight-balanced assignment across ``np`` ranks.
3. While imbalance > tolerance and iter < max_iterations:
     a. Pick the heaviest block ``b``.
     b. Pick a topology-safe split axis: longest extent first, but
        only an axis whose perpendicular faces don't break SHARED /
        DIRICHLET seams without paired splits on the partner block.
     c. Simulate the split at the cell-count midpoint of axis A.
        Record the cut on the original block label (or its already
        split parent), at the original-block-relative vertex index.
     d. Recompute weights and re-balance.
4. Group cuts on the same (label, axis) into one ``[[splits]]`` TOML
   entry with a sorted ``at = [...]`` list.

Key choices
-----------
- **Midpoint cut by cell count.**  Block weight tracks cells, so the
  midpoint along the cut axis is the cell-count midpoint, not vertex
  midpoint.  For uniform spacing they're the same; for stretched
  spacing the planner doesn't know the spacing law (it only has the
  ``.grd`` header).  Cell-count midpoint is the honest choice.
- **Co-splits across SHARED/DIRICHLET seams** are emitted as separate
  ``[[splits]]`` entries on each affected block.  ``shore split``
  validates compatibility.
- **k-axis is preferred** because for typical body-fitted topologies
  ``k_lo`` / ``k_hi`` are WALL / FREE (never SHARED), so k-cuts need
  no paired splits.  When it isn't safe (an exotic topology with
  k-coupled SHARED seams), the planner falls back to j or i with
  paired splits.

References
----------
Inspired by ``LoadBalance/src/Lib_Balance.F90`` and
``overset-exploded/src/lib/overset_exploded_library.F90``.  Both
those tools physically split block data while iterating; SHORE's
planner only schedules cuts and lets ``shore split`` execute them.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path

from shore.exceptions import ParameterError
from shore.io.grd import read_grd_metadata
from shore.io.proc_input import BalanceSummary

# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------


@dataclass
class _SimBlock:
    """One simulated block (leaf) in the planner's working state.

    Splits never overwrite the original block — they create chunks.
    Each chunk is one ``_SimBlock`` and tracks where it sits inside
    the original block via ``origin``.

    Attributes
    ----------
    label :
        Block label.  For the original input it's the adjacency-JSON
        label (or ``block_<idx>`` for grd-only planning).  For chunks,
        it's the original label — the planner emits cuts on the
        original label, so cuts stay anchored to one TOML entry per
        ``(label, axis)`` pair.
    shape :
        Cell counts ``[ni, nj, nk]`` for this chunk.
    parent :
        1-based index into the original ``.grd`` block list.  All
        chunks of an original block share the same ``parent``.
    origin :
        Cell-index offset of this chunk along each axis, relative to
        the **original** block.  ``origin["k"] = 12`` means the
        chunk starts at original-block vertex 12 on the k-axis.  An
        unsplit original block has ``origin = {"i": 0, "j": 0, "k": 0}``.
    """

    label: str
    shape: list[int]
    parent: int
    origin: dict[str, int]


@dataclass
class BalancePlan:
    """A plan to subdivide a block list for load balancing.

    Attributes
    ----------
    splits :
        List of ``(label, axis, at)`` tuples in the order they were
        decided.  The TOML output groups same-(label, axis) entries.
    iterations :
        Number of split iterations actually performed.
    converged :
        ``True`` if the final imbalance is within ``tolerance``.
    initial_summary :
        Balance summary before any split.
    final_summary :
        Balance summary after the last split.
    tolerance :
        The tolerance the planner targeted (fractional, e.g. 0.05 = 5%).
    """

    splits: list[tuple[str, str, int]]
    iterations: int
    converged: bool
    initial_summary: BalanceSummary
    final_summary: BalanceSummary
    tolerance: float

    def to_toml(self) -> str:
        """Render the plan as a ``shore split`` config TOML."""
        lines = [
            "version = 1",
            'kind    = "split"',
            "",
            f"# Generated by shore balance: target tolerance "
            f"{self.tolerance * 100:.1f}%, "
            f"final imbalance {self.final_summary.imbalance * 100:.1f}% "
            f"after {self.iterations} iteration(s).",
            "",
        ]
        # Group by (label, axis) and sort the cuts.
        grouped: dict[tuple[str, str], list[int]] = {}
        for label, axis, at in self.splits:
            grouped.setdefault((label, axis), []).append(at)
        for (label, axis), cuts in grouped.items():
            cuts_sorted = sorted(set(cuts))
            cuts_str = ", ".join(str(c) for c in cuts_sorted)
            lines.append("[[splits]]")
            lines.append(f'label = "{label}"')
            lines.append(f'axis  = "{axis}"')
            lines.append(f"at    = [{cuts_str}]")
            lines.append("")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Public planner
# ---------------------------------------------------------------------------


def plan_balance(
    grd_path: str | Path,
    *,
    np: int,
    tolerance: float = 0.05,
    max_iterations: int = 100,
    adjacency_path: str | Path | None = None,
) -> BalancePlan:
    """Plan the splits needed to bring a multi-block mesh's load imbalance
    below ``tolerance``.

    Reads the ``.grd``'s per-block cell counts and (optionally) the
    adjacency JSON's seam graph; iteratively schedules cuts on the
    heaviest block until imbalance ≤ tolerance or iter = max_iterations.

    Parameters
    ----------
    grd_path :
        Path to the input ``.grd``.  Block ordering and shapes are
        derived from this file.
    np :
        Total MPI rank count (>= 1).
    tolerance :
        Target fractional imbalance ``(max - min) / mean``.  Default
        ``0.05`` (5%).
    max_iterations :
        Maximum split iterations.  Each iteration adds one cut to one
        block (with paired co-splits when topology requires).  Default
        ``100``.
    adjacency_path :
        Optional path to a ``.adjacency.json`` sidecar.  When supplied,
        the planner refuses to cut along axes that would break
        SHARED / DIRICHLET seams without paired co-splits (and emits
        co-splits when an axis is chosen).  When ``None``, the planner
        uses longest-axis-first with no seam awareness — the resulting
        plan may then be rejected by ``shore split`` on a topology
        with k-coupled seams; a warning notes this.

    Returns
    -------
    BalancePlan
        ``plan.converged`` is True iff the final imbalance is within
        tolerance.  ``plan.splits`` lists the cuts in decision order.

    Raises
    ------
    ParameterError
        If ``np < 1``, ``tolerance <= 0``, or ``max_iterations < 1``.
    FileNotFoundError, ValueError
        Forwarded from :func:`shore.io.grd.read_grd_metadata` and the
        adjacency reader.
    """
    if np < 1:
        raise ParameterError(f"plan_balance: np must be >= 1, got {np}")
    if tolerance <= 0:
        raise ParameterError(f"plan_balance: tolerance must be > 0, got {tolerance}")
    if max_iterations < 1:
        raise ParameterError(f"plan_balance: max_iterations must be >= 1, got {max_iterations}")

    shapes = read_grd_metadata(grd_path)
    n_orig = len(shapes)

    if adjacency_path is None:
        warnings.warn(
            "plan_balance: no adjacency provided; falling back to longest-axis-first "
            "without seam awareness.  The resulting plan may be rejected by "
            "`shore split` on topologies with k-coupled SHARED seams.  Pass "
            "adjacency_path=<adj.json> for a topology-safe plan.",
            UserWarning,
            stacklevel=2,
        )
        seam_graph = {}
        # Default block labels — synthetic indices.
        labels = [f"block_{i + 1}" for i in range(n_orig)]
    else:
        seam_graph, labels = _load_adjacency(adjacency_path)
        if len(labels) != n_orig:
            raise ParameterError(
                f"plan_balance: .grd has {n_orig} blocks but adjacency JSON has "
                f"{len(labels)}; the two inputs disagree on block count."
            )

    # Initial simulated blocks, indexed by their position in the active list.
    sim: list[_SimBlock] = [
        _SimBlock(
            label=labels[i],
            shape=list(shapes[i]),  # cell counts
            parent=i + 1,
            origin={"i": 0, "j": 0, "k": 0},
        )
        for i in range(n_orig)
    ]

    initial_summary = _balance_summary(sim, np)
    splits: list[tuple[str, str, int]] = []
    # `shore split` v1 supports one axis per block per call.  Once we
    # commit to an axis for an original block, later iterations can
    # add more cuts on that same axis but cannot switch to another.
    locked_axis: dict[int, str] = {}  # parent (1-based block idx) -> axis

    if initial_summary.imbalance <= tolerance:
        # Already balanced — no plan needed.
        return BalancePlan(
            splits=splits,
            iterations=0,
            converged=True,
            initial_summary=initial_summary,
            final_summary=initial_summary,
            tolerance=tolerance,
        )

    # Map original block label → parent index for the seam-group walk.
    label_to_parent: dict[str, int] = {labels[i]: i + 1 for i in range(n_orig)}

    for iter_idx in range(1, max_iterations + 1):
        # Heaviest sim-leaf first.
        weights = [_weight(b) for b in sim]
        ordered = sorted(range(len(sim)), key=lambda i: -weights[i])
        progressed = False
        for b_idx in ordered:
            block = sim[b_idx]
            forced = locked_axis.get(block.parent)
            axis = _pick_split_axis(block, seam_graph, forced=forced)
            if axis is None:
                continue
            # Co-splits: walk the seam graph by *parent label* (not leaf
            # index).  All leaves of every parent in the group must be
            # cut at the same original-block vertex.
            group_parents = _seam_group_parents(block.label, axis, seam_graph, label_to_parent)
            if any(locked_axis.get(p) not in (None, axis) for p in group_parents):
                continue  # group conflicts with a previously-locked axis
            ax_idx = _axis_idx(axis)
            cell_mid_local = block.shape[ax_idx] // 2
            if cell_mid_local < 1:
                continue
            cut_vertex = block.origin[axis] + cell_mid_local

            # For every parent in the group, the cut at original-vertex
            # `cut_vertex` may already be present (idempotent), and may
            # exceed some leaves' ranges (skip those).  We require that
            # *some* leaf of each parent contains the cut.
            group_leaves = _leaves_containing_cut(sim, group_parents, axis, cut_vertex)
            covered_parents = {sim[idx].parent for idx in group_leaves}
            if not group_parents.issubset(covered_parents):
                # Not every parent has a leaf to receive the cut.  Try
                # a different leaf as the heaviest candidate.
                continue
            if not group_leaves:
                continue

            # Apply: for each leaf, split at its local midpoint along `axis`
            # such that the original-vertex matches `cut_vertex`.
            cut_applied_for_parent: set[int] = set()
            for leaf_idx in sorted(group_leaves, reverse=True):
                # reverse so later sim.append doesn't disturb earlier indices
                leaf = sim[leaf_idx]
                local_cut = cut_vertex - leaf.origin[axis]
                if not (0 < local_cut < leaf.shape[ax_idx]):
                    continue
                _apply_split(sim, leaf_idx, axis, local_cut)
                if leaf.parent not in cut_applied_for_parent:
                    splits.append((leaf.label, axis, cut_vertex))
                    cut_applied_for_parent.add(leaf.parent)
                progressed = True
            for p in group_parents:
                locked_axis[p] = axis
            if progressed:
                break
        if not progressed:
            break
        new_summary = _balance_summary(sim, np)
        if new_summary.imbalance <= tolerance:
            return BalancePlan(
                splits=splits,
                iterations=iter_idx,
                converged=True,
                initial_summary=initial_summary,
                final_summary=new_summary,
                tolerance=tolerance,
            )

    # Loop exited without converging — return the best-found plan.
    final_summary = _balance_summary(sim, np)
    return BalancePlan(
        splits=splits,
        iterations=len(splits),  # number of single-block cuts decided
        converged=False,
        initial_summary=initial_summary,
        final_summary=final_summary,
        tolerance=tolerance,
    )


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------


def _axis_idx(axis: str) -> int:
    return {"i": 0, "j": 1, "k": 2}[axis]


def _weight(b: _SimBlock) -> int:
    return b.shape[0] * b.shape[1] * b.shape[2]


def _greedy_assign(weights: list[int], np: int) -> list[int]:
    """Greedy descending weight-balanced assignment across `np` ranks."""
    n = len(weights)
    procs = [0] * n
    rank_load = [0] * np
    order = sorted(range(n), key=lambda i: -weights[i])
    for i in order:
        r = min(range(np), key=lambda r: rank_load[r])
        procs[i] = r
        rank_load[r] += weights[i]
    return procs


def _balance_summary(sim: list[_SimBlock], np: int) -> BalanceSummary:
    weights = [_weight(b) for b in sim]
    procs = _greedy_assign(weights, np)
    per_rank_blocks = [0] * np
    per_rank_cells = [0] * np
    for w, p in zip(weights, procs, strict=True):
        per_rank_blocks[p] += 1
        per_rank_cells[p] += w
    if np > 1 and per_rank_cells:
        mean = sum(per_rank_cells) / np
        imbalance = (max(per_rank_cells) - min(per_rank_cells)) / mean if mean > 0 else 0.0
    else:
        imbalance = 0.0
    return BalanceSummary(
        n_blocks=len(sim),
        np=np,
        per_rank_blocks=per_rank_blocks,
        per_rank_cells=per_rank_cells,
        imbalance=imbalance,
    )


def _pick_split_axis(
    block: _SimBlock,
    seam_graph: dict,
    *,
    forced: str | None = None,
) -> str | None:
    """Return the longest-extent topology-safe axis to split ``block``,
    or None if none is available.

    ``forced``: when set (because the block's parent is already locked
    to one axis by an earlier split), only that axis is considered.

    "Safe" means: the axis perpendicular face is not PERIODIC (the
    j-ring on cubed-sphere sub-blocks falls there) and has at least 2
    cells.  SHARED / DIRICHLET partners trigger co-splits, not
    rejection.  Without an adjacency graph the planner trusts the
    user — every axis with extent >= 2 is safe.
    """
    if forced is not None:
        if block.shape[_axis_idx(forced)] < 2:
            return None
        if not _axis_safe(block, forced, seam_graph):
            return None
        return forced
    candidates = []
    for axis in ("k", "j", "i"):
        if block.shape[_axis_idx(axis)] < 2:
            continue
        if not _axis_safe(block, axis, seam_graph):
            continue
        candidates.append(axis)
    if not candidates:
        return None

    # Two-level sort:
    #   1. Axes with no parallel SHARED / DIRICHLET partners win (they
    #      can be split without dragging a co-split group along —
    #      a single cut on one block, not N cuts on N blocks).
    #   2. Among ties, longest extent.
    def _seam_cost(axis: str) -> int:
        return _axis_parallel_partner_count(block, axis, seam_graph)

    candidates.sort(key=lambda a: (_seam_cost(a), -block.shape[_axis_idx(a)]))
    return candidates[0]


def _axis_parallel_partner_count(block: _SimBlock, axis: str, seam_graph: dict) -> int:
    """Number of SHARED / DIRICHLET partners on faces *parallel* to axis.

    A face is "parallel to axis" iff its perpendicular axis is not the
    split axis — equivalently, the face contains ``axis`` as one of
    its two along-axes.  Splitting the block along ``axis`` cuts
    these faces into pieces along ``axis``, so each partner block
    must also be split at the matching index for the seam to remain
    valid.
    """
    if not seam_graph:
        return 0
    info = seam_graph.get(block.label, {})
    return sum(
        1
        for face_name, face_info in info.items()
        if not _face_perpendicular_to(face_name, axis)
        and face_info.get("bc") in ("SHARED", "DIRICHLET")
    )


def _axis_safe(block: _SimBlock, axis: str, seam_graph: dict) -> bool:
    """Check whether splitting `block` along `axis` is topology-safe.

    Without an adjacency graph the planner trusts the user (every
    axis is safe).  With one, an axis is unsafe when:

    - any face *parallel* to the split axis (i.e. containing it as
      an along-axis) is PERIODIC.  Splitting that face would break
      the periodic loop on the partner block, which has no matching
      "first" / "last" cell.
    """
    if not seam_graph:
        return True
    info = seam_graph.get(block.label, {})
    for face_name, face_info in info.items():
        if _face_perpendicular_to(face_name, axis):
            continue
        if face_info.get("bc") == "PERIODIC":
            return False
    return True


def _face_perpendicular_to(face_name: str, axis: str) -> bool:
    """True iff face_name (e.g. 'i_lo') is perpendicular to axis ('i')."""
    return face_name.startswith(axis + "_")


def _seam_group_parents(
    start_label: str,
    axis: str,
    seam_graph: dict,
    label_to_parent: dict[str, int],
) -> set[int]:
    """Set of parent block indices that must be co-split with start_label
    along ``axis``.

    Walk SHARED / DIRICHLET partners on faces *parallel* to ``axis``
    (i.e. faces whose along-axes include ``axis``).  These are the
    seams that get cut along ``axis`` when the block is split, and
    each partner needs to mirror the cut.

    The walk is by *original block label*, so it returns parent
    indices.  An empty seam_graph yields {start parent} only.
    """
    if not seam_graph:
        return {label_to_parent[start_label]}

    visited_labels = {start_label}
    queue = [start_label]
    while queue:
        cur = queue.pop()
        info = seam_graph.get(cur, {})
        for face_name, face_info in info.items():
            if _face_perpendicular_to(face_name, axis):
                continue
            if face_info.get("bc") not in ("SHARED", "DIRICHLET"):
                continue
            partner = face_info.get("partner") or {}
            p_label = partner.get("label")
            if not p_label or p_label in visited_labels:
                continue
            if p_label not in label_to_parent:
                continue
            visited_labels.add(p_label)
            queue.append(p_label)

    return {label_to_parent[label] for label in visited_labels}


def _leaves_containing_cut(
    sim: list[_SimBlock],
    parents: set[int],
    axis: str,
    cut_vertex: int,
) -> list[int]:
    """Return sim-list indices of leaves of ``parents`` whose extent along
    ``axis`` strictly contains ``cut_vertex``.

    A leaf with ``origin[axis] = O`` and ``shape[axis] = N`` covers
    original vertices ``[O, O + N]``.  The cut at ``cut_vertex`` is
    interior iff ``O < cut_vertex < O + N``.
    """
    out = []
    ax_idx = _axis_idx(axis)
    for i, leaf in enumerate(sim):
        if leaf.parent not in parents:
            continue
        lo = leaf.origin[axis]
        hi = lo + leaf.shape[ax_idx]
        if lo < cut_vertex < hi:
            out.append(i)
    return out


def _seam_group(
    b_idx: int,
    axis: str,
    sim: list[_SimBlock],
    seam_graph: dict,
) -> list[int]:
    """Return the set of sim-list indices that must be co-split with b_idx
    along `axis` so every SHARED / DIRICHLET seam stays intact.

    When splitting block X along ``axis``, every face *parallel* to
    ``axis`` (i.e. containing ``axis`` as an along-axis) gets cut
    into pieces along ``axis``.  If that face is SHARED / DIRICHLET
    with a partner block, the partner's matching face must also be
    split at the same parametric index — and the partner's *other*
    parallel faces propagate the requirement transitively.

    Without an adjacency graph the group is just {b_idx}.
    """
    if not seam_graph:
        return [b_idx]
    label_to_idx = {b.label: i for i, b in enumerate(sim)}

    visited = {b_idx}
    queue = [b_idx]
    while queue:
        cur = queue.pop()
        cur_label = sim[cur].label
        info = seam_graph.get(cur_label, {})
        for face_name, face_info in info.items():
            # Only faces *parallel* to the split axis are affected.
            if _face_perpendicular_to(face_name, axis):
                continue
            if face_info.get("bc") not in ("SHARED", "DIRICHLET"):
                continue
            partner = face_info.get("partner")
            if not partner:
                continue
            p_label = partner.get("label")
            p_idx = label_to_idx.get(p_label)
            if p_idx is None or p_idx in visited:
                continue
            visited.add(p_idx)
            queue.append(p_idx)
    return sorted(visited)


def _apply_split(
    sim: list[_SimBlock],
    b_idx: int,
    axis: str,
    cell_mid: int,
) -> None:
    """Split sim[b_idx] along ``axis`` at local cell index ``cell_mid``.

    The block in place becomes the lower chunk (its shape shrinks
    along ``axis``; its origin is unchanged).  A new ``_SimBlock``
    for the upper chunk is appended to ``sim`` with origin shifted
    by ``cell_mid`` along ``axis``.

    Labels: chunks keep the **original** label.  ``shore split``
    works from the ``[[splits]]`` TOML which lists cuts per original
    block, so the planner's job is to record cuts against the
    original label, not synthetic chunk labels.
    """
    parent = sim[b_idx]
    ax_idx = _axis_idx(axis)
    full = parent.shape[ax_idx]

    # Lower chunk: in place — origin unchanged, shape[axis] = cell_mid.
    lower_shape = list(parent.shape)
    lower_shape[ax_idx] = cell_mid
    parent.shape = lower_shape

    # Upper chunk: new leaf at origin + cell_mid, length (full - cell_mid).
    upper_shape = list(parent.shape)
    upper_shape[ax_idx] = full - cell_mid
    upper_origin = dict(parent.origin)
    upper_origin[axis] = parent.origin[axis] + cell_mid
    sim.append(
        _SimBlock(
            label=parent.label,
            shape=upper_shape,
            parent=parent.parent,
            origin=upper_origin,
        )
    )


def _load_adjacency(path: str | Path) -> tuple[dict, list[str]]:
    """Read an adjacency JSON; return (seam_graph, labels in file order)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"plan_balance: adjacency file not found: {p}")
    data = json.loads(p.read_text(encoding="ascii", errors="replace"))
    if "blocks" not in data:
        raise ParameterError(f"plan_balance: {p} missing 'blocks' key")
    seam_graph: dict[str, dict] = {}
    labels: list[str] = []
    for blk in data["blocks"]:
        label = blk["label"]
        labels.append(label)
        seam_graph[label] = blk.get("faces", {})
    return seam_graph, labels
