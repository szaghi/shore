---
title: shore.balance
outline: [2, 3]
---

# `shore.balance`

Multi-way equal-partition split planner for Xall MPI load balancing.

Plans the splits that subdivide a multi-block assembly so its per-rank
load can be balanced — without mutating files — and emits a `split.toml`
that :func:`shore.split.split_blocks` executes.

Algorithm (replaces the old bisection planner; see GitHub issue #6)
-------------------------------------------------------------------
1. Read `.grd` cell counts.  Target chunk volume = total_cells /
   round(np \* blocks_per_rank).
2. For each block, choose per-axis chunk counts (:func:`_chunk_counts`)
   that drive the chunk volume toward the target while keeping every
   chunk roughly cubic and every axis `>= floor` cells, where
   `floor = min_coarse_cells * 2**(ngr-1)` (the multigrid coarse-cell
   floor) and cuts are snapped to multiples of `2**(ngr-1)`.
3. Propagate the candidate cuts across SHARED/DIRICHLET seams (exactly
   as :func:`shore.split.split_blocks` will), then floor-check the
   *propagated* chunk set and back off any cut that — combined with a
   neighbour's propagated cut — would create a sub-floor sliver.  Iterate
   to a fixpoint so the plan matches the split output block-for-block.
4. Group cuts per `(label, axis)` into one `[[splits]]` TOML entry.

Why this design
---------------
The old planner bisected the current-heaviest block repeatedly; with no
target block count it either stalled (heavy blocks left whole) or shredded
blocks into single-cell sheets (issue #6).  Equal multi-way partition with
a fixed target count avoids both.  Seam conformance is no longer the
planner's job to pre-pair: :func:`shore.split.split_blocks` propagates cuts
across seams automatically, so the planner only chooses per-block cuts and
verifies the floor survives propagation.

Best-effort: when the floor cannot be met (e.g. an input block already
smaller than `floor`, so no split can help), the planner warns, does what
it can, and returns `converged=False` with the achieved imbalance.

## `BalancePlan` <Badge type="tip" text="dataclass" /> {#BalancePlan}

```python
BalancePlan(
    splits: list[tuple[str, str, int]],
    iterations: int,
    converged: bool,
    initial_summary: BalanceSummary,
    final_summary: BalanceSummary,
    tolerance: float,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/balance.py#L89">source</a>

A plan to subdivide a block list for load balancing.

**Attributes**

- **`splits`** — List of `(label, axis, at)` tuples in the order they were decided.  The TOML output groups same-(label, axis) entries.
- **`iterations`** — Number of split iterations actually performed.
- **`converged`** — `True` if the final imbalance is within `tolerance`.
- **`initial_summary`** — Balance summary before any split.
- **`final_summary`** — Balance summary after the last split.
- **`tolerance`** — The tolerance the planner targeted (fractional, e.g. 0.05 = 5%).

### `BalancePlan.to_toml` <Badge type="info" text="method" /> {#BalancePlan-to_toml}

```python
to_toml() -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/balance.py#L116">source</a>

Render the plan as a `shore split` config TOML.

### `BalancePlan.to_bsplit_par` <Badge type="info" text="method" /> {#BalancePlan-to_bsplit_par}

```python
to_bsplit_par(
    *,
    ngr: int = 1,
    sol_root: str = 'sol',
    grd_root: str = 'cc',
    icc_root: str = 'cc',
    varn: int = 0,
    varo: int = 0,
    debug: int = 0,
) -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/balance.py#L142">source</a>

Render the plan as a `bsplit.par` for the external `bsplit` tool.

`bsplit` splits the framed **overset outputs** (`cc.0N` / `cc.0N.grd`)
*after* an overset run, remapping the donor pool by index — the
post-overset analogue of pre-overset :func:`to_toml` + `shore split`.
Use this when splitting before overset is impossible (e.g. the chimera
donor search degenerates on a split overlap block).

The plan must come from a **grd-only** :func:`plan_balance` (no
`adjacency_path`), whose synthetic `block_<n>` labels carry the
1-based grd block order — which is exactly `bsplit`'s block number.
Cut indices are fine-level vertex indices; `bsplit` rescales them per
multigrid level internally.

Format (see `bsplit/src/readpar.f90`):

```
<sol_root>
<grd_root>
<icc_root>
<ngr>
<varn>
<varo>
<debug>
<block_number> <direction 1|2|3> <fine_vertex_index>   (one per cut)
```

**Raises**

- `ParameterError` — If any split label is not a synthetic `block_<n>` label (a labeled plan cannot be mapped to bsplit block numbers without a number map).

## `plan_balance` <Badge type="info" text="function" /> {#plan_balance}

```python
plan_balance(
    grd_path: str | Path,
    *,
    np: int,
    blocks_per_rank: float = 1.0,
    ngr: int = 1,
    min_coarse_cells: int = 8,
    tolerance: float = 0.05,
    adjacency_path: str | Path | None = None,
    blocks: list | None = None,
) -> BalancePlan
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/balance.py#L203">source</a>

Plan a multi-way equal-partition split for MPI load balancing.

Replaces the old bisection-of-heaviest planner (issue #6), which
stalled (caps left whole) or over-fragmented (single-cell sheets).
This planner instead partitions **every** block into roughly equal,
roughly cubic chunks sized so the total lands near
`np * blocks_per_rank` blocks, then lets a greedy cell-weighted
assignment (and :func:`shore.io.proc_input`) place them.  ``shore
split`` auto-propagates the cuts across SHARED/DIRICHLET seams, so the
planner does **not** emit paired co-splits — it just chooses per-block
per-axis cuts.

Multigrid floor: every resulting chunk axis is kept ``&gt;=
min_coarse_cells \* 2\*\*(ngr - 1)`` fine cells and divisible by
`2**(ngr - 1)`, so the block survives `ngr` multigrid levels with
at least `min_coarse_cells` cells on the coarsest level.  When this
floor prevents reaching `tolerance` (e.g. a block too small to cut
further but still heavy), the planner does its best and returns
`converged=False` with the achieved imbalance — it never fails.

**Parameters**

- **`grd_path`** — Input `.grd` (block ordering + cell shapes).
- **`np`** — MPI rank count (&gt;= 1).
- **`blocks_per_rank`** — Target average blocks per rank.  `1.0` aims for ~`np` blocks; higher values give finer granularity / better balance at the cost of more blocks (and more overset donor-search work).
- **`ngr`** — Multigrid level count the mesh must support (&gt;= 1).
- **`min_coarse_cells`** — Minimum cells per axis at the coarsest (level-`ngr`) grid (&gt;= 1).
- **`tolerance`** — Target fractional imbalance `(max - min) / mean` for the converged flag.  The planner reports achieved imbalance regardless.
- **`adjacency_path`** — Optional `.adjacency.json`; used only for block labels (so the plan's `[[splits]]` entries name real blocks).  Seam conformance is handled at split time by cut propagation, not here.

**Returns**

- `BalancePlan`

## `chunk_axis_sizes_local` <Badge type="info" text="function" /> {#chunk_axis_sizes_local}

```python
chunk_axis_sizes_local(n_vertices: int, at: set[int]) -> list[int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/balance.py#L465">source</a>

Cell counts of each chunk when an axis of *n_vertices* is cut at *at*.
