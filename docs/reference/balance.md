# `shore balance`

Plan a **multi-way equal-partition** split for MPI load balancing. Reads a
`.grd` (and optionally an `.adjacency.json`), partitions every block into
roughly equal, roughly cubic chunks sized toward `np * blocks-per-rank`, and
emits either a `split.toml` for [`shore split`](./split) (pre-overset) or a
`bsplit.par` for the external `bsplit` tool (post-overset).

```bash
shore balance MESH.grd -o splits.toml --np N \
              [--blocks-per-rank 1.0] [--ngr 1] [--min-coarse-cells 8] \
              [--tolerance 0.05] [--adjacency MESH.adjacency.json] \
              [--format toml|bsplit-par]
```

## Why this exists

[`shore proc-input`](./proc-input)'s assignment is optimal **for the block list
it is given**. When the imbalance still exceeds the target, the only way to do
better is to subdivide blocks. `shore balance` plans that subdivision.

The planner replaces an earlier bisection-of-heaviest design that either stalled
(heavy blocks left whole) or over-fragmented (single-cell sheets) — see
[Design decisions → multi-way partition](/guide/design-decisions#multi-way-equal-partition-not-bisection-issue-6)
for the motivation and evidence ([issue #6](https://github.com/szaghi/shore/issues/6)).

The result is a *plan*, not an executed split: review it, run `shore split` (or
`bsplit`) to apply it, then `shore proc-input` on the larger block list.

## Arguments and options

| Argument / Option | Default | Description |
|---|---|---|
| `GRD` | *required* | Path to the input `.grd` (block ordering and cell counts). |
| `-o` / `--output` | *required* | Output path (`split.toml` or `bsplit.par`). |
| `--np INTEGER` | *required* | MPI rank count the plan targets (>= 1). |
| `--blocks-per-rank FLOAT` | `1.0` | Target average blocks per rank. Higher → finer granularity and better balance, at more blocks. |
| `--ngr INTEGER` | `1` | Multigrid levels the mesh must support. Every chunk axis is kept divisible by `2^(ngr-1)`. |
| `--min-coarse-cells INTEGER` | `8` | Minimum cells per axis on the coarsest (level-`ngr`) grid. Every chunk axis is kept `>= min-coarse-cells * 2^(ngr-1)` fine cells. |
| `--tolerance FLOAT` | `0.05` | Target fractional imbalance `(max - min) / mean` for the *converged* report. The planner reports achieved imbalance regardless and never fails on an unreachable tolerance. |
| `--adjacency PATH` | — | Path to a `.adjacency.json` sidecar. **Strongly recommended**: lets the planner floor-check against seam-*propagated* cuts (`shore split` mirrors cuts across SHARED / DIRICHLET seams). Without it the floor is only checked per block in isolation. |
| `--format {toml,bsplit-par}` | `toml` | Output dialect. `toml` for the pre-overset [`shore split`](./split); `bsplit-par` for the post-overset external `bsplit` (see [After-overset splitting](/guide/after-overset-bsplit)). `bsplit-par` requires a **grd-only** plan (do not pass `--adjacency`). |

## Algorithm

1. Read `.grd` cell counts. Target chunk volume = `total_cells / round(np * blocks-per-rank)`;
   target chunk **edge** = `target_volume^(1/3)` (keeps chunks roughly cubic).
2. For each block, choose per-axis chunk counts as a **pure function of that
   axis's extent** (not the 3-D volume), `round(extent / target_edge)`, capped so
   every chunk keeps `>= floor` cells where `floor = min-coarse-cells * 2^(ngr-1)`.
   Keying the count to the axis extent makes seam-linked blocks (which always have
   equal extent on the shared axis) choose **identical** cuts — conformal by
   construction.
   See [Design decisions → per-axis cut counts](/guide/design-decisions#per-axis-cut-counts-not-per-block-volume).
3. Propagate the candidate cuts across SHARED / DIRICHLET seams (exactly as
   `shore split` will), then floor-check the **propagated** chunk set and back off
   any cut that would create a sub-floor sliver. Iterate to a fixpoint so the plan
   matches the split output block-for-block.
4. Group cuts on the same `(label, axis)` into one `[[splits]]` TOML entry (or one
   line per cut in `bsplit.par`).

The planner does **not** emit paired co-splits — [`shore split`](./split)
propagates cuts across seams automatically. The planner only chooses per-block
cuts and verifies the floor survives propagation.

## Best-effort: the floor wins over the tolerance

The multigrid floor is a hard constraint; the tolerance is a target. When the
floor prevents reaching the tolerance — e.g. an input block already smaller than
the floor, so no split can help — the planner **warns, does its best, and
reports the achieved imbalance** rather than failing:

```
shore balance: 33 blocks (imbalance 235.3%) → 192 blocks (imbalance 1.94%) via 83 cut(s) across 29 block(s)
Note: achieved imbalance 1.94% > tolerance 5.0% ...   # only when not converged
```

This honesty is deliberate: some geometries cannot meet a requested tolerance
within the floor, and you need the real number.

## Workflow (pre-overset)

```bash
# 1. plan the splits (TOML)
shore balance wall.grd --np 16 --ngr 3 --min-coarse-cells 8 \
  --blocks-per-rank 8 --adjacency wall.adjacency.json -o split.toml

# 2. apply (cuts auto-propagate across seams)
shore split -c split.toml wall_*.geo --adjacency wall.adjacency.json

# 3. repack the larger block list
shore grd wall_split_*.geo -o wall.grd

# 4. final rank assignment
shore proc-input wall.grd wall.proc.input --np 16
```

## Workflow (post-overset, `bsplit-par`)

When pre-overset splitting breaks the chimera donor search, emit a `bsplit.par`
instead and split the overset outputs after the solver's donor search — see
[After-overset splitting](/guide/after-overset-bsplit):

```bash
shore balance raw.grd --np 16 --ngr 3 --min-coarse-cells 8 \
  --blocks-per-rank 8 --format bsplit-par -o bsplit.par
```

## Output schema (TOML)

The TOML follows [`shore split`](./split)'s `kind = "split"` schema:

```toml
version = 1
kind    = "split"

# Generated by shore balance: target tolerance 5.0%, final imbalance 1.9% after 83 iteration(s).

[[splits]]
label = "background_111"
axis  = "i"
at    = [60, 120]
# ... one [[splits]] per (label, axis)
```

## Output schema (`bsplit.par`)

Seven header lines then one cut per line — `<block_number> <direction 1|2|3>
<fine_index>`:

```
sol
cc
cc
3
0
0
0
20 1 60
20 1 120
20 2 60
...
```

The block number is the 1-based `.grd` block order (the grd-only plan's synthetic
`block_<n>` labels carry it). bsplit rescales the fine-level cut index per
multigrid level internally.

## Examples

```bash
# 16 ranks, 3 multigrid levels, ~8 blocks/rank, seam-aware floor check
shore balance wall.grd --np 16 --ngr 3 --min-coarse-cells 8 \
  --blocks-per-rank 8 --adjacency wall.adjacency.json -o split.toml

# finer granularity (better balance, more blocks)
shore balance wall.grd --np 16 --blocks-per-rank 16 --adjacency wall.adjacency.json -o split.toml

# post-overset: emit a bsplit.par (grd-only, no --adjacency)
shore balance raw.grd --np 16 --ngr 3 --format bsplit-par -o bsplit.par
```

## Python API

See [`shore.balance`](./api-balance) for `plan_balance(...)`, the `BalancePlan`
dataclass, and the `to_toml()` / `to_bsplit_par()` serializers.

## See also

- [`shore split`](./split) — applies a `split.toml` plan (auto-propagates seam cuts).
- [`shore proc-input`](./proc-input) — final cell-weighted rank assignment.
- [After-overset splitting](/guide/after-overset-bsplit) — the `bsplit` workflow.
- [Design decisions](/guide/design-decisions) — why the planner works this way.
