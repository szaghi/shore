# Split blocks and balance the MPI load

Xall gives each MPI rank whole blocks. A cubed sphere has 6 blocks of two
sizes, so on 10 ranks four ranks idle and the load is badly skewed. Split
the blocks into more, smaller chunks first, keeping enough cells per axis
for the multigrid levels you run.

The case is a cubed sphere with 33 × 64 × 33 nodes: four equatorial blocks
of 32 × 16 × 32 cells and two caps of 16 × 16 × 32, 81 920 cells in all.
The Python version is
[`examples/cookbook/split_balance.py`](https://github.com/szaghi/shore/blob/main/examples/cookbook/split_balance.py)
(about 2 s).

## Split blocks by hand

Write the cuts in a TOML file, one `[[splits]]` entry per block and axis,
with the cut positions as node indices:

```toml
version = 1
kind    = "split"

[[splits]]
label = "sub0"
axis  = "k"
at    = [16]
```

```bash
shore mesh sphere.stl --ni 33 --nj 64 --nk 33 --ds 1e-3 --theta-cap-deg 35 \
    -o sphere --adjacency-json sphere.adjacency.json
shore split -c manual.toml sphere_*.geo --adjacency sphere.adjacency.json -o manual
```

```text
✓ Wrote 12 blocks + manual.adjacency.json (1 split directive(s))
```

- One cut on `sub0` along k gave 12 blocks: the cut follows every seam
  that contains the k axis, round the equatorial ring and into both
  caps, so the chunks stay point-matched.
- Chunks are named `<label>__<axis><n>`, for example `sub0__k0`,
  `sub0__k1`; the output is `<stem>_<chunk>.geo` plus
  `<stem>.adjacency.json`.
- Always pass the `--adjacency` sidecar of the input. A `.geo` file holds
  only coordinates; without the sidecar every face is FREE, the seams are
  lost, and `cc.par` asks the overset pre-processor for donors across
  what was a 1-to-1 seam.
- In Python, `Project.split([(label, axis, [at, ...]), ...])` applies the
  same directives to a project's mesh.

## Plan a balanced split

`shore balance` reads the cell counts from a `.grd` and writes the split
TOML for you:

```bash
shore grd sphere_sub0.geo sphere_sub1.geo sphere_sub2.geo sphere_sub3.geo \
    sphere_cap_north.geo sphere_cap_south.geo -o sphere.grd
shore balance sphere.grd -o splits.toml --np 10 --ngr 3 --min-coarse-cells 4 \
    --adjacency sphere.adjacency.json
shore split -c splits.toml sphere_*.geo --adjacency sphere.adjacency.json -o split
shore grd split_*.geo -o split.grd
shore proc-input split.grd split.proc.input --np 10
```

```text
shore balance: 6 blocks (imbalance 200.0%) → 20 blocks (imbalance 0.00%) via 10 cut(s) across 6 block(s)
✓ Wrote 20 blocks + split.adjacency.json (10 split directive(s))
✓ Wrote split.grd (20 blocks)
...
  imbalance: 0.0%
```

<Figure name="load-balance" ext="svg" caption="Cells per rank before and after shore balance for a small case." />

- Imbalance is `(max - min) / mean` of the cells per rank. Before the
  split, 6 blocks on 10 ranks leave 4 ranks empty: 200 %.
- The planner cuts every block into roughly equal, roughly cubic chunks,
  aiming at `np × --blocks-per-rank` blocks (default 1.0), and chooses the
  same cut count on seam-linked axes so the result is conformal. The
  greedy cell-weighted assignment in `proc-input` then places them.
- `shore grd` packs blocks in the order you list them, and `cc.par` must
  use the same order. `shore split` writes the sidecar in sorted label
  order, which is the order the shell expands `split_*.geo` in.
- The planner never fails on the tolerance (`--tolerance`, default 5 %):
  it reports what it reached. Chunk granularity limits it. The same plan
  on 8 ranks gives 20 equal chunks on 8 ranks, 3 or 2 per rank, 40 %
  imbalance; choose `--np` and `--blocks-per-rank` so that the chunk
  count divides well.

The Python version builds the same sphere in a project, writes the `.grd`,
runs the same plan through `plan_balance` and applies it to the project:

<<< @/../examples/cookbook/split_balance.py#mesh{python}

<<< @/../examples/cookbook/split_balance.py#plan{python}

<<< @/../examples/cookbook/split_balance.py#apply{python}

```text
6 blocks, imbalance 200% -> 20 blocks, imbalance 0%
cuts: [('sub0', 'i', 16), ('sub0', 'k', 16), ('sub1', 'i', 16)] ...
20 blocks: sub0__i0__k0, sub0__i0__k1, sub0__i1__k0, sub0__i1__k1 ...
```

## Keep the multigrid floor

Xall coarsens each block `ngr - 1` times by 2. The planner keeps every
chunk axis at least `min_coarse_cells × 2^(ngr-1)` cells and divisible by
`2^(ngr-1)`:

| `--ngr` | `--min-coarse-cells` | Floor per axis | Effect here |
|---|---|---|---|
| 3 | 4 | 16 cells, multiple of 4 | 32-cell axes cut once, 16-cell axes not |
| 4 | 4 | 32 cells, multiple of 8 | no cut possible, warning |

```bash
shore balance sphere.grd -o t.toml --np 4 --ngr 4 --min-coarse-cells 4 \
    --adjacency sphere.adjacency.json
```

```text
UserWarning: plan_balance: 6 input block(s) already below the multigrid floor (>= 32 cells/axis,
divisible by 8) for ngr=4, min_coarse_cells=4 — they cannot be made conformal at this level count.
...
shore balance: 6 blocks (imbalance 40.0%) → 6 blocks (imbalance 40.00%) via 0 cut(s) across 0 block(s)
```

- Pass `--adjacency`: the floor is then checked against the cuts as they
  propagate through the seams, not only block by block.
- When the floor is not met the planner warns and returns its best plan;
  coarsen less (`--ngr`), lower `--min-coarse-cells`, or build the grid
  with more cells.

## See also

- [Load balancing](/concepts/load-balancing): seam propagation, per-axis
  cut counts, the floor, rank assignment.
- [Split after overset](/cookbook/after-overset-split): when the split must
  come after the donor search.
- [Multi-component Xall run](/cookbook/multi-component).
- CLI: [`shore split`](/reference/split), [`shore balance`](/reference/balance),
  [`shore grd`](/reference/grd), [`shore proc-input`](/reference/proc-input).
- API: [`plan_balance`](/reference/api/shore.balance#plan_balance),
  [`split_blocks`](/reference/api/shore.split#split_blocks).
