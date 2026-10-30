# Split after the overset pre-processor

Split before overset by default ([Split and balance](/cookbook/split-balance)).
Split after it when cutting the overlap region first breaks the donor
search: the overset pre-processor finds degenerate donor cells
(`Errore in calspig`, `sislin` / `tri2or` errors) only on the split
assembly. The fix runs the donor search once on the unsplit grid, then
splits the overset outputs with the external `bsplit` tool, which relabels
the donor pool by index instead of searching again. SHORE plans the cuts
and writes the `bsplit.par` that drives it.

[Load balancing → After-overset splitting](/concepts/load-balancing#after-overset-splitting-the-bsplit-workflow)
records the failure that made this necessary.

## Run the workflow

```bash
# 1. merge the UNSPLIT components into one assembly
shore grd-merge sphere.grd background.grd -o raw.grd
shore cc-par-merge sphere.cc.par background.cc.par -o raw.cc.par --grd raw.grd

# 2. run overset on the unsplit assembly (external Xall pre-processor)
#    -> cc.01/02/03 (+ .grd) and dist.0N, no degenerate donors
( cd run && overset )

# 3. plan the balance and write bsplit.par (grd-only plan: no --adjacency)
shore balance raw.grd --np 16 --ngr 3 --min-coarse-cells 4 \
      --blocks-per-rank 2 --format bsplit-par -o run/bsplit.par

# 4. split the overset outputs (external bsplit)
#    cc.0N (+ .grd) -> ccSPLIT.0N (+ .grd) at the balanced block count
( cd run && bsplit )

# 5. assign ranks on the split grid (cell-weighted)
shore proc-input run/ccSPLIT.01.grd run/ccSPLIT.proc.input --np 16
```

Steps 2 and 4 run Xall's own tools; SHORE does steps 1, 3 and 5. Step 3 on
an 11-block sphere-plus-background assembly:

```text
✓ Wrote bsplit.par (bsplit-par)
shore balance: 11 blocks (imbalance 200.0%) → 20 blocks (imbalance 0.00%) via 9 cut(s) across 9 block(s)
```

```text
sol
cc
cc
3
0
0
0
1 1 16
2 1 16
...
```

`bsplit.par` holds seven header lines (solution, grid and connectivity
roots `sol`, `cc`, `cc`; the multigrid level count; three flags at 0) and
one cut per line: `<grd block number> <direction 1|2|3> <fine node index>`.
`bsplit` rescales each fine-level index for every multigrid level, so the
same cuts serve `cc.01`, `cc.02` and `cc.03`.

## Notes

- `bsplit` identifies blocks by their 1-based position in the `.grd`, not
  by label. A plan made without `--adjacency` labels its blocks
  `block_<n>` with exactly that number, so `--format bsplit-par` refuses
  `--adjacency`:

  ```text
  Error: --format bsplit-par needs a grd-only plan; do not pass --adjacency (bsplit needs grd block numbers).
  ```

- Without the adjacency the multigrid floor is checked per block, not
  through the seams; the cut set is otherwise what the TOML plan would
  describe.
- In Python, `plan_balance(grd, np=..., ...).to_bsplit_par(ngr=...)` returns
  the same text.

## Overset output files (reference)

The files `bsplit` reads and writes come from Xall's tools, not from SHORE;
this summary is for orientation. All are Fortran unformatted binary, one
per multigrid level `N` (`cc.01` fine … `cc.03` coarsest). For a block of
`(ni, nj, nk)` interior cells at that level:

**`cc.0N.grd`**, framed node grid:

```
record: nbl                             (block count)
nbl ×   record: ni, nj, nk              (interior cells per block)
nbl ×   record: X[(ni+5)(nj+5)(nk+5)]   ┐ interior nodes plus 2 ghost
        record: Y[...]                  │ layers each side, i fastest
        record: Z[...]                  ┘
```

**`cc.0N`**, chimera connectivity:

```
record: nbl
nbl ×   record: ni, nj, nk
nbl ×   record: ICC[(ni+4)(nj+4)(nk+4)]   framed cell-centred int32 connectivity
record: prc                                (size of the donor pool)
record: bufrcc[1..prc]                     packed donor pool (real)
```

Each chimera cell's `ICC` is a 1-based pointer into `bufrcc`; there it finds
`[type, nd, then nd × (donor_block, i, j, k, weight)]`.

**`dist.0N`**: framed per-cell wall distance, framed as `ICC`.

`bsplit` builds an inverse map from every old framed cell to its new
`(chunk, local ijk)` and rewrites each donor reference through it: an index
relabelling, with no geometry and no new donor search. Cut faces become
adjacency boundaries; natural BCs are kept. Each chunk's grid is a copy of
its node range from the parent block. The split files load in the solver
as the unsplit ones would.

## See also

- [Split and balance](/cookbook/split-balance): the default, before overset.
- [Load balancing](/concepts/load-balancing).
- CLI: [`shore balance`](/reference/balance) (`--format bsplit-par`),
  [`shore proc-input`](/reference/proc-input),
  [`shore grd-merge`](/reference/grd-merge),
  [`shore cc-par-merge`](/reference/cc-par-merge).
- API: [`BalancePlan.to_bsplit_par`](/reference/api/shore.balance#BalancePlan-to_bsplit_par).
