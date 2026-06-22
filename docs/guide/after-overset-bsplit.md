# After-overset splitting (`bsplit` workflow)

The default load-balancing path splits the geometry **before** the overset
preprocessor runs ([`shore balance`](/reference/balance) →
[`shore split`](/reference/split)). For most cases that is the right approach.

But on some overset assemblies, splitting the chimera-overlap region *before*
overset breaks the solver's donor search — it produces degenerate donor
hexahedra (`Errore in calspig` / `tri2or`). When that happens, the fix is to
split **after** overset: run the donor search once on the intact, unsplit
overlap, then decompose the framed overset outputs for MPI with the external
`bsplit` tool. SHORE drives this by emitting a `bsplit.par` from a balance plan.

See [Design decisions → after-overset splitting](/guide/design-decisions#after-overset-splitting-the-bsplit-workflow)
for *why* this is necessary (with the experimental evidence).

::: tip When to use this
Use the **pre-overset** split by default. Reach for the **post-overset**
`bsplit` workflow only when a body-fitted block overlaps a background block and
splitting that overlap degenerates the donor search — the classic symptom is
`sislin`/`tri2or` errors that appear *only* when the assembly is split.
:::

## The idea in one line

Overset runs the donor search on the **unsplit** grid (clean); `bsplit` then
relabels the already-computed donor pool onto the split blocks by **pure index
remapping** — no donor re-search.

## Workflow

```bash
# 1. merge the UNSPLIT components into one assembly
shore grd-merge sphere.grd background.grd -o raw.grd
shore cc-par-merge sphere.cc.par background.cc.par -o raw.cc.par --grd raw.grd

# 2. run overset on the unsplit assembly (external solver preprocessor)
#    -> writes cc.01/02/03 (+ .grd) and dist.0N, with 0 degenerate donors
( cd run && overset )

# 3. plan the balance and emit a bsplit.par  (grd-only plan, NO --adjacency)
shore balance raw.grd --np 16 --ngr 3 --min-coarse-cells 8 \
      --blocks-per-rank 8 --format bsplit-par -o run/bsplit.par

# 4. split the overset outputs (external bsplit)
#    -> cc.0N (+ .grd) -> ccSPLIT.0N (+ .grd), now the balanced block count
( cd run && bsplit )

# 5. rank assignment on the split grid (cell-weighted)
shore proc-input run/ccSPLIT.01.grd run/ccSPLIT.proc.input --np 16
```

The case script `build_chimera_bsplit.sh` (in the shock-sphere example case)
automates steps 1, 3 and 5 and prints the exact commands for the manual overset
and bsplit steps.

## Why `--format bsplit-par` needs a grd-only plan

`bsplit` identifies blocks by their **1-based `.grd` block number**, not by a
label. A grd-only `shore balance` plan (no `--adjacency`) uses synthetic
`block_<n>` labels whose `n` is exactly that block number, so the mapping is
direct. Passing `--adjacency` with `--format bsplit-par` is rejected for this
reason. (The cut set is identical to what the TOML plan would describe; only the
rendering differs.)

`bsplit` rescales each fine-level cut index per multigrid level internally, so
the plan's fine-level cuts feed straight through to all of `cc.01/02/03`.

## What `bsplit` produces

For each multigrid level `N`, given inputs `cc.0N` and `cc.0N.grd`, `bsplit`
writes `ccSPLIT.0N` and `ccSPLIT.0N.grd` — the same files at the new (balanced)
block count.

## Overset output file formats (reference)

All files are Fortran-unformatted binary, per multigrid level `N`
(`cc.01` = fine … `cc.03` = coarsest). For a block with header dims
`(ni, nj, nk)` = **interior cell counts** at that level:

### `cc.0N.grd` — framed vertex grid
```
record: nbl                       (block count)
nbl ×   record: ni, nj, nk        (interior cell dims per block)
nbl ×   record: X[(ni+5)(nj+5)(nk+5)]   ┐ full FRAMED grid: interior vertices
        record: Y[...]                  │ plus 2 ghost layers each side,
        record: Z[...]                  ┘ i-fastest
```

### `cc.0N` — chimera connectivity
```
record: nbl
nbl ×   record: ni, nj, nk
nbl ×   record: ICC[(ni+4)(nj+4)(nk+4)]   framed cell-centred int32 connectivity
record: prc                                (size of the donor pool)
record: bufrcc[1..prc]                     packed donor pool (real)
```
Each chimera cell's `ICC` is a 1-based pointer into `bufrcc`; at that pointer:
`[type, nd, then nd × (donor_block, i, j, k, weight)]`. The donor `(block, i,
j, k)` references are what `bsplit` remaps onto the split block numbering.

### `dist.0N` — wall distance
Framed per-cell wall distance (geometry-derived), same framing as `ICC`.

## How `bsplit` remaps the donor pool

`bsplit` builds an **inverse map** (`bufrev`): for every *old* framed cell it
records which *new* `(chunk, local-ijk)` it became. It then relabels every donor
reference `(block, i, j, k)` in `bufrcc` through that map — a pure index
relabel, no geometry and no donor re-search. Split-induced internal boundaries
become adjacency BCs; natural BCs are preserved (promoted to edge/vertex codes at
cut corners). The grid (`.grd`) split is a straight copy of each chunk's framed
vertex range from its parent at the right offset.

The result references the new block numbering consistently — every donor pointer
in range, every donor block valid — so the split outputs load in the solver
exactly as the unsplit ones would, but decomposed for MPI.

## See also

- [`shore balance`](/reference/balance) — `--format bsplit-par`.
- [`shore proc-input`](/reference/proc-input) — final rank assignment.
- [Design decisions](/guide/design-decisions) — the failure that forced this.
