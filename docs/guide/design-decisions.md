# Design decisions

This page records *why* SHORE's load-balancing and splitting machinery works
the way it does. Each decision was forced by a concrete failure on a real
overset case (a cubed-sphere body inside a Cartesian background, the
shock-sphere case) — the motivations below are not hypothetical; they are the
lessons that shaped the current code.

[[toc]]

## Splitting auto-propagates cuts across seams

**Decision.** [`split_blocks`](/reference/api-split) cuts a block on multiple
axes in one pass, and any cut that subdivides a `SHARED` / `DIRICHLET` seam is
**automatically propagated to the partner block** (to a fixpoint). Callers never
author paired co-splits.

**Why.** The earlier splitter was single-axis-per-call and required the caller
to supply *matching* cuts on both sides of every shared seam by hand. On a
conformal multi-block topology that is both tedious and fragile: a cut on block
A's face implies a cut on block B at the mapped index (through the seam's
`AxisMap`, with the flip and the partner's extent), and getting it wrong
produces a non-conformal seam that the downstream solver rejects (`tri2or: due
piani coincidenti`). Worse, an early load-balance planner emitted per-block cuts
that were *not* paired, so the split silently dropped connectivity.

The fix moves the seam arithmetic into the splitter, where the `AxisMap` and
extents are already known. A cut is a property of an *axis equivalence class*
linked through seams; propagating to a fixpoint makes the whole class conformal
by construction.

```python
from shore.split import BlockSplit, split_blocks
# Cut sub0 on k only — the cut auto-propagates to sub1/sub3 (j-ring SHARED)
# and to both caps (DIRICHLET), so every block ends up k-split conformally.
new = split_blocks(blocks, [BlockSplit("sub0", "k", [40])])
```

Set `propagate_seams=False` only to reproduce the legacy "caller pre-pairs every
co-split" behaviour.

::: tip Consequence for the planner
Because the splitter guarantees conformance, the balance planner does **not**
need to emit paired co-splits — it just chooses per-block cuts. See below.
:::

## Per-axis cut counts (not per-block volume)

**Decision.** The balance planner chooses each axis's chunk count as a pure
function of **that axis's cell extent**, never the block's 3-D volume.

**Why.** Two blocks that share a seam on an axis always have *equal extent* on
that axis (a conformal tiling guarantees it). If the cut count were derived from
the 3-D volume, two neighbours with equal shared-axis extent but different
extents on the *other* axes would get **different** cut counts on the shared
axis — and when the splitter propagates them, the two interleaved cut sets
produce sub-floor sliver chunks (e.g. a block cut at `[30, 60]` receiving a
propagated `44` becomes `[30, 14, 16, 30]` — the `14` violates the multigrid
floor). Keying the count to the axis extent alone makes seam-linked blocks
choose **identical** cuts, so propagation is idempotent and no slivers form.

This is the single change that took the level-2 test case from 150 % imbalance
(slivers backed off) to under 1 % (clean, aligned cuts).

## Multi-way equal partition, not bisection (issue #6)

**Decision.** [`plan_balance`](/reference/api-balance) partitions every block
into roughly equal, roughly cubic chunks sized so the total lands near
`np * blocks_per_rank`, in a single pass.

**Why.** The original planner did *iterated bisection of the current-heaviest
block*. On a 33-block, ~31 M-cell assembly targeting 40 ranks it had two failure
modes ([issue #6](https://github.com/szaghi/shore/issues/6)):

| approach | result | usable? |
|---|---|---|
| bisection, under tolerance | leaves the heavy polar caps **whole**, stalls at ~35 % imbalance | no |
| bisection, forced to converge | shreds wake blocks into **single-cell sheets** (2388 blocks at 0.1 %) | no — halo/MPI overhead dominates |
| **multi-way equal partition** | ~400 well-shaped bricks at a few % | **yes** |

Bisection-of-heaviest has no notion of a *target block count* and no multi-way
partition, so it either under-cuts or over-fragments. A fixed target with
equal-edge chunks avoids both. Cuts are seam-aligned (above), and the result is
floor-checked against seam-*propagated* cuts so the plan matches the split output
block-for-block.

```bash
shore balance assembly.grd --np 16 --ngr 3 --min-coarse-cells 8 \
  --blocks-per-rank 8 --adjacency assembly.adjacency.json -o split.toml
```

## The multigrid cell floor (`--ngr` / `--min-coarse-cells`)

**Decision.** Every chunk axis is kept `>= min_coarse_cells * 2^(ngr-1)` fine
cells and divisible by `2^(ngr-1)`.

**Why.** The solver coarsens each block by `2^(ngr-1)` for its multigrid
levels, and its overset stencils reach ±2 *coarse* cells. A chunk that is too
thin collapses to too few coarse cells → a degenerate donor (`tri2or`). The floor
is the planner's hard constraint; it is **best-effort** — if an input block is
already below the floor (so no split can help), the planner warns and proceeds
rather than failing, reporting the achieved imbalance. This honesty matters:
some geometries genuinely cannot meet a requested tolerance within the floor, and
the user needs the real number, not a crash.

## Per-face natural-BC codes

**Decision.** A `FREE` [`Face`](/reference/api-mesh) can carry an explicit
`bc_code` (a positive Xall natural-BC code in `[1, 19]`), overriding the
per-block `free_face_bc` default in [cc.par](/reference/api-cc-par).

**Why.** A single block often needs *heterogeneous* outer BCs — e.g. inflow on
its upstream face and extrapolation on the others. The old model collapsed all
`FREE` faces of a block to one code, which cannot express that. `bc_code` lets
one block carry inflow on one face and extrapolation on the rest.

::: warning Natural BC codes are POSITIVE
Overset stores the natural-BC code **verbatim** and rejects any face `icc < 1`
(`Errore in calspig: icc sulle facce deve essere > 0`). Use `9` for
free/far-field extrapolation, `3`/`4` for inflow/outflow — **not** a negative
code like `-9`. This corrected a long-standing wrong example in the docs.
:::

```python
# inflow on the -X (i_lo) face; the block's other FREE faces stay default (9)
blk._set_face_bc("i_lo", FaceBC.FREE, bc_code=3)
```

## Cell-weighted rank assignment with local search (issue #7)

**Decision.** [`proc-input`](/reference/api-proc-input) weights blocks by
interior cells `(ni-1)(nj-1)(nk-1)` and refines the greedy LPT seed with a
move/swap local search (pinned blocks held fixed).

**Why.** [Issue #7](https://github.com/szaghi/shore/issues/7) reported node
weighting; investigation showed the `.grd` metadata already carries cell counts,
so the weight was correct — but the **single greedy LPT pass had no
refinement**, leaving a few percent of imbalance a local search recovers (the
issue's case: ~5.6 % → ~1.6 %). The refinement is deterministic (stable
tie-breaking) so results are reproducible, and explicitly-pinned blocks (`--meta`)
are never moved.

## After-overset splitting (the `bsplit` workflow)

**Decision.** When splitting *before* overset breaks the chimera donor search,
run overset on the **unsplit** assembly and split its outputs *afterward* with
the external `bsplit` tool. SHORE drives this by emitting a `bsplit.par` from a
balance plan: `shore balance --format bsplit-par`.

**Why — the failure that forced it.** On the cubed-sphere-in-background case, the
chimera-overlap block (the central plateau the body sits inside) **cannot be
split before overset**. The evidence, isolated experimentally:

| configuration | degenerate donors (`sislin`) |
|---|---|
| unsplit assembly | **0** |
| sphere split, background whole | 0 |
| background split, sphere whole | ~968 |
| any split of the overlap (pre- or post-tiled) | ~968 |

Splitting the background — by *any* means, including re-tiling the plateau at
mesh-generation time — produces degenerate donor hexahedra. The fringe-cell set
is **identical** between split and unsplit (verified: 100 % overlap), so it is
not a geometry or fringe-placement problem: overset's donor *search* returns
worse donors (sliver near-wall cubed-sphere cells, point outside the hexahedron)
when more blocks compete in the overlap. This is a property of the solver's
donor search on a curved, near-wall-stretched donor block — not something any
mesh decomposition or pre-split strategy can avoid, because *overset needs the
overlap region to be one block*, which is incompatible with load-balancing it.

**The resolution** sidesteps the conflict entirely:

1. Run overset on the unsplit assembly → the donor search runs **once on the
   intact overlap** → zero degenerate donors.
2. `bsplit` decomposes the framed overset outputs (`cc.0N`, `cc.0N.grd`) for MPI
   by **pure index remapping** of the donor pool (its `bufrev` inverse map
   relabels every donor `(block, i, j, k)` → `(chunk, local-ijk)`) — **no donor
   re-search**.

The cut indices come from the same `plan_balance`, rendered to bsplit's `.par`
dialect instead of a `split.toml`. The two renderings describe the identical cut
set; only the consumer differs.

```bash
# 1. merge the UNSPLIT components, run overset (clean donor search)
shore grd-merge sphere.grd background.grd -o raw.grd
shore cc-par-merge sphere.cc.par background.cc.par -o raw.cc.par --grd raw.grd
overset                                   # writes cc.0N (+ .grd), 0 degenerate donors

# 2. plan the balance and emit a bsplit.par (grd-only plan)
shore balance raw.grd --np 16 --ngr 3 --min-coarse-cells 8 \
  --blocks-per-rank 8 --format bsplit-par -o bsplit.par

# 3. split the overset outputs (external bsplit), then assign ranks
bsplit                                    # cc.0N -> ccSPLIT.0N (+ .grd)
shore proc-input ccSPLIT.01.grd ccSPLIT.proc.input --np 16
```

See [After-overset splitting](/guide/after-overset-bsplit) for the full
workflow and the file-format details.

::: tip When to use which split
**Pre-overset** ([`shore split`](/reference/split)) is the default — it splits
the geometry before the solver ever runs, and the solver's donor search handles
the (now smaller) blocks. Use **post-overset** (`bsplit`) only when the donor
search degenerates on a split overlap, as above.
:::
