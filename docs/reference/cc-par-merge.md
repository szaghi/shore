# `shore cc-par-merge`

Concatenate multiple `cc.par` files into one for an Xall multi-mesh run.

```bash
shore cc-par-merge INPUT1.cc.par INPUT2.cc.par ... -o OUTPUT.cc.par --grd MESH.grd
```

## Motivation

Several independent Chimera components in one Xall run (e.g. a
near-body wall + a far-field background) are generated as separate
`.grd` + `cc.par` pairs. To submit them as a single Xall job, the two
`cc.par` files have to be combined into one whose blocks list is the
concatenation of the inputs, with patch block-indices and `family`
pointers shifted so they still refer to the right targets.

`shore cc-par-merge` does the merge:

```bash
shore cc-par-merge background.cc.par wall.cc.par \
                   -o assembly.cc.par --grd assembly.grd
```

The matching binary-grid concatenation is [`shore grd-merge`](./grd-merge);
the two commands must be invoked with the same input ordering. See
[Coordinated invocation](#coordinated-invocation) below.

## What gets shifted, what doesn't

For inputs in caller-specified order `[A, B, C, ...]`:

- **Block list**: simple concatenation. Block 1 of `B` becomes
  block `len(A.blocks) + 1` in the merged file, and so on.
- **Patch `block` field**: each patch in `B` has its 1-based block
  index shifted by `len(A.blocks)`; same recipe per subsequent input.
- **Patch `family` field**: `cc.par` connection patches reference
  their partner by 1-based patch index. Since each input is
  internally complete (no SHARED across the merge boundary), every
  family pointer stays *inside* its own input — but the absolute
  index changes. The merger shifts `family` by `len(A.patches)` for
  patches in `B`, etc.
- **Patch extents** (Imin/Imax/Jmin/Jmax/Kmin/Kmax): unchanged.
- **Boxes**: concatenated; each box's associated-block index is
  shifted by the cumulative block count of the preceding inputs,
  mirroring the patch-`block` shift. Box vertex coordinates pass
  through unchanged.

Symmetric pairing survives the merge: every connection patch with
`bc ≥ 100` continues to point at a partner patch whose `family`
points back at it.

## Header reconciliation

Header fields other than the `.grd` filename are taken from `inputs[0]`
verbatim (output basename, ghost-cells flag, multigrid levels,
boundary-layer thickness, numerical-beach width). The `.grd` filename
on header line 1 is overridden by the `--grd` flag — that's how the
merged file knows which `.grd` it pairs with.

If subsequent inputs disagree on any header field, a `UserWarning` is
emitted but the merge proceeds with `inputs[0]`'s value. Surface drift
loudly without blocking.

## Block-comment auto-prefixing

Each block line in the merged `cc.par` has its trailing comment
prefixed with the input file's stem so blocks from different
topologies are distinguishable when reading the merged file:

```
        1         0         1 ! 1 : background.center
        2         0         1 ! 6 : wall.sub0
        2         0         1 ! 10 : wall.cap_north
```

`cc.par` treats comments as advisory only — Xall indexes blocks by
their 1-based position — so this is purely for human readability.
A trailing `.cc` in the input filename (typical for `wall.cc.par`)
is stripped from the prefix.

## Coordinated invocation

`cc-par-merge` and [`grd-merge`](./grd-merge) must be called with the
**same** input ordering — `cc.par` patch block-indices are 1-based
positions into the `.grd` block list. A typical end-to-end pair:

```bash
shore grd-merge    background.grd    wall.grd    -o assembly.grd
shore cc-par-merge background.cc.par wall.cc.par -o assembly.cc.par \
                   --grd assembly.grd
```

The example at `examples/05-chimera-assembly/cc_par_assembly.py` builds the wall and
background meshes from scratch, runs both merges, and prints the
equivalent CLI commands. See the
[cc.par + boxes pipeline](/guide/cc-par-pipeline) walkthrough for the
end-to-end story.

## Python API

```python
from shore.io.cc_par import merge_cc_par

merge_cc_par(
    ["background.cc.par", "wall.cc.par"],
    "assembly.cc.par",
    grd_filename="assembly.grd",
)
```

Raises `ValueError` for fewer than 2 inputs, `FileNotFoundError` for
a missing input.

## See also

- [`shore cc-par`](./cc-par) — produce one `cc.par` per Chimera
  component.
- [`shore grd-merge`](./grd-merge) — companion `.grd` merger;
  invoke as a pair.
- [Marker boxes](./boxes) — boxes section of `cc.par`, also handled
  by the merger.
