---
title: shore.io.grd
outline: [2, 3]
---

# `shore.io.grd`

Xall binary `.grd` writer.

Format
------
Fortran unformatted sequential binary (each `write` statement wraps its
payload in 4-byte record-length markers):

  Record 1   : nb                          - int32, number of blocks
  Records 2..nb+1 : ni_n, nj_n, nk_n      - 3x int32 per block,
                                             stored as *cell* counts
                                             (= node count - 1)
  For each block n = 1..nb:
    Record nb+1+3(n-1)+1 : X flat          - ni_n\*nj_n\*nk_n float64,
                                             i-fastest (Fortran column-major)
    Record nb+1+3(n-1)+2 : Y flat
    Record nb+1+3(n-1)+3 : Z flat

Reference implementation: `geogrd/src/geogrd.f90` in
`Xnavis-projects/geogrd`.

Note on node vs. cell counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
`geogrd` reads `.geo` files whose first line is the *node* count
`(ni, nj, nk)` (the SHORE convention) and immediately subtracts 1
before writing to `.grd` (`ni = ni - 1`).  SHORE `(ni, nj, nk, 3)`
arrays therefore map directly: pass the array as-is and this writer
performs the subtraction.

## `write_grd` <Badge type="info" text="function" /> {#write_grd}

```python
write_grd(path: str | Path, geo_paths: list[str | Path]) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/grd.py#L38">source</a>

Write multiple structured blocks to an Xall binary `.grd` file.

Reads each block from an ASCII `.geo` file and packs all blocks into
one Fortran unformatted binary `.grd`.  The input list determines both
the block count and block ordering, which must match the ordering used
when generating the companion `cc.par`.

**Parameters**

- **`path`** — Output file path.  `.grd` extension is appended if absent.
- **`geo_paths`** — Ordered list of `.geo` file paths, one per block.

**Returns**

- `Path` — The path actually written.

## `read_grd_metadata` <Badge type="info" text="function" /> {#read_grd_metadata}

```python
read_grd_metadata(path: str | Path) -> list[tuple[int, int, int]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/grd.py#L90">source</a>

Return the per-block `(ni, nj, nk)` cell counts from a `.grd` file.

Parses the `.grd` header (block count + per-block shape records)
without decoding any coordinate payload.  Cheap — costs O(n_blocks)
bytes of I/O regardless of mesh size.

The returned tuples are **cell counts** (the form stored in the
`.grd` header), matching the values written by
:func:`write_grd`.  Block ordering is the file order — block 1 is
the first tuple, block 2 the second, etc.

**Parameters**

- **`path`** — Path to a `.grd` file.

**Returns**

- `list[tuple[int, int, int]]` — One `(ni, nj, nk)` per block, in file order.

**Raises**

- `FileNotFoundError` — If `path` does not exist.
- `ValueError` — If the file is too short, the header markers are malformed, or the block count is non-positive.

## `merge_grd` <Badge type="info" text="function" /> {#merge_grd}

```python
merge_grd(inputs: list[str | Path], out: str | Path) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/grd.py#L177">source</a>

Concatenate multiple `.grd` files into one.

Used when an Xall run combines several independent Chimera
components (e.g. wall + background): the merged `.grd` lists all
blocks of *inputs[0]* first, then *inputs[1]*, etc., with the header
block count rewritten to the sum.  Block payloads are passed
through byte-for-byte; no node values are touched.

The companion :func:`shore.io.cc_par.merge_cc_par` shifts patch
block-indices and family pointers to match this concatenation.
Invoke them as a pair, with the **same** input ordering.

**Parameters**

- **`inputs`** — Two-or-more input `.grd` paths, in the desired output order.
- **`out`** — Output path; `.grd` extension is appended if absent.

**Returns**

- `Path` — The path actually written.

**Raises**

- `ValueError` — If `inputs` has fewer than 2 entries, or any input is empty / malformed.
