---
title: shore.io.plot3d
outline: [2, 3]
---

# `shore.io.plot3d`

Plot3D export of structured block lists (ADR 0035).

Multi-block, whole (no IBLANK), 3-D grid file, Fortran unformatted
sequential, double precision, little endian (the layout ParaView's and
VTK's PLOT3D readers and most CFD tools read):

  Record 1 : nblocks                       - int32
  Record 2 : (ni, nj, nk) per block        - 3 \* nblocks int32, \*node\* counts
  Record 2 + n : X, Y, Z of block n        - 3 \* ni \* nj \* nk float64,
                                             each i-fastest

Plot3D carries no block names, connectivity or boundary conditions; CGNS
(:mod:`shore.io.cgns`) does.

## `write_plot3d` <Badge type="info" text="function" /> {#write_plot3d}

```python
write_plot3d(blocks: Sequence[HexBlock], path: str | Path) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/plot3d.py#L30">source</a>

Write *blocks* as one multi-block Plot3D grid file at *path*.
