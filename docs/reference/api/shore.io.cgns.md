---
title: shore.io.cgns
outline: [2, 3]
---

# `shore.io.cgns`

CGNS export of structured block lists (ADR 0035).

One CGNS Base per component, one structured Zone per block, written in the
CGNS/HDF5 file mapping with `h5py` (optional dependency, `shore[cgns]`).

- `GridCoordinates`: `CoordinateX/Y/Z`, double precision.
- `ZoneGridConnectivity`: one `GridConnectivity1to1_t` per connected
  face (SHARED, DIRICHLET, PERIODIC), with `PointRange`,
  `PointRangeDonor` and `Transform` derived from the face's
  :class:`shore.mesh.face.AxisMap`. CGNS 1-to-1 connections are point
  matched: a face whose nodes do not coincide with its partner's is refused.
- `ZoneBC`: one `BC_t` per other face: WALL is `BCWall`, SYMMETRY is
  `BCSymmetryPlane` and FREE (an overset fringe or a natural boundary,
  decided downstream) is `BCTypeUserDefined` with the descriptor
  `shore:FREE`.

Indices are 1-based and Fortran-ordered as the SIDS require; an array of
SIDS dimensions `(d1, d2, ...)` is stored as an HDF5 dataset of shape
`(..., d2, d1)`. The O-grid's periodic j seam is not a face of the block
graph (as in `cc.par`), so it is not written either.

## `write_cgns` <Badge type="info" text="function" /> {#write_cgns}

```python
write_cgns(
    components: Mapping[str, Sequence[HexBlock]],
    path: str | Path,
    *,
    tolerance: float = 1e-09,
) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/cgns.py#L48">source</a>

Write *components* (`{base name: blocks}`) as one CGNS/HDF5 file.

**Parameters**

- **`components`** — Base name -&gt; its blocks; the bases are written in this order and the zones in block order. Connections must stay inside a base.
- **`path`** — Output file (conventionally `.cgns`).
- **`tolerance`** — Largest distance between a connected face's node and its partner's, relative to the base's bounding-box diagonal.

**Raises**

- `ParameterError` — A base or zone name longer than 32 characters or repeated, a connection to a block outside its base, or a connection whose nodes do not coincide.
