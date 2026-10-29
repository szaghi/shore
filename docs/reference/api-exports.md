# shore.io.cgns and shore.io.plot3d

Exports of structured block lists ([ADR 0035](../adr/0035-cgns-plot3d-export.md)).
The CLI front ends are [`shore cgns`](./cgns) and [`shore plot3d`](./plot3d);
on a project, `Project.write_cgns(path, components=())` and
`Project.write_plot3d(path, components=())` are the journaled output
commands `WriteCgns` and `WritePlot3d`.

## `write_cgns`

```python
from shore.io.cgns import write_cgns

write_cgns(components: Mapping[str, Sequence[HexBlock]], path, *, tolerance=1e-9) -> Path
```

One CGNS/HDF5 file (needs `h5py`, the `[cgns]` extra). `components` maps
a base name to its blocks; bases and zones are written in order.
Connections must stay inside a base.

| Parameter | Description |
|---|---|
| `components` | Base name → its blocks. |
| `path` | Output file (conventionally `.cgns`). |
| `tolerance` | Largest distance between a connected face's node and its partner's, relative to the base's bounding-box diagonal. |

Raises `ParameterError` for:
- a base or zone name longer than 32 characters, repeated, or containing `/`;
- a connection to a block outside its base;
- a connection whose nodes do not coincide.

`transform_vector(face, partner, along)` returns the CGNS `Transform` of a
1-to-1 face connection from its `AxisMap.along`.

## `write_plot3d`

```python
from shore.io.plot3d import write_plot3d

write_plot3d(blocks: Sequence[HexBlock], path) -> Path
```

One multi-block, whole, 3-D, double-precision Plot3D grid file (Fortran
unformatted, little endian). See [`shore plot3d`](./plot3d) for the record
layout.
