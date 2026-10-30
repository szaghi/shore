# Rebuild a block face on a CAD surface

Use this when a block face must lie on one face of a STEP, IGES or BREP model:
you keep the face's four boundary rows and SHORE fills the interior on the
surface, by transfinite interpolation and elliptic smoothing in its (u, v)
parameters.

The script needs the `cad` extra (`pip install "shore-mesh[cad]"`) and runs in
about 3 s.

## Find the surface under a block face

<<< @/../examples/cookbook/cad_face.py#find_surface{python}

```text
3 surfaces; the face under the probe is surface 0
```

- A CAD geometry lists its faces as SURFACE entities:
  `geometry.entities(EntityKind.SURFACE)`. The `surface=` argument below is the
  index into that tuple.
- Projecting a point you know lies on the wanted face (`BODY` projection)
  returns the index of the face it lands on in `entity_index`. The
  [Geometry](/concepts/geometry) page describes the entity model.

## Fill the face interior on the surface

<<< @/../examples/cookbook/cad_face.py#fill_face{python}

```text
largest distance from the CAD face: 2.2e-16
largest node move: 2.4e-05
other nodes unchanged: True
p.surface_grid('annulus', 'i_lo', surface=0)
```

- The four boundary rows are projected to (u, v) once and kept bit for bit;
  the interior is TFI in (u, v), then Winslow smoothing with control functions
  that keep the boundary clustering. Every interior node is evaluated on the
  surface, so it lies on the CAD face to round-off.
- On this cylinder (a developable surface whose (u, v) is a scaled isometry)
  the primitive's own clustered face is already the exact answer: the command
  reproduces it to 2.4e-5. On a curved, non-developable face it is the step
  that puts the face on the surface with a smooth interior.
- `smoothing="none"` keeps the plain TFI in (u, v).
- Only the face moves; the rest of the block is unchanged. If the face is
  SHARED or DIRICHLET with a neighbour, the neighbour's copy is written through
  ([ADR 0040](/adr/0040-seams-kept-closed)).
- It is a journaled command: undo restores the nodes, and a saved project
  replays it bit for bit. With several geometries, pass `geometry="name"`.

## Limits

The command raises `ParameterError` and changes nothing when:

- the geometry has no parametric surfaces (an STL body);
- `surface` is out of range, or the block or face name is unknown;
- a boundary node is off the surface by more than 1e-6 of the model diagonal:
  put the boundary on the face first (for example with
  [associate and snap](/cookbook/associate));
- the face boundary crosses the periodic seam of the surface's
  parametrisation (here the cylinder seam, which the sector avoids), or
  collapses to a point at a pole or apex. Split the block at the seam instead.

## See also

- [Surface grids](/concepts/surface-grids): the TFI and elliptic smoothing in (u, v)
- [Load STL and CAD](/cookbook/load-geometry)
- [`shore.surface.tfi`](/reference/api/shore.surface.tfi),
  [`shore.project.project`](/reference/api/shore.project.project)
- [ADR 0007](/adr/0007-surface-grids)
