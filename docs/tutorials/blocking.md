---
title: Blocking a wing–body
---

# Blocking a wing–body

This tutorial builds a multi-block grid round a half wing–body from an empty
box: you split a blocking skeleton, cut a hole where the body is, tie the
hole's faces to the CAD body and its base edges to the wing–fuselage
junction, wrap the hole in an O-grid, and premesh it with exact projection,
wall-normal lines and a 2e-3 first cell. Then you change the CAD model's
span and fuselage radius and replay the journal to get the new grid. It is
the phase-7 acceptance case ([ADR 0051](/adr/0051-wing-body-acceptance)),
and every number below is checked by `tests/spec/test_wing_body.py`.

::: info You need
- SHORE with the `cad` extra (build123d, for STEP) and the `cgns` extra
  (h5py): `pip install 'shore-mesh[cad,cgns]'`.
- About 10 s for the script
  [`examples/tutorials/blocking.py`](https://github.com/szaghi/shore/blob/main/examples/tutorials/blocking.py)
  (8 s on a laptop, both models included).
:::

<Figure name="hero-wing-body" caption="The premeshed half wing–body: wall faces coloured per block, a slice of the volume grid on the symmetry plane." />

## 1. Build the CAD model

The model is a half model, cut at z = 0, which becomes a symmetry plane: a
capsule fuselage along x (radius 0.5, cylinder length 5) and a NACA 0012
wing with a blunt trailing edge (chord 1 from x = 0) along +z to the tip at
z = 3, fused into one solid and written as STEP.

<<< @/../examples/tutorials/blocking.py#model{python}

<Figure name="cad-entities" caption="The CAD wing–body: faces coloured per face, edges as dark curves." />

Blocking needs exact geometry: SHORE reads the STEP file's faces, edges and
vertices (`OcctGeometry`), and can project onto them exactly. A fused STL
would do for the walls, but not for the junction curve.

## 2. Start the blocking

A blocking is a skeleton of hex blocks whose edges carry node counts
([ADR 0042](/adr/0042-blocking-skeleton-premesh)). You start with one block
over the far field: 53 × 33 × 25 nodes on [-6, 7] × [-4, 4] × [0, 6]. The
z = 0 face lies on the symmetry plane.

<<< @/../examples/tutorials/blocking.py#box{python}

```text
CAD model: 8 faces, 19 edges
```

## 3. Split it round the body

`split_blocking(block, axis, fraction)` cuts a block across an axis and
carries the cut through every block sharing the cut edges, so the sheet
runs through the whole blocking. The helper `cut` finds the block holding a
point on the cut and turns the coordinate into that block's fraction. The
cuts bound the fuselage (x = ±3, |y| < 0.8, z < 0.8), the wing root
(-0.3 < x < 1.3, |y| < 0.3) and the wing beyond its tip (z < 3.3).

<<< @/../examples/tutorials/blocking.py#splits{python}

```text
blocks after the splits: 75
```

A cut lands on the node nearest the fraction (the box has 4 nodes per unit
length, so x = -0.3 lands on x = -0.25), and the blocks on either side keep
the nodes they hold. In the GUI you cut through a
picked edge node instead ([The same in the GUI](/tutorials/gui)).

## 4. Find the hole and its walls

The hole is the 10 blocks the body occupies: the 3 × 3 blocks round the
fuselage and the block round the wing. Their faces that bound the hole are
the faces owned by only one of them; the faces on z = 0 are left out, since
they lie on the symmetry plane, not on the body.

<<< @/../examples/tutorials/blocking.py#hole{python}

```text
hole: 10 blocks, 25 wall faces
```

<Figure name="blocking-skeleton" caption="The blocking skeleton before the premesh: edges as lines, block-centre markers; the hole round the body." />

## 5. Associate the walls and the junction

An association says where a blocking entity belongs on the geometry; the
premesh projects onto it ([ADR 0003](/adr/0003-association),
[ADR 0039](/adr/0039-association-editing)). The 25 wall faces go to the
body as a whole (`kind="body"`): each node lands on the nearest point of any
CAD face.

The 4 edges at the base of the wing's block run round the wing root, on the
fuselage. Associated with the body only, their nodes would land anywhere
near the junction, and the grid would cut the corner. Associated with the
junction curves (the CAD edges lying on both the fuselage and the wing:
upper, lower and trailing-edge base), they land on the junction itself.
Without these 4 associations no wall node reaches the junction: all stay
more than 1e-3 away (spec case A5).

<<< @/../examples/tutorials/blocking.py#associate{python}

```text
junction curves: [6, 7, 8]
```

An edge without an association of its own inherits that of its face; the
base edges have their own, so the premesh puts them on the curves and the
faces round them on the body. The journal stores
the curve indices together with fingerprints, so an edited model's edges
are found again even if their indices change.

## 6. Delete the hole and wrap it in an O-grid

Deleting the 10 blocks leaves a hole whose faces are now owned by the
blocks round it. Those faces still carry the associations (they are the same
skeleton faces). A blocking face mapped straight onto a curved wall gives
skewed cells at the hole's corners, so `ogrid` inserts a layer of blocks on
the hole's 25 faces ([ADR 0043](/adr/0043-ogrid-insertion)): every vertex of a layered
face gets an inner copy moved away from the face by `(1 - factor)` times half
its edges leaving the face (`factor=0.7`: 15 % of the neighbouring block's
depth), and `nodes=17` sets the node count across the layer.

<<< @/../examples/tutorials/blocking.py#ogrid{python}

```text
blocks after the O-grid: 90
```

The new layer blocks take the walls; the old blocks keep the outer part.

## 7. Premesh, orthogonal lines, wall layer

Three commands make the grid:

<<< @/../examples/tutorials/blocking.py#premesh{python}

- **`premesh(exact=True)`** meshes every edge and face once, projects them
  onto their associations and copies them into every block that shares them,
  so seams match bit for bit; block interiors are filled by transfinite
  interpolation. `exact=True` is the tier-2 projection onto the CAD
  surfaces and curves. The tessellation-based tier 1 left wall nodes on the
  junction up to 4e-3 off the body, because the CAD edge's tessellation and
  the faces' tessellation do not coincide; exact projection puts them within
  1e-9.
- **`orthogonalize(0.5)`** reshapes the lines leaving each wall into cubic
  curves that start along the wall normal
  ([ADR 0048](/adr/0048-orthogonal-wall-lines)). Each line becomes a cubic
  Hermite curve whose start tangent is the normal times `strength` times the
  line's chord; the nodes keep their arc-length fractions. Strength 1 is
  refused here (the command refuses to invert a valid block): at the concave
  junction the normal is the bisector of the fuselage and wing normals, and
  the wing's side layers would invert. 0.5 brings the wall angle from 24.7°
  to 9.0° median (95th percentile 68.6° to 25.3°, ADR 0051).
- **`wall_layer(2e-3)`** clusters every wall line geometrically so its first
  cell is exactly 2e-3 long ([ADR 0050](/adr/0050-wall-layer)).

```text
premeshed in 3 s
```

<Figure name="blocking-premesh" caption="The premeshed blocking: the wall grid on the body and the grid on the symmetry plane." />

<Figure name="wall-layer" caption="The wall layer at the wing–body junction: first cell 2e-3, lines leaving the wall along its normal." />

## 8. Check the grid

`report` measures what the acceptance asks for: the smallest scaled
Jacobian, how far the wall nodes are from the CAD body (exact projection),
the first cell lengths, and the angle between each wall node's first grid
segment and the CAD normal.

<<< @/../examples/tutorials/blocking.py#quality{python}

```text
original: 90 blocks, min scaled Jacobian 0.173, wall off the body 9.5e-10, first cell 0.002-0.002, wall angle median 14.3 deg
```

The median here covers every wall node off the symmetry plane. The trailing
edge, the tip and the junction have no single normal, and their nodes keep
angles near 45°; without them the median is 9.0° (ADR 0051).

## 9. Write the solver files

<<< @/../examples/tutorials/blocking.py#write{python}

The journal is saved before the writes, so it holds the 46 commands that
build the grid and nothing else. `write_cgns` writes one zone per block with
1-to-1 connections for the seams and `BCWall` for the walls. The premesh
makes unassociated boundary faces FREE, so the faces on z = 0 and the far
field are written as `BCTypeUserDefined` (`shore:FREE`): set the symmetry
condition in your solver set-up. `write_xall` writes the `.geo` files, `mesh.grd`,
`cc.par` and `proc.input` under `xall/`.

## 10. Edit the model and replay

The edit: span 2.7 instead of 3 and fuselage radius 0.55 instead of 0.5.
The journal reads `wb.step`, so writing the edited file under the same name
in `edited/` and replaying with `input_dir="edited"` runs the same 46
commands on it.

<<< @/../examples/tutorials/blocking.py#replay{python}

```text
edited: 90 blocks, min scaled Jacobian 0.135, wall off the body 2.5e-09, first cell 0.002-0.002, wall angle median 15.1 deg
total 7 s
```

The blocking itself has not moved: the same box, cuts and hole. The fatter
fuselage moves the junction, and the junction edges are found again by their
fingerprints, not by index; the premesh then projects onto the new faces and
curves. The edit must stay inside the hole: a span beyond 3.3 would pierce
the wing's block, and the blocking would need new cuts.

## Limits

From ADR 0051:

- a half model; full-span and multi-component models are open;
- edits must stay inside the blocking's hole;
- concave junctions need `orthogonalize` below strength 1;
- sharp edges (trailing edge, tip) keep wall angles near 45°;
- wall layers are algebraic: hyperbolic marching of a blocking's walls is
  deferred.

## What you learned

- A blocking starts as a box and is shaped by splits, deletions and O-grids;
  each is a journaled command.
- Faces are associated with the body, edges with the curves they must
  follow; the junction needs its own association.
- `premesh(exact=True)` projects onto the CAD entities exactly and keeps
  every seam closed; `orthogonalize` and `wall_layer` shape the near-wall
  grid.
- The journal replays on an edited CAD model and gives a valid grid on the
  new geometry.

## Next

- [The same in the GUI](/tutorials/gui): this blocking built with the mouse.
- [Blocking](/concepts/blocking) and [Association](/concepts/association):
  how it works.
- [Blocking recipes](/cookbook/blocking): vertex moves, partial O-grids.
- [Check grid quality](/cookbook/quality),
  [Export CGNS and Plot3D](/cookbook/exports).
