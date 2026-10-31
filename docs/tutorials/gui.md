---
title: The same in the GUI
---

# The same in the GUI

This tutorial builds the half wing–body grid of
[Blocking a wing–body](/tutorials/blocking) with the mouse: forms for the
geometry, the box, the splits and the grid commands; the Blocking dock to
select faces, edges and blocks; the Association dock to pick the body and
the junction curves. Every action is a project command and a line of the
journal, so the session you save replays headless and gives the same grid
bit for bit. `tests_gui/spec/test_gui_wing_body.py` drives exactly these
steps through the window ([ADR 0051](/adr/0051-wing-body-acceptance)).

::: info You need
- SHORE with the `gui` and `cad` extras: `pip install 'shore-mesh[gui,cad]'`.
- The STEP model `wb.step`: run step 1 of
  [`examples/tutorials/blocking.py`](https://github.com/szaghi/shore/blob/main/examples/tutorials/blocking.py),
  or the whole script.
- A GPU driver with Metal, Direct3D 12 or Vulkan. Under WSL2 wgpu finds
  only a CPU adapter and the view is slow: run the GUI with native Windows
  Python.
- The commands themselves take seconds; the premesh about 3 s.
:::

## 1. Open the window

```bash
shore gui            # or: shore-gui
shore gui proj       # open a saved project (add --no-replay to skip the rebuild)
```

<Figure name="gui-window" caption="The SHORE window with the wing–body grid loaded." />

| Area | What you use it for here |
|---|---|
| **Model** / **Geometry** (left) | blocks and their faces with BC and association; the CAD model's faces, edges and vertices |
| Viewport (centre) | orbit (left button), pan (right), zoom (wheel); **double click** picks |
| **Command** (right) | a form per command; **Run** executes it |
| **Blocking** (right) | select skeleton edges, faces, blocks; split, O-grid, delete, associate |
| **Association** (right) | the picked geometry entity |
| **Console** (bottom) | **Journal**: one line per command; **Log**: what ran, and errors |
| Toolbar | New, Open…, Save, Save as…, Undo, Redo, Fit, Quality |

All controls are listed in the [GUI reference](/reference/gui).

## 2. Load the geometry

In the Command dock pick **Load geometry**. Set *STL / STEP / IGES / BREP*
to `wb.step` (the folder button browses) and *name* to `wb`, then **Run**.
The model appears in the viewport, and the Geometry tab lists its faces,
edges and vertices. The Journal tab shows:

```python
p.load_geometry('wb.step', name='wb')
```

## 3. Create and split the blocking

Pick **Create blocking**: *ni* `53`, *nj* `33`, *nk* `25`, *vmin*
`(-6, -4, 0)`, *vmax* `(7, 4, 6)`; **Run**. The skeleton is drawn as orange
lines with their nodes. Then pick **Symmetry plane**: *point* `(0, 0, 0)`,
*normal* `(0, 0, 1)`; **Run**. The faces on z = 0 will be premeshed as
`SYMMETRY` faces and stay on the plane through the orthogonal lines
([ADR 0053](/adr/0053-blocking-symmetry-planes)).

Make the 10 cuts of the headless tutorial: x = −3, −0.3, 1.3, 3;
y = ±0.8, ±0.3; z = 0.8, 3.3. Two ways:

- **Split here:** in the Blocking dock set *select* to `edge`, double click
  an orange edge at the node you want (the dock shows `block edge node n of
  N`, the status bar the node's `xyz`), then **Split here**. The cut runs
  through the picked node across the whole sheet.
- **Split blocking form:** *block*, *axis* and *fraction*, as
  `p.split_blocking` takes them. The spec uses this form, with the fraction
  of the cut coordinate along the block that holds it (the headless
  tutorial's `cut` helper computes it).

Either way, each cut is one `p.split_blocking(...)` line, and the result is
75 blocks.

## 4. Associate the hole's faces with the body

The hole is the 9 blocks round the fuselage and the block round the wing
(x in [−0.3, 1.3], |y| < 0.3, z in [0.8, 3.3]). Its wall faces are the 25
faces of those blocks that face the rest of the blocking, less those on
z = 0.

<Figure name="gui-blocking" caption="The Blocking dock with a face selection on the skeleton." />

1. In the Blocking dock set *select* to `face`. A purple marker sits at the
   centre of every face, drawn over everything.
2. Double click the first wall face's marker, then **shift**-double click
   the other 24. Shift adds a marker, or removes it if it is already
   selected; *selection* lists them as `b37.k_hi, ...`.
3. Select the geometry: double click the model. On a CAD model a double
   click picks a face (`wb: surface 3` in the Association dock); shift adds
   the other faces of the model, all 8 of them.
4. **Associate with geometry** runs one `p.associate_blocking(...)` per
   selected face.

<Figure name="gui-association" caption="The Association dock with a geometry entity selected." />

::: tip The body or its faces
The headless tutorial associates the walls with the body as a whole
(`kind='body'`). A double click on a CAD model selects faces, not the body,
so step 3 associates with all 8 faces (`'surface', (0, ..., 7)`): the
premesh projects onto the nearest of them, and on this model gives the same
quality figures as the body (checked headlessly). To get exactly the
headless journal, use the **Associate blocking entity** form with *kind*
`body` for each face instead. The spec sets the dock's geometry to the
body directly.
:::

## 5. Associate the wing root with the junction

1. Set *select* to `edge` and **Clear selection**.
2. Double click the 4 edges at the bottom of the wing's block (z = 0.8,
   round the wing root), shift for the 2nd to 4th.
3. Double click the junction curves on the model: the CAD edges where the
   wing meets the fuselage (upper, lower and the trailing-edge base). Shift
   adds them; the Association dock shows `wb: curve 6, 7, 8`.
4. **Associate with geometry.**

These 4 associations put the wall grid's corner on the junction; without
them it cuts the corner (see the
[headless tutorial](/tutorials/blocking#_5-associate-the-walls-and-the-junction)).

## 6. Delete the hole and insert the O-grid

1. Set *select* to `block`, **Clear selection**, and select the 10 hole
   blocks' markers (shift adds). **Delete blocks**: one
   `p.delete_blocks(...)` line; 65 blocks remain.
2. Set *select* to `face`, **Clear selection**, and select the 25 faces
   round the hole: they are now faces of the blocks outside it, and keep
   their associations.
3. Set *O-grid factor* to `0.7` and *O-grid nodes* to `17`, then
   **O-grid**. With no block selected the set is every block, and the
   selected faces are the layered ones. The blocking has 90 blocks.

## 7. Premesh, orthogonal lines, wall layer

Three forms, each **Run** in turn (they run on a worker thread with a
progress bar):

| Form | Field | Value |
|---|---|---|
| **Premesh blocking** | *exact* | checked |
| **Orthogonal wall lines** | *strength* | `0.5` |
| **Wall layer** | *first_step* | `2e-3` |

The blocks appear over the skeleton. **Quality** in the toolbar colours each
boundary quad by the scaled Jacobian of its cell (red ≤ 0, white, blue = 1);
the Display dock adds slices through any block. The headless tutorial
explains each choice: exact projection, strength 0.5 at the concave
junction, the 2e-3 first cell.

## 8. Read the journal and save

The Journal tab now holds the whole session, the same kind of script the
headless tutorial writes:

```python
p.associate_blocking('b37', 'i_hi_k_lo', 'wb', 'curve', (6, 7, 8))
p.delete_blocks(('b6', 'b11', 'b16', 'b7', 'b12', 'b17', 'b8', 'b13', 'b18', 'b37'))
p.ogrid((), faces=(('b5', 'i_hi'), ...), factor=0.7, nodes=17)
p.premesh(exact=True)
p.orthogonalize(strength=0.5)
p.wall_layer(0.002)
```

Undo and Redo (toolbar, or the usual shortcuts) walk this history. Display
settings (visibility, slices, Quality) are not commands and are not
journaled. **Save** writes a project directory (`project.toml`,
`journal.py`, `inputs/` with a copy of `wb.step`, the grids); **Write CGNS**
and **Write Xall run** forms write the solver files.

## 9. Replay headless

The saved project replays without the GUI:

```bash
shore run proj/journal.py --save-project proj2
shore cgns proj -o wing_body.cgns
```

or `Project.load("proj")` in Python, which replays the journal and requires
every grid to equal the saved one bit for bit. The GUI spec checks exactly
this (case I1), and also replays the session's journal on the edited model
(span 2.7, radius 0.55) to a valid grid on the new body (case I2), as step
10 of the headless tutorial does.

## What you learned

- The Command dock's forms run any project command; long ones run on a
  worker thread.
- The Blocking dock selects skeleton edges, faces and blocks (double click,
  shift adds) and splits, O-grids, deletes and associates them.
- The Association dock holds the geometry entity to associate with: faces
  and curves of a CAD model.
- The Journal tab is the session as a script; a saved session replays
  headless, bit for bit.

## Next

- [GUI reference](/reference/gui): every dock and control.
- [Blocking](/concepts/blocking) and [Association](/concepts/association).
- [Projects and replay](/tutorials/projects): saving, loading and replaying
  in depth.
