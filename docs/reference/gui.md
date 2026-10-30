# shore gui

The SHORE desktop GUI ([ADR 0006](../adr/0006-gui-stack.md),
[ADR 0037](../adr/0037-gui-viewer-command-shell.md)): a viewer and command
shell over a project. Every action is a project command, shown as a line of
the journal, so a session saved from the GUI replays headless
(`shore run`, `Project.load`).

Needs the `gui` extra: `pip install 'shore-mesh[gui]'` (PySide6-Essentials,
pygfx, wgpu, rendercanvas; about 90 MB). Rendering uses the GPU through
Metal, Direct3D 12 or Vulkan.

::: warning WSL2
Under WSL2, wgpu finds only a CPU adapter, and orbiting a 10⁷-node mesh
takes about 160 ms per frame. Run the GUI with native Windows Python.
:::

## Synopsis

```
shore gui [OPTIONS] [PROJECT_DIR]
shore-gui [PROJECT_DIR] [--no-replay]
```

| Argument / option | Description |
|---|---|
| `PROJECT_DIR` | A saved project to open (`Project.save`, `--project`). |
| `--no-replay` | Open it from its saved grids without replaying its journal ([ADR 0033](../adr/0033-load-without-replay.md)). |
| `--smoke` | Build the window, show it briefly, print `shore gui smoke: N blocks, journal M lines` and exit. |

## The window

| Area | Content |
|---|---|
| Left: **Model** / **Geometry** | Model: components (the current one marked), their blocks (label and `ni x nj x nk`), each block's six faces with their BC and association. Geometry: each loaded geometry, with its triangles (STL) or its faces, edges and vertices (CAD). |
| Centre: viewport | Every block of every component (its six boundary faces and twelve edges), every geometry (triangles; CAD edges and vertices drawn over the faces) and the slice planes. Orbit with the left button, pan with the right, zoom with the wheel. **Double click** picks: on a block or slice, the status bar shows `component/block`, the face (or `k=2` for a slice), the cell `(i, j, k)` and the float64 centre `xyz`; on a geometry, its name, the entity (`surface 3`, `curve 7`, `point 2`, or the STL `body`) and the picked point. |
| Right: **Command** | A form per command: load geometry, use component, build topology, build primitive, march, snap to geometry, the blocking commands (create, split, O-grid, delete blocks, move a vertex, associate an entity, premesh), orthogonal wall lines, write the Xall run, CGNS or Plot3D. **Run** executes it. Builds, the march and the Xall run go on a worker thread, with a progress bar; **Cancel** stops a march (the mesh is left as before). |
| Right: **Display** | Show or hide blocks, geometries, slices and the blocking. Slices: pick a block, an axis and a node index, **Add slice**; the list shows them, **Remove slice** drops the selected one. Display settings are not commands: they do not enter the journal and are not saved ([ADR 0038](../adr/0038-gui-geometry-and-slices.md)). |
| Right: **Association** | **Select block** (face / edge / vertex): a double click on a block selects the picked face, its nearest edge or its nearest corner. A double click on a geometry selects its entity; with shift it adds another of the same geometry and kind (several surfaces, a curve chain). **Associate** links the block entity with the geometry entities (`p.associate(...)` in the journal; the block's component is made current first if needed); **Clear association** removes it. Then run **Snap to geometry** to move the nodes. The model tree's third column shows each face's association, and associated edges and vertices are listed after the faces ([ADR 0039](../adr/0039-association-editing.md)). |
| Right: **Spacing** | Acts on an edge selected in the Association dock's *edge* mode: shows its axis, node count, length, and first and last cell. Pick a law (uniform, tanh, tanh2, pinned) and its parameters; **Redistribute** runs `p.redistribute(...)` ([ADR 0041](../adr/0041-redistribute-edge.md)). |
| Right: **Blocking** | The blocking's skeleton is drawn as orange lines over the blocks, each edge with its nodes ([ADR 0042](../adr/0042-blocking-skeleton-premesh.md)). **Select** edge / face / block / vertex: in face, block and vertex mode a purple marker sits at each face's or block's centre or at each vertex, drawn over everything ([ADR 0044](../adr/0044-blocking-selection.md)). A double click on an edge or marker selects it; with shift it is added, or removed if already selected; the selection survives mode changes. **Split here** cuts through the picked edge node (`p.split_blocking`). **O-grid** (factor, nodes) runs `p.ogrid` with the selected blocks as the set (all blocks when none is selected) and the selected faces as the layered faces: select a hole's faces for a layer around a body. **Delete blocks** deletes the selected blocks; in vertex mode, **Move vertex** moves the one selected vertex to the typed `x, y, z` and **Move vertex to geometry pick** to the last double-clicked geometry point (`p.move_vertex`, onto its association if it has one, [ADR 0045](../adr/0045-move-blocking-vertex.md)); **Associate with geometry** associates each selected face or edge with the geometry selected in the Association dock. The Display dock's *Blocking* box hides skeleton and markers. |
| Bottom: **Console** | **Journal**: the project's journal, one line per command. **Log**: what ran, and errors. |
| Toolbar | New, Open…, Save, Save as…, Undo, Redo, Fit, **Quality** (colours each boundary face and slice quad by the scaled Jacobian of its cell, red ≤ 0 → white → blue = 1). |

Structured arguments are typed as Python literals. For example, *Build
topology* takes `{'ni': 33, 'nj': 64}` for the ray-cast topologies, and
`{'ni_body': 60, 'ni_wake': 15, 'n_stations': 4, 'wake_length': 10.0}` for
the C-grid. *Build primitive* takes the primitive's parameters.

## Precision

Positions stay in float64 on the CPU. The GPU receives them relative to a
floating origin that follows the view, and they are re-uploaded when the
view moves far ("rebase"). Wall cells of 10⁻⁶ at 10³ from the origin
therefore display and pick exactly, and picked coordinates come from the
float64 nodes.

## Python

```python
from shore.gui.session import Session          # Qt-free: Session.run("march"), undo, open, save
from shore.gui.window import MainWindow, main  # the window; main(project=None, replay=True)
```
