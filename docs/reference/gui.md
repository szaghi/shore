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
| Left: **Model** | Components (the current one marked), their blocks (label and `ni x nj x nk`), each block's six faces with their BC. |
| Centre: viewport | Every block of every component: its six boundary faces and twelve edges. Orbit with the left button, pan with the right, zoom with the wheel. **Double click** picks: the status bar shows `component/block`, the face, the cell `(i, j, k)` and the float64 centre `xyz`. |
| Right: **Command** | A form per command: load geometry, use component, build topology, build primitive, march, snap to geometry, write the Xall run, CGNS or Plot3D. **Run** executes it. Builds, the march and the Xall run go on a worker thread, with a progress bar; **Cancel** stops a march (the mesh is left as before). |
| Bottom: **Console** | **Journal**: the project's journal, one line per command. **Log**: what ran, and errors. |
| Toolbar | New, Open…, Save, Save as…, Undo, Redo, Fit, **Quality** (colours each boundary face by the scaled Jacobian of its cell, red ≤ 0 → white → blue = 1). |

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
