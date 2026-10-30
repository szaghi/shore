---
title: shore.gui
outline: [2, 3]
---

# `shore.gui`

SHORE's desktop GUI (ADR 0006, 0037): `shore gui` / `shore-gui`.

Needs the `gui` extra (PySide6-Essentials, pygfx, wgpu, rendercanvas).
`shore.gui.session` and `shore.gui.forms` are Qt-free; the window lives
in `shore.gui.window` and the renderer behind `shore.gui.viewport`.

**Modules**

- [`shore.gui.forms`](./shore.gui.forms) — Command forms of the GUI
- [`shore.gui.session`](./shore.gui.session) — The GUI's owner of a project
- [`shore.gui.viewport`](./shore.gui.viewport) — The GUI's renderer interface
- [`shore.gui.window`](./shore.gui.window) — SHORE's main window

## `main` <Badge type="info" text="function" /> {#main}

```python
main(argv: list[str] | None = None) -> int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/__init__.py#L13">source</a>

`shore-gui [PROJECT_DIR] [--no-replay] [--smoke]`.
