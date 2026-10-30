"""Screenshots of the SHORE GUI (``shore.gui.window.MainWindow``), off-screen.

The window runs under ``QT_QPA_PLATFORM=offscreen`` with an off-screen
rendercanvas for the 3-D view. The window is grabbed with ``QWidget.grab()``
and the view's own rendering (``PygfxViewport.draw_now``) is painted into the
central widget's rectangle, so the picture is the window as a user sees it.
The phase-7 wing-body blocking (ADR 0051) is built through the window as
``tests_gui/spec/test_gui_wing_body.py`` builds it. Needs ``shore[gui,cad]``.
"""

from __future__ import annotations

import functools
import os
import tempfile
from pathlib import Path

import numpy as np
from theme import OUT, webp

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")
HOLE = [(x, y, 0.4) for x in (-1.5, 0.5, 2.0) for y in (-0.55, 0.0, 0.55)] + [(0.5, 0.0, 2.0)]
BASE_EDGES = ("j_lo_k_lo", "j_hi_k_lo", "i_lo_k_lo", "i_hi_k_lo")
SIZE = (1400, 900)
# the 3-D view's background gradient (top, bottom) per mode; light is the GUI's own
BACKGROUND = {"light": ("#f4f5f7", "#d9dce1"), "dark": ("#30363d", "#1c2025")}


def _app():
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    app.setStyle("Fusion")
    return app


def _palette(mode: str):
    """Fusion's own palette (light) or a dark one."""
    from PySide6.QtGui import QColor, QPalette

    app = _app()
    if mode == "light":
        return app.style().standardPalette()
    p = QPalette()
    roles = QPalette.ColorRole
    for role, colour in (
        (roles.Window, "#2b2f35"),
        (roles.WindowText, "#e3e8ee"),
        (roles.Base, "#1f2328"),
        (roles.AlternateBase, "#2b2f35"),
        (roles.ToolTipBase, "#2b2f35"),
        (roles.ToolTipText, "#e3e8ee"),
        (roles.Text, "#e3e8ee"),
        (roles.Button, "#353a41"),
        (roles.ButtonText, "#e3e8ee"),
        (roles.BrightText, "#ff6b6b"),
        (roles.Link, "#4aa3c2"),
        (roles.Highlight, "#4aa3c2"),
        (roles.HighlightedText, "#101418"),
        (roles.PlaceholderText, "#7d8894"),
        (roles.Mid, "#3a4048"),
        (roles.Dark, "#181b1f"),
        (roles.Light, "#4a515a"),
        (roles.Shadow, "#0d0f12"),
    ):
        p.setColor(role, QColor(colour))
    for role in (roles.WindowText, roles.Text, roles.ButtonText):
        p.setColor(QPalette.ColorGroup.Disabled, role, QColor("#7d8894"))
    return p


def _block_at(blocking, point) -> str:
    for label in blocking.labels():
        pts = blocking.points[list(blocking.block(label))]
        if np.all(pts.min(axis=0) < point) and np.all(point < pts.max(axis=0)):
            return label
    raise AssertionError(f"no block holds {point}")


class Driver:
    """One window and the wing-body blocking built through it."""

    def __init__(self) -> None:
        from PySide6.QtCore import Qt
        from PySide6.QtWidgets import QDockWidget
        from renders_cad import wing_body

        from shore.gui.window import MainWindow

        _app()
        self.model = wing_body(Path(tempfile.mkdtemp(prefix="shore-gui-")) / "wb.step")
        w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
        self.w = w
        # the right-hand docks as tabs, so each panel has its full height

        docks = {d.objectName(): d for d in w.findChildren(QDockWidget)}
        right = [docks[n] for n in ("command", "blocking", "association", "display", "spacing")]
        for d in right[1:]:
            w.tabifyDockWidget(right[0], d)
        w.resizeDocks([docks["model"]], [300], Qt.Orientation.Horizontal)
        w.resizeDocks([right[0]], [330], Qt.Orientation.Horizontal)
        w.resizeDocks([docks["console"]], [150], Qt.Orientation.Vertical)
        self.docks = docks
        w.resize(*SIZE)
        w.show()
        self._events()
        self.walls: dict = {}
        self.hole: list[str] = []
        self.junction: list[int] = []
        self._build()

    # -- the build of tests_gui/spec/test_gui_wing_body.py ---------------------------------

    def _ok(self) -> None:
        self.w.wait_idle()
        if self.w.last_error is not None:
            raise RuntimeError(self.w.log.toPlainText())

    def _build(self) -> None:
        from shore.geometry.base import EntityKind
        from shore.gui.session import SkeletonKey

        w = self.w
        here = Path.cwd()
        os.chdir(self.model.parent)  # the journal then shows 'wb.step', not a temp path
        try:
            w.run_form("load_geometry", {"path": self.model.name, "name": "wb"})
            self._ok()
        finally:
            os.chdir(here)
        w.run_form("create_blocking", {"ni": "53", "nj": "33", "nk": "25",
                                       "vmin": "(-6, -4, 0)", "vmax": "(7, 4, 6)"})  # fmt: skip
        self._ok()
        corner = np.array([-5.99, -3.99, 0.01])
        for axis, cuts in (("i", (-3.0, -0.3, 1.3, 3.0)), ("j", (-0.8, -0.3, 0.3, 0.8)),
                           ("k", (0.8, 3.3))):  # fmt: skip
            dim = "ijk".index(axis)
            for c in cuts:
                point = corner.copy()
                point[dim] = c
                blocking = w.session.project.blocking
                label = _block_at(blocking, point)
                pts = blocking.points[list(blocking.block(label))]
                a, z = pts.min(axis=0)[dim], pts.max(axis=0)[dim]
                fraction = repr(float((c - a) / (z - a)))
                w.run_form("split_blocking", {"block": label, "axis": axis, "fraction": fraction})
                self._ok()
        blocking = w.session.project.blocking
        hole = [_block_at(blocking, np.array(c)) for c in HOLE]
        owners: dict = {}
        for label in hole:
            for f in FACES:
                owners.setdefault(blocking.entity(label, f)[0], []).append((label, f))
        on_plane = {k for k in owners if np.allclose(blocking.points[list(k[1:])][:, 2], 0.0)}
        walls = {k: v[0] for k, v in owners.items() if len(v) == 1 and k not in on_plane}
        w.blocking_mode.setCurrentText("face")
        w.clear_blocking_selection()
        for label, f in walls.values():
            w.select_skeleton(SkeletonKey("main", label, f), add=True)
        w.selected_geometry = ("wb", "body", [])
        w.associate_blocking_selection()
        self._ok()
        g = w.session.project.geometries["wb"].geometry
        curves = g.display().curves
        junction = [k for k in range(len(g.entities(EntityKind.CURVE)))
                    if np.allclose(np.hypot(curves[k][:, 1], curves[k][:, 2]), 0.5, atol=1e-3)
                    and np.all(np.abs(curves[k][:, 1]) < 0.07)]  # fmt: skip
        w.blocking_mode.setCurrentText("edge")
        w.clear_blocking_selection()
        for e in BASE_EDGES:
            w.select_skeleton(SkeletonKey("main", hole[-1], e), add=True)
        w.selected_geometry = ("wb", "curve", junction)
        w.associate_blocking_selection()
        self._ok()
        w.blocking_mode.setCurrentText("block")
        w.clear_blocking_selection()
        for label in hole:
            w.select_skeleton(SkeletonKey("main", label, ""), add=True)
        w.delete_selection()
        self._ok()
        w.blocking_mode.setCurrentText("face")
        w.clear_blocking_selection()
        blocking = w.session.project.blocking
        for k, _ in w.session.skeleton_markers("face"):
            if blocking.entity(k.block, k.edge)[0] in walls:
                w.select_skeleton(k, add=True)
        w.ogrid_factor.setText("0.7")
        w.ogrid_nodes.setText("17")
        w._ogrid_from_panel()
        self._ok()
        w.run_form("premesh", {"exact": True})
        self._ok()
        w.run_form("orthogonalize", {"strength": "0.5"})
        self._ok()
        w.run_form("wall_layer", {"first_step": "2e-3"})
        self._ok()
        w.clear_blocking_selection()
        self.walls, self.hole, self.junction = walls, hole, junction

    # -- view and screenshot ------------------------------------------------------------

    def _events(self) -> None:
        app = _app()
        for _ in range(5):
            app.processEvents()

    def mode(self, mode: str) -> None:
        import pygfx as gfx

        _app().setPalette(_palette(mode))
        scene = self.w.viewport.scene
        for child in list(scene.children):
            if isinstance(child, gfx.Background):
                scene.remove(child)
        scene.add(gfx.Background.from_color(*BACKGROUND[mode]))
        self._events()

    def look(self, focal, distance: float, view_dir) -> None:
        self.w.viewport.look_at_point(np.asarray(focal, float), distance, np.asarray(view_dir))

    def shot(self, name: str, mode: str) -> Path:
        from PySide6.QtGui import QImage, QPainter, QTextCursor

        w = self.w
        w.console.setCurrentWidget(w.journal)
        w.journal.moveCursor(QTextCursor.MoveOperation.End)
        self._events()
        rect = w.centralWidget().geometry()
        canvas = w.canvas
        canvas.set_logical_size(rect.width(), rect.height())
        w.viewport.camera.aspect = rect.width() / rect.height()
        img = np.ascontiguousarray(w.viewport.draw_now()[..., :4]).astype(np.uint8)
        h, wd = img.shape[:2]
        view = QImage(img.data, wd, h, 4 * wd, QImage.Format.Format_RGBA8888).copy()
        frame = w.grab().toImage()
        painter = QPainter(frame)
        painter.drawImage(rect, view)
        painter.end()
        out = OUT / mode / f"{name}.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        frame.save(str(out))
        out = webp(out)
        return out

    def raise_dock(self, name: str) -> None:
        self.docks[name].raise_()
        self._events()


@functools.cache
def driver() -> Driver:
    return Driver()


def _wall_slices(d: Driver) -> None:
    """Show, as Display-dock slices, every wall face of the grid and the
    symmetry-plane face of the blocks on the wall."""
    from shore.gui.session import BlockKey

    w = d.w
    if w.slices:
        return
    for b in w.session.project.mesh.blocks:
        walls = [f.name for f in b.faces() if f.bc.name == "WALL"]
        for name in walls:
            axis = name[0]
            index = 0 if name.endswith("lo") else b.nodes.shape["ijk".index(axis)] - 1
            w.add_slice(BlockKey("main", b.label), axis, index)
        if walls and np.allclose(b.nodes[:, :, 0, 2], 0.0, atol=0.05):
            w.add_slice(BlockKey("main", b.label), "k", 0)


def gui_window(mode: str) -> None:
    """The final grid: wall and symmetry-plane slices coloured by the scaled Jacobian."""
    d = driver()
    w = d.w
    d.mode(mode)
    w.clear_blocking_selection()
    w.viewport._highlight.visible = False
    w.set_visible("geometry", False)
    w.set_visible("skeleton", False)
    w.set_visible("block", False)
    _wall_slices(d)
    w.quality_action.setChecked(True)
    w.refresh_quality()
    w.command.setCurrentIndex(w.command.findData("wall_layer"))
    w._editors["first_step"].setText("2e-3")
    d.raise_dock("display")
    d.look((0.4, 0.0, 0.9), 7.0, (0.45, 0.8, -0.45))
    d.shot("gui-window", mode)
    w.quality_action.setChecked(False)
    w.refresh_quality()


def gui_blocking(mode: str) -> None:
    """The blocking before the O-grid, the hole's faces selected in the dock."""
    from shore.gui.session import SkeletonKey

    d = driver()
    w = d.w
    d.mode(mode)
    s = w.session
    for _ in range(4):  # wall_layer, orthogonalize, premesh, ogrid
        s.undo()
    try:
        w.set_visible("slice", False)
        w.set_visible("geometry", True)
        w.set_visible("skeleton", True)
        w.blocking_mode.setCurrentText("face")
        w.clear_blocking_selection()
        blocking = s.project.blocking
        outlines = []
        for k, _ in s.skeleton_markers("face"):
            key = blocking.entity(k.block, k.edge)[0]
            if key in d.walls:
                w.select_skeleton(SkeletonKey("main", k.block, k.edge), add=True)
                pts = blocking.points[[key[1], key[2], key[4], key[3], key[1]]]
                outlines.append(np.vstack([pts, np.full((1, 3), np.nan)]))
        w.ogrid_factor.setText("0.7")
        w.ogrid_nodes.setText("17")
        d.raise_dock("blocking")
        d.look((0.4, 0.0, 1.2), 11.5, (0.45, 0.8, -0.45))
        w.viewport.highlight_polyline(np.concatenate(outlines))
        d.shot("gui-blocking", mode)
    finally:
        w.viewport._highlight.visible = False
        w.clear_blocking_selection()
        for _ in range(4):
            s.redo()
        w.set_visible("slice", True)


def gui_association(mode: str) -> None:
    """A wall face of the premeshed grid and a CAD face selected in the Association dock."""
    from shore.gui.session import BlockKey

    d = driver()
    w = d.w
    d.mode(mode)
    w.clear_blocking_selection()
    w.set_visible("skeleton", False)
    w.set_visible("block", False)
    w.set_visible("slice", False)
    w.set_visible("geometry", True)
    w.selection_mode.setCurrentText("face")
    project = w.session.project
    # the wall block on the wing's upper surface near the root
    target = None
    for b in project.mesh.blocks:
        for f in b.faces():
            if f.bc.name != "WALL":
                continue
            ax = "ijk".index(f.name[0])
            g = np.take(b.nodes, 0 if f.name.endswith("lo") else -1, axis=ax)
            c = g.reshape(-1, 3).mean(axis=0)
            if c[1] > 0.02 and 1.0 < c[2] < 3.0 and 0.0 < c[0] < 1.0:
                target = (b, f.name, g, c)
    b, face, g, c = target
    w._select_block(BlockKey("main", b.label), face, c)
    geometry = project.geometries["wb"].geometry
    proj = geometry.project(c[None], _body(), exact=True)
    surface = _surface_of(geometry, proj.points[0])
    w._select_geometry("wb", "surface", surface, False)
    outline = np.vstack([g[0], g[1:, -1], g[-1, ::-1][1:], g[::-1, 0][1:]])
    d.raise_dock("association")
    d.look((0.6, 0.0, 1.3), 7.0, (0.5, 0.75, -0.45))
    w.viewport.highlight_polyline(outline)
    d.shot("gui-association", mode)
    w.viewport._highlight.visible = False
    w.set_visible("slice", True)


def _body():
    from shore.geometry.base import BODY

    return BODY


def _surface_of(geometry, point: np.ndarray) -> int:
    """The CAD face nearest *point*."""
    from shore.geometry.base import EntityKind, EntityRef

    d = [
        float(geometry.project(point[None], EntityRef(EntityKind.SURFACE, k)).distances[0])
        for k in range(len(geometry.entities(EntityKind.SURFACE)))
    ]
    return int(np.argmin(d))


FIGURES = {
    "gui-window": gui_window,
    "gui-blocking": gui_blocking,
    "gui-association": gui_association,
}
