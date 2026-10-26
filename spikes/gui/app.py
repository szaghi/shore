"""Interactive GUI spike (ADR 0006, items 1-4 on screen).

A PySide6 main window with an ICEM/Pointwise-style layout: blocks tree (left),
pygfx viewport (centre), properties (right), console (bottom). Toolbar actions
load a synthetic assembly, run a 5 s orbit benchmark at uncapped frame rate,
and run a 30 s numpy "march" in a worker thread while UI latency is sampled.
Double-click picks a boundary cell. "Zoom to wall" animates a user-style zoom
onto a wall cell; the viewport's rebase trigger keeps it exact, and every
rebase is logged to the console.

The default assembly is the hard case: block extent ~285, ``ds_wall / extent =
1e-9``, wall-cell aspect ratio 1e4, at a 1e3 world offset (one origin per block
fails there; see README).

Usage::

    python app.py --nodes 10000000 [--no-vsync] [--origin-mode rebase]
    QT_QPA_PLATFORM=offscreen python app.py --smoke   # CI / headless check
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np
from interaction import wall_cell, zoom_to
from latency import HEARTBEAT_MS, Worker
from PySide6.QtCore import Qt, QThreadPool, QTimer
from PySide6.QtGui import QAction, QColor, QFont, QGuiApplication, QPalette
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDockWidget,
    QFormLayout,
    QLabel,
    QMainWindow,
    QPlainTextEdit,
    QTreeWidget,
    QTreeWidgetItem,
    QWidget,
)
from rendercanvas.qt import QRenderWidget  # must follow the PySide6 import
from synthetic import make_assembly
from viewport import PygfxViewport

NODE_CHOICES = (100_000, 1_000_000, 10_000_000, 20_000_000)
ORBIT_SECONDS = 5.0
MARCH_SECONDS = 30.0


def _palette(dark: bool) -> QPalette:
    """A restrained Fusion palette for dark or light mode."""
    p = QPalette()
    if dark:
        base, alt, text, window, button = "#1e1f22", "#26282c", "#e6e6e6", "#2b2d31", "#33363b"
    else:
        base, alt, text, window, button = "#ffffff", "#f3f4f6", "#1d1d1f", "#eceef1", "#f7f7f8"
    accent = "#3d7eff"
    for role, color in (
        (QPalette.ColorRole.Window, window),
        (QPalette.ColorRole.WindowText, text),
        (QPalette.ColorRole.Base, base),
        (QPalette.ColorRole.AlternateBase, alt),
        (QPalette.ColorRole.Text, text),
        (QPalette.ColorRole.Button, button),
        (QPalette.ColorRole.ButtonText, text),
        (QPalette.ColorRole.Highlight, accent),
        (QPalette.ColorRole.HighlightedText, "#ffffff"),
        (QPalette.ColorRole.ToolTipBase, base),
        (QPalette.ColorRole.ToolTipText, text),
    ):
        p.setColor(role, QColor(color))
    return p


class MainWindow(QMainWindow):
    """Spike main window."""

    def __init__(self, args: argparse.Namespace) -> None:
        super().__init__()
        self.args = args
        self.setWindowTitle("SHORE GUI spike")
        self.resize(1600, 960)

        self.canvas = QRenderWidget(self, update_mode="ondemand", vsync=not args.no_vsync)
        self.setCentralWidget(self.canvas)
        self.vp = PygfxViewport(
            self.canvas, origin_mode=args.origin_mode, rebase_factor=args.rebase_factor
        )
        self.vp.on_rebase = self._on_rebase
        self.blocks: list = []
        self.vp.renderer.add_event_handler(self._on_double_click, "double_click")

        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Block", "ni x nj x nk"])
        self._dock("Model", self.tree, Qt.DockWidgetArea.LeftDockWidgetArea)

        props = QWidget()
        form = QFormLayout(props)
        self.prop = {k: QLabel("-") for k in ("block", "face", "cell", "xyz", "pick ms")}
        for k, w in self.prop.items():
            w.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            form.addRow(k, w)
        self._dock("Properties", props, Qt.DockWidgetArea.RightDockWidgetArea)

        self.console = QPlainTextEdit(readOnly=True)
        self.console.setFont(QFont("monospace"))
        self._dock("Console", self.console, Qt.DockWidgetArea.BottomDockWidgetArea)

        tb = self.addToolBar("Main")
        self.nodes_box = QComboBox()
        for n in NODE_CHOICES:
            self.nodes_box.addItem(f"{n:.0e} nodes", n)
        self.nodes_box.setCurrentIndex(
            min(range(len(NODE_CHOICES)), key=lambda i: abs(NODE_CHOICES[i] - args.nodes))
        )
        tb.addWidget(self.nodes_box)
        for text, slot in (
            ("Load", self.load),
            ("Fit", self.vp.show_all),
            ("Orbit bench", self.orbit_bench),
            ("Fake march", self.fake_march),
            ("Zoom to wall", self.zoom_to_wall),
        ):
            act = QAction(text, self)
            act.triggered.connect(slot)
            tb.addAction(act)

        self.fps_timer = QTimer(self, interval=1000)
        self.fps_timer.timeout.connect(self._update_status)
        self.fps_timer.start()

        hints = QGuiApplication.styleHints()
        hints.colorSchemeChanged.connect(lambda _s: self.apply_theme())
        self.apply_theme()

    # -- helpers -----------------------------------------------------------------

    def _dock(self, title: str, widget: QWidget, area: Qt.DockWidgetArea) -> None:
        dock = QDockWidget(title, self)
        dock.setWidget(widget)
        dock.setObjectName(title)
        self.addDockWidget(area, dock)

    def log(self, line: str) -> None:
        self.console.appendPlainText(line)
        print(line, flush=True)

    def apply_theme(self) -> None:
        dark = QGuiApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark
        QApplication.setPalette(_palette(dark))
        if dark:
            self.vp.set_background("#2b2e35", "#15161a")
        else:
            self.vp.set_background("#f4f5f7", "#d3d7de")

    def _update_status(self) -> None:
        s = self.vp.frame_stats()
        self.statusBar().showMessage(
            f"{self.vp.n_nodes:,} nodes | frames {s.n} | mean {s.mean_ms:.1f} ms"
            f" | p95 {s.p95_ms:.1f} ms | rebases {self.vp.rebase_count}"
            f" (last {self.vp.last_rebase_ms:.1f} ms)"
            f" | adapter {self.vp.renderer.device.adapter.summary}"
        )

    # -- actions -----------------------------------------------------------------

    def load(self) -> None:
        n = int(self.nodes_box.currentData())
        t0 = time.perf_counter()
        a = self.args
        blocks = make_assembly(
            n_nodes=n,
            n_phi=5,
            n_theta=4,
            radius=a.radius,
            thickness=2.0 * a.radius,
            ds_wall=a.ds_over_extent * 2.85 * a.radius,
            aspect_ratio=a.aspect_ratio,
        )
        self.blocks = blocks
        t_build = time.perf_counter() - t0
        self.vp.clear()
        self.tree.clear()
        t0 = time.perf_counter()
        for b in blocks:
            self.vp.add_block(b)
            QTreeWidgetItem(self.tree, [b.name, "x".join(map(str, b.nodes.shape[:3]))])
        self.vp.show_all()
        self.canvas.force_draw()
        self.log(
            f"load(nodes={self.vp.n_nodes:,}) build {t_build:.2f} s,"
            f" upload+first frame {time.perf_counter() - t0:.2f} s"
        )

    def orbit_bench(self) -> None:
        target = self.vp.target.copy()
        eye0 = np.asarray(self.vp.camera.local.position, dtype=np.float64) - target
        t_start = time.perf_counter()

        def rotate() -> None:
            ang = 2.0 * np.pi * (time.perf_counter() - t_start) / ORBIT_SECONDS
            c, s = np.cos(ang), np.sin(ang)
            rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
            self.vp.camera.local.position = target + rot @ eye0
            self.vp.camera.look_at(target)

        self.vp.reset_stats()
        self.vp.before_draw = rotate
        self.canvas.set_update_mode("continuous", max_fps=1000)

        def stop() -> None:
            self.canvas.set_update_mode("ondemand", max_fps=30)
            self.vp.before_draw = None
            s = self.vp.frame_stats()
            ok = s.p95_ms <= 16.0
            self.log(
                f"orbit_bench: {s.n} frames, mean {s.mean_ms:.2f} ms, p95 {s.p95_ms:.2f} ms,"
                f" max {s.max_ms:.2f} ms -> {'PASS' if ok else 'FAIL'} (budget 16 ms)"
            )

        QTimer.singleShot(int(ORBIT_SECONDS * 1000), stop)

    def fake_march(self, seconds: float = MARCH_SECONDS) -> None:
        late: list[float] = []
        last = [time.perf_counter()]

        def beat() -> None:
            now = time.perf_counter()
            late.append(max(0.0, 1e3 * (now - last[0]) - HEARTBEAT_MS))
            last[0] = now

        self._beat = QTimer(self, interval=HEARTBEAT_MS)
        self._beat.timeout.connect(beat)

        def done(n: int) -> None:
            self._beat.stop()
            a = np.asarray(late)
            ok = a.max() <= 100.0
            self.log(
                f"fake_march: {n} sweeps, UI lateness p99 {np.percentile(a, 99):.1f} ms,"
                f" max {a.max():.1f} ms -> {'PASS' if ok else 'FAIL'} (budget 100 ms)"
            )

        worker = Worker("numpy", seconds)
        worker.signals.done.connect(done)
        self._beat.start()
        QThreadPool.globalInstance().start(worker)
        self.log(f"fake_march(seconds={seconds}) started; orbit the view to feel latency")

    def zoom_to_wall(self, steps: int = 90, on_done: object = None) -> None:
        """Animate a user-style zoom onto a wall cell, then pick it."""
        if not self.blocks:
            self.log("zoom_to_wall: load an assembly first")
            return
        cell = wall_cell(self.blocks[len(self.blocks) // 2])
        start = self.vp.rebase_count
        frames = zoom_to(self.vp, cell, d_start=3.0 * self.args.radius, n_steps=steps)
        self.log(
            f"zoom_to_wall(block={cell.block!r}) cell height {cell.height:.3e}, {steps} frames"
        )
        self._zoom_timer = QTimer(self, interval=16)

        def step() -> None:
            if next(frames, None) is not None:
                self.canvas.force_draw()
                return
            self._zoom_timer.stop()
            res = self.vp.pick(*self.vp.project(cell.centre))
            err = None if res is None else float(np.linalg.norm(res.xyz - cell.centre))
            ok = err is not None and err < 0.5 * cell.height
            if res is not None:
                self._show_pick(res)
            self.log(
                f"zoom_to_wall: rebases {self.vp.rebase_count - start}, pick error"
                f" {err if err is None else f'{err:.2e}'} -> {'PASS' if ok else 'FAIL'}"
            )
            if callable(on_done):
                on_done()

        self._zoom_timer.timeout.connect(step)
        self._zoom_timer.start()

    def _on_rebase(self, info: dict) -> None:
        o = np.array2string(info["origin"], precision=6, separator=", ")
        self.log(
            f"rebase(origin={o}) drift {info['drift_views']:.0f} views,"
            f" view width {info['view_width']:.2e}, {info['ms']:.1f} ms"
        )

    def _show_pick(self, res: object) -> None:
        self.prop["block"].setText(res.block)
        self.prop["face"].setText(res.face)
        self.prop["cell"].setText(str(res.cell))
        self.prop["xyz"].setText(np.array2string(res.xyz, precision=9))
        self.prop["pick ms"].setText(f"{res.elapsed_ms:.1f}")

    def _on_double_click(self, event: object) -> None:
        res = self.vp.pick(event.x, event.y)
        if res is None:
            return
        self._show_pick(res)
        ok = "PASS" if res.elapsed_ms <= 50.0 else "FAIL"
        self.log(f"pick(block={res.block!r}, cell={res.cell}) {res.elapsed_ms:.1f} ms -> {ok}")


def _smoke(win: MainWindow, app: QApplication) -> None:
    """Load, orbit, pick, short march, zoom to wall, quit: every code path."""
    win.nodes_box.setCurrentIndex(0)
    win.load()
    win.orbit_bench()

    def pick_centre() -> None:
        w, h = win.vp.renderer.logical_size
        res = win.vp.pick(w / 2, h / 2)
        win.log(f"smoke pick at centre: {None if res is None else (res.block, res.cell)}")
        win.fake_march(seconds=3.0)
        QTimer.singleShot(
            4000, lambda: win.zoom_to_wall(on_done=lambda: QTimer.singleShot(500, app.quit))
        )

    QTimer.singleShot(int(ORBIT_SECONDS * 1000) + 500, pick_centre)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--nodes", type=int, default=10_000_000)
    ap.add_argument(
        "--origin-mode",
        choices=("absolute", "block", "floating", "tiled", "rebase"),
        default="rebase",
    )
    ap.add_argument("--rebase-factor", type=float, default=100.0, help="trigger, view widths")
    ap.add_argument("--radius", type=float, default=100.0, help="assembly wall radius")
    ap.add_argument("--ds-over-extent", type=float, default=1e-9, help="wall spacing / extent")
    ap.add_argument("--aspect-ratio", type=float, default=1e4, help="wall-cell aspect ratio")
    ap.add_argument("--no-vsync", action="store_true", help="uncap frame rate for benchmarks")
    ap.add_argument("--smoke", action="store_true", help="scripted run, then quit")
    args = ap.parse_args()

    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow(args)
    win.show()
    if args.smoke:
        QTimer.singleShot(0, lambda: _smoke(win, app))
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
