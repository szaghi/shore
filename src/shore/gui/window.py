"""SHORE's main window (ADR 0006, 0037): model tree, viewport, command forms, console.

Every action goes through :class:`~shore.gui.session.Session`, so it is a
Project command and a journal line (ADR 0004). Long commands run on a worker
thread with the ADR 0004 progress and cancel hooks; the window refreshes from
the session when they finish.

Layout (ADR 0006): model tree on the left, the viewport in the centre,
command forms on the right, the journal and the log at the bottom.
"""

from __future__ import annotations

import time
import traceback
from collections.abc import Callable
from typing import Any

from PySide6.QtCore import QObject, QRunnable, Qt, QThreadPool, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QTabWidget,
    QToolBar,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from shore.gui.forms import FORMS, Field, Form, form
from shore.gui.session import BlockKey, Session
from shore.gui.viewport import PickResult, boundary_values
from shore.project.commands import CancelToken

AskPath = Callable[[str, str], str | None]  # (kind: "file" | "dir" | "save", title) -> path


def _ask_path(kind: str, title: str) -> str | None:
    if kind == "dir":
        path = QFileDialog.getExistingDirectory(None, title)
    elif kind == "save":
        path, _ = QFileDialog.getSaveFileName(None, title)
    else:
        path, _ = QFileDialog.getOpenFileName(None, title)
    return path or None


class _Signals(QObject):
    progress = Signal(float, str)
    done = Signal(object)  # None on success, else the exception


class _Task(QRunnable):
    def __init__(self, work: Callable[[Callable[[float, str], None]], None]) -> None:
        super().__init__()
        self.work = work
        self.signals = _Signals()

    def run(self) -> None:
        try:
            self.work(self.signals.progress.emit)
        except BaseException as exc:  # reported to the UI thread, not lost
            self.signals.done.emit(exc)
        else:
            self.signals.done.emit(None)


class MainWindow(QMainWindow):
    """The SHORE GUI over one :class:`Session`."""

    def __init__(
        self,
        session: Session | None = None,
        *,
        ask_path: AskPath = _ask_path,
        offscreen_viewport: bool = False,
    ) -> None:
        super().__init__()
        from shore.gui.viewport.pygfx_backend import PygfxViewport, create_canvas

        self.session = session if session is not None else Session()
        self.ask_path = ask_path
        self.setWindowTitle("SHORE")
        self.resize(1400, 900)
        self._pool = QThreadPool()
        self._pool.setMaxThreadCount(1)  # commands run one at a time, in order
        self._busy = False
        self._cancel: CancelToken | None = None
        self.last_error: BaseException | None = None
        self.last_pick: PickResult | None = None
        self._shown: list[BlockKey] = []

        self.canvas = create_canvas(self, offscreen=offscreen_viewport)
        self.viewport = PygfxViewport(self.canvas)
        self.viewport.on_pick(self._picked)
        if offscreen_viewport:
            self.setCentralWidget(QLabel("offscreen viewport"))
        else:
            self.setCentralWidget(self.canvas)

        self._build_tree()
        self._build_forms()
        self._build_console()
        self._build_actions()
        self._progress = QProgressBar()
        self._progress.setRange(0, 1000)
        self._progress.setVisible(False)
        self.statusBar().addPermanentWidget(self._progress)
        self.session.subscribe(self.refresh)
        self.refresh()

    # -- panels ----------------------------------------------------------------

    def _build_tree(self) -> None:
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Model", "Detail"])
        dock = QDockWidget("Model", self)
        dock.setObjectName("model")
        dock.setWidget(self.tree)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)

    def _build_forms(self) -> None:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.command = QComboBox()
        for f in FORMS:
            self.command.addItem(f.title, f.method)
        self.command.currentIndexChanged.connect(lambda _: self._show_form())
        layout.addWidget(self.command)
        self._form_host = QWidget()
        self._form_layout = QFormLayout(self._form_host)
        layout.addWidget(self._form_host)
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._run_current_form)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.cancel)
        self.cancel_button.setEnabled(False)
        layout.addWidget(self.run_button)
        layout.addWidget(self.cancel_button)
        layout.addStretch(1)
        dock = QDockWidget("Command", self)
        dock.setObjectName("command")
        dock.setWidget(panel)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        self._editors: dict[str, Any] = {}
        self._show_form()

    def _show_form(self) -> None:
        while self._form_layout.rowCount():
            self._form_layout.removeRow(0)
        self._editors = {}
        f = form(self.command.currentData())
        for field in f.fields:
            editor = self._editor(field)
            self._editors[field.name] = editor
            self._form_layout.addRow(field.label or field.name, editor)

    def _editor(self, field: Field) -> Any:
        if field.kind == "bool":
            box = QCheckBox()
            box.setChecked(bool(field.default))
            return box
        if field.kind == "choice":
            combo = QComboBox()
            combo.addItems(list(field.choices))
            combo.setCurrentText(str(field.default))
            return combo
        edit = QLineEdit(
            ""
            if field.default is None
            else repr(field.default)
            if field.kind == "literal"
            else str(field.default)
        )
        if field.kind in ("file", "dir"):
            kind = "dir" if field.kind == "dir" else "file"
            action = edit.addAction(
                self.style().standardIcon(self.style().StandardPixmap.SP_DirOpenIcon),
                QLineEdit.ActionPosition.TrailingPosition,
            )
            action.triggered.connect(lambda: self._browse(edit, kind, field.label or field.name))
        return edit

    def _browse(self, edit: QLineEdit, kind: str, title: str) -> None:
        path = self.ask_path(kind, title)
        if path:
            edit.setText(path)

    def _form_values(self) -> dict[str, Any]:
        values: dict[str, Any] = {}
        for name, editor in self._editors.items():
            if isinstance(editor, QCheckBox):
                values[name] = editor.isChecked()
            elif isinstance(editor, QComboBox):
                values[name] = editor.currentText()
            else:
                values[name] = editor.text()
        return values

    def _build_console(self) -> None:
        self.console = QTabWidget()
        self.journal = QPlainTextEdit()
        self.journal.setReadOnly(True)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.console.addTab(self.journal, "Journal")
        self.console.addTab(self.log, "Log")
        dock = QDockWidget("Console", self)
        dock.setObjectName("console")
        dock.setWidget(self.console)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, dock)

    def _build_actions(self) -> None:
        bar = QToolBar("Main")
        bar.setObjectName("main")
        self.addToolBar(bar)

        def action(text: str, slot: Callable[[], None], shortcut: Any = None) -> QAction:
            a = QAction(text, self)
            a.triggered.connect(slot)
            if shortcut is not None:
                a.setShortcut(QKeySequence(shortcut))
            bar.addAction(a)
            return a

        self.new_action = action("New", self.new_project, QKeySequence.StandardKey.New)
        self.open_action = action("Open…", self.open_project, QKeySequence.StandardKey.Open)
        self.save_action = action("Save", self.save_project, QKeySequence.StandardKey.Save)
        self.save_as_action = action("Save as…", lambda: self.save_project(ask=True))
        bar.addSeparator()
        self.undo_action = action("Undo", self.undo, QKeySequence.StandardKey.Undo)
        self.redo_action = action("Redo", self.redo, QKeySequence.StandardKey.Redo)
        bar.addSeparator()
        self.fit_action = action("Fit", self.viewport.show_all)
        self.quality_action = action("Quality", self.refresh_quality)
        self.quality_action.setCheckable(True)

    # -- refresh -------------------------------------------------------------------

    def refresh(self) -> None:
        """Bring the tree, the journal and the viewport in line with the session."""
        self.tree.clear()
        for name, current, blocks in self.session.tree():
            item = QTreeWidgetItem(
                [name + (" (current)" if current else ""), f"{len(blocks)} blocks"]
            )
            for label, shape, faces in blocks:
                child = QTreeWidgetItem([label, "x".join(map(str, shape))])
                for face, bc in faces.items():
                    child.addChild(QTreeWidgetItem([face, bc]))
                item.addChild(child)
            self.tree.addTopLevelItem(item)
            item.setExpanded(True)
        self.journal.setPlainText("\n".join(self.session.journal_lines()))
        blocks = self.session.blocks()
        keys = [key for key, _ in blocks]
        # node values change in place (a march fills the layers): always re-upload
        self.viewport.set_blocks([(key, b.nodes) for key, b in blocks])
        if keys != self._shown:  # frame the view when the set of blocks changes
            self._shown = keys
            self.viewport.show_all()
        if self.quality_action.isChecked():
            self.refresh_quality()
        self._update_actions()

    def refresh_quality(self) -> None:
        """Colour boundary quads by the scaled Jacobian of their cells, or clear."""
        if not self.quality_action.isChecked():
            self.viewport.set_face_scalars(None)
            return
        from shore.quality import scaled_jacobian

        values = {
            key: boundary_values(scaled_jacobian(b.nodes)) for key, b in self.session.blocks()
        }
        self.viewport.set_face_scalars(values, clim=(-1.0, 1.0))

    def _update_actions(self) -> None:
        idle = not self._busy
        for a in (
            self.new_action,
            self.open_action,
            self.save_action,
            self.save_as_action,
            self.undo_action,
            self.redo_action,
        ):
            a.setEnabled(idle)
        self.run_button.setEnabled(idle)
        self.cancel_button.setEnabled(self._busy and self._cancel is not None)

    # -- actions -------------------------------------------------------------------

    def _report(self, exc: BaseException | None) -> None:
        self.last_error = exc
        if exc is None:
            return
        self.log.appendPlainText(f"Error: {exc}")
        self.log.appendPlainText("".join(traceback.format_exception_only(type(exc), exc)).strip())
        self.statusBar().showMessage(f"Error: {exc}", 10000)
        self.console.setCurrentWidget(self.log)

    def _guard(self, work: Callable[[], None]) -> None:
        try:
            work()
        except Exception as exc:
            self._report(exc)
        else:
            self.last_error = None

    def run_form(self, method: str, values: dict[str, Any] | None = None) -> None:
        """Run the form of *method* with *values* (field name -> text or value)."""
        f = form(method)
        try:
            args, kwargs = f.arguments(values or {})
        except Exception as exc:
            self._report(exc)
            return
        self._submit(f, args, kwargs)

    def _run_current_form(self) -> None:
        self.run_form(self.command.currentData(), self._form_values())

    def _submit(self, f: Form, args: list[Any], kwargs: dict[str, Any]) -> None:
        if self._busy:
            self._report(RuntimeError("a command is already running"))
            return
        self.log.appendPlainText(f"run {f.method}")
        if not f.long:
            self._guard(lambda: self.session.run(f.method, *args, **kwargs))
            return
        self._busy = True
        self._cancel = CancelToken() if f.method == "march" else None
        cancel = self._cancel
        self._progress.setValue(0)
        self._progress.setVisible(True)
        self._update_actions()
        self._held_listeners = self.session._listeners
        self.session._listeners = []  # refresh once, on the UI thread, when done

        def work(progress: Callable[[float, str], None]) -> None:
            self.session.run(f.method, *args, progress=progress, cancel=cancel, **kwargs)

        task = _Task(work)
        # bound methods of this window: Qt queues the calls to the UI thread
        task.signals.progress.connect(self._on_progress)
        task.signals.done.connect(self._task_done)
        self._task = task  # keep the signals object alive until done
        self._pool.start(task)

    def _task_done(self, exc: object) -> None:
        self.session._listeners = self._held_listeners
        self._busy = False
        self._cancel = None
        self._task = None
        self._progress.setVisible(False)
        self._report(exc if isinstance(exc, BaseException) else None)
        self.refresh()

    def _on_progress(self, fraction: float, message: str) -> None:
        self._progress.setValue(int(1000 * fraction))
        self.statusBar().showMessage(message)

    def cancel(self) -> None:
        if self._cancel is not None:
            self._cancel.cancel()

    def wait_idle(self, timeout: float = 600.0) -> None:
        """Process events until no command runs (tests and scripted use)."""
        end = time.monotonic() + timeout
        app = QApplication.instance()
        while self._busy and time.monotonic() < end:
            app.processEvents()
            time.sleep(0.005)
        app.processEvents()
        if self._busy:
            raise TimeoutError("command still running")

    def undo(self) -> None:
        self._guard(self.session.undo)

    def redo(self) -> None:
        self._guard(self.session.redo)

    def new_project(self) -> None:
        self._guard(self.session.new)

    def open_project(self, directory: str | None = None, *, replay: bool = True) -> None:
        path = directory or self.ask_path("dir", "Open SHORE project")
        if path:
            self._guard(lambda: self.session.open(path, replay=replay))

    def save_project(self, directory: str | None = None, *, ask: bool = False) -> None:
        path = directory
        if path is None and (ask or self.session.directory is None):
            path = self.ask_path("dir", "Save SHORE project")
            if not path:
                return
        self._guard(lambda: self.session.save(path))
        if self.last_error is None:
            self.statusBar().showMessage(f"Saved {self.session.directory}", 5000)

    def _picked(self, result: PickResult | None) -> None:
        self.last_pick = result
        if result is None:
            self.statusBar().clearMessage()
            return
        key = result.key
        x, y, z = (float(v) for v in result.xyz)
        self.statusBar().showMessage(
            f"{key.component}/{key.label}  {result.face}  cell {result.cell}  "
            f"xyz ({x:.9g}, {y:.9g}, {z:.9g})"
        )


def main(project: str | None = None, *, replay: bool = True, smoke: bool = False) -> int:
    """Run the GUI; *smoke* builds the window, shows it briefly and exits."""
    app = QApplication.instance() or QApplication([])
    app.setStyle("Fusion")
    window = MainWindow()
    if project:
        window.open_project(project, replay=replay)
        if window.last_error is not None:
            print(f"Error: {window.last_error}")
            return 1
    window.show()
    if smoke:
        for _ in range(20):
            app.processEvents()
            time.sleep(0.01)
        n = len(window.session.blocks())
        print(f"shore gui smoke: {n} blocks, journal {len(window.session.journal_lines())} lines")
        window.close()
        return 0
    return int(app.exec())
