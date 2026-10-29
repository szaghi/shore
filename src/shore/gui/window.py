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

import numpy as np
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
    QListWidget,
    QMainWindow,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTabWidget,
    QToolBar,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from shore.gui.forms import FORMS, Field, Form, form
from shore.gui.session import BlockKey, Session, SkeletonKey
from shore.gui.viewport import (
    GeometryPick,
    PickResult,
    SkeletonPick,
    SliceKey,
    boundary_values,
    slice_values,
)
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
        self._build_display()
        self._build_association()
        self._build_spacing()
        self._build_blocking()
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
        self.tree.setHeaderLabels(["Model", "Detail", "Association"])
        self.geometry_tree = QTreeWidget()
        self.geometry_tree.setHeaderLabels(["Geometry", "Entities"])
        tabs = QTabWidget()
        tabs.addTab(self.tree, "Model")
        tabs.addTab(self.geometry_tree, "Geometry")
        dock = QDockWidget("Model", self)
        dock.setObjectName("model")
        dock.setWidget(tabs)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)

    def _build_association(self) -> None:
        """Pick a block entity and geometry entities, then associate them (ADR 0039)."""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.selection_mode = QComboBox()
        self.selection_mode.addItems(["face", "edge", "vertex"])
        form = QFormLayout()
        form.addRow("select block", self.selection_mode)
        self.block_selection_label = QLabel("-")
        self.geometry_selection_label = QLabel("-")
        form.addRow("block entity", self.block_selection_label)
        form.addRow("geometry", self.geometry_selection_label)
        layout.addLayout(form)
        layout.addWidget(QLabel("Double click a block, then geometry (shift adds)."))
        self.associate_button = QPushButton("Associate")
        self.associate_button.clicked.connect(self.associate_selection)
        self.dissociate_button = QPushButton("Clear association")
        self.dissociate_button.clicked.connect(self.clear_association)
        layout.addWidget(self.associate_button)
        layout.addWidget(self.dissociate_button)
        layout.addStretch(1)
        dock = QDockWidget("Association", self)
        dock.setObjectName("association")
        dock.setWidget(panel)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        self.selected_block: tuple[BlockKey, str] | None = None
        self.selected_geometry: tuple[str, str, list[int]] | None = None

    def _build_spacing(self) -> None:
        """Redistribute the nodes along a picked block edge (ADR 0041)."""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.addWidget(QLabel("Select a block edge (Association: select edge)."))
        self.spacing_edge_label = QLabel("-")
        form = QFormLayout()
        form.addRow("edge", self.spacing_edge_label)
        self.spacing_law = QComboBox()
        self.spacing_law.addItems(["uniform", "tanh", "tanh2", "pinned"])
        self.spacing_beta = QLineEdit("3.0")
        self.spacing_cluster = QLineEdit("")
        self.spacing_first = QLineEdit("")
        self.spacing_last = QLineEdit("")
        self.spacing_interior = QComboBox()
        self.spacing_interior.addItems(["uniform", "tanh", "tanh2"])
        form.addRow("law", self.spacing_law)
        form.addRow("beta", self.spacing_beta)
        form.addRow("cluster at", self.spacing_cluster)
        form.addRow("first step", self.spacing_first)
        form.addRow("last step", self.spacing_last)
        form.addRow("interior", self.spacing_interior)
        layout.addLayout(form)
        button = QPushButton("Redistribute")
        button.clicked.connect(self._redistribute_from_panel)
        layout.addWidget(button)
        layout.addStretch(1)
        dock = QDockWidget("Spacing", self)
        dock.setObjectName("spacing")
        dock.setWidget(panel)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)

    def _build_blocking(self) -> None:
        """Split the blocking at a picked skeleton node (ADR 0042)."""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.addWidget(QLabel("Double click a blocking edge (orange) at a node."))
        self.skeleton_selection_label = QLabel("-")
        form = QFormLayout()
        form.addRow("edge", self.skeleton_selection_label)
        layout.addLayout(form)
        button = QPushButton("Split here")
        button.clicked.connect(self.split_selection)
        layout.addWidget(button)
        layout.addStretch(1)
        dock = QDockWidget("Blocking", self)
        dock.setObjectName("blocking")
        dock.setWidget(panel)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        self.selected_skeleton: SkeletonPick | None = None
        self._skeleton_box: tuple[float, ...] | None = None

    def split_selection(self) -> None:
        """Run ``Project.split_blocking`` through the selected skeleton node:
        the picked edge's block, its axis and the node's fraction along it
        (making the edge's component current first if needed)."""
        from shore.mesh.blocking import EDGES

        pick = self.selected_skeleton
        if pick is None:
            self._report(ValueError("pick a blocking edge first"))
            return
        key: SkeletonKey = pick.key
        axis = EDGES[key.edge][0]
        fraction = pick.node / (pick.nodes - 1)

        def work() -> None:
            if self.session.project.current_component != key.component:
                self.session.run("use_component", key.component)
            self.session.run("split_blocking", key.block, axis, fraction)

        self._guard(work)
        if self.last_error is None:
            self.selected_skeleton = None
            self.skeleton_selection_label.setText("-")

    def edge_info(self) -> dict[str, Any] | None:
        """The selected edge: block, edge, axis, nodes, length, first and last cell."""
        from shore.gui.viewport import EDGE_NAMES, edge_nodes
        from shore.mesh.redistribute import edge_axis

        if self.selected_block is None or self.selected_block[1] not in EDGE_NAMES:
            return None
        key, edge = self.selected_block
        block = self.session.block(key)
        if block is None:
            return None
        pts = edge_nodes(block.nodes, edge)
        cells = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        return {
            "key": key,
            "edge": edge,
            "axis": edge_axis(edge),
            "nodes": len(pts),
            "length": float(cells.sum()),
            "first": float(cells[0]),
            "last": float(cells[-1]),
        }

    def _show_edge_info(self) -> None:
        info = self.edge_info()
        if info is None:
            self.spacing_edge_label.setText("-")
            return
        key = info["key"]
        self.spacing_edge_label.setText(
            f"{key.component}/{key.label} {info['edge']} ({info['axis']}, {info['nodes']} nodes, "
            f"length {info['length']:.6g}, first {info['first']:.3g}, last {info['last']:.3g})"
        )

    def redistribute_selection(self, law: str = "uniform", **params: Any) -> None:
        """Run ``Project.redistribute`` on the selected edge (making its
        component current first if needed)."""
        info = self.edge_info()
        if info is None:
            self._report(ValueError("select a block edge first (selection mode: edge)"))
            return
        key = info["key"]

        def work() -> None:
            if self.session.project.current_component != key.component:
                self.session.run("use_component", key.component)
            self.session.run("redistribute", key.label, info["edge"], law, **params)

        self._guard(work)
        self._show_edge_info()

    def _redistribute_from_panel(self) -> None:
        params: dict[str, Any] = {}
        try:
            params["beta"] = float(self.spacing_beta.text())
            for name, edit in (
                ("cluster_at", self.spacing_cluster),
                ("first_step", self.spacing_first),
                ("last_step", self.spacing_last),
            ):
                if edit.text().strip():
                    params[name] = float(edit.text())
        except ValueError as exc:
            self._report(exc)
            return
        params["interior"] = self.spacing_interior.currentText()
        self.redistribute_selection(self.spacing_law.currentText(), **params)

    def _build_display(self) -> None:
        """Visibility toggles and slice planes: view state, not commands (ADR 0038)."""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.visibility: dict[str, QCheckBox] = {}
        for kind, label in (
            ("block", "Blocks"),
            ("geometry", "Geometry"),
            ("slice", "Slices"),
            ("skeleton", "Blocking"),
        ):
            box = QCheckBox(label)
            box.setChecked(True)
            box.toggled.connect(lambda on, k=kind: self.set_visible(k, on))
            layout.addWidget(box)
            self.visibility[kind] = box
        form = QFormLayout()
        self.slice_block = QComboBox()
        self.slice_axis = QComboBox()
        self.slice_axis.addItems(["i", "j", "k"])
        self.slice_index = QSpinBox()
        self.slice_index.setRange(0, 0)
        self.slice_block.currentIndexChanged.connect(lambda _: self._slice_range())
        self.slice_axis.currentIndexChanged.connect(lambda _: self._slice_range())
        form.addRow("block", self.slice_block)
        form.addRow("axis", self.slice_axis)
        form.addRow("index", self.slice_index)
        layout.addLayout(form)
        add = QPushButton("Add slice")
        add.clicked.connect(self._add_slice_from_panel)
        layout.addWidget(add)
        self.slice_list = QListWidget()
        layout.addWidget(self.slice_list)
        remove = QPushButton("Remove slice")
        remove.clicked.connect(self._remove_selected_slice)
        layout.addWidget(remove)
        layout.addStretch(1)
        dock = QDockWidget("Display", self)
        dock.setObjectName("display")
        dock.setWidget(panel)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        self.slices: list[SliceKey] = []
        self._geometry_shown: list[tuple[str, int]] = []

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
                linked = self.session.associations(BlockKey(name, label))
                for face, bc in faces.items():
                    child.addChild(QTreeWidgetItem([face, bc, linked.pop(face, "")]))
                for entity, text in linked.items():  # associated edges and vertices
                    child.addChild(QTreeWidgetItem([entity, "", text]))
                item.addChild(child)
            self.tree.addTopLevelItem(item)
            item.setExpanded(True)
        self.journal.setPlainText("\n".join(self.session.journal_lines()))
        blocks = self.session.blocks()
        keys = [key for key, _ in blocks]
        # node values change in place (a march fills the layers): always re-upload
        self.viewport.set_blocks([(key, b.nodes) for key, b in blocks])
        geometries = list(self.session.project.geometries.items())
        shown = [(name, id(body)) for name, body in geometries]
        if shown != self._geometry_shown:
            self._geometry_shown = shown
            self.viewport.set_geometries([(n, b.geometry.display()) for n, b in geometries])
            self.geometry_tree.clear()
            for name, body in geometries:
                d = body.geometry.display()
                if d.curves or len(d.vertices):
                    n_faces = int(d.triangle_entity.max()) + 1 if len(d.triangle_entity) else 0
                    text = f"{n_faces} faces, {len(d.curves)} edges, {len(d.vertices)} vertices"
                else:
                    text = f"{len(d.triangles)} triangles"
                self.geometry_tree.addTopLevelItem(QTreeWidgetItem([name, text]))
        skeleton = self.session.skeleton()
        self.viewport.set_skeleton(skeleton)
        box = None
        if skeleton:
            pts = np.concatenate([pts for _, pts in skeleton])
            box = (*pts.min(axis=0), *pts.max(axis=0))
        # frame the view when the set of blocks or the blocking's extent changes
        if keys != self._shown or (box is not None and box != self._skeleton_box):
            self._shown = keys
            self.viewport.show_all()
        self._skeleton_box = box
        self._refresh_slices()
        if self.quality_action.isChecked():
            self.refresh_quality()
        self._update_actions()

    # -- view state (ADR 0038): not commands, never journaled ---------------------------

    def set_visible(self, kind: str, visible: bool) -> None:
        """Show or hide all blocks, geometries or slices."""
        self.viewport.show(kind, visible)
        box = self.visibility[kind]
        if box.isChecked() != visible:
            box.setChecked(visible)

    def add_slice(self, block: BlockKey, axis: str, index: int) -> SliceKey:
        """Show the constant-*axis* plane *index* of *block*."""
        node = self.session.block(block)
        if node is None:
            raise KeyError(f"no block {block}")
        n = node.nodes.shape["ijk".index(axis)]
        if not 0 <= index < n:
            raise IndexError(f"{axis} = {index} outside 0..{n - 1}")
        key = SliceKey(block, axis, index)
        if key not in self.slices:
            self.slices.append(key)
        self._refresh_slices()
        if self.quality_action.isChecked():
            self.refresh_quality()
        return key

    def remove_slice(self, key: SliceKey) -> None:
        if key in self.slices:
            self.slices.remove(key)
        self._refresh_slices()

    def _refresh_slices(self) -> None:
        """Keep the slices whose block still exists and holds the index."""
        kept = []
        for key in self.slices:
            block = self.session.block(key.block)
            if block is not None and key.index < block.nodes.shape["ijk".index(key.axis)]:
                kept.append(key)
        self.slices = kept
        self.viewport.set_slices([(k, self.session.block(k.block).nodes) for k in kept])
        self.slice_list.clear()
        for k in kept:
            self.slice_list.addItem(f"{k.block.component}/{k.block.label}  {k.axis} = {k.index}")
        current = self.slice_block.currentData()
        self.slice_block.blockSignals(True)
        self.slice_block.clear()
        for bkey, _ in self.session.blocks():
            self.slice_block.addItem(f"{bkey.component}/{bkey.label}", bkey)
        if current is not None:
            i = self.slice_block.findData(current)
            self.slice_block.setCurrentIndex(max(i, 0))
        self.slice_block.blockSignals(False)
        self._slice_range()

    def _slice_range(self) -> None:
        block = (
            self.session.block(self.slice_block.currentData()) if self.slice_block.count() else None
        )
        n = (
            block.nodes.shape["ijk".index(self.slice_axis.currentText())]
            if block is not None
            else 1
        )
        self.slice_index.setRange(0, n - 1)

    def _add_slice_from_panel(self) -> None:
        if self.slice_block.count():
            self.add_slice(
                self.slice_block.currentData(),
                self.slice_axis.currentText(),
                self.slice_index.value(),
            )

    def _remove_selected_slice(self) -> None:
        row = self.slice_list.currentRow()
        if 0 <= row < len(self.slices):
            self.remove_slice(self.slices[row])

    def refresh_quality(self) -> None:
        """Colour boundary quads by the scaled Jacobian of their cells, or clear."""
        if not self.quality_action.isChecked():
            self.viewport.set_face_scalars(None)
            return
        from shore.quality import scaled_jacobian

        metrics = {key: scaled_jacobian(b.nodes) for key, b in self.session.blocks()}
        values: dict[Any, np.ndarray] = {key: boundary_values(m) for key, m in metrics.items()}
        for sk in self.slices:
            if sk.block in metrics:
                values[sk] = slice_values(metrics[sk.block], sk.axis, sk.index)
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

    def pick_at(self, x: float, y: float, *, add: bool = False) -> Any:
        """Pick at pixel (x, y) as a double click does (shift: *add*)."""
        result = self.viewport.pick(x, y)
        self._picked(result, add)
        return result

    def _picked(
        self, result: PickResult | GeometryPick | SkeletonPick | None, add: bool = False
    ) -> None:
        self.last_pick = result
        if result is None:
            self.statusBar().clearMessage()
            return
        x, y, z = (float(v) for v in result.xyz)
        where = f"xyz ({x:.9g}, {y:.9g}, {z:.9g})"
        if isinstance(result, SkeletonPick):
            k = result.key
            self.selected_skeleton = result
            text = f"{k.component}/{k.block} {k.edge}  node {result.node} of {result.nodes}"
            self.skeleton_selection_label.setText(text)
            self.statusBar().showMessage(f"blocking {text}  {where}")
            return
        if isinstance(result, GeometryPick):
            e = result.entity
            self._select_geometry(result.geometry, e.kind.value, e.index, add)
            self.statusBar().showMessage(f"{result.geometry}  {e.kind.value} {e.index}  {where}")
            return
        if isinstance(result.key, SliceKey):
            key = result.key.block
        else:
            key = result.key
            self._select_block(key, result.face, result.xyz)
        self.statusBar().showMessage(
            f"{key.component}/{key.label}  {result.face}  cell {result.cell}  {where}"
        )

    # -- association (ADR 0039) ------------------------------------------------------

    def _select_block(self, key: BlockKey, face: str, xyz: np.ndarray) -> None:
        """Select the face, or its nearest edge / corner, per the selection mode."""
        from shore.gui.viewport import (
            VERTEX_IJK,
            edge_nodes,
            nearest_edge,
            nearest_vertex,
        )

        block = self.session.block(key)
        if block is None:
            return
        mode = self.selection_mode.currentText()
        if mode == "edge":
            entity = nearest_edge(block.nodes, face, xyz)
            self.viewport.highlight_polyline(edge_nodes(block.nodes, entity))
        elif mode == "vertex":
            v = nearest_vertex(block.nodes, face, xyz)
            entity = f"v{v}"
            self.viewport.highlight_polyline(block.nodes[VERTEX_IJK[v]])
        else:
            entity = face
        self.selected_block = (key, entity)
        self.block_selection_label.setText(f"{key.component}/{key.label} {entity}")
        self._show_edge_info()

    def _select_geometry(self, name: str, kind: str, index: int, add: bool) -> None:
        current = self.selected_geometry
        if add and current is not None and current[:2] == (name, kind):
            indices = current[2] + ([index] if index not in current[2] else [])
        else:
            indices = [index]
        self.selected_geometry = (name, kind, indices)
        text = name if kind == "body" else f"{name}: {kind} {', '.join(map(str, indices))}"
        self.geometry_selection_label.setText(text)

    def associate_selection(self) -> None:
        """Associate the selected block entity with the selected geometry entities."""
        if self.selected_block is None or self.selected_geometry is None:
            self._report(ValueError("select a block entity and geometry first"))
            return
        key, entity = self.selected_block
        name, kind, indices = self.selected_geometry
        self._run_on(key, entity, name, kind, () if kind == "body" else tuple(indices))

    def clear_association(self) -> None:
        if self.selected_block is None:
            self._report(ValueError("select a block entity first"))
            return
        key, entity = self.selected_block
        self._run_on(key, entity, "body", "none", ())

    def _run_on(self, key: BlockKey, *args: Any) -> None:
        def work() -> None:
            if self.session.project.current_component != key.component:
                self.session.run("use_component", key.component)
            self.session.run("associate", key.label, *args)

        self._guard(work)


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
