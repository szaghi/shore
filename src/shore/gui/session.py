"""The GUI's owner of a project (ADR 0037): every action is a Project command.

``Session`` holds no Qt object, so the GUI logic is tested headless. Widgets
call its methods only; each state change goes through ``Project`` (and so
into its journal, ADR 0004), and listeners are told what changed.

Long commands (march, builds) are run by the window on a worker thread, with
the ADR 0004 progress and cancel hooks; ``Session.run`` itself is
synchronous.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from shore.mesh.block import HexBlock
from shore.project import Project
from shore.project.commands import CancelToken, ProgressFn
from shore.project.project import JOURNAL_HEADER

Listener = Callable[[], None]


@dataclass(frozen=True)
class BlockKey:
    """A block of a component, as the viewport and the model tree name it."""

    component: str
    label: str


@dataclass(frozen=True)
class SkeletonKey:
    """A blocking entity of a component (ADR 0042, 0044): ``edge`` is an edge
    or face name of blocking block ``block``, or ``""`` for the block itself."""

    component: str
    block: str
    edge: str


class Session:
    """One project, its saved location and the GUI's change listeners."""

    def __init__(self, project: Project | None = None, directory: Path | None = None) -> None:
        self.project = project if project is not None else Project()
        self.directory = directory  # where the project was opened from / saved to
        self._listeners: list[Listener] = []
        #: surfaces shown but not part of the project (``view_file`` on a .npz)
        self.overlays: dict[str, Any] = {}

    # -- listeners -------------------------------------------------------------

    def subscribe(self, listener: Listener) -> None:
        """Call *listener* after every change of the project."""
        self._listeners.append(listener)

    def _changed(self) -> None:
        for listener in list(self._listeners):
            listener()

    # -- actions ---------------------------------------------------------------

    def run(
        self,
        method: str,
        *args: Any,
        progress: ProgressFn | None = None,
        cancel: CancelToken | None = None,
        **kwargs: Any,
    ) -> None:
        """Call ``Project.<method>``: one command, one journal line.

        The project is unchanged if the command raises (the command layer
        restores a failed or cancelled march, and a command that fails while
        building leaves no state behind).
        """
        call = getattr(self.project, method, None)
        if call is None or method.startswith("_") or not callable(call):
            raise AttributeError(f"Project has no command {method!r}")
        if method == "march":
            kwargs.update(progress=progress, cancel=cancel)
        try:
            call(*args, **kwargs)
        finally:
            self._changed()

    def run_cli(self, command: Sequence[str]) -> str:
        """Run a CLI command (``["mesh", "body.stl", "--nk", "6"]``) on this
        project (ADR 0062) and return its output.

        Its Project commands land in this project's history, each one
        undoable. A command that fails (exit code not 0) is rolled back:
        the project is as it was, and the error carries the command's output.
        """
        from typer.testing import CliRunner

        from shore.cli import app, using_project

        before = len(self.project.history)
        try:
            with using_project(self.project):
                result = CliRunner().invoke(app, [str(a) for a in command])
            if result.exit_code != 0:
                self.project.rollback(before)
                detail = result.output.strip()
                if not detail and result.exception is not None:
                    detail = f"{type(result.exception).__name__}: {result.exception}"
                raise RuntimeError(f"shore {' '.join(map(str, command[:2]))}: {detail}")
            return result.output
        finally:
            self._changed()

    def run_journal(self, path: str | Path) -> None:
        """Replace the project with the one a journal builds (``shore run``)."""
        from shore.project import run_journal

        new = run_journal(path)
        old, self.project, self.directory = self.project, new, None
        old.close()
        self._changed()

    def view_file(self, path: str | Path) -> None:
        """Show a file as ``shore view`` does: a ``.geo`` file loads as the
        current component's blocks (``load_blocks``, one journal line); a
        ``.npz`` surface (key ``"k0"``) is shown only, outside the project."""
        import numpy as np

        p = Path(path)
        if p.suffix == ".geo":
            self.run("load_blocks", [str(p)])
            return
        if p.suffix != ".npz":
            raise ValueError(f"view_file: a .geo or .npz file is needed, got {p.name}")
        with np.load(p) as data:
            nodes = np.asarray(data["k0"], dtype=np.float64)
        self.overlays[p.stem] = surface_display(nodes)
        self._changed()

    def undo(self) -> None:
        try:
            self.project.undo()
        finally:
            self._changed()

    def redo(self) -> None:
        try:
            self.project.redo()
        finally:
            self._changed()

    def open(self, directory: str | Path, *, replay: bool = True) -> None:
        """Replace the project with a saved one (ADR 0016; ADR 0033 without replay)."""
        old = self.project
        self.project = Project.load(directory, replay=replay)
        self.directory = Path(directory)
        old.close()
        self._changed()

    def new(self) -> None:
        self.project.close()
        self.project = Project()
        self.directory = None
        self._changed()

    def save(self, directory: str | Path | None = None) -> Path:
        target = Path(directory) if directory is not None else self.directory
        if target is None:
            raise ValueError("save: no directory given and the project was never saved")
        out = self.project.save(target)
        self.directory = target
        self._changed()
        return out

    # -- views -----------------------------------------------------------------

    def journal_lines(self) -> list[str]:
        """The journal's command lines (without the header)."""
        return self.project.journal()[len(JOURNAL_HEADER) :].splitlines()

    def blocks(self) -> list[tuple[BlockKey, HexBlock]]:
        """Every block of every component with a mesh, components in project order."""
        out = []
        for name, mesh in self.project.components.items():
            if mesh is None:
                continue
            out.extend((BlockKey(name, b.label), b) for b in mesh.blocks)
        return out

    def skeleton(self) -> list[tuple[SkeletonKey, Any]]:
        """Every blocking edge of every component, with its nodes on the
        straight line between its skeleton vertices (``Blocking.polylines``)."""
        out = []
        for name, blocking in self.project.blockings.items():
            if blocking is not None:
                out.extend(
                    (SkeletonKey(name, block, edge), pts)
                    for block, edge, pts in blocking.polylines()
                )
        return out

    def skeleton_markers(self, kind: str) -> list[tuple[SkeletonKey, Any]]:
        """A marker per blocking face (``kind="face"``: ``SkeletonKey(component,
        block, face)``, the face named by its first block) or block (``"block"``:
        ``SkeletonKey(component, block, "")``), at its skeleton centre."""
        out = []
        for name, blocking in self.project.blockings.items():
            if blocking is not None:
                out.extend(
                    (SkeletonKey(name, block, face), xyz)
                    for block, face, xyz in blocking.centres(kind)
                )
        return out

    def tree(
        self,
    ) -> list[tuple[str, bool, list[tuple[str, tuple[int, int, int], dict[str, str]]]]]:
        """``(component, is current, [(block label, shape, {face: BC name})])``."""
        current = self.project.current_component
        out = []
        for name, mesh in self.project.components.items():
            blocks = []
            for b in mesh.blocks if mesh is not None else []:
                shape = tuple(int(n) for n in b.nodes.shape[:3])
                blocks.append((b.label, shape, {f.name: f.bc.name for f in b.faces()}))
            out.append((name, name == current, blocks))
        return out

    def associations(self, key: BlockKey) -> dict[str, str]:
        """``{entity: label}`` for the associated faces, edges and vertices
        (``"v0"``..) of a block; labels like ``"body"`` or ``"part: surface 2, 5"``."""
        block = self.block(key)
        if block is None:
            return {}
        out = {
            f.name: self.association_label(f.association) for f in block.faces() if f.association
        }
        out.update(
            {e.name: self.association_label(e.association) for e in block.edges() if e.association}
        )
        out.update(
            {
                f"v{v.index}": self.association_label(v.association)
                for v in block.vertices()
                if v.association
            }
        )
        return out

    def association_label(self, association: Any) -> str:
        """``"geometry"`` for a body, else ``"geometry: kind i, j"`` (indices
        re-found from the association's fingerprints in the current geometry)."""
        from shore.geometry.base import EntityKind
        from shore.geometry.keys import match_keys

        if association.kind is EntityKind.BODY:
            return association.geometry
        body = self.project.geometries.get(association.geometry)
        if body is None:
            return f"{association.geometry}: {association.kind.value} (geometry missing)"
        report = match_keys(list(association.keys), body.geometry)
        if not report.all_matched:
            return f"{association.geometry}: {association.kind.value} (unresolved)"
        indices = ", ".join(str(m.ref.index) for m in report.matches)
        return f"{association.geometry}: {association.kind.value} {indices}"

    def block(self, key: BlockKey) -> HexBlock | None:
        mesh = self.project.components.get(key.component)
        if mesh is None:
            return None
        return next((b for b in mesh.blocks if b.label == key.label), None)


def surface_display(nodes: Any) -> Any:
    """A structured surface ``(ni, nj, 3)`` as a ``geometry.base.Display``:
    two triangles per quad, one body entity."""
    import numpy as np

    from shore.geometry.base import Display, EntityKind

    ni, nj = nodes.shape[:2]
    idx = np.arange(ni * nj).reshape(ni, nj)
    a, b, c, d = idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]
    tris = np.concatenate(
        [np.stack([a, b, c], -1).reshape(-1, 3), np.stack([a, c, d], -1).reshape(-1, 3)]
    )
    return Display(
        points=nodes.reshape(-1, 3),
        triangles=tris,
        triangle_entity=np.zeros(len(tris), dtype=np.int64),
        surface_kind=EntityKind.BODY,
    )
