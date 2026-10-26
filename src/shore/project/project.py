"""Project: the stateful document the CLI, the Python API and the GUI drive.

Every change goes through a ``Command`` (``shore.project.commands``); the
public methods below build one and execute it. The project keeps the applied
commands as its history, so it can undo and redo them and render them as a
*journal*: a plain Python script of the same method calls that
``shore run`` replays (ADR 0004).

On disk a project is a directory::

    project.toml      metadata (format version, shore version, labels)
    journal.py        the history; loading a project replays it
    adjacency.json    face graph and associations (when a mesh exists)
    grids/<n>.npz     node arrays, written for inspection and checked on load

``Project.load`` replays ``journal.py`` and then requires every grid to be
bit-identical to the saved one, so a load is also a determinism check.
"""

from __future__ import annotations

import runpy
import tomllib
from pathlib import Path
from typing import Any

import numpy as np

from shore.exceptions import MeshInputError, ParameterError
from shore.mesh.body import Body
from shore.mesh.mesh import Mesh
from shore.project.commands import (
    BuildState,
    BuildTopology,
    CancelToken,
    Command,
    Context,
    LoadGeometry,
    March,
    ProgressFn,
    Snap,
    Undo,
    UndoCache,
    WriteAdjacency,
    WriteGeo,
)

PROJECT_FORMAT = 1
JOURNAL_HEADER = (
    "# SHORE journal: replay with `shore run <this file>`.\n"
    "from shore.project import Project\n\n"
    "p = Project()\n"
)


class Project:
    """Meshing project: body, mesh and the command history that made them."""

    def __init__(self, undo_root: str | Path | None = None) -> None:
        self.body: Body | None = None
        self.mesh: Mesh | None = None
        self.build: BuildState | None = None
        self.undo_cache = UndoCache(Path(undo_root) if undo_root else None)
        self._done: list[tuple[Command, Undo | None]] = []
        self._redo: list[Command] = []

    # -- commands ------------------------------------------------------------

    def load_geometry(
        self, path: str | Path, anchor_strategy: str = "auto", skip_hygiene: bool = False
    ) -> None:
        """Load the body from STL or STEP/IGES/BREP."""
        self.execute(LoadGeometry(str(path), anchor_strategy, skip_hygiene))

    def build_topology(
        self,
        topology: str,
        *,
        nk: int,
        ds: float = 1e-3,
        growth: float = 1.15,
        k_law: str = "geometric",
        k_beta: float = 3.0,
        k_thickness: float | None = None,
        options: dict[str, Any] | None = None,
        **params: Any,
    ) -> None:
        """Build the blocks of *topology* (see ``shore.pipeline.build_topology``)."""
        self.execute(
            BuildTopology(
                topology, nk, ds, growth, k_law, k_beta, k_thickness, dict(options or {}), params
            )
        )

    def march(
        self,
        smooth_strength: float = 0.2,
        smooth_iters: int = 2,
        blend_normals_k: int = 0,
        *,
        progress: ProgressFn | None = None,
        cancel: CancelToken | None = None,
    ) -> None:
        """Extrude the blocks from the wall; *progress* / *cancel* are hooks."""
        self.execute(
            March(smooth_strength, smooth_iters, blend_normals_k), progress=progress, cancel=cancel
        )

    def snap(self, exact: bool = False) -> None:
        """Project associated nodes onto their geometry."""
        self.execute(Snap(exact))

    def write_geo(self, stem: str | Path) -> None:
        """Write the ``.geo`` files (``shore mesh`` naming)."""
        self.execute(WriteGeo(str(stem)))

    def write_adjacency(self, path: str | Path) -> None:
        """Write the ``.adjacency.json`` sidecar."""
        self.execute(WriteAdjacency(str(path)))

    # -- execution, undo, redo -------------------------------------------------

    def execute(
        self,
        command: Command,
        *,
        progress: ProgressFn | None = None,
        cancel: CancelToken | None = None,
    ) -> None:
        """Apply *command*, record it, and clear the redo stack."""
        undo = command.apply(self, Context(progress, cancel))
        self._done.append((command, undo))
        self._redo.clear()

    def undo(self) -> Command:
        """Revert the last state-changing command (and any outputs after it)."""
        while self._done:
            command, undo = self._done.pop()
            if command.mutates:
                assert undo is not None
                undo.restore()
                self._redo.append(command)
                return command
        raise ParameterError("Nothing to undo.")

    def redo(self) -> Command:
        """Re-apply the last undone command."""
        if not self._redo:
            raise ParameterError("Nothing to redo.")
        command = self._redo.pop()
        undo = command.apply(self, Context())
        self._done.append((command, undo))
        return command

    @property
    def history(self) -> list[Command]:
        """Applied commands, oldest first (undone ones excluded)."""
        return [c for c, _ in self._done]

    def require_mesh(self) -> tuple[Mesh, BuildState]:
        if self.mesh is None or self.build is None:
            raise ParameterError("No mesh: call build_topology first.")
        return self.mesh, self.build

    # -- journal and persistence -------------------------------------------------

    def journal(self) -> str:
        """The history as a replayable Python script."""
        return JOURNAL_HEADER + "".join(c.to_source() + "\n" for c in self.history)

    def save_journal(self, path: str | Path) -> Path:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(self.journal())
        return out

    def save(self, directory: str | Path) -> Path:
        """Write the project directory (see module docstring)."""
        from shore import __version__
        from shore.io.adjacency import write_adjacency_json

        d = Path(directory)
        (d / "grids").mkdir(parents=True, exist_ok=True)
        self.save_journal(d / "journal.py")
        labels: list[str] = []
        if self.mesh is not None:
            write_adjacency_json(self.mesh.blocks, d / "adjacency.json")
            for i, block in enumerate(self.mesh.blocks):
                np.savez(d / "grids" / f"{i:04d}.npz", nodes=block.nodes)
                labels.append(block.label)
        meta = [
            f"format = {PROJECT_FORMAT}",
            f'shore_version = "{__version__}"',
            f"commands = {len(self.history)}",
            "labels = [" + ", ".join(f'"{lbl}"' for lbl in labels) + "]",
        ]
        (d / "project.toml").write_text("\n".join(meta) + "\n")
        return d

    @classmethod
    def load(cls, directory: str | Path) -> Project:
        """Replay a project directory's journal and verify its saved grids."""
        d = Path(directory)
        try:
            meta = tomllib.loads((d / "project.toml").read_text())
        except (OSError, tomllib.TOMLDecodeError) as err:
            raise MeshInputError(f"{d}: not a SHORE project ({err})") from err
        if meta.get("format") != PROJECT_FORMAT:
            raise MeshInputError(f"{d}: unsupported project format {meta.get('format')!r}")
        project = run_journal(d / "journal.py")
        labels = meta.get("labels", [])
        blocks = project.mesh.blocks if project.mesh is not None else []
        if [b.label for b in blocks] != labels:
            raise MeshInputError(f"{d}: replayed blocks {[b.label for b in blocks]} != {labels}")
        for i, block in enumerate(blocks):
            with np.load(d / "grids" / f"{i:04d}.npz") as saved:
                if not np.array_equal(saved["nodes"], block.nodes):
                    raise MeshInputError(
                        f"{d}: replay of block '{block.label}' differs from the saved grid"
                        " (journal not deterministic, or inputs changed)"
                    )
        return project

    def close(self) -> None:
        """Delete the on-disk undo snapshots."""
        self.undo_cache.close()

    def __del__(self) -> None:
        self.undo_cache.close()


def run_journal(path: str | Path) -> Project:
    """Execute a journal script and return the ``Project`` it built (``p``).

    A journal is Python: run only journals you trust.
    """
    namespace = runpy.run_path(str(path), run_name="__shore_journal__")
    project = namespace.get("p")
    if not isinstance(project, Project):
        raise MeshInputError(f"{path}: journal does not define a Project named 'p'.")
    return project
