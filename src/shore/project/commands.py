"""Commands: the only way to change a ``Project`` (ADR 0004).

A command is a frozen dataclass. ``apply`` changes the project and returns an
``Undo`` (or ``None`` for output-only commands such as ``WriteGeo``);
``to_source`` renders the command as the one line of public Python API that
re-issues it, which is what the journal records.

Node-array changes (march, snap) are undone from snapshots that the
``UndoCache`` keeps on disk, not in RAM: a 1e7-node grid is 240 MB.
"""

from __future__ import annotations

import dataclasses
import shutil
import tempfile
from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from shore.exceptions import ParameterError, ShoreError
from shore.mesh.body import Body
from shore.mesh.mesh import Mesh
from shore.pipeline import build_topology, k_spacing

if TYPE_CHECKING:
    from shore.project.project import Project

ProgressFn = Callable[[float, str], None]


class Cancelled(ShoreError):
    """Raised inside a long command when its ``CancelToken`` was set."""


class CancelToken:
    """Cooperative cancellation flag, checked once per layer / sweep."""

    def __init__(self) -> None:
        self.cancelled = False

    def cancel(self) -> None:
        self.cancelled = True

    def check(self) -> None:
        if self.cancelled:
            raise Cancelled("command cancelled")


@dataclass
class Context:
    """Per-execution hooks: progress reporting and cancellation."""

    progress: ProgressFn | None = None
    cancel: CancelToken | None = None

    def report(self, fraction: float, message: str = "") -> None:
        if self.progress is not None:
            self.progress(fraction, message)

    def check(self) -> None:
        if self.cancel is not None:
            self.cancel.check()


@dataclass
class Undo:
    """How to revert one applied command."""

    restore: Callable[[], None]


class UndoCache:
    """On-disk store of node-array snapshots, deleted with the project."""

    def __init__(self, root: Path | None = None) -> None:
        self._dir = Path(tempfile.mkdtemp(prefix="shore-undo-", dir=root))
        self._n = 0

    def snapshot(self, mesh: Mesh) -> Callable[[], None]:
        """Save every block's nodes; return a function restoring them in place."""
        path = self._dir / f"{self._n:06d}.npz"
        self._n += 1
        np.savez(path, **{f"b{i}": b.nodes for i, b in enumerate(mesh.blocks)})
        blocks = list(mesh.blocks)

        def restore() -> None:
            with np.load(path) as data:
                for i, b in enumerate(blocks):
                    b.nodes[...] = data[f"b{i}"]  # in place: keeps seam views aliased

        return restore

    def close(self) -> None:
        shutil.rmtree(self._dir, ignore_errors=True)


class Command(ABC):
    """Base of all commands. Subclasses are frozen dataclasses."""

    method: ClassVar[str]  # the Project method that issues this command
    positional: ClassVar[tuple[str, ...]] = ()  # fields rendered positionally
    mutates: ClassVar[bool] = True  # False: output only, never undone

    @abstractmethod
    def apply(self, project: Project, ctx: Context) -> Undo | None: ...

    def to_source(self) -> str:
        """``p.<method>(...)`` with positional fields and non-default keywords."""
        args = [repr(getattr(self, name)) for name in self.positional]
        for f in dataclasses.fields(self):  # type: ignore[arg-type]
            if f.name in self.positional:
                continue
            value = getattr(self, f.name)
            default = f.default if f.default is not dataclasses.MISSING else None
            if f.default_factory is not dataclasses.MISSING:  # type: ignore[misc]
                default = f.default_factory()  # type: ignore[misc]
            if value != default:
                args.append(f"{f.name}={value!r}")
        return f"p.{self.method}({', '.join(args)})"


@dataclass(frozen=True)
class LoadGeometry(Command):
    """Load the body (STL or STEP/IGES/BREP). Discards any existing mesh."""

    method: ClassVar[str] = "load_geometry"
    positional: ClassVar[tuple[str, ...]] = ("path",)

    path: str
    anchor_strategy: str = "auto"
    skip_hygiene: bool = False

    def apply(self, project: Project, ctx: Context) -> Undo:
        prev = (project.body, project.mesh, project.build)
        project.body = Body.from_file(
            self.path, anchor_strategy=self.anchor_strategy, skip_hygiene=self.skip_hygiene
        )
        project.mesh = None
        project.build = None

        def restore() -> None:
            project.body, project.mesh, project.build = prev

        return Undo(restore)


@dataclass(frozen=True)
class BuildTopology(Command):
    """Build the topology's blocks on the body (``shore mesh`` stage 2).

    The k-spacing fields define the wall-normal law used both by the
    topology (seam pinning) and by ``March``. *options* are the topology
    constructor's arguments (``ch``, ``oh``); *params* go to its ``build``.
    """

    method: ClassVar[str] = "build_topology"
    positional: ClassVar[tuple[str, ...]] = ("topology",)

    topology: str
    nk: int
    ds: float = 1e-3
    growth: float = 1.15
    k_law: str = "geometric"
    k_beta: float = 3.0
    k_thickness: float | None = None
    options: dict[str, Any] = field(default_factory=dict)
    params: dict[str, Any] = field(default_factory=dict)

    def to_source(self) -> str:
        # *params* read better as keywords: p.build_topology('ogrid', nk=6, ni=10)
        head = Command.to_source(dataclasses.replace(self, params={}))
        extra = ", ".join(f"{k}={v!r}" for k, v in self.params.items())
        return head if not extra else f"{head[:-1]}, {extra})"

    def apply(self, project: Project, ctx: Context) -> Undo:
        if project.body is None:
            raise ParameterError("build_topology needs a body: call load_geometry first.")
        spacing = k_spacing(
            self.k_law, ds=self.ds, growth=self.growth, beta=self.k_beta, thickness=self.k_thickness
        )
        topo, blocks = build_topology(
            project.body,
            self.topology,
            nk=self.nk,
            spacing=spacing,
            options=self.options,
            **self.params,
        )
        prev = (project.mesh, project.build)
        project.mesh = Mesh(body=project.body, blocks=blocks, topology=topo)
        project.build = BuildState(self.topology, self.nk, self.ds, self.growth, spacing)

        def restore() -> None:
            project.mesh, project.build = prev

        return Undo(restore)


@dataclass
class BuildState:
    """What later commands need from the last ``BuildTopology``."""

    topology: str
    nk: int
    ds: float
    growth: float
    spacing: Any  # Spacing1D


@dataclass(frozen=True)
class March(Command):
    """Extrude all blocks from the wall (``shore mesh`` stage 3)."""

    method: ClassVar[str] = "march"

    smooth_strength: float = 0.2
    smooth_iters: int = 2
    blend_normals_k: int = 0
    smooth_mode: str = "limited"

    def apply(self, project: Project, ctx: Context) -> Undo:
        mesh, build = project.require_mesh()
        restore = project.undo_cache.snapshot(mesh)
        n_layers = build.nk - 1

        def layer(k: int, jmin: float) -> None:  # k: the layer just accepted, 1-based
            ctx.check()
            ctx.report(k / n_layers, f"layer {k}/{n_layers}, min J {jmin:.2e}")

        try:
            mesh.march(
                ds=build.ds,
                growth=build.growth,
                smooth_strength=self.smooth_strength,
                smooth_iters=self.smooth_iters,
                steps=build.spacing.build_steps(build.nk),
                layer_callback=layer,
                blend_normals_k=self.blend_normals_k,
                smooth_mode=self.smooth_mode,
            )
        except BaseException:
            restore()  # a failed or cancelled march leaves no half-extruded grid
            raise
        return Undo(restore)


@dataclass(frozen=True)
class Snap(Command):
    """Project associated nodes onto their geometry (ADR 0003)."""

    method: ClassVar[str] = "snap"

    exact: bool = False

    def apply(self, project: Project, ctx: Context) -> Undo:
        mesh, _ = project.require_mesh()
        restore = project.undo_cache.snapshot(mesh)
        mesh.snap(exact=self.exact)
        return Undo(restore)


@dataclass(frozen=True)
class WriteGeo(Command):
    """Write ``.geo`` files as ``shore mesh`` does (output only)."""

    method: ClassVar[str] = "write_geo"
    positional: ClassVar[tuple[str, ...]] = ("stem",)
    mutates: ClassVar[bool] = False

    stem: str

    def apply(self, project: Project, ctx: Context) -> None:
        mesh, build = project.require_mesh()
        stem = Path(self.stem)
        if build.topology != "cubed_sphere" and len(mesh.blocks) == 1:
            # single-block O-grid / 1-block C-grid: one <stem>.geo, no label
            from shore.io.geo import write_geo

            write_geo(
                stem if stem.suffix == ".geo" else stem.with_suffix(".geo"), mesh.blocks[0].nodes
            )
        else:
            mesh.write_geo(stem)


@dataclass(frozen=True)
class WriteAdjacency(Command):
    """Write the ``.adjacency.json`` sidecar (output only)."""

    method: ClassVar[str] = "write_adjacency"
    positional: ClassVar[tuple[str, ...]] = ("path",)
    mutates: ClassVar[bool] = False

    path: str

    def apply(self, project: Project, ctx: Context) -> None:
        from shore.io.adjacency import write_adjacency_json

        mesh, _ = project.require_mesh()
        out = Path(self.path)
        out.parent.mkdir(parents=True, exist_ok=True)
        write_adjacency_json(mesh.blocks, out)
