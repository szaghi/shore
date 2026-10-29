"""Project: the stateful document the CLI, the Python API and the GUI drive.

Every change goes through a ``Command`` (``shore.project.commands``); the
public methods below build one and execute it. The project keeps the applied
commands as its history, so it can undo and redo them and render them as a
*journal*: a plain Python script of the same method calls that
``shore run`` replays (ADR 0004).

A project holds named geometries and named components (ADR 0016): every
mesh command acts on the current component (``use_component``; the default
is ``"main"``), and ``write_xall`` writes the Xall files of all of them.

On disk a project is a directory (format 3; 1 and 2 still load)::

    project.toml                 metadata: format, shore version, components
                                 and their block labels, input hashes
    journal.py                   the state-changing history, reading inputs/
    inputs/                      copies of every file the journal reads
    adjacency.json               the current component's face graph
    grids/<component>/<n>.npz    node arrays, checked on load
    grids/<component>/adjacency.json  each component's face graph (ADR 0033)

``Project.load`` checks the input copies against their hashes, replays
``journal.py`` with the inputs resolved against the project directory
(wherever it now is), and requires every grid to be bit-identical to the
saved one, so a load is also a determinism check. Format-1 projects (one
mesh, inputs at their journaled paths) still load. ``load(replay=False)``
(ADR 0033) skips the rebuild: the journal is recorded (geometries loaded,
components selected), and each component's blocks are read from its grids and
rewired from its face graph.
"""

from __future__ import annotations

import dataclasses
import hashlib
import runpy
import shutil
import tomllib
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from shore.exceptions import MeshInputError, ParameterError
from shore.mesh.blocking import Blocking
from shore.mesh.body import Body
from shore.mesh.mesh import BLEND_ZERO, Mesh
from shore.project.commands import (
    Associate,
    AssociateBlocking,
    BuildCollar,
    BuildPrimitive,
    BuildState,
    BuildTipCap,
    BuildTopology,
    CancelToken,
    Command,
    Context,
    CreateBlocking,
    DeleteBlocks,
    FillFace,
    LoadBlocks,
    LoadGeometry,
    March,
    Ogrid,
    Premesh,
    ProgressFn,
    Redistribute,
    Snap,
    Split,
    SplitBlocking,
    Undo,
    UndoCache,
    UseComponent,
    WriteAdjacency,
    WriteCcPar,
    WriteCgns,
    WriteGeo,
    WriteGrd,
    WritePlot3d,
    WriteProcInput,
    WriteXall,
    input_root,
    legacy_collar_volume,
    recording,
)

PROJECT_FORMAT = 3
DEFAULT_COMPONENT = "main"
JOURNAL_HEADER = (
    "# SHORE journal: replay with `shore run <this file>`.\n"
    "from shore.project import Project\n\n"
    "p = Project()\n"
)


@dataclasses.dataclass
class _Component:
    mesh: Mesh | None = None
    build: BuildState | None = None
    blocking: Blocking | None = None


class Project:
    """Meshing project: geometries, components and the command history that made them."""

    def __init__(self, undo_root: str | Path | None = None) -> None:
        self._geometries: dict[str, Body] = {}
        self._digests: dict[str, str] = {}  # geometry name -> SHA-256 of its file
        self._components: dict[str, _Component] = {DEFAULT_COMPONENT: _Component()}
        self._current = DEFAULT_COMPONENT
        self._input_root = input_root.get()  # set while a saved project replays
        self.last_result: Any = None  # what the last command reported (paths, summaries)
        self.undo_cache = UndoCache(Path(undo_root) if undo_root else None)
        self._done: list[tuple[Command, Undo | None]] = []
        self._built_on: dict[str, tuple[str, ...]] = {}  # recorded (not replayed) meshes
        self._redo: list[Command] = []

    # -- state ---------------------------------------------------------------

    @property
    def mesh(self) -> Mesh | None:
        """The current component's mesh."""
        return self._components[self._current].mesh

    @mesh.setter
    def mesh(self, value: Mesh | None) -> None:
        self._components[self._current].mesh = value

    @property
    def build(self) -> BuildState | None:
        """What later commands need from the command that made the current mesh."""
        return self._components[self._current].build

    @build.setter
    def build(self, value: BuildState | None) -> None:
        self._components[self._current].build = value

    @property
    def blocking(self) -> Blocking | None:
        """The current component's blocking (ADR 0042)."""
        return self._components[self._current].blocking

    @blocking.setter
    def blocking(self, value: Blocking | None) -> None:
        self._components[self._current].blocking = value

    @property
    def blockings(self) -> dict[str, Blocking | None]:
        """Each component's blocking (``None`` without one), in creation order."""
        return {name: c.blocking for name, c in self._components.items()}

    def require_blocking(self, command: str) -> Blocking:
        if self.blocking is None:
            raise ParameterError(f"{command}: no blocking: call create_blocking first.")
        return self.blocking

    @property
    def body(self) -> Body | None:
        """The geometry named ``"body"`` (the default of every command)."""
        return self._geometries.get("body")

    @property
    def geometries(self) -> Mapping[str, Body]:
        """Loaded geometries by name (read-only)."""
        return MappingProxyType(self._geometries)

    @property
    def components(self) -> dict[str, Mesh | None]:
        """Each component's mesh (``None`` while empty), in creation order."""
        return {name: c.mesh for name, c in self._components.items()}

    @property
    def current_component(self) -> str:
        return self._current

    def require_geometry(self, name: str, command: str) -> Body:
        if name not in self._geometries:
            if name == "body" and not self._geometries:
                raise ParameterError(f"{command} needs a body: call load_geometry first.")
            raise ParameterError(
                f"{command}: no geometry named {name!r}; loaded: {sorted(self._geometries)}."
            )
        return self._geometries[name]

    def geometry_digest(self, name: str) -> str | None:
        """SHA-256 of the file geometry *name* was loaded from (None if not loaded)."""
        return self._digests.get(name)

    def set_geometry(self, name: str, body: Body | None, digest: str | None = None) -> None:
        """Register (or, with ``None``, remove) geometry *name*; for commands."""
        if body is None:
            self._geometries.pop(name, None)
            self._digests.pop(name, None)
        else:
            self._geometries[name] = body
            if digest is not None:
                self._digests[name] = digest

    def discard_built_on(self, name: str) -> dict[str, tuple[Mesh | None, BuildState | None]]:
        """Empty the components whose mesh was built on geometry *name*; return
        their previous state for ``restore_components``. For commands."""
        old: dict[str, tuple[Mesh | None, BuildState | None]] = {}
        for key, c in self._components.items():
            built_on = c.build.geometries if c.build is not None else ()
            # as before ADR 0016, reloading "body" also drops the current
            # component's geometry-free mesh (primitives, loaded blocks)
            legacy = name == "body" and key == self._current and not built_on
            if name in built_on or (legacy and c.mesh is not None):
                old[key] = (c.mesh, c.build)
                c.mesh, c.build = None, None
        return old

    def restore_components(self, old: dict[str, tuple[Mesh | None, BuildState | None]]) -> None:
        for key, (mesh, build) in old.items():
            if key in self._components:
                self._components[key].mesh, self._components[key].build = mesh, build

    def select_component(self, name: str) -> bool:
        """Make *name* current, creating it if new; True if it was created. For commands."""
        created = name not in self._components
        if created:
            self._components[name] = _Component()
        self._current = name
        return created

    def drop_component(self, name: str) -> None:
        """Remove component *name* (not the current one). For commands."""
        if name != self._current:
            self._components.pop(name, None)

    # -- commands ------------------------------------------------------------

    def load_geometry(
        self,
        path: str | Path,
        anchor_strategy: str = "auto",
        skip_hygiene: bool = False,
        *,
        name: str = "body",
    ) -> None:
        """Load geometry *name* from STL or STEP/IGES/BREP (ADR 0016)."""
        self.execute(LoadGeometry(str(path), anchor_strategy, skip_hygiene, name))

    def use_component(self, name: str) -> None:
        """Make component *name* current, creating it empty if new (ADR 0016)."""
        self.execute(UseComponent(name))

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
        geometry: str = "body",
        **params: Any,
    ) -> None:
        """Build the blocks of *topology* on *geometry* (see
        ``shore.pipeline.build_topology``)."""
        self.execute(
            BuildTopology(
                topology,
                nk,
                ds,
                growth,
                k_law,
                k_beta,
                k_thickness,
                dict(options or {}),
                params,
                geometry,
            )
        )

    def march(
        self,
        smooth_strength: float = 0.2,
        smooth_iters: int = 2,
        blend_normals_k: int | None = None,
        smooth_mode: str = "limited",
        march_scheme: str = "algebraic",
        *,
        progress: ProgressFn | None = None,
        cancel: CancelToken | None = None,
    ) -> None:
        """Extrude the blocks from the wall; *progress* / *cancel* are hooks.

        ``blend_normals_k`` is removed (ADR 0034): 0 warns, anything else is
        refused when the command runs (a journal loaded without replay keeps it).
        """
        if blend_normals_k == 0:
            warnings.warn(BLEND_ZERO, DeprecationWarning, stacklevel=2)
        self.execute(
            March(smooth_strength, smooth_iters, blend_normals_k or 0, smooth_mode, march_scheme),
            progress=progress,
            cancel=cancel,
        )

    def snap(self, exact: bool = False) -> None:
        """Project associated nodes onto their geometry (any project geometry, ADR 0039)."""
        self.execute(Snap(exact))

    def redistribute(
        self,
        block: str,
        edge: str,
        law: str = "uniform",
        *,
        beta: float = 3.0,
        cluster_at: float | None = None,
        first_step: float | None = None,
        last_step: float | None = None,
        interior: str = "uniform",
    ) -> None:
        """Apply a spacing law along *edge* of *block* to the lines of its axis,
        through the seams (``Redistribute``, ADR 0041); then ``snap`` if needed."""
        self.execute(
            Redistribute(block, edge, law, beta, cluster_at, first_step, last_step, interior)
        )

    def associate(
        self,
        block: str,
        entity: str,
        geometry: str = "body",
        kind: str = "body",
        indices: Sequence[int] = (),
    ) -> None:
        """Associate face, edge or vertex *entity* of *block* (current component)
        with entities *indices* of kind *kind* of geometry *geometry*, or clear
        it (``kind="none"``); see ``Associate`` (ADR 0003, 0039)."""
        self.execute(Associate(block, entity, geometry, kind, tuple(int(i) for i in indices)))

    def create_blocking(
        self,
        ni: int,
        nj: int,
        nk: int,
        *,
        geometry: str = "body",
        margin: float = 0.0,
        vmin: Sequence[float] | None = None,
        vmax: Sequence[float] | None = None,
    ) -> None:
        """Start the blocking with one block of *ni* x *nj* x *nk* nodes around
        *geometry* (bounding box grown by *margin* x its extent per side) or
        between *vmin* and *vmax* (``CreateBlocking``, ADR 0042)."""
        corners = [None if c is None else tuple(float(x) for x in c) for c in (vmin, vmax)]
        self.execute(CreateBlocking(int(ni), int(nj), int(nk), geometry, float(margin), *corners))

    def split_blocking(self, block: str, axis: str, fraction: float) -> None:
        """Cut *block* across *axis* at *fraction*, through its sheet (``SplitBlocking``)."""
        self.execute(SplitBlocking(block, axis, float(fraction)))

    def delete_blocks(self, blocks: Sequence[str]) -> None:
        """Remove blocking blocks *blocks* (``DeleteBlocks``)."""
        self.execute(DeleteBlocks(tuple(blocks)))

    def ogrid(
        self,
        blocks: Sequence[str] = (),
        *,
        faces: Sequence[Sequence[str]] = (),
        exclude: Sequence[Sequence[str]] = (),
        factor: float = 0.5,
        nodes: int = 5,
    ) -> None:
        """Insert an O-grid on the faces of blocking blocks *blocks* (all when
        empty): *faces*, or the set's boundary less *exclude* (``(block, face)``
        pairs), sized by *factor*, *nodes* across (``Ogrid``, ADR 0043)."""

        def pairs(items: Sequence[Sequence[str]]) -> tuple[tuple[str, str], ...]:
            return tuple((str(b), str(f)) for b, f in items)

        self.execute(Ogrid(tuple(blocks), pairs(faces), pairs(exclude), float(factor), int(nodes)))

    def associate_blocking(
        self,
        block: str,
        entity: str,
        geometry: str = "body",
        kind: str = "body",
        indices: Sequence[int] = (),
    ) -> None:
        """Associate face, edge or vertex *entity* of blocking block *block* with
        geometry, as ``associate`` does for mesh blocks (``AssociateBlocking``)."""
        self.execute(
            AssociateBlocking(block, entity, geometry, kind, tuple(int(i) for i in indices))
        )

    def premesh(self, exact: bool = False) -> None:
        """Make the component's mesh from its blocking (``Premesh``, ADR 0042)."""
        self.execute(Premesh(exact))

    def surface_grid(
        self,
        block: str,
        face: str,
        *,
        surface: int,
        smoothing: str = "elliptic",
        geometry: str = "body",
    ) -> None:
        """Rebuild *face* of *block* from its boundary on SURFACE entity *surface*
        of *geometry* (TFI + elliptic smoothing in (u, v), ADR 0007)."""
        self.execute(FillFace(block, face, surface, smoothing, geometry))

    def build_primitive(self, kind: str, **params: Any) -> None:
        """Build a parametric primitive (``box``, ``annulus``, ``flat_caps``,
        ``h_grid_rect``, ``h_grid_tube``); see ``shore.pipeline.build_primitive``."""
        self.execute(BuildPrimitive(kind, params))

    def build_collar(
        self,
        other: str | Path,
        *,
        seed: tuple[float, float, float],
        n_loop: int,
        ds_surface: float,
        growth_surface: float,
        m_a: int,
        m_b: int,
        ds_wall: float,
        growth_wall: float,
        nk: int,
        label: str = "collar",
        plane: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None,
        geometry: str = "body",
        loop_law: str = "uniform",
        loop_beta: float = 3.0,
        loop_cluster_at: float | None = None,
        loop_curvature: float = 0.0,
        split: bool = False,
        wake_length: float | None = None,
        volume: str | None = None,
    ) -> None:
        """Collar block around the junction of *geometry* and *other*, a
        geometry name or a file path (ADR 0010, 0016); with *plane* the
        half-collar of a half model (ADR 0013); ``loop_*`` distribute the
        nodes along the junction curve (ADR 0020); *split* cuts the collar
        at sharp corners, with wakes of *wake_length* (ADR 0024); *volume*
        the outer-boundary construction (ADR 0025). Replaces the current
        component's mesh. See ``BuildCollar``."""
        self.execute(
            BuildCollar(
                str(other),
                tuple(float(c) for c in seed),  # type: ignore[arg-type]
                n_loop,
                ds_surface,
                growth_surface,
                m_a,
                m_b,
                ds_wall,
                growth_wall,
                nk,
                label,
                None if plane is None else _plane_tuple(plane),
                geometry,
                loop_law,
                float(loop_beta),
                None if loop_cluster_at is None else float(loop_cluster_at),
                float(loop_curvature),
                bool(split),
                None if wake_length is None else float(wake_length),
                None if volume is None else str(volume),
            )
        )

    def build_tip_cap(self, body_geo: str | Path, **options: Any) -> None:
        """Wing-tip cap from a C-grid body ``.geo``; options as ``BuildTipCap``."""
        self.execute(BuildTipCap(str(body_geo), **options))

    def load_blocks(
        self, geo_paths: list[str | Path], adjacency: list[str | Path] | tuple = ()
    ) -> None:
        """Replace the mesh with blocks read from ``.geo`` files (+ sidecars)."""
        self.execute(LoadBlocks(tuple(map(str, geo_paths)), tuple(map(str, adjacency))))

    def split(self, splits: list[tuple[str, str, list[int] | tuple[int, ...]]]) -> None:
        """Split blocks by ``(label, axis, at)`` directives."""
        self.execute(
            Split(tuple((str(lb), str(ax), tuple(int(i) for i in at)) for lb, ax, at in splits))
        )

    def write_geo(self, stem: str | Path, per_block: bool = False) -> None:
        """Write the ``.geo`` files (see ``WriteGeo`` for the naming rules)."""
        self.execute(WriteGeo(str(stem), per_block))

    def write_grd(self, path: str | Path, geo_paths: list[str | Path]) -> None:
        """Pack ``.geo`` files into an Xall ``.grd``."""
        self.execute(WriteGrd(str(path), tuple(map(str, geo_paths))))

    def write_cc_par(self, adjacency: str | Path, path: str | Path, **options: Any) -> None:
        """Write an Xall ``cc.par``; options as ``WriteCcPar``."""
        if "boxes" in options:
            options["boxes"] = tuple(options["boxes"] or ())
        if "multigrid_levels" in options:
            options["multigrid_levels"] = tuple(options["multigrid_levels"])
        self.execute(WriteCcPar(str(adjacency), str(path), **options))

    def write_proc_input(self, grd: str | Path, path: str | Path, **options: Any) -> None:
        """Write an Xall ``proc.input``; options as ``WriteProcInput``."""
        self.execute(WriteProcInput(str(grd), str(path), **options))

    def write_adjacency(self, path: str | Path) -> None:
        """Write the ``.adjacency.json`` sidecar."""
        self.execute(WriteAdjacency(str(path)))

    def write_cgns(self, path: str | Path, components: Sequence[str] = ()) -> None:
        """Write a CGNS file of *components* (all with a mesh when empty): a base
        per component, a zone per block, 1-to-1 connections and BCs from the
        face graph (ADR 0035). Needs ``h5py`` (``shore-mesh[cgns]``)."""
        self.execute(WriteCgns(str(path), tuple(components)))

    def write_plot3d(self, path: str | Path, components: Sequence[str] = ()) -> None:
        """Write the blocks of *components* (all with a mesh when empty) as one
        multi-block Plot3D grid file (ADR 0035)."""
        self.execute(WritePlot3d(str(path), tuple(components)))

    def write_xall(
        self,
        directory: str | Path,
        *,
        grd: str = "mesh.grd",
        np: int = 1,
        components: dict[str, dict[str, Any]] | None = None,
        output_basename: str = "cc.",
        multigrid_levels: tuple[int, int] = (4, 1),
        boundary_layer_thickness: float = -0.1,
        numerical_beach: float = -1.0,
    ) -> None:
        """Write the Xall files of every component, merged; per-component and
        run-wide options as ``WriteXall`` (ADR 0016, 0017)."""
        self.execute(
            WriteXall(
                str(directory),
                grd,
                np,
                {name: dict(options) for name, options in (components or {}).items()},
                output_basename,
                tuple(multigrid_levels),  # type: ignore[arg-type]
                boundary_layer_thickness,
                numerical_beach,
            )
        )

    # -- execution, undo, redo -------------------------------------------------

    def execute(
        self,
        command: Command,
        *,
        progress: ProgressFn | None = None,
        cancel: CancelToken | None = None,
    ) -> None:
        """Apply *command*, record it, and clear the redo stack.

        While a project loads without replay (ADR 0033) only geometries and
        component switches are applied; the other commands are recorded."""
        if recording.get() and not isinstance(command, LoadGeometry | UseComponent):
            if command.mutates:
                names = (getattr(command, f, None) for f in ("geometry", "other"))
                built_on = tuple(v for v in names if v in self._geometries)
                if built_on:
                    self._built_on[self._current] = built_on
            self._done.append((command, None))
            return
        undo = command.apply(self, Context(progress, cancel))
        self._done.append((command, undo))
        self._redo.clear()

    def undo(self) -> Command:
        """Revert the last state-changing command (and any outputs after it)."""
        while self._done:
            command, undo = self._done.pop()
            if command.mutates:
                if undo is None:
                    self._done.append((command, undo))
                    raise ParameterError(
                        "the history before this point was loaded without replay "
                        "(Project.load(replay=False)): it cannot be undone."
                    )
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
        """Write the project directory, inputs included (see module docstring)."""
        from shore import __version__
        from shore.io.adjacency import write_adjacency_json

        d = Path(directory)
        (d / "grids").mkdir(parents=True, exist_ok=True)
        # stage the copies: the project may be saved back where its inputs live
        staging = d / ".inputs.new"
        shutil.rmtree(staging, ignore_errors=True)
        commands, copies = self._portable_history(staging)
        shutil.rmtree(d / "inputs", ignore_errors=True)
        if staging.exists():
            staging.rename(d / "inputs")
        (d / "journal.py").write_text(
            JOURNAL_HEADER + "".join(c.to_source() + "\n" for c in commands)
        )
        if self.mesh is not None:
            write_adjacency_json(self.mesh.blocks, d / "adjacency.json")
        meta = [
            f"format = {PROJECT_FORMAT}",
            f'shore_version = "{__version__}"',
            f"commands = {len(commands)}",
            "",
            "[components]",
        ]
        for name, component in self._components.items():
            labels = []
            if component.mesh is not None:
                (d / "grids" / name).mkdir(parents=True, exist_ok=True)
                for i, block in enumerate(component.mesh.blocks):
                    np.savez(d / "grids" / name / f"{i:04d}.npz", nodes=block.nodes)
                    labels.append(block.label)
                # the component's face graph: a load without replay rewires from it
                write_adjacency_json(component.mesh.blocks, d / "grids" / name / "adjacency.json")
            meta.append(f"{_toml_key(name)} = [" + ", ".join(_toml_str(x) for x in labels) + "]")
        meta += ["", "[inputs]"]
        meta += [f"{_toml_key(rel)} = {_toml_str(digest)}" for rel, digest in copies.items()]
        (d / "project.toml").write_text("\n".join(meta) + "\n")
        return d

    def _portable_history(self, inputs: Path) -> tuple[list[Command], dict[str, str]]:
        """The state-changing history with every input file copied into *inputs*.

        Returns the rewritten commands and ``{relative copy path: sha256}``.
        A ``build_collar`` whose *other* names a loaded geometry has no file.
        """
        commands: list[Command] = []
        copies: dict[str, str] = {}
        by_source: dict[Path, str] = {}
        names: set[str] = set()
        for command in self.history:
            if not command.mutates:
                continue  # outputs: replaying them would write where the journal says
            changes: dict[str, Any] = {}
            for field in command.inputs:
                if isinstance(command, BuildCollar) and command.other in names:
                    continue
                value = getattr(command, field)
                paths = value if isinstance(value, tuple) else (value,)
                moved = tuple(self._copy_input(p, inputs, by_source, copies) for p in paths)
                changes[field] = moved if isinstance(value, tuple) else moved[0]
            if isinstance(command, LoadGeometry):
                names.add(command.name)
            commands.append(dataclasses.replace(command, **changes) if changes else command)  # type: ignore[type-var]
        return commands, copies

    def _copy_input(
        self, path: str, inputs: Path, by_source: dict[Path, str], copies: dict[str, str]
    ) -> str:
        root = self._input_root
        source = Path(path) if root is None or Path(path).is_absolute() else root / path
        source = source.resolve()
        if source not in by_source:
            if not source.is_file():
                raise MeshInputError(f"save: input file {source} no longer exists.")
            name = f"{len(by_source):02d}_{_original_name(source)}"
            inputs.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, inputs / name)
            rel = f"inputs/{name}"
            by_source[source] = rel
            copies[rel] = _sha256(inputs / name)
        return by_source[source]

    @classmethod
    def load(cls, directory: str | Path, *, replay: bool = True) -> Project:
        """Replay a project directory's journal and verify its saved grids.

        With ``replay=False`` (ADR 0033) the journal is only recorded:
        geometries are loaded and components selected, but no mesh is rebuilt;
        each component's blocks come from its saved grids and face graph. The
        loaded history cannot be undone, and a loaded mesh cannot be marched.
        """
        d = Path(directory)
        try:
            meta = tomllib.loads((d / "project.toml").read_text())
        except (OSError, tomllib.TOMLDecodeError) as err:
            raise MeshInputError(f"{d}: not a SHORE project ({err})") from err
        version = meta.get("format")
        if version not in (1, 2, PROJECT_FORMAT):
            raise MeshInputError(f"{d}: unsupported project format {version!r}")
        if not replay and not _saved_face_graphs(d, meta):
            raise MeshInputError(
                f"{d}: loading without replay needs each component's face graph "
                "(grids/<component>/adjacency.json, saved since ADR 0033): load the project "
                "with replay once and save it again."
            )
        # before format 3 collars defaulted to the "stations" volume (ADR 0026)
        legacy = legacy_collar_volume.set(version < 3)
        try:
            return cls._load(d, meta, version) if replay else cls._load_saved(d, meta)
        finally:
            legacy_collar_volume.reset(legacy)

    @classmethod
    def _load_saved(cls, d: Path, meta: dict[str, Any]) -> Project:
        """``load(replay=False)``: the journal recorded, the meshes read back."""
        from shore.io.adjacency import read_adjacency_json, wire_blocks_from_adjacency
        from shore.mesh.block import HexBlock

        cls._check_inputs(d, meta)
        token = recording.set(True)
        try:
            project = run_journal(d / "journal.py", input_dir=d)
        finally:
            recording.reset(token)
        saved = meta.get("components", {})
        if list(saved) != list(project._components):
            raise MeshInputError(
                f"{d}: recorded components {list(project._components)} != {list(saved)}"
            )
        for name, labels in saved.items():
            component = project._components[name]
            if not labels:
                continue
            folder = d / "grids" / name
            blocks = []
            for i, label in enumerate(labels):
                with np.load(folder / f"{i:04d}.npz") as grid:
                    blocks.append(HexBlock(label=label, nodes=grid["nodes"]))
            adjacency = folder / "adjacency.json"
            payload = read_adjacency_json(adjacency)
            wire_blocks_from_adjacency(blocks, payload, source=str(adjacency))
            component.mesh = Mesh(body=None, blocks=blocks, topology=None)  # type: ignore[arg-type]
            component.build = BuildState(
                "blocks",
                int(blocks[0].nodes.shape[2]),
                marchable=False,
                geometries=project._built_on.get(name, ()),
            )
        return project

    @staticmethod
    def _check_inputs(d: Path, meta: dict[str, Any]) -> None:
        for rel, digest in meta.get("inputs", {}).items():
            copy = d / rel
            if not copy.is_file() or _sha256(copy) != digest:
                raise MeshInputError(f"{d}: input {rel} is missing or changed since the save.")

    @classmethod
    def _load(cls, d: Path, meta: dict[str, Any], version: int) -> Project:
        if version == 1:
            project = run_journal(d / "journal.py")
            _check_grids(d, d / "grids", project.mesh, meta.get("labels", []))
            return project
        cls._check_inputs(d, meta)
        project = run_journal(d / "journal.py", input_dir=d)
        saved = meta.get("components", {})
        if list(saved) != list(project._components):
            raise MeshInputError(
                f"{d}: replayed components {list(project._components)} != {list(saved)}"
            )
        for name, labels in saved.items():
            _check_grids(d, d / "grids" / name, project._components[name].mesh, labels)
        return project

    def close(self) -> None:
        """Delete the on-disk undo snapshots."""
        self.undo_cache.close()

    def __del__(self) -> None:
        self.undo_cache.close()


def run_journal(path: str | Path, input_dir: str | Path | None = None) -> Project:
    """Execute a journal script and return the ``Project`` it built (``p``).

    With *input_dir*, relative input paths in the journal resolve against it
    (a saved project's directory). A journal is Python: run only journals you
    trust.
    """
    token = input_root.set(None if input_dir is None else Path(input_dir).resolve())
    try:
        namespace = runpy.run_path(str(path), run_name="__shore_journal__")
    finally:
        input_root.reset(token)
    project = namespace.get("p")
    if not isinstance(project, Project):
        raise MeshInputError(f"{path}: journal does not define a Project named 'p'.")
    return project


def _plane_tuple(plane: Any) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """``(point, normal)`` as plain float tuples, so the journal line replays."""
    try:
        point, normal = plane
        return tuple(float(c) for c in point), tuple(float(c) for c in normal)  # type: ignore[return-value]
    except (TypeError, ValueError):
        raise ParameterError("plane must be (point, normal), two 3-vectors.") from None


def _original_name(source: Path) -> str:
    """The file name of an input, without the ``NN_`` index of a saved project's copy
    (so re-saving a loaded project does not stack indices: ``00_00_plate.stl``)."""
    name = source.name
    is_copy = source.parent.name == "inputs" and (source.parent.parent / "project.toml").is_file()
    if is_copy and len(name) > 3 and name[:2].isdigit() and name[2] == "_":
        return name[3:]
    return name


def _saved_face_graphs(d: Path, meta: dict[str, Any]) -> bool:
    """Whether every non-empty component of a saved project has its face graph."""
    if meta.get("format") != PROJECT_FORMAT:
        return False
    return all(
        (d / "grids" / name / "adjacency.json").is_file()
        for name, labels in meta.get("components", {}).items()
        if labels
    )


def _check_grids(d: Path, grids: Path, mesh: Mesh | None, labels: list[str]) -> None:
    blocks = mesh.blocks if mesh is not None else []
    if [b.label for b in blocks] != labels:
        raise MeshInputError(f"{d}: replayed blocks {[b.label for b in blocks]} != {labels}")
    for i, block in enumerate(blocks):
        with np.load(grids / f"{i:04d}.npz") as saved:
            if not np.array_equal(saved["nodes"], block.nodes):
                raise MeshInputError(
                    f"{d}: replay of block '{block.label}' differs from the saved grid"
                    " (journal not deterministic, or inputs changed)"
                )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _toml_str(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _toml_key(value: str) -> str:
    return _toml_str(value)
