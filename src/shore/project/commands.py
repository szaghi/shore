"""Commands: the only way to change a ``Project`` (ADR 0004).

A command is a frozen dataclass. ``apply`` changes the project and returns an
``Undo`` (or ``None`` for output-only commands such as ``WriteGeo``);
``to_source`` renders the command as the one line of public Python API that
re-issues it, which is what the journal records.

Node-array changes (march, snap) are undone from snapshots that the
``UndoCache`` keeps on disk, not in RAM: a 1e7-node grid is 240 MB.

A project holds named geometries and named components (ADR 0016); mesh
commands act on the current component. Commands that read files name those
fields in ``inputs``: ``Project.save`` copies the files into the project, and
a replay resolves relative input paths against ``input_root`` (the saved
project's directory).
"""

from __future__ import annotations

import contextvars
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
from shore.mesh.mesh import BLEND_REMOVED, Mesh
from shore.pipeline import (
    apply_splits,
    build_primitive,
    build_tip_cap,
    build_topology,
    k_spacing,
    load_blocks,
)

if TYPE_CHECKING:
    from shore.project.project import Project

ProgressFn = Callable[[float, str], None]


#: Directory that relative input paths resolve against while a saved
#: project's journal replays (ADR 0016); unset: the working directory.
input_root: contextvars.ContextVar[Path | None] = contextvars.ContextVar(
    "shore_input_root", default=None
)


#: Set while a project saved before format 3 replays: its collars were built when
#: ``CollarParams.volume`` defaulted to ``"stations"`` (ADR 0026), so a ``BuildCollar``
#: without a volume gets that one, and the saved grids still match.
legacy_collar_volume: contextvars.ContextVar[bool] = contextvars.ContextVar(
    "shore_legacy_collar_volume", default=False
)


#: Set while a project loads without replay (ADR 0033): ``Project.execute``
#: records the commands instead of applying them (geometries and component
#: switches excepted).
recording: contextvars.ContextVar[bool] = contextvars.ContextVar("shore_recording", default=False)


def resolve_input(path: str) -> str:
    """*path*, relative to ``input_root`` when that is set and *path* is relative."""
    root = input_root.get()
    if root is None or Path(path).is_absolute():
        return path
    return str(root / path)


def _file_digest(path: str) -> str:
    """SHA-256 of a file (geometry identity for reuse, ADR 0018)."""
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check_name(kind: str, name: object) -> str:
    if not isinstance(name, str) or not name.strip():
        raise ParameterError(f"{kind} name must be a non-empty string, got {name!r}.")
    return name


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
    inputs: ClassVar[tuple[str, ...]] = ()  # fields holding input file paths (str or tuple)

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
    """Load a geometry (STL or STEP/IGES/BREP) under *name*.

    Reloading a name replaces that geometry and discards the meshes of the
    components built on it (ADR 0016).
    """

    method: ClassVar[str] = "load_geometry"
    positional: ClassVar[tuple[str, ...]] = ("path",)
    inputs: ClassVar[tuple[str, ...]] = ("path",)

    path: str
    anchor_strategy: str = "auto"
    skip_hygiene: bool = False
    name: str = "body"

    def apply(self, project: Project, ctx: Context) -> Undo:
        name = _check_name("geometry", self.name)
        body = Body.from_file(
            resolve_input(self.path),
            anchor_strategy=self.anchor_strategy,
            skip_hygiene=self.skip_hygiene,
        )
        prev = (project.geometries.get(name), project.geometry_digest(name))
        discarded = project.discard_built_on(name)
        project.set_geometry(name, body, _file_digest(resolve_input(self.path)))

        def restore() -> None:
            project.set_geometry(name, *prev)
            project.restore_components(discarded)

        return Undo(restore)


@dataclass(frozen=True)
class UseComponent(Command):
    """Make *name* the current component, creating it empty if new (ADR 0016)."""

    method: ClassVar[str] = "use_component"
    positional: ClassVar[tuple[str, ...]] = ("name",)

    name: str

    def apply(self, project: Project, ctx: Context) -> Undo:
        name = _check_name("component", self.name)
        prev = project.current_component
        created = project.select_component(name)

        def restore() -> None:
            project.select_component(prev)
            if created:
                project.drop_component(name)

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
    geometry: str = "body"

    def to_source(self) -> str:
        # *params* read better as keywords: p.build_topology('ogrid', nk=6, ni=10)
        return _params_source(self, self.params)

    def apply(self, project: Project, ctx: Context) -> Undo:
        body = project.require_geometry(self.geometry, "build_topology")
        spacing = k_spacing(
            self.k_law, ds=self.ds, growth=self.growth, beta=self.k_beta, thickness=self.k_thickness
        )
        topo, blocks = build_topology(
            body,
            self.topology,
            nk=self.nk,
            spacing=spacing,
            options=self.options,
            **self.params,
        )
        prev = (project.mesh, project.build)
        project.mesh = Mesh(body=body, blocks=blocks, topology=topo)
        project.build = BuildState(
            self.topology, self.nk, self.ds, self.growth, spacing, geometries=(self.geometry,)
        )

        def restore() -> None:
            project.mesh, project.build = prev

        return Undo(restore)


@dataclass
class BuildState:
    """What later commands need from the command that made the mesh."""

    topology: str  # topology name, "primitive:<kind>", "tip_cap" or "blocks"
    nk: int
    ds: float = 0.0
    growth: float = 0.0
    spacing: Any = None  # Spacing1D, for March
    writer: Callable[[Path], None] | None = None  # own .geo naming (primitives)
    marchable: bool = True
    geometries: tuple[str, ...] = ()  # the named geometries the mesh was built on


def _params_source(command: Command, params: dict[str, Any]) -> str:
    """``to_source`` with *params* rendered as trailing keywords."""
    head = Command.to_source(dataclasses.replace(command, params={}))  # type: ignore[type-var]
    extra = ", ".join(f"{k}={v!r}" for k, v in params.items())
    return head if not extra else f"{head[:-1]}, {extra})"


@dataclass(frozen=True)
class March(Command):
    """Extrude all blocks from the wall (``shore mesh`` stage 3)."""

    method: ClassVar[str] = "march"

    smooth_strength: float = 0.2
    smooth_iters: int = 2
    blend_normals_k: int = 0  # removed (ADR 0034): kept so old journals load without replay
    smooth_mode: str = "limited"
    march_scheme: str = "algebraic"

    def apply(self, project: Project, ctx: Context) -> Undo:
        if self.blend_normals_k:
            raise ParameterError(BLEND_REMOVED)
        mesh, build = project.require_mesh()
        if not build.marchable:
            raise ParameterError(
                f"march: the current mesh ({build.topology}) cannot be marched; "
                "march before splitting, and primitives / caps are built volumetric."
            )
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
                smooth_mode=self.smooth_mode,
                march_scheme=self.march_scheme,
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
class FillFace(Command):
    """Rebuild one block face from its boundary on a CAD surface (ADR 0007).

    The face's four boundary rows stay as they are; the interior is filled by
    ``shore.surface.tfi.surface_grid`` on SURFACE entity number *surface* of
    the loaded geometry.
    """

    method: ClassVar[str] = "surface_grid"
    positional: ClassVar[tuple[str, ...]] = ("block", "face")

    block: str
    face: str
    surface: int
    smoothing: str = "elliptic"
    geometry: str = "body"

    def apply(self, project: Project, ctx: Context) -> Undo:
        from shore.geometry import EntityKind
        from shore.mesh.face import FACE_ALONG_AXES, _face_slice
        from shore.surface.tfi import surface_grid

        mesh, _ = project.require_mesh()
        body = project.require_geometry(self.geometry, "surface_grid")
        blocks = {b.label: b for b in mesh.blocks}
        if self.block not in blocks:
            raise ParameterError(f"Unknown block {self.block!r}; blocks: {sorted(blocks)}.")
        if self.face not in FACE_ALONG_AXES:
            raise ParameterError(f"Unknown face {self.face!r}; faces: {sorted(FACE_ALONG_AXES)}.")
        geometry = body.geometry
        surfaces = geometry.entities(EntityKind.SURFACE)
        if not 0 <= self.surface < len(surfaces):
            raise ParameterError(
                f"The geometry has {len(surfaces)} parametric surface(s); "
                f"surface={self.surface} does not exist (an STL body has none)."
            )
        view = _face_slice(blocks[self.block], self.face)  # type: ignore[arg-type]
        grid = surface_grid(geometry, surfaces[self.surface], view, smoothing=self.smoothing)
        restore = project.undo_cache.snapshot(mesh)
        view[...] = grid.points
        return Undo(restore)


@dataclass(frozen=True)
class WriteGeo(Command):
    """Write ``.geo`` files (output only).

    Naming: a primitive's own rule (``shore primitive``), unless
    *per_block*; otherwise ``shore mesh``'s rule: one ``<stem>.geo`` for a
    single non-cubed-sphere block, else ``<stem>_<label>.geo`` per block
    (always the latter with *per_block*, as ``shore split`` and
    ``shore cap-tip`` write).
    """

    method: ClassVar[str] = "write_geo"
    positional: ClassVar[tuple[str, ...]] = ("stem",)
    mutates: ClassVar[bool] = False

    stem: str
    per_block: bool = False

    def apply(self, project: Project, ctx: Context) -> None:
        from shore.io.geo import write_geo

        mesh, build = project.require_mesh()
        stem = Path(self.stem)
        stem.parent.mkdir(parents=True, exist_ok=True)
        if build.writer is not None and not self.per_block:
            build.writer(stem)
        elif self.per_block or build.topology == "cubed_sphere" or len(mesh.blocks) != 1:
            for block in mesh.blocks:
                write_geo(stem.parent / f"{stem.name}_{block.label}.geo", block.nodes)
        else:
            single = stem if stem.suffix == ".geo" else stem.with_suffix(".geo")
            write_geo(single, mesh.blocks[0].nodes)


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


@dataclass(frozen=True)
class BuildPrimitive(Command):
    """Build a parametric primitive (``shore primitive``); replaces the mesh."""

    method: ClassVar[str] = "build_primitive"
    positional: ClassVar[tuple[str, ...]] = ("kind",)

    kind: str
    params: dict[str, Any] = field(default_factory=dict)

    def to_source(self) -> str:
        return _params_source(self, self.params)

    def apply(self, project: Project, ctx: Context) -> Undo:
        prim = build_primitive(self.kind, **self.params)
        prev = (project.mesh, project.build)
        project.mesh = Mesh(body=None, blocks=prim.blocks, topology=None)  # type: ignore[arg-type]
        nk = int(prim.blocks[0].nodes.shape[2])
        project.build = BuildState(
            f"primitive:{self.kind}", nk, writer=prim.write_geo, marchable=False
        )

        def restore() -> None:
            project.mesh, project.build = prev

        return Undo(restore)


@dataclass(frozen=True)
class BuildTipCap(Command):
    """Wing-tip cap from a C-grid body ``.geo`` (``shore cap-tip``); replaces the mesh."""

    method: ClassVar[str] = "build_tip_cap"
    positional: ClassVar[tuple[str, ...]] = ("body_geo",)
    inputs: ClassVar[tuple[str, ...]] = ("body_geo",)

    body_geo: str
    tip: str = "jlo"
    le_frac: float = 0.05
    te_frac: float = 0.95
    span_axis: str = "z"
    tip_inset: float | None = None
    nk: int | None = None
    n_cap_j: int = 12
    cap_j_spacing: str = "tanh"
    cap_j_beta: float = 4.0
    auto_round_te: bool = True
    te_thickness: float | None = None
    te_thickness_frac: float | None = 0.002

    def apply(self, project: Project, ctx: Context) -> Undo:
        kwargs = {f.name: getattr(self, f.name) for f in dataclasses.fields(self)}
        mesh, resolved = build_tip_cap(resolve_input(kwargs.pop("body_geo")), **kwargs)
        prev = (project.mesh, project.build, project.last_result)
        project.mesh = mesh
        project.build = BuildState("tip_cap", int(mesh.blocks[0].nodes.shape[2]), marchable=False)
        project.last_result = resolved

        def restore() -> None:
            project.mesh, project.build, project.last_result = prev

        return Undo(restore)


@dataclass(frozen=True)
class BuildCollar(Command):
    """Collar block around the junction of the loaded body and *other* (ADR 0010).

    Body A is the project's geometry, body B is loaded from the file *other*
    (its path is what the journal records). Surface steps are
    ``ds_surface * growth_surface**j`` (*m_a* on A, *m_b* on B), wall-normal
    steps ``ds_wall * growth_wall**k`` (*nk* of them). With *plane*
    (``(point, normal)``) the half-collar of a half model (ADR 0013). The
    nodes along the junction curve follow ``Spacing1D(loop_law,
    beta=loop_beta, cluster_at=loop_cluster_at)`` and *loop_curvature*
    (ADR 0020). With *split*, ``split_collar`` (ADR 0024): the loop is cut at
    sharp corners, with wakes of *wake_length* (default: the collar's width
    on A), six blocks per piece (closed loops only). *volume* is
    ``CollarParams.volume`` (ADR 0025): ``"stations"`` or ``"smoothed"``; ``None``,
    the default (``"smoothed"``; ``"stations"`` when a project saved before
    format 3 replays, ADR 0026).
    Replaces the mesh with the collar's block(s), a separate Chimera
    component.
    """

    method: ClassVar[str] = "build_collar"
    positional: ClassVar[tuple[str, ...]] = ("other",)
    inputs: ClassVar[tuple[str, ...]] = ("other",)  # when *other* is a file, not a geometry

    other: str
    seed: tuple[float, float, float]
    n_loop: int
    ds_surface: float
    growth_surface: float
    m_a: int
    m_b: int
    ds_wall: float
    growth_wall: float
    nk: int
    label: str = "collar"
    plane: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None
    geometry: str = "body"
    loop_law: str = "uniform"
    loop_beta: float = 3.0
    loop_cluster_at: float | None = None
    loop_curvature: float = 0.0
    split: bool = False
    wake_length: float | None = None
    volume: str | None = None

    def apply(self, project: Project, ctx: Context) -> Undo:
        from shore.collar import CollarParams, build_collar, split_collar
        from shore.geometry import load_geometry
        from shore.mesh.spacing import Spacing1D

        if self.loop_law not in ("uniform", "tanh", "tanh2"):
            raise ParameterError(
                f"loop_law must be 'uniform', 'tanh' or 'tanh2', got {self.loop_law!r}."
            )
        loop_spacing = (
            None
            if self.loop_law == "uniform"
            else Spacing1D(self.loop_law, beta=self.loop_beta, cluster_at=self.loop_cluster_at)  # type: ignore[arg-type]
        )
        body = project.require_geometry(self.geometry, "build_collar")
        if self.other in project.geometries:
            other = project.geometries[self.other].geometry
            built_on: tuple[str, ...] = (self.geometry, self.other)
        else:
            other = load_geometry(resolve_input(self.other))
            built_on = (self.geometry,)
        for name in ("m_a", "m_b", "nk"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ParameterError(f"{name} must be an integer >= 1, got {value!r}.")
        if not self.split and self.wake_length is not None:
            raise ParameterError("wake_length applies to split collars (split=True) only.")
        if self.split and self.plane is not None:
            raise ParameterError("split collars are closed: not with a symmetry plane.")
        common = {
            "n_loop": self.n_loop,
            "steps_a": self.ds_surface * self.growth_surface ** np.arange(self.m_a),
            "steps_b": self.ds_surface * self.growth_surface ** np.arange(self.m_b),
            "steps_k": self.ds_wall * self.growth_wall ** np.arange(self.nk),
            "label": self.label,
            "loop_spacing": loop_spacing,
            "loop_curvature": self.loop_curvature,
            "params": CollarParams(volume=self._volume()),
        }
        seed = np.asarray(self.seed, dtype=np.float64)
        if self.split:
            blocks = split_collar(
                body.geometry, other, seed, wake_length=self.wake_length, **common
            )
        else:
            blocks = [build_collar(body.geometry, other, seed, plane=self.plane, **common)]
        prev = (project.mesh, project.build)
        # no body association: k_lo lies on two bodies
        project.mesh = Mesh(body=None, blocks=blocks, topology=None)  # type: ignore[arg-type]
        project.build = BuildState("collar", self.nk + 1, marchable=False, geometries=built_on)

        def restore() -> None:
            project.mesh, project.build = prev

        return Undo(restore)

    def _volume(self) -> str:
        from shore.collar import DEFAULT

        if self.volume is not None:
            return self.volume
        return "stations" if legacy_collar_volume.get() else DEFAULT.volume


@dataclass(frozen=True)
class LoadBlocks(Command):
    """Blocks from ``.geo`` files, wired from adjacency sidecars; replaces the mesh."""

    method: ClassVar[str] = "load_blocks"
    positional: ClassVar[tuple[str, ...]] = ("geo_paths",)
    inputs: ClassVar[tuple[str, ...]] = ("geo_paths", "adjacency")

    geo_paths: tuple[str, ...]
    adjacency: tuple[str, ...] = ()

    def apply(self, project: Project, ctx: Context) -> Undo:
        blocks = load_blocks(
            tuple(map(resolve_input, self.geo_paths)), tuple(map(resolve_input, self.adjacency))
        )
        prev = (project.mesh, project.build)
        project.mesh = Mesh(body=None, blocks=blocks, topology=None)  # type: ignore[arg-type]
        project.build = BuildState("blocks", int(blocks[0].nodes.shape[2]), marchable=False)

        def restore() -> None:
            project.mesh, project.build = prev

        return Undo(restore)


@dataclass(frozen=True)
class Split(Command):
    """Split blocks by ``(label, axis, at)`` directives (``shore split``)."""

    method: ClassVar[str] = "split"
    positional: ClassVar[tuple[str, ...]] = ("splits",)

    splits: tuple[tuple[str, str, tuple[int, ...]], ...]

    def apply(self, project: Project, ctx: Context) -> Undo:
        mesh, build = project.require_mesh()
        blocks = apply_splits(mesh.blocks, list(self.splits))
        prev = (project.mesh, project.build)
        project.mesh = Mesh(body=mesh.body, blocks=blocks, topology=mesh.topology)
        project.build = dataclasses.replace(build, writer=None, marchable=False)

        def restore() -> None:
            project.mesh, project.build = prev

        return Undo(restore)


@dataclass(frozen=True)
class WriteGrd(Command):
    """Pack ``.geo`` files into an Xall ``.grd`` (``shore grd``, output only).

    Packs from the ``.geo`` files, not from memory, so the ``.grd`` carries
    the ``.geo`` print precision exactly as the file pipeline does.
    """

    method: ClassVar[str] = "write_grd"
    positional: ClassVar[tuple[str, ...]] = ("path", "geo_paths")
    mutates: ClassVar[bool] = False

    path: str
    geo_paths: tuple[str, ...]

    def apply(self, project: Project, ctx: Context) -> None:
        from shore.io.grd import write_grd

        project.last_result = write_grd(self.path, list(self.geo_paths))


@dataclass(frozen=True)
class WriteCcPar(Command):
    """Xall ``cc.par`` from an adjacency sidecar (``shore cc-par``, output only).

    *meta* maps block labels to ``BlockMeta`` overrides (``level``,
    ``group``, ``priority``, ``comment``, ``free_face_bc``); *boxes* is a
    list of marker-box JSON entries.
    """

    method: ClassVar[str] = "write_cc_par"
    positional: ClassVar[tuple[str, ...]] = ("adjacency", "path")
    mutates: ClassVar[bool] = False

    adjacency: str
    path: str
    grd_filename: str
    output_basename: str = "cc."
    free_face_bc: int = 40
    multigrid_levels: tuple[int, int] = (4, 1)
    boundary_layer_thickness: float = -0.1
    numerical_beach: float = -1.0
    meta: dict[str, dict[str, Any]] = field(default_factory=dict)
    boxes: tuple[dict[str, Any], ...] = ()

    def apply(self, project: Project, ctx: Context) -> None:
        from shore.io.adjacency import read_adjacency_json
        from shore.io.boxes import box_from_json
        from shore.io.cc_par import BlockMeta, write_cc_par

        payload = read_adjacency_json(self.adjacency)
        block_meta = {}
        for blk in payload["blocks"]:
            label = blk["label"]
            ov = self.meta.get(label, {})
            block_meta[label] = BlockMeta(
                level=ov.get("level", 1),
                group=ov.get("group", 0),
                priority=ov.get("priority", 1),
                comment=ov.get("comment", label),
                free_face_bc=ov.get("free_face_bc", self.free_face_bc),
            )
        boxes = [box_from_json(entry) for entry in self.boxes] if self.boxes else None
        project.last_result = write_cc_par(
            payload,
            self.path,
            grd_filename=self.grd_filename,
            output_basename=self.output_basename,
            multigrid_levels=self.multigrid_levels,
            boundary_layer_thickness=self.boundary_layer_thickness,
            numerical_beach=self.numerical_beach,
            block_metadata=block_meta,
            boxes=boxes,
        )


@dataclass(frozen=True)
class WriteXall(Command):
    """The Xall files of the whole run, every component (ADR 0016, 0017; output only).

    For each component with a mesh, in creation order, into
    ``directory/<component>/``: the ``.geo`` files (one per block), the
    adjacency sidecar, ``<component>.grd`` and ``cc.par``, written by the
    same commands as ``write_geo`` / ``write_adjacency`` / ``write_grd`` /
    ``write_cc_par``. Then, in *directory*, *grd* and ``cc.par`` merged from
    them (``merge_grd``, ``merge_cc_par``) and one ``proc.input`` balanced
    over the merged ``.grd`` on *np* ranks.

    The ``cc.par`` header options apply to every component (``merge_cc_par``
    requires them to agree). *components* maps a component to its options:
    ``free_face_bc`` / ``meta`` / ``boxes`` for its ``cc.par`` (as
    ``write_cc_par``), ``body`` for all its blocks in ``proc.input``, and
    ``proc_meta`` (block label -> ``proc`` / ``group`` / ``body``) per block.
    All of it is checked before anything is written.
    """

    method: ClassVar[str] = "write_xall"
    positional: ClassVar[tuple[str, ...]] = ("directory",)
    mutates: ClassVar[bool] = False
    component_keys: ClassVar[frozenset[str]] = frozenset(
        {"free_face_bc", "meta", "boxes", "body", "proc_meta"}
    )

    directory: str
    grd: str = "mesh.grd"
    np: int = 1
    components: dict[str, dict[str, Any]] = field(default_factory=dict)
    output_basename: str = "cc."
    multigrid_levels: tuple[int, int] = (4, 1)
    boundary_layer_thickness: float = -0.1
    numerical_beach: float = -1.0

    def apply(self, project: Project, ctx: Context) -> None:
        from shore.io.cc_par import merge_cc_par
        from shore.io.grd import merge_grd
        from shore.io.proc_input import write_proc_input

        meshes = project.components
        self._check(meshes)
        names = [name for name, mesh in meshes.items() if mesh is not None]
        if not names:
            raise ParameterError("write_xall: no component has a mesh.")
        out = Path(self.directory)
        grds, cc_pars = [], []
        assignments: dict[str, dict[int, int]] = {"proc": {}, "group": {}, "body": {}}
        offset = 0
        current = project.current_component
        try:
            for name in names:
                options = self.components.get(name, {})
                project.select_component(name)
                d = out / name
                WriteGeo(str(d / name), per_block=True).apply(project, ctx)
                WriteAdjacency(str(d / f"{name}.adjacency.json")).apply(project, ctx)
                blocks = project.require_mesh()[0].blocks
                geos = tuple(str(d / f"{name}_{b.label}.geo") for b in blocks)
                WriteGrd(str(d / f"{name}.grd"), geos).apply(project, ctx)
                WriteCcPar(
                    str(d / f"{name}.adjacency.json"),
                    str(d / "cc.par"),
                    grd_filename=f"{name}.grd",
                    output_basename=self.output_basename,
                    free_face_bc=options.get("free_face_bc", WriteCcPar.free_face_bc),
                    multigrid_levels=tuple(self.multigrid_levels),  # type: ignore[arg-type]
                    boundary_layer_thickness=self.boundary_layer_thickness,
                    numerical_beach=self.numerical_beach,
                    meta=dict(options.get("meta", {})),
                    boxes=tuple(options.get("boxes", ())),
                ).apply(project, ctx)
                for index, block in enumerate(blocks, start=offset + 1):
                    if "body" in options:
                        assignments["body"][index] = int(options["body"])
                    for key, value in options.get("proc_meta", {}).get(block.label, {}).items():
                        assignments[key][index] = int(value)
                offset += len(blocks)
                grds.append(d / f"{name}.grd")
                cc_pars.append(d / "cc.par")
        finally:
            project.select_component(current)
        merge_grd(grds, out / self.grd)
        merge_cc_par(cc_pars, out / "cc.par", grd_filename=self.grd)
        project.last_result = write_proc_input(
            out / self.grd,
            out / "proc.input",
            np=self.np,
            proc_assignment=assignments["proc"] or None,
            group_assignment=assignments["group"] or None,
            body_assignment=assignments["body"] or None,
        )

    def _check(self, meshes: dict[str, Mesh | None]) -> None:
        for name, options in self.components.items():
            if meshes.get(name) is None:
                raise ParameterError(
                    f"write_xall: no component {name!r} with a mesh; "
                    f"components: {[n for n, m in meshes.items() if m is not None]}."
                )
            if not isinstance(options, dict):
                raise ParameterError(f"write_xall: options of {name!r} must be a mapping.")
            unknown = set(options) - self.component_keys
            if unknown:
                raise ParameterError(
                    f"write_xall: unknown option(s) {sorted(unknown)} for {name!r}; "
                    f"known: {sorted(self.component_keys)}."
                )
            labels = {b.label for b in meshes[name].blocks}  # type: ignore[union-attr]
            for key in ("meta", "proc_meta"):
                missing = set(options.get(key, {})) - labels
                if missing:
                    raise ParameterError(
                        f"write_xall: {key} of {name!r} names unknown block(s) {sorted(missing)}."
                    )
            for label, values in options.get("proc_meta", {}).items():
                bad = set(values) - {"proc", "group", "body"}
                if bad:
                    raise ParameterError(
                        f"write_xall: proc_meta of {name!r}.{label} has unknown "
                        f"key(s) {sorted(bad)}."
                    )


@dataclass(frozen=True)
class WriteProcInput(Command):
    """Xall ``proc.input`` from a ``.grd`` (``shore proc-input``, output only).

    *meta* maps 1-based block indices (as strings) to ``proc`` / ``group``
    / ``body`` overrides.
    """

    method: ClassVar[str] = "write_proc_input"
    positional: ClassVar[tuple[str, ...]] = ("grd", "path")
    mutates: ClassVar[bool] = False

    grd: str
    path: str
    np: int = 1
    meta: dict[str, dict[str, int]] = field(default_factory=dict)
    header_comment: str = "generated by SHORE"

    def apply(self, project: Project, ctx: Context) -> None:
        from shore.io.proc_input import write_proc_input

        assignments: dict[str, dict[int, int]] = {"proc": {}, "group": {}, "body": {}}
        for key, values in self.meta.items():
            try:
                index = int(key)
            except (TypeError, ValueError) as err:
                raise ParameterError(f"proc-input meta key {key!r} is not a block index") from err
            if not isinstance(values, dict):
                raise ParameterError(f"proc-input meta for block {index} must be a mapping")
            for name in assignments:
                if name in values:
                    assignments[name][index] = int(values[name])
        project.last_result = write_proc_input(
            self.grd,
            self.path,
            np=self.np,
            proc_assignment=assignments["proc"] or None,
            group_assignment=assignments["group"] or None,
            body_assignment=assignments["body"] or None,
            header_comment=self.header_comment,
        )
