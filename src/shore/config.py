"""TOML-driven configuration for the ``shore mesh`` pipeline.

A config file is a *single-body* meshing recipe.  Every key has a default
that mirrors the corresponding ``shore mesh`` CLI flag.  CLI flags
override config-file values; config-file values override built-in
defaults.

Schema (single body, ``version = 1``)
-------------------------------------

::

    version = 1                              # required, integer

    input  = "fuselage.stl"                  # required, path to STL
    output = "fuselage.geo"                  # optional output stem

    [grid]
    ni = 60                                  # latitude  resolution
    nj = 80                                  # longitude resolution (multiple of 4)
    nk = 40                                  # wall-normal layers
    theta_cap_deg = 30.0                     # polar exclusion (degrees)

    [march]
    ds            = 1.0e-4                   # first-layer thickness
    growth        = 1.15                     # geometric growth ratio
    smooth        = 0.2                      # Laplacian blend strength
    smooth_iters  = 2                        # Laplacian sweeps per layer
    smooth_mode   = "limited"                # limited | legacy (unlimited, folds thin layers)
    march_scheme  = "algebraic"              # algebraic | hyperbolic (Steger-Chaussee)
    blend_normals_k = 0                      # deprecated (ADR 0019); keep 0

    [surface]
    i_spacing = "uniform"                    # uniform | tanh | tanh2
    i_beta    = 3.0
    j_spacing = "uniform"                    # uniform | tanh2
    j_beta    = 3.0

    [volume]
    k_spacing   = "geometric"                # geometric | tanh | tanh2
    k_beta      = 3.0
    k_thickness = 5.0                        # required for tanh / tanh2

Unknown keys raise ``ParameterError``: the schema is closed so that
typos surface immediately rather than being silently ignored.
"""

from __future__ import annotations

import sys
import tomllib
from dataclasses import dataclass, field, fields
from pathlib import Path

from shore.exceptions import ParameterError

CONFIG_VERSION = 1

# Closed schema: keys allowed in each TOML table.  A key not listed
# here triggers ``ParameterError`` at load time.
_TOP_KEYS = frozenset({"version", "kind", "input", "output", "grid", "march", "surface", "volume"})
_GRID_KEYS = frozenset({"ni", "nj", "nk", "theta_cap_deg", "theta_cap", "topology"})
_MESH_TOPOLOGIES = frozenset({"cubed_sphere", "ogrid"})
_MARCH_KEYS = frozenset(
    {"ds", "growth", "smooth", "smooth_iters", "smooth_mode", "march_scheme", "blend_normals_k"}
)
_SURFACE_KEYS = frozenset({"i_spacing", "i_beta", "j_spacing", "j_beta"})
_VOLUME_KEYS = frozenset({"k_spacing", "k_beta", "k_thickness"})

# Primitive schema (kind = "primitive")
_PRIM_TOP_KEYS = frozenset({"version", "kind", "output", "geometry", "grid", "spacing"})
_GEOMETRY_KEYS = frozenset(
    {
        "type",
        "vmin",
        "vmax",  # box
        "axis",
        "center",
        "z0",
        "z1",  # annulus / flat_caps
        "r_in",
        "r_out",  # annulus
        "theta_start",
        "theta_end",  # annulus sector
        "inner_frac",  # flat_caps
    }
)
_PRIM_GRID_KEYS = frozenset({"ni", "nj", "nk"})
_PRIM_SPACING_KEYS = frozenset({"i_law", "j_law", "k_law", "i_beta", "j_beta", "k_beta"})
_PRIM_TYPES = frozenset({"box", "annulus", "flat_caps"})

# Split schema (kind = "split")
_SPLIT_TOP_KEYS = frozenset({"version", "kind", "output", "splits"})
_SPLIT_KEYS = frozenset({"label", "axis", "at"})


@dataclass
class GridConfig:
    ni: int = 40
    nj: int = 60
    nk: int = 30
    theta_cap: float | None = None  # radians; resolved from theta_cap_deg if given
    topology: str = "cubed_sphere"  # 'cubed_sphere' (6 blocks) or 'ogrid' (1 block)


@dataclass
class MarchConfig:
    ds: float = 1e-3
    growth: float = 1.15
    smooth: float = 0.2
    smooth_iters: int = 2
    smooth_mode: str = "limited"
    march_scheme: str = "algebraic"
    blend_normals_k: int = 0


@dataclass
class SurfaceConfig:
    i_spacing: str = "uniform"
    i_beta: float = 3.0
    j_spacing: str = "uniform"
    j_beta: float = 3.0


@dataclass
class VolumeConfig:
    k_spacing: str = "geometric"
    k_beta: float = 3.0
    k_thickness: float | None = None


@dataclass
class GeometryConfig:
    """Geometric inputs for a primitive.

    The ``type`` field selects which primitive: ``"box"``, ``"annulus"``,
    or ``"flat_caps"``.  Only the keys relevant to that primitive are
    consulted; the rest are ignored *but still validated* against the
    closed schema (typos surface).
    """

    type: str = "box"
    # box
    vmin: tuple[float, float, float] | None = None
    vmax: tuple[float, float, float] | None = None
    # annulus / flat_caps shared
    axis: str | tuple[float, float, float] = "z"
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    z0: float | None = None
    z1: float | None = None
    # annulus
    r_in: float | None = None
    r_out: float | None = None
    theta_start: float = 0.0
    theta_end: float | None = None
    # flat_caps: ``r_out`` / ``z0`` / ``z1`` from above; inner_frac
    # controls the central square's corner radius as a fraction of
    # r_out (default 0.5 in FlatCapsCylinder.build when None).
    inner_frac: float | None = None


@dataclass
class PrimSpacingConfig:
    """Per-axis spacing laws for primitives.

    Box / annulus / flat_caps share the same three axes (i, j, k for
    box and annulus; i, c, k for flat_caps where ``c`` maps to the
    central-square / arc direction — written as ``j`` in the config
    for visual consistency with the other primitives).
    """

    i_law: str = "uniform"
    j_law: str = "uniform"
    k_law: str = "uniform"
    i_beta: float = 3.0
    j_beta: float = 3.0
    k_beta: float = 3.0


@dataclass
class PrimGridConfig:
    """Node counts for primitives."""

    ni: int = 40
    nj: int = 40
    nk: int = 30


@dataclass
class PrimitiveConfig:
    """Resolved configuration for a single ``shore primitive`` invocation."""

    output: Path | None = None
    geometry: GeometryConfig = field(default_factory=GeometryConfig)
    grid: PrimGridConfig = field(default_factory=PrimGridConfig)
    spacing: PrimSpacingConfig = field(default_factory=PrimSpacingConfig)


@dataclass
class SplitDirective:
    """One split directive parsed from a TOML split config.

    Mirrors :class:`shore.split.BlockSplit`; kept separate so the
    config layer doesn't import the split module (avoids a circular
    import for users who only want config validation).
    """

    label: str
    axis: str
    at: tuple[int, ...]


@dataclass
class SplitConfig:
    """Resolved configuration for a single ``shore split`` invocation."""

    output: Path | None = None
    splits: list[SplitDirective] = field(default_factory=list)


@dataclass
class MeshConfig:
    """Resolved configuration for a single ``shore mesh`` invocation.

    All sections are present with defaults; ``input`` is required.
    """

    input: Path | None = None
    output: Path | None = None
    grid: GridConfig = field(default_factory=GridConfig)
    march: MarchConfig = field(default_factory=MarchConfig)
    surface: SurfaceConfig = field(default_factory=SurfaceConfig)
    volume: VolumeConfig = field(default_factory=VolumeConfig)


def load_config(path: str | Path) -> MeshConfig:
    """Load and validate a SHORE TOML config.

    Parameters
    ----------
    path :
        Path to the ``.toml`` file.

    Returns
    -------
    MeshConfig
        Fully populated config; missing keys take their dataclass defaults.

    Raises
    ------
    ParameterError
        On file-not-found, malformed TOML, missing required keys,
        unknown keys, or version mismatch.
    """
    p = Path(path)
    if not p.exists():
        raise ParameterError(f"Config file not found: {p}")

    try:
        with p.open("rb") as f:
            raw = tomllib.load(f)
    except tomllib.TOMLDecodeError as exc:
        raise ParameterError(f"Invalid TOML in {p}: {exc}") from exc

    _reject_unknown_keys(raw, _TOP_KEYS, where=str(p))

    version = raw.get("version")
    if version is None:
        raise ParameterError(f"{p}: missing required key 'version' (current: {CONFIG_VERSION})")
    if version != CONFIG_VERSION:
        raise ParameterError(
            f"{p}: unsupported config version {version!r} (this build supports {CONFIG_VERSION})"
        )

    kind = raw.get("kind", "mesh")
    if kind != "mesh":
        raise ParameterError(f"{p}: kind={kind!r} is not a mesh config (use load_primitive_config)")

    if "input" not in raw:
        raise ParameterError(f"{p}: missing required key 'input'")

    cfg = MeshConfig()
    cfg.input = _resolve_path(raw["input"], base=p.parent)
    if "output" in raw:
        cfg.output = _resolve_path(raw["output"], base=p.parent)

    if "grid" in raw:
        _reject_unknown_keys(raw["grid"], _GRID_KEYS, where=f"{p}:[grid]")
        g = raw["grid"]
        if "theta_cap" in g and "theta_cap_deg" in g:
            raise ParameterError(
                f"{p}:[grid]: 'theta_cap' and 'theta_cap_deg' are mutually exclusive"
            )
        topology = g.get("topology", cfg.grid.topology)
        if topology not in _MESH_TOPOLOGIES:
            raise ParameterError(
                f"{p}:[grid]: unknown topology {topology!r}; "
                f"expected one of {sorted(_MESH_TOPOLOGIES)!r}"
            )
        cfg.grid = GridConfig(
            ni=g.get("ni", cfg.grid.ni),
            nj=g.get("nj", cfg.grid.nj),
            nk=g.get("nk", cfg.grid.nk),
            theta_cap=_resolve_theta(g),
            topology=topology,
        )

    if "march" in raw:
        _reject_unknown_keys(raw["march"], _MARCH_KEYS, where=f"{p}:[march]")
        m = raw["march"]
        cfg.march = MarchConfig(
            ds=m.get("ds", cfg.march.ds),
            growth=m.get("growth", cfg.march.growth),
            smooth=m.get("smooth", cfg.march.smooth),
            smooth_iters=m.get("smooth_iters", cfg.march.smooth_iters),
            smooth_mode=m.get("smooth_mode", cfg.march.smooth_mode),
            march_scheme=m.get("march_scheme", cfg.march.march_scheme),
            blend_normals_k=m.get("blend_normals_k", cfg.march.blend_normals_k),
        )
        if cfg.march.smooth_mode not in ("limited", "legacy"):
            raise ParameterError(
                f"{p}:[march] smooth_mode must be 'limited' or 'legacy',"
                f" got {cfg.march.smooth_mode!r}"
            )

    if "surface" in raw:
        _reject_unknown_keys(raw["surface"], _SURFACE_KEYS, where=f"{p}:[surface]")
        s = raw["surface"]
        cfg.surface = SurfaceConfig(
            i_spacing=s.get("i_spacing", cfg.surface.i_spacing),
            i_beta=s.get("i_beta", cfg.surface.i_beta),
            j_spacing=s.get("j_spacing", cfg.surface.j_spacing),
            j_beta=s.get("j_beta", cfg.surface.j_beta),
        )

    if "volume" in raw:
        _reject_unknown_keys(raw["volume"], _VOLUME_KEYS, where=f"{p}:[volume]")
        v = raw["volume"]
        cfg.volume = VolumeConfig(
            k_spacing=v.get("k_spacing", cfg.volume.k_spacing),
            k_beta=v.get("k_beta", cfg.volume.k_beta),
            k_thickness=v.get("k_thickness", cfg.volume.k_thickness),
        )

    return cfg


def _reject_unknown_keys(d: dict, allowed: frozenset[str], where: str) -> None:
    extra = set(d) - allowed
    if extra:
        raise ParameterError(
            f"{where}: unknown key(s) {sorted(extra)!r}; allowed: {sorted(allowed)!r}"
        )


def _resolve_path(value: object, base: Path) -> Path:
    if not isinstance(value, str):
        raise ParameterError(f"path must be a string, got {type(value).__name__}: {value!r}")
    pth = Path(value)
    return pth if pth.is_absolute() else (base / pth).resolve()


def _resolve_theta(grid: dict) -> float | None:
    if "theta_cap" in grid:
        return float(grid["theta_cap"])
    if "theta_cap_deg" in grid:
        import math

        return math.radians(float(grid["theta_cap_deg"]))
    return None


def dump_template(stream=None) -> None:
    """Write an annotated single-body TOML template to *stream* (default: stdout).

    The template uses the dataclass defaults so it stays in sync with
    the code automatically.
    """
    if stream is None:
        stream = sys.stdout

    g, m, s, v = GridConfig(), MarchConfig(), SurfaceConfig(), VolumeConfig()

    def emit(line: str = "") -> None:
        print(line, file=stream)

    emit(f"# SHORE meshing config — version {CONFIG_VERSION}")
    emit("# CLI flags override these values; remove a key to fall back to its built-in default.")
    emit("")
    emit(f"version = {CONFIG_VERSION}")
    emit("")
    emit('input  = "body.stl"      # required: path to the input STL (relative to this file)')
    emit('# output = "body.geo"    # optional: output stem for .geo files')
    emit("")
    emit("[grid]")
    for f in fields(g):
        if f.name == "theta_cap":
            emit("# theta_cap_deg = 30.0     # polar exclusion in degrees")
            emit("# theta_cap     = 0.5236   # alternative: in radians (mutually exclusive)")
            continue
        emit(f"{f.name} = {getattr(g, f.name)!r}")
    emit("")
    emit("[march]")
    for f in fields(m):
        emit(f"{f.name} = {getattr(m, f.name)!r}")
    emit("")
    emit("[surface]")
    for f in fields(s):
        emit(f"{f.name} = {getattr(s, f.name)!r}")
    emit("")
    emit("[volume]")
    for f in fields(v):
        if f.name == "k_thickness":
            emit("# k_thickness = 5.0       # required when k_spacing is 'tanh' or 'tanh2'")
            continue
        emit(f"{f.name} = {getattr(v, f.name)!r}")


# ---------------------------------------------------------------------------
# Primitive config loader
# ---------------------------------------------------------------------------


def _coerce_axis(value: object) -> str | tuple[float, float, float]:
    """Validate an axis spec coming out of TOML.

    TOML allows:
    - a string ``"x"`` / ``"y"`` / ``"z"``
    - a string ``"1.0 0.0 0.0"`` (whitespace-separated triple)
    - an array ``[1.0, 0.0, 0.0]``

    The string forms are passed through verbatim; arrays are returned
    as 3-tuples for ``resolve_axis`` to consume.
    """
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple)) and len(value) == 3:
        try:
            return (float(value[0]), float(value[1]), float(value[2]))
        except (TypeError, ValueError) as exc:
            raise ParameterError(f"axis: cannot parse {value!r} as a 3-vector") from exc
    raise ParameterError(f"axis: expected 'x'/'y'/'z' or 3-component array, got {value!r}")


def _coerce_3vec(value: object, name: str) -> tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ParameterError(f"{name}: expected 3-component array, got {value!r}")
    try:
        return (float(value[0]), float(value[1]), float(value[2]))
    except (TypeError, ValueError) as exc:
        raise ParameterError(f"{name}: components must be numeric, got {value!r}") from exc


def load_primitive_config(path: str | Path) -> PrimitiveConfig:
    """Load and validate a SHORE TOML primitive config.

    A primitive config is identified by ``kind = "primitive"`` at the
    top level.  Mesh and primitive configs share the same versioning
    but otherwise have disjoint schemas.

    Parameters
    ----------
    path :
        Path to the ``.toml`` file.

    Returns
    -------
    PrimitiveConfig
        Fully populated config; missing keys take their dataclass defaults.

    Raises
    ------
    ParameterError
        On file-not-found, malformed TOML, missing required keys,
        unknown keys, version mismatch, or unknown geometry type.
    """
    p = Path(path)
    if not p.exists():
        raise ParameterError(f"Config file not found: {p}")

    try:
        with p.open("rb") as f:
            raw = tomllib.load(f)
    except tomllib.TOMLDecodeError as exc:
        raise ParameterError(f"Invalid TOML in {p}: {exc}") from exc

    _reject_unknown_keys(raw, _PRIM_TOP_KEYS, where=str(p))

    version = raw.get("version")
    if version is None:
        raise ParameterError(f"{p}: missing required key 'version' (current: {CONFIG_VERSION})")
    if version != CONFIG_VERSION:
        raise ParameterError(
            f"{p}: unsupported config version {version!r} (this build supports {CONFIG_VERSION})"
        )

    kind = raw.get("kind")
    if kind != "primitive":
        raise ParameterError(
            f'{p}: primitive configs require kind = "primitive" (got kind={kind!r})'
        )

    cfg = PrimitiveConfig()
    if "output" in raw:
        cfg.output = _resolve_path(raw["output"], base=p.parent)

    if "geometry" in raw:
        _reject_unknown_keys(raw["geometry"], _GEOMETRY_KEYS, where=f"{p}:[geometry]")
        g = raw["geometry"]
        gtype = g.get("type", cfg.geometry.type)
        if gtype not in _PRIM_TYPES:
            raise ParameterError(
                f"{p}:[geometry]: unknown type {gtype!r}; expected one of {sorted(_PRIM_TYPES)!r}"
            )
        cfg.geometry = GeometryConfig(
            type=gtype,
            vmin=_coerce_3vec(g["vmin"], "vmin") if "vmin" in g else None,
            vmax=_coerce_3vec(g["vmax"], "vmax") if "vmax" in g else None,
            axis=_coerce_axis(g["axis"]) if "axis" in g else cfg.geometry.axis,
            center=_coerce_3vec(g["center"], "center") if "center" in g else cfg.geometry.center,
            z0=g.get("z0"),
            z1=g.get("z1"),
            r_in=g.get("r_in"),
            r_out=g.get("r_out"),
            theta_start=g.get("theta_start", cfg.geometry.theta_start),
            theta_end=g.get("theta_end"),
            inner_frac=g.get("inner_frac"),
        )

    if "grid" in raw:
        _reject_unknown_keys(raw["grid"], _PRIM_GRID_KEYS, where=f"{p}:[grid]")
        gg = raw["grid"]
        cfg.grid = PrimGridConfig(
            ni=gg.get("ni", cfg.grid.ni),
            nj=gg.get("nj", cfg.grid.nj),
            nk=gg.get("nk", cfg.grid.nk),
        )

    if "spacing" in raw:
        _reject_unknown_keys(raw["spacing"], _PRIM_SPACING_KEYS, where=f"{p}:[spacing]")
        sp = raw["spacing"]
        cfg.spacing = PrimSpacingConfig(
            i_law=sp.get("i_law", cfg.spacing.i_law),
            j_law=sp.get("j_law", cfg.spacing.j_law),
            k_law=sp.get("k_law", cfg.spacing.k_law),
            i_beta=sp.get("i_beta", cfg.spacing.i_beta),
            j_beta=sp.get("j_beta", cfg.spacing.j_beta),
            k_beta=sp.get("k_beta", cfg.spacing.k_beta),
        )

    return cfg


def dump_primitive_template(stream=None, type_: str = "box") -> None:
    """Write an annotated primitive TOML template to *stream*.

    *type_* is one of ``"box"``, ``"annulus"``, ``"flat_caps"`` and
    determines which geometry block is emitted.
    """
    if stream is None:
        stream = sys.stdout
    if type_ not in _PRIM_TYPES:
        raise ParameterError(
            f"dump_primitive_template: unknown type {type_!r}; "
            f"expected one of {sorted(_PRIM_TYPES)!r}"
        )

    grid, sp = PrimGridConfig(), PrimSpacingConfig()

    def emit(line: str = "") -> None:
        print(line, file=stream)

    emit(f"# SHORE primitive config — version {CONFIG_VERSION}")
    emit("# CLI flags override these values; remove a key to fall back to its built-in default.")
    emit("")
    emit(f"version = {CONFIG_VERSION}")
    emit('kind    = "primitive"')
    emit("")
    emit('# output = "background.geo"   # optional: where to write the .geo file(s)')
    emit("")
    emit("[geometry]")
    emit(f'type = "{type_}"')
    if type_ == "box":
        emit("vmin = [0.0, 0.0, 0.0]")
        emit("vmax = [10.0, 5.0, 3.0]")
    elif type_ == "annulus":
        emit("r_in   = 0.5")
        emit("r_out  = 2.0")
        emit("z0     = 0.0")
        emit("z1     = 4.0")
        emit('axis   = "z"          # "x" / "y" / "z" or [vx, vy, vz]')
        emit("center = [0.0, 0.0, 0.0]")
        emit("# theta_start = 0.0          # partial-sector start (radians)")
        emit("# theta_end   = 3.141593     # partial-sector end; omit for full annulus")
    else:  # flat_caps
        emit("r_out      = 2.0")
        emit("z0         = 0.0")
        emit("z1         = 4.0")
        emit('axis       = "z"      # "x" / "y" / "z" or [vx, vy, vz]')
        emit("center     = [0.0, 0.0, 0.0]")
        emit("inner_frac = 0.5      # central square corner radius / r_out (in (0, 1))")
    emit("")
    emit("[grid]")
    for f in fields(grid):
        emit(f"{f.name} = {getattr(grid, f.name)!r}")
    emit("")
    emit("[spacing]")
    for f in fields(sp):
        emit(f"{f.name} = {getattr(sp, f.name)!r}")


# ---------------------------------------------------------------------------
# Split config loader
# ---------------------------------------------------------------------------


def load_split_config(path: str | Path) -> SplitConfig:
    """Load and validate a SHORE TOML split config.

    A split config is identified by ``kind = "split"`` at the top
    level.  Schema (single example)::

        version = 1
        kind    = "split"

        # Optional output stem; the CLI supplies the rest.
        # output = "split_out/mesh"

        [[splits]]
        label = "sub0"
        axis  = "k"
        at    = [4]

        [[splits]]
        label = "background_center"
        axis  = "i"
        at    = [8, 16, 24]

    Returns
    -------
    SplitConfig
        Fully populated config.

    Raises
    ------
    ParameterError
        On file-not-found, malformed TOML, missing required keys,
        unknown keys, or version / kind mismatch.
    """
    p = Path(path)
    if not p.exists():
        raise ParameterError(f"Config file not found: {p}")

    try:
        with p.open("rb") as f:
            raw = tomllib.load(f)
    except tomllib.TOMLDecodeError as exc:
        raise ParameterError(f"Invalid TOML in {p}: {exc}") from exc

    _reject_unknown_keys(raw, _SPLIT_TOP_KEYS, where=str(p))

    version = raw.get("version")
    if version is None:
        raise ParameterError(f"{p}: missing required key 'version' (current: {CONFIG_VERSION})")
    if version != CONFIG_VERSION:
        raise ParameterError(
            f"{p}: unsupported config version {version!r} (this build supports {CONFIG_VERSION})"
        )

    kind = raw.get("kind")
    if kind != "split":
        raise ParameterError(f'{p}: split configs require kind = "split" (got kind={kind!r})')

    cfg = SplitConfig()
    if "output" in raw:
        cfg.output = _resolve_path(raw["output"], base=p.parent)

    raw_splits = raw.get("splits")
    if not isinstance(raw_splits, list) or not raw_splits:
        raise ParameterError(f"{p}: a split config must define at least one [[splits]] entry")

    for i, entry in enumerate(raw_splits):
        if not isinstance(entry, dict):
            raise ParameterError(f"{p}:[[splits]][{i}]: must be a TOML table")
        _reject_unknown_keys(entry, _SPLIT_KEYS, where=f"{p}:[[splits]][{i}]")
        for required in ("label", "axis", "at"):
            if required not in entry:
                raise ParameterError(f"{p}:[[splits]][{i}]: missing required key {required!r}")
        label = entry["label"]
        axis = entry["axis"]
        at = entry["at"]
        if not isinstance(label, str):
            raise ParameterError(f"{p}:[[splits]][{i}]: 'label' must be a string, got {label!r}")
        if axis not in ("i", "j", "k"):
            raise ParameterError(
                f"{p}:[[splits]][{i}]: 'axis' must be 'i' / 'j' / 'k', got {axis!r}"
            )
        if not isinstance(at, list) or not at:
            raise ParameterError(f"{p}:[[splits]][{i}]: 'at' must be a non-empty list, got {at!r}")
        try:
            cuts = tuple(int(x) for x in at)
        except (TypeError, ValueError) as exc:
            raise ParameterError(
                f"{p}:[[splits]][{i}]: 'at' entries must be integers, got {at!r}"
            ) from exc
        cfg.splits.append(SplitDirective(label=label, axis=axis, at=cuts))

    return cfg
