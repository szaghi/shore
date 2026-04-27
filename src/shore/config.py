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
    blend_normals_k = 0                      # STL-normal blend window (0 = off)

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
from dataclasses import dataclass, field, fields
from pathlib import Path

import tomllib

from shore.exceptions import ParameterError

CONFIG_VERSION = 1

# Closed schema: keys allowed in each TOML table.  A key not listed
# here triggers ``ParameterError`` at load time.
_TOP_KEYS = frozenset({"version", "input", "output", "grid", "march", "surface", "volume"})
_GRID_KEYS = frozenset({"ni", "nj", "nk", "theta_cap_deg", "theta_cap"})
_MARCH_KEYS = frozenset({"ds", "growth", "smooth", "smooth_iters", "blend_normals_k"})
_SURFACE_KEYS = frozenset({"i_spacing", "i_beta", "j_spacing", "j_beta"})
_VOLUME_KEYS = frozenset({"k_spacing", "k_beta", "k_thickness"})


@dataclass
class GridConfig:
    ni: int = 40
    nj: int = 60
    nk: int = 30
    theta_cap: float | None = None  # radians; resolved from theta_cap_deg if given


@dataclass
class MarchConfig:
    ds: float = 1e-3
    growth: float = 1.15
    smooth: float = 0.2
    smooth_iters: int = 2
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
        cfg.grid = GridConfig(
            ni=g.get("ni", cfg.grid.ni),
            nj=g.get("nj", cfg.grid.nj),
            nk=g.get("nk", cfg.grid.nk),
            theta_cap=_resolve_theta(g),
        )

    if "march" in raw:
        _reject_unknown_keys(raw["march"], _MARCH_KEYS, where=f"{p}:[march]")
        m = raw["march"]
        cfg.march = MarchConfig(
            ds=m.get("ds", cfg.march.ds),
            growth=m.get("growth", cfg.march.growth),
            smooth=m.get("smooth", cfg.march.smooth),
            smooth_iters=m.get("smooth_iters", cfg.march.smooth_iters),
            blend_normals_k=m.get("blend_normals_k", cfg.march.blend_normals_k),
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
