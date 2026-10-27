"""SHORE command-line interface."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import numpy as np
import typer
from rich.console import Console
from rich.markup import escape as _esc
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from shore.config import (
    MeshConfig,
    PrimitiveConfig,
    dump_template,
    load_config,
    load_primitive_config,
    load_split_config,
)
from shore.mesh import (
    Body,
    Spacing1D,
    check_geo,
    export_geo,
    export_surface,
    info_geo,
    info_surface,
    view_geo_file,
    view_surface_file,
)
from shore.project import Project

app = typer.Typer(
    name="shore",
    help="SHORE - Structured Hyperbolic Overlapping Remeshing Engine.",
    add_completion=False,
    no_args_is_help=True,
)

_console = Console()
_err = Console(stderr=True)


def _version_callback(value: bool) -> None:
    if value:
        from shore import __version__

        _console.print(f"shore {__version__}")
        raise typer.Exit()


@app.callback()
def _main(
    version: Annotated[
        bool | None,
        typer.Option(
            "--version",
            help="Show the SHORE version and exit.",
            callback=_version_callback,
            is_eager=True,
        ),
    ] = None,
) -> None:
    """SHORE root callback — only handles --version for now."""


# ---------------------------------------------------------------------------
# Shared option aliases
# ---------------------------------------------------------------------------

# SHORE convention: --ni, --nj, --nk are *node* counts (the array shape
# of the (ni, nj, nk, 3) coordinate block).  A block with --nk 30 has
# 30 wall-normal layers including the wall (k=0), i.e. 29 cells extruded
# from the wall.  The .grd binary header is the only place where cell
# counts (= node count - 1) appear; the conversion is automatic.
_NI = Annotated[int, typer.Option(help="Latitude node count (i-axis of the surface).")]
_NJ = Annotated[
    int,
    typer.Option(
        help="Longitude node count (j-axis of the surface).  For "
        "'cubed_sphere' this is the full equator-ring node count and must "
        "be a multiple of 4: it is split into 4 sub-blocks of nj//4 cells "
        "each (nj//4 + 1 nodes per sub-block), and the polar caps are also "
        "nj//4 cells per side.  E.g. --nj 128 -> 32 cells per sub-block.  "
        "For 'ogrid' it is the single block's j-node count directly."
    ),
]
_NK = Annotated[int, typer.Option(help="Wall-normal node count (includes the wall at k=0).")]
_DS = Annotated[float, typer.Option(help="First-layer thickness (same units as the STL).")]
_GROWTH = Annotated[float, typer.Option(help="Geometric growth ratio between successive layers.")]
_SMOOTH = Annotated[float, typer.Option(help="Laplacian smooth strength per layer [0-1].")]
_SITERS = Annotated[int, typer.Option(help="Laplacian sweeps per layer.")]
_JOURNAL = Annotated[
    Path | None,
    typer.Option("--journal", help="Also write the run as a replayable journal (`shore run`)."),
]
_PROJECT = Annotated[
    Path | None,
    typer.Option(
        "--project",
        help="Add the result to the project in this directory (created if missing; ADR 0018).",
    ),
]
_COMPONENT = Annotated[
    str | None,
    typer.Option("--component", help="Component of --project to build (default: see command)."),
]


def _open_project(directory: Path | None, component: str | None, default: str) -> Project:
    """The project in *directory* (or a new one) with component selected.

    Without *directory*: a throwaway project, as before ADR 0018.
    """
    from shore.exceptions import MeshInputError, ParameterError

    if directory is None:
        if component is not None:
            raise ParameterError("--component needs --project.")
        return Project()
    if (directory / "project.toml").exists():
        proj = Project.load(directory)
    elif directory.exists() and any(directory.iterdir()):
        raise MeshInputError(f"{directory} exists and is not a SHORE project.")
    else:
        proj = Project()
    name = component or default
    if proj.current_component != name:
        proj.use_component(name)
    return proj


def _ensure_geometry(
    proj: Project, path: Path, name: str | None, *, skip_hygiene: bool = False
) -> str:
    """Load *path* as geometry *name* (default: its stem) unless already loaded.

    The same name loaded from a file with the same SHA-256 is reused, so the
    components built on it survive. A default (stem) name already taken by
    a different file is an error: replacing it would discard its components.
    """
    from shore.exceptions import ParameterError
    from shore.project.commands import _file_digest

    key = name or path.stem
    if key in proj.geometries:
        if proj.geometry_digest(key) == _file_digest(str(path)):
            return key
        if name is None:
            raise ParameterError(
                f"geometry {key!r} is already loaded from another file; "
                "give this one an explicit name (--name-a / --name-b / --geometry-name)."
            )
    proj.load_geometry(path, name=key, skip_hygiene=skip_hygiene)
    return key


def _save_project(proj: Project, directory: Path | None) -> None:
    if directory is not None:
        proj.save(directory)
        _console.print(
            f"[green]✓[/green] Saved project {_esc(str(directory))} "
            f"(component '{_esc(proj.current_component)}')"
        )


_BLEND_K = Annotated[
    int,
    typer.Option(
        "--blend-normals-k",
        help="Deprecated (ADR 0019): layers over which to blend body normals into the "
        "marching direction (0 = off). No measured benefit; refused on C-grids.",
    ),
]
_VERBOSE = Annotated[bool, typer.Option("--verbose", "-v", help="Show rich progress bar.")]

_THETA_CAP = Annotated[
    float | None,
    typer.Option(
        "--theta-cap", help="Polar exclusion angle in radians [default 0.05*pi (~9 deg)]."
    ),
]
_THETA_CAP_DEG = Annotated[
    float | None,
    typer.Option(
        "--theta-cap-deg",
        help="Polar exclusion angle in degrees. Mutually exclusive with --theta-cap.",
    ),
]
_CAP_PULL = Annotated[
    float | None,
    typer.Option(
        "--cap-pull",
        help="[DEPRECATED, Phase 4] Ignored - cap interior is now built by "
        "Shirley + stereographic projection with C1 seam continuity.",
        hidden=True,
    ),
]
_CAP_SQUARENESS = Annotated[
    float | None,
    typer.Option(
        "--cap-squareness",
        help="[DEPRECATED] Ignored - cap k=0 is now built by gnomonic "
        "flat-square projection, which produces square cells by construction.",
        hidden=True,
    ),
]

_SURFACE_I_SPACING = Annotated[
    str,
    typer.Option(
        "--surface-i-spacing",
        help="Latitude spacing: 'uniform' (default), 'tanh', 'tanh2'.",
    ),
]
_SURFACE_I_BETA = Annotated[
    float,
    typer.Option(
        "--surface-i-beta",
        help="Clustering strength for tanh/tanh2 i-spacing (default 3.0).",
    ),
]
_SURFACE_J_SPACING = Annotated[
    str,
    typer.Option(
        "--surface-j-spacing",
        help="Longitude spacing: 'uniform' (default), 'tanh2'.",
    ),
]
_SURFACE_J_BETA = Annotated[
    float,
    typer.Option(
        "--surface-j-beta",
        help="Clustering strength for tanh2 j-spacing (default 3.0).",
    ),
]
_J_CLUSTER_LON = Annotated[
    float | None,
    typer.Option(
        "--j-cluster-lon",
        help="Longitude (radians) to cluster toward with tanh2 j-spacing.",
    ),
]
_VOLUME_K_SPACING = Annotated[
    str,
    typer.Option(
        "--volume-k-spacing",
        help="Wall-normal spacing: 'geometric' (default), 'tanh', 'tanh2'.",
    ),
]
_VOLUME_K_BETA = Annotated[
    float,
    typer.Option(
        "--volume-k-beta",
        help="Clustering strength for tanh/tanh2 k-spacing (default 3.0).",
    ),
]
_VOLUME_K_THICKNESS = Annotated[
    float | None,
    typer.Option(
        "--volume-k-thickness",
        help="Total wall-normal thickness. Required for tanh/tanh2 k-spacing.",
    ),
]

# C-grid topology flags (only used when --topology cgrid).
_NI_BODY = Annotated[
    int | None,
    typer.Option(
        "--ni-body",
        help="[cgrid] Per-station body wrap node count (upper-TE → LE → lower-TE).",
    ),
]
_NI_WAKE = Annotated[
    int | None,
    typer.Option(
        "--ni-wake",
        help="[cgrid] Per-branch wake node count (TE node included).",
    ),
]
_N_STATIONS = Annotated[
    int | None,
    typer.Option(
        "--n-stations",
        help="[cgrid] Number of span-axis stations.",
    ),
]
_WAKE_LENGTH = Annotated[
    float | None,
    typer.Option(
        "--wake-length",
        help="[cgrid] Length of each wake branch from the TE, in body units.",
    ),
]
_WAKE_GROWTH = Annotated[
    float,
    typer.Option(
        "--wake-growth",
        help="[cgrid] Geometric growth ratio for wake i-spacing (default 1.15).",
    ),
]
_SPAN_AXIS = Annotated[
    str,
    typer.Option(
        "--span-axis",
        help="[cgrid] Span axis: 'x', 'y', or 'z' (default 'z').",
    ),
]
_TE_DIHEDRAL_DEG = Annotated[
    float,
    typer.Option(
        "--te-dihedral-deg",
        help="[cgrid] Sharpness threshold (degrees) for TE detection (default 30.0).",
    ),
]
_CGRID_BLOCKS = Annotated[
    int,
    typer.Option(
        "--cgrid-blocks",
        help=(
            "[cgrid] Block layout: 3 (default, Xall-correct sharp-TE: "
            "upper / body / lower with SHARED seams), 4 (Xall-correct "
            "blunt-TE: 3 + a wake-base block W bridging the TE gap), "
            "or 1 (legacy single-block, kept for opt-in compatibility — "
            "Xall cannot apply BCs cleanly on its wake-cut overlap)."
        ),
    ),
]
_CAP_TIPS = Annotated[
    str,
    typer.Option(
        "--cap-tips",
        help=(
            "[cgrid] Conformal wing-tip caps.  'none' (default), "
            "'jlo' (j=0 station only), 'jhi' (j=n_stations-1 only), "
            "or 'both'.  When set, the C-grid mesh is augmented with "
            "5 cap blocks per requested tip + 1 cap_wake_base per tip "
            "(total +6 blocks per tip, since the C-grid auto-promotes "
            "to 4-block).  Sharp-TE bodies are auto-rounded over their "
            "full span by --cap-te-thickness-frac (default 0.2% chord) "
            "so the cap-body interface is bit-conformal."
        ),
    ),
]
_CAP_TE_THICKNESS = Annotated[
    float | None,
    typer.Option(
        "--cap-te-thickness",
        help=(
            "[cgrid + caps] Absolute TE thickness for body-rounding "
            "(world units).  Mutually exclusive with "
            "--cap-te-thickness-frac."
        ),
    ),
]
_CAP_TE_THICKNESS_FRAC = Annotated[
    float,
    typer.Option(
        "--cap-te-thickness-frac",
        help=(
            "[cgrid + caps] TE thickness as a fraction of chord "
            "(default 0.002 = 0.2%).  Used only on sharp-TE bodies."
        ),
    ),
]
_CAP_LE_FRAC = Annotated[
    float,
    typer.Option("--cap-le-frac", help="[cgrid + caps] LE-cut chord fraction (default 0.05)."),
]
_CAP_TE_FRAC = Annotated[
    float,
    typer.Option("--cap-te-frac", help="[cgrid + caps] TE-cut chord fraction (default 0.95)."),
]
_CAP_TIP_INSET = Annotated[
    float | None,
    typer.Option(
        "--cap-tip-inset",
        help=("[cgrid + caps] Spanwise wall thickness of the cap.  Default = 1e-3 * chord."),
    ),
]
_CAP_N_CAP_J = Annotated[
    int,
    typer.Option(
        "--cap-n-cap-j",
        help="[cgrid + caps] Cap j-axis (airfoil-thickness) node count (default 12).",
    ),
]
_CAP_J_SPACING = Annotated[
    str,
    typer.Option(
        "--cap-j-spacing-cgrid",
        help=("[cgrid + caps] Cap j-axis spacing law: uniform / tanh / tanh2 (default tanh)."),
    ),
]
_CAP_J_BETA = Annotated[
    float,
    typer.Option(
        "--cap-j-beta-cgrid",
        help="[cgrid + caps] Clustering strength for tanh / tanh2 (default 4.0).",
    ),
]
_CAP_NK = Annotated[
    int | None,
    typer.Option(
        "--cap-nk",
        help=("[cgrid + caps] Cap k-axis (spanwise) node count.  Default = matches --nk."),
    ),
]


# CH topology flags (only used when --topology ch).
_CH_CHORD_MIN = Annotated[
    float | None,
    typer.Option(
        "--ch-chord-min",
        help="[ch] Far-field box minimum chord coordinate (world units).",
    ),
]
_CH_CHORD_MAX = Annotated[
    float | None,
    typer.Option(
        "--ch-chord-max",
        help="[ch] Far-field box maximum chord coordinate (world units).",
    ),
]
_CH_NORMAL_MIN = Annotated[
    float | None,
    typer.Option(
        "--ch-normal-min",
        help="[ch] Far-field box minimum normal coordinate (world units).",
    ),
]
_CH_NORMAL_MAX = Annotated[
    float | None,
    typer.Option(
        "--ch-normal-max",
        help="[ch] Far-field box maximum normal coordinate (world units).",
    ),
]
_NK_CHANNEL = Annotated[
    int,
    typer.Option(
        "--nk-channel",
        help="[ch/oh] Total channel-block layers including shared seam (default 5).",
    ),
]
_K_LAW_CHANNEL = Annotated[
    str,
    typer.Option(
        "--k-law-channel",
        help="[ch/oh] Channel k-axis spacing law: uniform / tanh / tanh2 (default uniform).",
    ),
]
_K_BETA_CHANNEL = Annotated[
    float,
    typer.Option(
        "--k-beta-channel",
        help="[ch/oh] Clustering strength for tanh / tanh2 channel spacing (default 3.0).",
    ),
]
_C1_SEAM_CHANNEL = Annotated[
    bool,
    typer.Option(
        "--c1-seam-channel/--no-c1-seam-channel",
        help=(
            "[ch/oh] Pin the channel's first step at every (i, j) to the wrap's "
            "last marching step — bit-exact C1 spacing match across the "
            "wrap-channel seam (default on)."
        ),
    ),
]

# CH frame/outflow knobs not previously exposed.  freestream_world is a
# 3-vector (exposed as three scalars); i_split_body_* and ni_outflow are ints.
_CH_FREESTREAM_X = Annotated[
    float | None,
    typer.Option("--ch-freestream-x", help="[ch] Freestream direction x component (world)."),
]
_CH_FREESTREAM_Y = Annotated[
    float | None,
    typer.Option("--ch-freestream-y", help="[ch] Freestream direction y component (world)."),
]
_CH_FREESTREAM_Z = Annotated[
    float | None,
    typer.Option("--ch-freestream-z", help="[ch] Freestream direction z component (world)."),
]
_CH_I_SPLIT_UPPER = Annotated[
    int | None,
    typer.Option(
        "--ch-i-split-body-upper",
        help="[ch] Body-wrap i index where the upper channel segment splits (auto if unset).",
    ),
]
_CH_I_SPLIT_LOWER = Annotated[
    int | None,
    typer.Option(
        "--ch-i-split-body-lower",
        help="[ch] Body-wrap i index where the lower channel segment splits (auto if unset).",
    ),
]
_CH_NI_OUTFLOW = Annotated[
    int,
    typer.Option("--ch-ni-outflow", help="[ch] Outflow sub-block i resolution (default 8)."),
]

# OH far-field box extents (axis-aligned, world coords).  OHTopology takes
# box_vmin / box_vmax as 3-vectors; the CLI exposes them as six scalars to
# match the scalar-flag house style used by the CH frame.
_OH_BOX_XMIN = Annotated[
    float | None,
    typer.Option("--oh-box-xmin", help="[oh] Far-field box minimum x (world units)."),
]
_OH_BOX_XMAX = Annotated[
    float | None,
    typer.Option("--oh-box-xmax", help="[oh] Far-field box maximum x (world units)."),
]
_OH_BOX_YMIN = Annotated[
    float | None,
    typer.Option("--oh-box-ymin", help="[oh] Far-field box minimum y (world units)."),
]
_OH_BOX_YMAX = Annotated[
    float | None,
    typer.Option("--oh-box-ymax", help="[oh] Far-field box maximum y (world units)."),
]
_OH_BOX_ZMIN = Annotated[
    float | None,
    typer.Option("--oh-box-zmin", help="[oh] Far-field box minimum z (world units)."),
]
_OH_BOX_ZMAX = Annotated[
    float | None,
    typer.Option("--oh-box-zmax", help="[oh] Far-field box maximum z (world units)."),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _progress() -> Progress:
    return Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]{task.description}"),
        BarColumn(bar_width=32),
        MofNCompleteColumn(),
        TextColumn("[dim]{task.fields[info]}"),
        TimeElapsedColumn(),
        console=_console,
        transient=False,
    )


def _default_output(input_path: Path, suffix: str) -> Path:
    return input_path.parent / (input_path.stem + suffix)


def _resolve_theta(
    theta_cap: float | None,
    theta_cap_deg: float | None,
) -> float | None:
    if theta_cap is not None and theta_cap_deg is not None:
        _err.print("[red]Error:[/red] supply at most one of --theta-cap and --theta-cap-deg.")
        raise typer.Exit(1) from None
    if theta_cap_deg is not None:
        return np.deg2rad(theta_cap_deg)
    return theta_cap


def _is_commandline(ctx, name: str) -> bool:  # type: ignore[no-untyped-def]
    """True iff parameter *name* was supplied on the command line.

    The source is read from the *passed* context rather than
    ``click.get_current_context()``: typer >= 0.26 vendors its own copy of
    click (``typer._click``) with a separate context stack, so the global
    click context is empty during a typer command and the lookup would
    silently return ``None``.  The comparison is by enum *name* (not identity
    against ``click.core.ParameterSource``) for the same reason — typer's
    vendored ``ParameterSource`` is a distinct enum object from the
    top-level click one, so identity comparison is ``False`` under vendoring.
    """
    if ctx is None:
        return False
    try:
        src = ctx.get_parameter_source(name)
    except Exception:
        return False
    return src is not None and getattr(src, "name", None) == "COMMANDLINE"


def _make_picker(cfg, ctx):  # type: ignore[no-untyped-def]
    """Build a (name, cli_value, cfg_value) -> resolved-value function.

    Precedence: explicit CLI flag wins; otherwise config (when present and
    not None); otherwise the CLI default.  *ctx* is the typer/click context
    of the running command (see :func:`_is_commandline`).
    """

    def pick(name: str, cli_value, cfg_value):  # type: ignore[no-untyped-def]
        if cfg is None:
            return cli_value
        if _is_commandline(ctx, name):
            return cli_value
        if cfg_value is not None:
            return cfg_value
        return cli_value

    return pick


def _build_spacing_k(
    volume_k_spacing: str,
    volume_k_beta: float,
    volume_k_thickness: float | None,
    nk: int,
    ds: float,
    growth: float,
) -> Spacing1D:
    if volume_k_spacing in {"tanh", "tanh2"} and volume_k_thickness is None:
        _err.print(
            "[red]Error:[/red] --volume-k-thickness is required when "
            f"--volume-k-spacing={volume_k_spacing}."
        )
        raise typer.Exit(1) from None
    return Spacing1D(
        law=volume_k_spacing,
        ds=ds,
        growth=growth,
        beta=volume_k_beta,
        total_thickness=volume_k_thickness,
    )


# ---------------------------------------------------------------------------
# shore mesh
# ---------------------------------------------------------------------------


@app.command()
def mesh(
    ctx: typer.Context,
    input: Annotated[
        Path | None,
        typer.Argument(
            help="Input STL or CAD file (STEP/IGES/BREP, needs shore-mesh[cad]);"
            " optional when --config supplies it."
        ),
    ] = None,
    config: Annotated[
        Path | None,
        typer.Option(
            "-c",
            "--config",
            help="TOML config file (CLI flags take precedence over its values).",
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option("-o", "--output", help="Output stem for .geo files."),
    ] = None,
    ni: _NI = 40,
    nj: _NJ = 60,
    nk: _NK = 30,
    ds: _DS = 1e-3,
    growth: _GROWTH = 1.15,
    smooth: _SMOOTH = 0.2,
    smooth_iters: _SITERS = 2,
    smooth_mode: Annotated[
        str,
        typer.Option(
            "--smooth-mode",
            help=(
                "limited (default): limit the smoothing motion where it would fold or "
                "over-skew a layer; legacy: unlimited (reproduces pre-0.2 grids, folds "
                "thin first layers)."
            ),
        ),
    ] = "limited",
    march_scheme: Annotated[
        str,
        typer.Option(
            "--march-scheme",
            help=(
                "algebraic (default): normal extrusion + smoothing; hyperbolic: "
                "Steger-Chaussee marching (orthogonal k-lines, prescribed cell volumes)."
            ),
        ),
    ] = "algebraic",
    theta_cap: _THETA_CAP = None,
    theta_cap_deg: _THETA_CAP_DEG = None,
    cap_pull: _CAP_PULL = None,
    cap_squareness: _CAP_SQUARENESS = None,
    surface_i_spacing: _SURFACE_I_SPACING = "uniform",
    surface_i_beta: _SURFACE_I_BETA = 3.0,
    surface_j_spacing: _SURFACE_J_SPACING = "uniform",
    surface_j_beta: _SURFACE_J_BETA = 3.0,
    j_cluster_lon: _J_CLUSTER_LON = None,
    volume_k_spacing: _VOLUME_K_SPACING = "geometric",
    volume_k_beta: _VOLUME_K_BETA = 3.0,
    volume_k_thickness: _VOLUME_K_THICKNESS = None,
    blend_normals_k: _BLEND_K = 0,
    topology: Annotated[
        str,
        typer.Option(
            "--topology",
            help="Topology: 'cubed_sphere' (default), 'ogrid', 'cgrid' "
            "(C-grid for sharp-TE bodies — see --ni-body / --ni-wake / "
            "--n-stations / --wake-length), or 'ch' (C-grid wrap + "
            "H-channel far-field — adds --ch-chord-min/max / "
            "--ch-normal-min/max / --nk-channel).",
        ),
    ] = "cubed_sphere",
    ni_body: _NI_BODY = None,
    ni_wake: _NI_WAKE = None,
    n_stations: _N_STATIONS = None,
    wake_length: _WAKE_LENGTH = None,
    wake_growth: _WAKE_GROWTH = 1.15,
    span_axis: _SPAN_AXIS = "z",
    te_dihedral_deg: _TE_DIHEDRAL_DEG = 30.0,
    cgrid_blocks: _CGRID_BLOCKS = 3,
    cap_tips: _CAP_TIPS = "none",
    cap_te_thickness: _CAP_TE_THICKNESS = None,
    cap_te_thickness_frac: _CAP_TE_THICKNESS_FRAC = 0.002,
    cap_le_frac: _CAP_LE_FRAC = 0.05,
    cap_te_frac: _CAP_TE_FRAC = 0.95,
    cap_tip_inset: _CAP_TIP_INSET = None,
    cap_n_cap_j: _CAP_N_CAP_J = 12,
    cap_j_spacing_cgrid: _CAP_J_SPACING = "tanh",
    cap_j_beta_cgrid: _CAP_J_BETA = 4.0,
    cap_nk: _CAP_NK = None,
    ch_chord_min: _CH_CHORD_MIN = None,
    ch_chord_max: _CH_CHORD_MAX = None,
    ch_normal_min: _CH_NORMAL_MIN = None,
    ch_normal_max: _CH_NORMAL_MAX = None,
    nk_channel: _NK_CHANNEL = 5,
    k_law_channel: _K_LAW_CHANNEL = "uniform",
    k_beta_channel: _K_BETA_CHANNEL = 3.0,
    c1_seam_channel: _C1_SEAM_CHANNEL = True,
    ch_freestream_x: _CH_FREESTREAM_X = None,
    ch_freestream_y: _CH_FREESTREAM_Y = None,
    ch_freestream_z: _CH_FREESTREAM_Z = None,
    ch_i_split_body_upper: _CH_I_SPLIT_UPPER = None,
    ch_i_split_body_lower: _CH_I_SPLIT_LOWER = None,
    ch_ni_outflow: _CH_NI_OUTFLOW = 8,
    oh_box_xmin: _OH_BOX_XMIN = None,
    oh_box_xmax: _OH_BOX_XMAX = None,
    oh_box_ymin: _OH_BOX_YMIN = None,
    oh_box_ymax: _OH_BOX_YMAX = None,
    oh_box_zmin: _OH_BOX_ZMIN = None,
    oh_box_zmax: _OH_BOX_ZMAX = None,
    adjacency_json: Annotated[
        Path | None,
        typer.Option(
            "--adjacency-json",
            help=(
                "Also write a .adjacency.json sidecar describing the SHARED "
                "/ DIRICHLET face graph.  Required by `shore cc-par`."
            ),
        ),
    ] = None,
    journal: _JOURNAL = None,
    project: _PROJECT = None,
    component: _COMPONENT = None,
    geometry_name: Annotated[
        str | None,
        typer.Option("--geometry-name", help="Geometry name in --project (default: file stem)."),
    ] = None,
    verbose: _VERBOSE = False,
) -> None:
    """Full pipeline: STL -> structured volume -> geo files.

    Topologies (``--topology``): ``cubed_sphere`` (default), ``ogrid``,
    ``cgrid``, ``ch``, ``oh``.

    - ``cubed_sphere`` (default): 6-block decomposition (4 equatorial
      sub-blocks + 2 polar caps); writes 6 ``<stem>_<label>.geo`` files.
    - ``ogrid``: single-block lat-lon grid, j-periodic; writes one
      ``<stem>.geo`` file.

    Node vs cell counts: ``--ni``/``--nj``/``--nk`` are NODE counts; a
    block has (nodes - 1) cells per axis.  The one exception is ``--nj``
    under ``cubed_sphere``: it is the full equator-ring node count and is
    consumed as ``nj // 4``, so each of the 4 sub-blocks (and the polar
    caps) gets ``nj // 4`` cells, i.e. ``nj // 4 + 1`` nodes.  To land on
    a target cell count C per axis: ``--ni C+1``, ``--nk C+1``, and
    ``--nj 4*C`` (cubed_sphere) or ``--nj C+1`` (ogrid).  Example for
    power-of-two cells: ``--ni 65 --nj 128 --nk 65`` yields 64x32x64-cell
    sub-blocks and 32x32x64-cell caps.

    Parameter precedence (highest first): explicit CLI flag -> config file -> built-in default.
    """
    cfg: MeshConfig | None = None
    if config is not None:
        try:
            cfg = load_config(config)
        except Exception as exc:
            _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
            raise typer.Exit(1) from None

    # Resolve input STL: CLI positional > config > error.
    if input is None and cfg is not None:
        input = cfg.input
    if input is None:
        _err.print("[red]Error:[/red] no input STL given (pass as positional arg or via --config).")
        raise typer.Exit(1) from None

    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    # CLI > config > default for every scalar option.  ``_pick`` returns
    # the CLI value when the user typed the flag explicitly (click's
    # parameter source is COMMANDLINE), otherwise the config value
    # (when a config file was loaded), otherwise the CLI default.
    _pick = _make_picker(cfg, ctx)
    output = _pick("output", output, cfg.output if cfg else None)
    ni = _pick("ni", ni, cfg.grid.ni if cfg else None)
    nj = _pick("nj", nj, cfg.grid.nj if cfg else None)
    nk = _pick("nk", nk, cfg.grid.nk if cfg else None)
    ds = _pick("ds", ds, cfg.march.ds if cfg else None)
    growth = _pick("growth", growth, cfg.march.growth if cfg else None)
    smooth = _pick("smooth", smooth, cfg.march.smooth if cfg else None)
    smooth_iters = _pick("smooth_iters", smooth_iters, cfg.march.smooth_iters if cfg else None)
    smooth_mode = _pick("smooth_mode", smooth_mode, cfg.march.smooth_mode if cfg else None)
    march_scheme = _pick("march_scheme", march_scheme, cfg.march.march_scheme if cfg else None)
    if march_scheme not in ("algebraic", "hyperbolic"):
        _err.print(
            "[red]Error:[/red] --march-scheme must be 'algebraic' or 'hyperbolic', "
            f"got {march_scheme!r}."
        )
        raise typer.Exit(1)
    if smooth_mode not in ("limited", "legacy"):
        _err.print(
            f"[red]Error:[/red] --smooth-mode must be 'limited' or 'legacy', got {smooth_mode!r}."
        )
        raise typer.Exit(1)
    blend_normals_k = _pick(
        "blend_normals_k", blend_normals_k, cfg.march.blend_normals_k if cfg else None
    )
    if blend_normals_k > 0:
        _err.print(
            "[yellow]Warning:[/yellow] --blend-normals-k (and [march] blend_normals_k) is "
            "deprecated (ADR 0019): it showed no benefit over the grid normals, is refused on "
            "C-grid topologies, and will be removed. Use 0."
        )
    surface_i_spacing = _pick(
        "surface_i_spacing", surface_i_spacing, cfg.surface.i_spacing if cfg else None
    )
    surface_i_beta = _pick("surface_i_beta", surface_i_beta, cfg.surface.i_beta if cfg else None)
    surface_j_spacing = _pick(
        "surface_j_spacing", surface_j_spacing, cfg.surface.j_spacing if cfg else None
    )
    surface_j_beta = _pick("surface_j_beta", surface_j_beta, cfg.surface.j_beta if cfg else None)
    volume_k_spacing = _pick(
        "volume_k_spacing", volume_k_spacing, cfg.volume.k_spacing if cfg else None
    )
    volume_k_beta = _pick("volume_k_beta", volume_k_beta, cfg.volume.k_beta if cfg else None)
    volume_k_thickness = _pick(
        "volume_k_thickness", volume_k_thickness, cfg.volume.k_thickness if cfg else None
    )
    topology = _pick("topology", topology, cfg.grid.topology if cfg else None)
    if topology not in ("cubed_sphere", "ogrid", "cgrid", "ch", "oh"):
        _err.print(
            f"[red]Error:[/red] unknown --topology {topology!r}; "
            "expected 'cubed_sphere', 'ogrid', 'cgrid', 'ch', or 'oh'."
        )
        raise typer.Exit(1) from None

    if topology in ("cgrid", "ch"):
        for flag_name, val in (
            ("--ni-body", ni_body),
            ("--ni-wake", ni_wake),
            ("--n-stations", n_stations),
            ("--wake-length", wake_length),
        ):
            if val is None:
                _err.print(f"[red]Error:[/red] --topology {topology} requires {flag_name}.")
                raise typer.Exit(1) from None
        if span_axis not in ("x", "y", "z"):
            _err.print(
                f"[red]Error:[/red] --span-axis must be 'x', 'y', or 'z'; got {span_axis!r}."
            )
            raise typer.Exit(1) from None
        if topology == "cgrid" and cgrid_blocks not in (1, 3, 4):
            _err.print(
                f"[red]Error:[/red] --cgrid-blocks must be 1, 3, or 4; "
                f"got {cgrid_blocks}.  Use 3 (default, sharp-TE), "
                "4 (blunt-TE), or 1 (legacy single-block)."
            )
            raise typer.Exit(1) from None

    if topology == "ch":
        for flag_name, val in (
            ("--ch-chord-min", ch_chord_min),
            ("--ch-chord-max", ch_chord_max),
            ("--ch-normal-min", ch_normal_min),
            ("--ch-normal-max", ch_normal_max),
        ):
            if val is None:
                _err.print(f"[red]Error:[/red] --topology ch requires {flag_name}.")
                raise typer.Exit(1) from None

    if topology == "oh":
        for flag_name, val in (
            ("--oh-box-xmin", oh_box_xmin),
            ("--oh-box-xmax", oh_box_xmax),
            ("--oh-box-ymin", oh_box_ymin),
            ("--oh-box-ymax", oh_box_ymax),
            ("--oh-box-zmin", oh_box_zmin),
            ("--oh-box-zmax", oh_box_zmax),
        ):
            if val is None:
                _err.print(f"[red]Error:[/red] --topology oh requires {flag_name}.")
                raise typer.Exit(1) from None
        for axis, lo, hi in (
            ("x", oh_box_xmin, oh_box_xmax),
            ("y", oh_box_ymin, oh_box_ymax),
            ("z", oh_box_zmin, oh_box_zmax),
        ):
            if lo >= hi:
                _err.print(
                    f"[red]Error:[/red] --oh-box-{axis}min must be strictly "
                    f"less than --oh-box-{axis}max (got {lo} >= {hi})."
                )
                raise typer.Exit(1) from None
        # OH pins theta_cap = pi/4 internally (cube-face planarity); a
        # user-supplied --theta-cap / --theta-cap-deg is ignored.
        if _is_commandline(ctx, "theta_cap") or _is_commandline(ctx, "theta_cap_deg"):
            import warnings as _warn

            _warn.warn(
                "--theta-cap is ignored for --topology oh: the OH wrap pins "
                "theta_cap to pi/4 so each sub-block maps onto one cube face.",
                UserWarning,
                stacklevel=2,
            )

    # theta_cap is special: CLI exposes two flags (--theta-cap and
    # --theta-cap-deg) but the config has a single resolved radian
    # value.  Apply the same picker logic to whichever CLI flag was set.
    if cfg is not None and cfg.grid.theta_cap is not None:
        cli_set_rad = _is_commandline(ctx, "theta_cap")
        cli_set_deg = _is_commandline(ctx, "theta_cap_deg")
        if not (cli_set_rad or cli_set_deg):
            theta_cap = cfg.grid.theta_cap
            theta_cap_deg = None

    if cap_pull is not None:
        import warnings as _warn

        _warn.warn(
            "--cap-pull is deprecated (Phase 4): cap interior is now built by "
            "Shirley concentric + stereographic projection and the flag is ignored.",
            DeprecationWarning,
            stacklevel=2,
        )

    if cap_squareness is not None:
        import warnings as _warn

        _warn.warn(
            "--cap-squareness is deprecated and ignored: cap k=0 is now built "
            "by gnomonic flat-square projection.",
            DeprecationWarning,
            stacklevel=2,
        )

    if j_cluster_lon is not None:
        import warnings as _warn

        _warn.warn(
            "--j-cluster-lon is deprecated and ignored: gnomonic cap projection "
            "is symmetric in the 4 quadrants; use --surface-j-spacing tanh2 to "
            "cluster nodes near the cap-corner longitudes.",
            DeprecationWarning,
            stacklevel=2,
        )

    theta = _resolve_theta(theta_cap, theta_cap_deg)
    # validate --volume-k-* up front (prints the CLI error and exits); the
    # BuildTopology command builds the same Spacing1D from the same values
    _build_spacing_k(volume_k_spacing, volume_k_beta, volume_k_thickness, nk, ds, growth)
    out_stem = output or _default_output(input, "")

    # --- load body ---
    # C-grid topology doesn't ray-cast through the body — it slices the
    # triangulation directly with mesh_plane (for CAD input: the merged
    # tessellation).  Skip the trimesh.repair pass to avoid the optional
    # networkx dependency.
    _skip_hygiene = topology in ("cgrid", "ch")
    try:
        proj = _open_project(project, component, "main")
        if project is None:  # a throwaway project: the legacy default name
            proj.load_geometry(input, skip_hygiene=_skip_hygiene)
            geometry = "body"
        else:
            geometry = _ensure_geometry(proj, input, geometry_name, skip_hygiene=_skip_hygiene)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    if verbose:
        _console.print(
            f"[green]✓[/green] Loaded {_esc(str(input))}  "
            f"({len(proj.geometries[geometry].mesh.faces)} triangles)"
        )

    # --- build topology + march ---
    _theta_topo = theta if theta is not None else 0.05 * np.pi
    options: dict[str, object] = {}
    if topology == "ch":
        # Freestream is optional: pass it only if any component was given.
        _fs = (ch_freestream_x, ch_freestream_y, ch_freestream_z)
        _freestream = (
            tuple(float(v) if v is not None else 0.0 for v in _fs)
            if any(v is not None for v in _fs)
            else None
        )
        options = {
            "chord_min": ch_chord_min,
            "chord_max": ch_chord_max,
            "normal_min": ch_normal_min,
            "normal_max": ch_normal_max,
            "nk_channel": nk_channel,
            "freestream_world": _freestream,
            "span_axis": span_axis,
            "k_law_channel": k_law_channel,
            "k_beta_channel": k_beta_channel,
            "i_split_body_upper": ch_i_split_body_upper,
            "i_split_body_lower": ch_i_split_body_lower,
            "ni_outflow": ch_ni_outflow,
            "c1_seam": c1_seam_channel,
        }
    elif topology == "oh":
        options = {
            "box_vmin": (oh_box_xmin, oh_box_ymin, oh_box_zmin),
            "box_vmax": (oh_box_xmax, oh_box_ymax, oh_box_zmax),
            "nk_channel": nk_channel,
            "k_law_channel": k_law_channel,
            "k_beta_channel": k_beta_channel,
            "c1_seam": c1_seam_channel,
        }
    params: dict[str, object]
    if topology == "cgrid":
        if cap_tips not in ("none", "jlo", "jhi", "both"):
            _err.print(
                f"[red]Error:[/red] --cap-tips must be 'none' / 'jlo' / "
                f"'jhi' / 'both', got {cap_tips!r}."
            )
            raise typer.Exit(1) from None
        params = {
            "ni_body": ni_body,
            "ni_wake": ni_wake,
            "n_stations": n_stations,
            "wake_length": wake_length,
            "wake_growth": wake_growth,
            "span_axis": span_axis,
            "te_dihedral_deg": te_dihedral_deg,
            "blocks": cgrid_blocks,
            "cap_tips": cap_tips,
            "cap_te_thickness": cap_te_thickness,
            "cap_te_thickness_frac": cap_te_thickness_frac,
            "cap_le_frac": cap_le_frac,
            "cap_te_frac": cap_te_frac,
            "cap_tip_inset": cap_tip_inset,
            "cap_n_cap_j": cap_n_cap_j,
            "cap_j_spacing": cap_j_spacing_cgrid,
            "cap_j_beta": cap_j_beta_cgrid,
            "cap_nk": cap_nk,
        }
    elif topology == "ch":
        params = {
            "ni_body": ni_body,
            "ni_wake": ni_wake,
            "n_stations": n_stations,
            "wake_length": wake_length,
            "wake_growth": wake_growth,
            "te_dihedral_deg": te_dihedral_deg,
        }
    else:
        params = {
            "ni": ni,
            "nj": nj,
            "theta_cap": _theta_topo,
            "i_spacing": surface_i_spacing,
            "i_beta": surface_i_beta,
            "j_spacing": surface_j_spacing,
            "j_beta": surface_j_beta,
        }
    try:
        proj.build_topology(
            topology,
            nk=nk,
            ds=ds,
            growth=growth,
            k_law=volume_k_spacing,
            k_beta=volume_k_beta,
            k_thickness=volume_k_thickness,
            options=options,
            geometry=geometry,
            **params,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    try:
        if verbose:
            with _progress() as prog:
                t = prog.add_task("Marching      ", total=nk - 1, info="")

                def _cb(fraction: float, message: str) -> None:
                    prog.update(t, advance=1, info=message)

                proj.march(
                    smooth, smooth_iters, blend_normals_k, smooth_mode, march_scheme, progress=_cb
                )
        else:
            proj.march(smooth, smooth_iters, blend_normals_k, smooth_mode, march_scheme)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    # --- write ---
    m = proj.mesh
    assert m is not None
    proj.write_geo(out_stem)
    if verbose:
        if topology == "cubed_sphere":
            nc = nj // 4 + 1
            _console.print(
                f"[green]✓[/green] Wrote 6 .geo files  (sub: {ni}x{nc}x{nk}, caps: {nc}x{nc}x{nk})"
            )
        elif len(m.blocks) == 1:
            shape = m.blocks[0].nodes.shape
            single = out_stem if out_stem.suffix == ".geo" else out_stem.with_suffix(".geo")
            _console.print(f"[green]✓[/green] Wrote {single}  ({shape[0]}x{shape[1]}x{shape[2]})")
        else:
            n = len(m.blocks)
            labels = ", ".join(b.label for b in m.blocks)
            _console.print(f"[green]✓[/green] Wrote {n} .geo files ({labels})")

    if adjacency_json is not None:
        proj.write_adjacency(adjacency_json)
        if verbose:
            _console.print(f"[green]✓[/green] Wrote {adjacency_json}")

    if journal is not None:
        proj.save_journal(journal)
        if verbose:
            _console.print(f"[green]✓[/green] Wrote {journal}  (replay: shore run {journal})")
    _save_project(proj, project)


# ---------------------------------------------------------------------------
# shore run
# ---------------------------------------------------------------------------


@app.command()
def run(
    journal: Annotated[
        Path,
        typer.Argument(help="Journal written by `shore mesh --journal` or Project.save_journal."),
    ],
    save_project: Annotated[
        Path | None,
        typer.Option("--save-project", help="Also save the resulting project directory."),
    ] = None,
) -> None:
    """Replay a SHORE journal.

    A journal is a Python script of Project method calls: replaying it
    re-runs the same commands and rewrites the same outputs, bit for bit.
    It is executed as Python, so run only journals you trust.
    """
    from shore.project import run_journal

    if not journal.exists():
        _err.print(f"[red]Error:[/red] journal not found: {_esc(str(journal))}")
        raise typer.Exit(1)
    try:
        proj = run_journal(journal)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    if save_project is not None:
        proj.save(save_project)
    _console.print(
        f"[green]✓[/green] Replayed {len(proj.history)} commands from {_esc(str(journal))}"
        + (f"; project saved to {_esc(str(save_project))}" if save_project else "")
    )


# ---------------------------------------------------------------------------
# shore collar
# ---------------------------------------------------------------------------


@app.command()
def collar(
    body_a: Annotated[Path, typer.Argument(help="Body A (STL or CAD): the collar's wall side A.")],
    body_b: Annotated[Path, typer.Argument(help="Body B (STL or CAD), meeting body A.")],
    seed: Annotated[
        tuple[float, float, float],
        typer.Option("--seed", help="A point near the junction (selects the curve)."),
    ],
    n_loop: Annotated[int, typer.Option("--n-loop", help="Nodes along the junction curve.")],
    ds_surface: Annotated[float, typer.Option("--ds-surface", help="First surface step.")],
    growth_surface: Annotated[float, typer.Option("--growth-surface", help="Surface growth.")],
    m_a: Annotated[int, typer.Option("--m-a", help="Surface layers on body A.")],
    m_b: Annotated[int, typer.Option("--m-b", help="Surface layers on body B.")],
    ds_wall: Annotated[float, typer.Option("--ds-wall", help="First wall-normal step.")],
    growth_wall: Annotated[float, typer.Option("--growth-wall", help="Wall-normal growth.")],
    nk: Annotated[int, typer.Option("--nk", help="Wall-normal layers.")],
    label: Annotated[str, typer.Option("--label", help="Block label.")] = "collar",
    plane: Annotated[
        tuple[float, float, float, float, float, float] | None,
        typer.Option(
            "--plane", help="PX PY PZ NX NY NZ: symmetry plane of a half model (half-collar)."
        ),
    ] = None,
    name_a: Annotated[
        str | None, typer.Option("--name-a", help="Geometry name of body A (default: stem).")
    ] = None,
    name_b: Annotated[
        str | None, typer.Option("--name-b", help="Geometry name of body B (default: stem).")
    ] = None,
    output: Annotated[
        Path | None, typer.Option("--output", "-o", help="Output stem for the .geo file.")
    ] = None,
    adjacency_json: Annotated[
        Path | None, typer.Option("--adjacency-json", help="Also write a .adjacency.json.")
    ] = None,
    project: _PROJECT = None,
    component: _COMPONENT = None,
) -> None:
    """Build the collar block around the junction of two bodies (ADR 0010, 0013).

    With --plane, the half-collar of a half model on that symmetry plane.
    The result goes to the .geo file(s) (-o), to component --component
    (default 'collar') of --project, or both.
    """
    if output is None and project is None:
        _err.print("[red]Error:[/red] nothing to write: give -o and/or --project.")
        raise typer.Exit(1)
    try:
        proj = _open_project(project, component, "collar")
        a = _ensure_geometry(proj, body_a, name_a)
        b = _ensure_geometry(proj, body_b, name_b)
        proj.build_collar(
            b,
            geometry=a,
            seed=seed,
            n_loop=n_loop,
            ds_surface=ds_surface,
            growth_surface=growth_surface,
            m_a=m_a,
            m_b=m_b,
            ds_wall=ds_wall,
            growth_wall=growth_wall,
            nk=nk,
            label=label,
            plane=None if plane is None else (plane[:3], plane[3:]),
        )
        if output is not None:
            proj.write_geo(output)
        if adjacency_json is not None:
            proj.write_adjacency(adjacency_json)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    assert proj.mesh is not None
    shape = "x".join(map(str, proj.mesh.blocks[0].nodes.shape[:3]))
    _console.print(f"[green]✓[/green] Collar '{_esc(label)}' {shape}")
    _save_project(proj, project)


# ---------------------------------------------------------------------------
# shore xall
# ---------------------------------------------------------------------------

_XALL_RUN_KEYS = {
    "grd", "np", "output_basename", "multigrid_levels", "boundary_layer_thickness",
    "numerical_beach",
}  # fmt: skip
_XALL_COMPONENT_KEYS = {"free_face_bc", "body", "meta", "proc_meta", "boxes"}


def _xall_options(path: Path) -> dict[str, object]:
    """``write_xall`` keyword arguments from an options TOML file (ADR 0018)."""
    import json
    import tomllib

    from shore.exceptions import ParameterError

    data = tomllib.loads(path.read_text())
    unknown = set(data) - {"run", "components"}
    if unknown:
        raise ParameterError(
            f"{path}: unknown table(s) {sorted(unknown)}; use [run] and [components]."
        )
    run = dict(data.get("run", {}))
    unknown = set(run) - _XALL_RUN_KEYS
    if unknown:
        raise ParameterError(f"{path}: unknown [run] key(s) {sorted(unknown)}.")
    if "multigrid_levels" in run:
        run["multigrid_levels"] = tuple(run["multigrid_levels"])
    components: dict[str, dict[str, object]] = {}
    for name, options in data.get("components", {}).items():
        unknown = set(options) - _XALL_COMPONENT_KEYS
        if unknown:
            raise ParameterError(
                f"{path}: unknown key(s) {sorted(unknown)} in [components.{name}]."
            )
        options = dict(options)
        if "boxes" in options:
            options["boxes"] = json.loads((path.parent / options["boxes"]).read_text())
        components[name] = options
    return {**run, "components": components}


@app.command(name="xall")
def xall(
    project_dir: Annotated[Path, typer.Argument(help="A saved SHORE project directory.")],
    out: Annotated[Path, typer.Option("--out", help="Directory for the Xall run files.")],
    options: Annotated[
        Path | None,
        typer.Option("--options", help="TOML file: [run] and [components.<name>] options."),
    ] = None,
) -> None:
    """Write the Xall run of every component of a project (ADR 0016-0018).

    Per component: .geo files, adjacency sidecar, .grd and cc.par under
    OUT/<component>/; then the merged .grd, cc.par and one proc.input in OUT.
    """
    from shore.project import Project as _Project

    try:
        kwargs = _xall_options(options) if options is not None else {}
        proj = _Project.load(project_dir)
        proj.write_xall(out, **kwargs)  # type: ignore[arg-type]
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    names = [n for n, m in proj.components.items() if m is not None]
    _console.print(f"[green]✓[/green] Wrote the Xall run of {names} to {_esc(str(out))}")


# ---------------------------------------------------------------------------
# shore config
# ---------------------------------------------------------------------------


@app.command(name="config")
def config_cmd(
    action: Annotated[
        str,
        typer.Argument(
            help="One of: 'template' (emit a default config), 'validate' (check a config).",
        ),
    ],
    path: Annotated[
        Path | None,
        typer.Argument(help="Path to validate (required for 'validate'); ignored for 'template'."),
    ] = None,
) -> None:
    """Inspect SHORE TOML configs.

    Examples
    --------

        shore config template > shore.toml
        shore config validate shore.toml
    """
    if action == "template":
        dump_template()
        return

    if action == "validate":
        if path is None:
            _err.print("[red]Error:[/red] 'validate' requires a config-file path.")
            raise typer.Exit(1) from None
        try:
            cfg = load_config(path)
        except Exception as exc:
            _err.print(f"[red]✗[/red] {_esc(str(exc))}")
            raise typer.Exit(1) from None
        _console.print(f"[green]✓[/green] {path} parses cleanly")
        _console.print(f"  input         = {cfg.input}")
        _console.print(
            f"  grid          ni={cfg.grid.ni} nj={cfg.grid.nj} nk={cfg.grid.nk} "
            f"topology={cfg.grid.topology!r}"
        )
        _console.print(
            f"  march         ds={cfg.march.ds} growth={cfg.march.growth} "
            f"smooth={cfg.march.smooth}/{cfg.march.smooth_iters}"
        )
        _console.print(
            f"  surface       i={cfg.surface.i_spacing}({cfg.surface.i_beta}) "
            f"j={cfg.surface.j_spacing}({cfg.surface.j_beta})"
        )
        _console.print(
            f"  volume        k={cfg.volume.k_spacing}({cfg.volume.k_beta}) "
            f"thickness={cfg.volume.k_thickness}"
        )
        return

    _err.print(f"[red]Error:[/red] unknown action {action!r}; expected 'template' or 'validate'.")
    raise typer.Exit(1) from None


# ---------------------------------------------------------------------------
# shore project
# ---------------------------------------------------------------------------


@app.command()
def project(
    input: Annotated[
        Path,
        typer.Argument(help="Input STL or CAD file (STEP/IGES/BREP, needs shore-mesh[cad])."),
    ],
    output: Annotated[
        Path | None,
        typer.Option("-o", "--output", help="Output .npz file."),
    ] = None,
    ni: _NI = 40,
    nj: _NJ = 60,
    theta_cap: _THETA_CAP = None,
    theta_cap_deg: _THETA_CAP_DEG = None,
    surface_i_spacing: _SURFACE_I_SPACING = "uniform",
    surface_i_beta: _SURFACE_I_BETA = 3.0,
    surface_j_spacing: _SURFACE_J_SPACING = "uniform",
    surface_j_beta: _SURFACE_J_BETA = 3.0,
    j_cluster_lon: _J_CLUSTER_LON = None,
    verbose: _VERBOSE = False,
) -> None:
    """Project a structured sphere grid onto an STL surface (k=0 layer only)."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    theta = _resolve_theta(theta_cap, theta_cap_deg)
    out = output or _default_output(input, "_k0.npz")

    if verbose:
        with _progress() as prog:
            t = prog.add_task("Loading body  ", total=1, info="")
            body = Body.from_file(input)
            prog.update(t, advance=1, info=f"{len(body.mesh.faces)} triangles")

        with _progress() as prog:
            t = prog.add_task("Projecting    ", total=1, info="")
            k0 = body.project(
                ni,
                nj,
                theta_cap=theta,
                i_spacing=surface_i_spacing,
                i_beta=surface_i_beta,
                j_spacing=surface_j_spacing,
                j_beta=surface_j_beta,
                j_cluster_lon=j_cluster_lon,
            )
            prog.update(t, advance=1, info=f"({ni}x{nj}) surface layer")
    else:
        body = Body.from_file(input)
        k0 = body.project(
            ni,
            nj,
            theta_cap=theta,
            i_spacing=surface_i_spacing,
            i_beta=surface_i_beta,
            j_spacing=surface_j_spacing,
            j_beta=surface_j_beta,
            j_cluster_lon=j_cluster_lon,
        )

    Body.save_surface(k0, out)

    if verbose:
        _console.print(f"[green]✓[/green] Wrote {out}  ({ni}x{nj}x3)")


# ---------------------------------------------------------------------------
# shore info
# ---------------------------------------------------------------------------


@app.command()
def info(
    input: Annotated[Path, typer.Argument(help="Input .geo or .npz file.")],
) -> None:
    """Print metadata of a .geo block or .npz surface file."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    suffix = input.suffix.lower()

    if suffix == ".geo":
        data = info_geo(input)
        ni, nj, nk = data["shape"]
        table = Table(title=str(input), show_header=True, header_style="bold cyan")
        table.add_column("Property", style="dim")
        table.add_column("Value")
        table.add_row("Dimensions", f"ni={ni}  nj={nj}  nk={nk}")
        table.add_row("Points", f"{data['n_points']:,}")
        table.add_row("Cells", f"{data['n_cells']:,}")
        lo, hi = data["xyz_min"], data["xyz_max"]
        for ax, lo_val, hi_val in zip("xyz", lo, hi, strict=True):
            table.add_row(f"{ax} bounds", f"[{lo_val:.6g},  {hi_val:.6g}]")
        _console.print(table)

    elif suffix == ".npz":
        data = info_surface(input)
        ni, nj = data["shape"]
        table = Table(title=str(input), show_header=True, header_style="bold cyan")
        table.add_column("Property", style="dim")
        table.add_column("Value")
        table.add_row("Dimensions", f"ni={ni}  nj={nj}")
        table.add_row("Points", f"{data['n_points']:,}")
        lo, hi = data["xyz_min"], data["xyz_max"]
        for ax, lo_val, hi_val in zip("xyz", lo, hi, strict=True):
            table.add_row(f"{ax} bounds", f"[{lo_val:.6g},  {hi_val:.6g}]")
        _console.print(table)

    else:
        _err.print(f"[red]Error:[/red] unsupported extension '{suffix}'. Expected .geo or .npz.")
        raise typer.Exit(1) from None


# ---------------------------------------------------------------------------
# shore check
# ---------------------------------------------------------------------------


@app.command()
def check(
    input: Annotated[Path, typer.Argument(help="Input .geo file.")],
    warn_threshold: Annotated[
        float,
        typer.Option(help="Warn if jmin drops below this value."),
    ] = 1e-6,
) -> None:
    """Validate cell Jacobians layer-by-layer in a .geo file.

    Exits with code 1 if any layer has a negative Jacobian.
    """
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    result = check_geo(input, warn_threshold=warn_threshold)

    table = Table(title=f"Jacobian check - {input}", show_header=True, header_style="bold cyan")
    table.add_column("Layer", justify="right", style="dim")
    table.add_column("jmin", justify="right")
    table.add_column("Status")

    for k, jmin in enumerate(result["jmins"]):
        if jmin <= 0.0:
            status = "[bold red]INVERTED[/bold red]"
        elif jmin < warn_threshold:
            status = "[yellow]warn[/yellow]"
        else:
            status = "[green]ok[/green]"
        table.add_row(str(k + 1), f"{jmin:.4e}", status)

    _console.print(table)

    if not result["ok"]:
        _err.print("[red]✗[/red] Inverted cells detected.")
        raise typer.Exit(1) from None

    _console.print("[green]✓[/green] All Jacobians positive.")


# ---------------------------------------------------------------------------
# shore export
# ---------------------------------------------------------------------------


@app.command()
def export(
    inputs: Annotated[
        list[Path],
        typer.Argument(
            help=(
                "One or more input files.  A single .geo / .npz writes a "
                ".vts (or a .vtm with one block when --output ends in .vtm); "
                "multiple .geo files are bundled into a single .vtm "
                "MultiBlock."
            ),
        ),
    ],
    output: Annotated[
        Path | None,
        typer.Option("-o", "--output", help="Output file (.vts for single input, .vtm for multi)."),
    ] = None,
    no_jacobian: Annotated[
        bool,
        typer.Option("--no-jacobian", help="Skip Jacobian cell data (geo only)."),
    ] = False,
) -> None:
    """Export .geo / .npz files to VTK XML for ParaView.

    - A single ``.geo`` or ``.npz`` writes one ``.vts``.
    - Two or more ``.geo`` files are bundled into one ``.vtm`` MultiBlock
      (block labels are derived from each input's filename stem).
    """
    if not inputs:
        _err.print("[red]Error:[/red] no input files given.")
        raise typer.Exit(1) from None

    for inp in inputs:
        if not inp.exists():
            _err.print(f"[red]Error:[/red] file not found: {inp}")
            raise typer.Exit(1) from None

    # Single-input dispatch (preserved from the original command).
    # Exception: if the user explicitly asks for a .vtm output, fall
    # through to the multi-input path with one block — useful for
    # single-block topologies (O-grid, C-grid) whose authors want a
    # MultiBlock for consistency with multi-block topologies.
    if len(inputs) == 1 and not (output is not None and output.suffix.lower() == ".vtm"):
        input = inputs[0]
        suffix = input.suffix.lower()
        try:
            if suffix == ".geo":
                out = export_geo(input, output, with_jacobian=not no_jacobian)
            elif suffix == ".npz":
                out = export_surface(input, output)
            else:
                _err.print(
                    f"[red]Error:[/red] unsupported extension '{suffix}'. Expected .geo or .npz."
                )
                raise typer.Exit(1) from None
        except ImportError as exc:
            _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
            raise typer.Exit(1) from None
        _console.print(f"[green]✓[/green] Wrote {out}")
        return

    # Multi-input: every file must be .geo; bundle into one .vtm.
    bad = [p for p in inputs if p.suffix.lower() != ".geo"]
    if bad:
        _err.print(
            "[red]Error:[/red] multi-input export accepts .geo only; got: "
            + ", ".join(str(p) for p in bad)
        )
        raise typer.Exit(1) from None

    out_path = output if output is not None else inputs[0].with_suffix(".vtm")
    if out_path.suffix.lower() != ".vtm":
        out_path = out_path.with_suffix(".vtm")

    try:
        from shore.io.geo import read_geo
        from shore.io.vtk import write_caps_vtm
    except ImportError as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    blocks = [read_geo(p) for p in inputs]
    labels = [p.stem for p in inputs]
    try:
        write_caps_vtm(
            out_path,
            blocks,
            with_jacobian=not no_jacobian,
            labels=labels,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    _console.print(f"[green]✓[/green] Wrote {out_path} ({len(blocks)} blocks: {', '.join(labels)})")


# ---------------------------------------------------------------------------
# shore view
# ---------------------------------------------------------------------------


@app.command()
def view(
    input: Annotated[Path, typer.Argument(help="Input .geo or .npz file.")],
    colorby: Annotated[
        str,
        typer.Option(help="Colour field: 'Jacobian' (geo) or 'radius' (npz) or 'none'."),
    ] = "",
    no_edges: Annotated[bool, typer.Option("--no-edges", help="Hide cell edges.")] = False,
    opacity: Annotated[float, typer.Option(help="Grid opacity [0-1] (geo only).")] = 0.6,
) -> None:
    """Open an interactive PyVista window to inspect a .geo or .npz file."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    suffix = input.suffix.lower()

    try:
        if suffix == ".geo":
            cb = colorby if colorby else "Jacobian"
            view_geo_file(input, colorby=cb, show_edges=not no_edges, opacity=opacity)
        elif suffix == ".npz":
            cb = colorby if colorby else "radius"
            view_surface_file(input, colorby=cb, show_edges=not no_edges)
        else:
            _err.print(
                f"[red]Error:[/red] unsupported extension '{suffix}'. Expected .geo or .npz."
            )
            raise typer.Exit(1) from None
    except ImportError as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None


# ---------------------------------------------------------------------------
# shore primitive
# ---------------------------------------------------------------------------


primitive_app = typer.Typer(
    name="primitive",
    help="Parametric structured-grid primitives (box, annulus, flat-caps).",
    no_args_is_help=True,
)
app.add_typer(primitive_app)


_AXIS_HELP = (
    "Cylinder axis: 'x', 'y', 'z' (case-insensitive), or three "
    "components like '1.0 0.0 0.0' for an arbitrary direction."
)


def _resolve_axis_cli(spec: str | None) -> str | tuple[float, float, float] | None:
    """CLI helper: pass-through string axis; the primitive constructor
    handles both the label and the whitespace-separated triple form."""
    return spec


def _emit_primitive_writeback(
    proj: Project,
    journal: Path | None,
    output: Path | None,
    default_name: str,
    adjacency_json: Path | None = None,
) -> None:
    """Write the primitive's ``.geo`` files (and sidecar, journal) via *proj*.

    When ``adjacency_json`` is set, also emit a ``.adjacency.json``
    sidecar describing the SHARED / DIRICHLET face graph (for
    ``shore cc-par`` consumption).  Single-block primitives (box,
    annulus) emit a sidecar with one block and no connection patches.
    """
    out = output if output is not None else Path(f"{default_name}.geo")
    proj.write_geo(out)
    assert proj.mesh is not None
    if len(proj.mesh.blocks) > 1:
        stem = out.with_suffix("") if out.suffix == ".geo" else out
        for blk in proj.mesh.blocks:
            _console.print(f"[green]wrote[/green] {stem}_{blk.label}.geo")
    else:
        target = out if out.suffix == ".geo" else out.with_suffix(".geo")
        _console.print(f"[green]wrote[/green] {target}")

    if adjacency_json is not None:
        proj.write_adjacency(adjacency_json)
        _console.print(f"[green]wrote[/green] {adjacency_json}")
    if journal is not None:
        proj.save_journal(journal)
        _console.print(f"[green]wrote[/green] {journal}")


def _load_prim_config(path: Path | None) -> PrimitiveConfig | None:
    if path is None:
        return None
    try:
        return load_primitive_config(path)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None


def _pick_prim(cfg: PrimitiveConfig | None, ctx):  # type: ignore[no-untyped-def]
    """Same precedence picker as ``shore mesh``: explicit CLI > config > default."""

    def pick(name: str, cli_value, cfg_value):  # type: ignore[no-untyped-def]
        if cfg is None:
            return cli_value
        if _is_commandline(ctx, name):
            return cli_value
        if cfg_value is not None:
            return cfg_value
        return cli_value

    return pick


# ----- shore primitive box ------------------------------------------------


@primitive_app.command("box")
def primitive_box(
    ctx: typer.Context,
    config: Annotated[
        Path | None,
        typer.Option(
            "-c",
            "--config",
            help='TOML primitive config (kind = "primitive", type = "box").',
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option("-o", "--output", help="Output .geo file."),
    ] = None,
    vmin: Annotated[
        tuple[float, float, float],
        typer.Option("--min", help="Minimum corner of the box (x y z)."),
    ] = (0.0, 0.0, 0.0),
    vmax: Annotated[
        tuple[float, float, float],
        typer.Option("--max", help="Maximum corner of the box (x y z)."),
    ] = (1.0, 1.0, 1.0),
    ni: Annotated[int, typer.Option(help="Node count along x (i-axis).")] = 40,
    nj: Annotated[int, typer.Option(help="Node count along y (j-axis).")] = 40,
    nk: Annotated[int, typer.Option(help="Node count along z (k-axis).")] = 30,
    i_law: Annotated[
        str, typer.Option("--i-spacing", help="Spacing along i: uniform | tanh | tanh2.")
    ] = "uniform",
    j_law: Annotated[
        str, typer.Option("--j-spacing", help="Spacing along j: uniform | tanh | tanh2.")
    ] = "uniform",
    k_law: Annotated[
        str, typer.Option("--k-spacing", help="Spacing along k: uniform | tanh | tanh2.")
    ] = "uniform",
    i_beta: Annotated[float, typer.Option("--i-beta")] = 3.0,
    j_beta: Annotated[float, typer.Option("--j-beta")] = 3.0,
    k_beta: Annotated[float, typer.Option("--k-beta")] = 3.0,
    adjacency_json: Annotated[
        Path | None,
        typer.Option(
            "--adjacency-json",
            help="Also write a .adjacency.json sidecar (for shore cc-par).",
        ),
    ] = None,
    journal: _JOURNAL = None,
    project: _PROJECT = None,
    component: _COMPONENT = None,
) -> None:
    """Build a Cartesian box primitive (single block)."""
    cfg = _load_prim_config(config)
    if cfg is not None and cfg.geometry.type != "box":
        _err.print(
            f"[red]Error:[/red] config geometry.type='{cfg.geometry.type}' does not match 'box'."
        )
        raise typer.Exit(1) from None

    pick = _pick_prim(cfg, ctx)
    g = cfg.geometry if cfg else None
    grid = cfg.grid if cfg else None
    sp = cfg.spacing if cfg else None

    vmin = pick("vmin", vmin, g.vmin if g else None)
    vmax = pick("vmax", vmax, g.vmax if g else None)
    ni = pick("ni", ni, grid.ni if grid else None)
    nj = pick("nj", nj, grid.nj if grid else None)
    nk = pick("nk", nk, grid.nk if grid else None)
    i_law = pick("i_law", i_law, sp.i_law if sp else None)
    j_law = pick("j_law", j_law, sp.j_law if sp else None)
    k_law = pick("k_law", k_law, sp.k_law if sp else None)
    i_beta = pick("i_beta", i_beta, sp.i_beta if sp else None)
    j_beta = pick("j_beta", j_beta, sp.j_beta if sp else None)
    k_beta = pick("k_beta", k_beta, sp.k_beta if sp else None)
    output = pick("output", output, cfg.output if cfg else None)

    try:
        proj = _open_project(project, component, "main")
        proj.build_primitive(
            "box",
            vmin=tuple(vmin),
            vmax=tuple(vmax),
            ni=ni,
            nj=nj,
            nk=nk,
            i_law=i_law,
            j_law=j_law,
            k_law=k_law,
            i_beta=i_beta,
            j_beta=j_beta,
            k_beta=k_beta,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    _emit_primitive_writeback(
        proj, journal, output, default_name="box", adjacency_json=adjacency_json
    )
    _save_project(proj, project)


# ----- shore primitive annulus -------------------------------------------


@primitive_app.command("annulus")
def primitive_annulus(
    ctx: typer.Context,
    config: Annotated[
        Path | None,
        typer.Option(
            "-c",
            "--config",
            help='TOML primitive config (kind = "primitive", type = "annulus").',
        ),
    ] = None,
    output: Annotated[Path | None, typer.Option("-o", "--output", help="Output .geo file.")] = None,
    r_in: Annotated[float, typer.Option("--r-in", help="Inner radius.")] = 0.5,
    r_out: Annotated[float, typer.Option("--r-out", help="Outer radius.")] = 2.0,
    z0: Annotated[float, typer.Option("--z0", help="Axial start (along --axis).")] = 0.0,
    z1: Annotated[float, typer.Option("--z1", help="Axial end (along --axis).")] = 1.0,
    axis: Annotated[str, typer.Option("--axis", help=_AXIS_HELP)] = "z",
    center: Annotated[
        tuple[float, float, float],
        typer.Option("--center", help="World-space anchor for the cylinder axis."),
    ] = (0.0, 0.0, 0.0),
    theta_start: Annotated[
        float, typer.Option("--theta-start", help="Sector start in radians (default 0).")
    ] = 0.0,
    theta_end: Annotated[
        float | None,
        typer.Option(
            "--theta-end",
            help="Sector end in radians; omit for a full annulus.",
        ),
    ] = None,
    ni: Annotated[int, typer.Option(help="Radial node count (i-axis).")] = 20,
    nj: Annotated[int, typer.Option(help="Azimuthal node count (j-axis).")] = 64,
    nk: Annotated[int, typer.Option(help="Axial node count (k-axis).")] = 30,
    i_law: Annotated[str, typer.Option("--i-spacing")] = "uniform",
    j_law: Annotated[str, typer.Option("--j-spacing")] = "uniform",
    k_law: Annotated[str, typer.Option("--k-spacing")] = "uniform",
    i_beta: Annotated[float, typer.Option("--i-beta")] = 3.0,
    j_beta: Annotated[float, typer.Option("--j-beta")] = 3.0,
    k_beta: Annotated[float, typer.Option("--k-beta")] = 3.0,
    adjacency_json: Annotated[
        Path | None,
        typer.Option(
            "--adjacency-json",
            help="Also write a .adjacency.json sidecar (for shore cc-par).",
        ),
    ] = None,
    journal: _JOURNAL = None,
    project: _PROJECT = None,
    component: _COMPONENT = None,
) -> None:
    """Build an annular-cylinder primitive (single block, optionally periodic-j)."""
    cfg = _load_prim_config(config)
    if cfg is not None and cfg.geometry.type != "annulus":
        _err.print(
            f"[red]Error:[/red] config geometry.type='{cfg.geometry.type}' "
            "does not match 'annulus'."
        )
        raise typer.Exit(1) from None

    pick = _pick_prim(cfg, ctx)
    g = cfg.geometry if cfg else None
    grid = cfg.grid if cfg else None
    sp = cfg.spacing if cfg else None

    r_in = pick("r_in", r_in, g.r_in if g else None)
    r_out = pick("r_out", r_out, g.r_out if g else None)
    z0 = pick("z0", z0, g.z0 if g else None)
    z1 = pick("z1", z1, g.z1 if g else None)
    axis = pick("axis", axis, g.axis if g else None)
    center = pick("center", center, g.center if g else None)
    theta_start = pick("theta_start", theta_start, g.theta_start if g else None)
    theta_end = pick("theta_end", theta_end, g.theta_end if g else None)
    ni = pick("ni", ni, grid.ni if grid else None)
    nj = pick("nj", nj, grid.nj if grid else None)
    nk = pick("nk", nk, grid.nk if grid else None)
    i_law = pick("i_law", i_law, sp.i_law if sp else None)
    j_law = pick("j_law", j_law, sp.j_law if sp else None)
    k_law = pick("k_law", k_law, sp.k_law if sp else None)
    i_beta = pick("i_beta", i_beta, sp.i_beta if sp else None)
    j_beta = pick("j_beta", j_beta, sp.j_beta if sp else None)
    k_beta = pick("k_beta", k_beta, sp.k_beta if sp else None)
    output = pick("output", output, cfg.output if cfg else None)

    try:
        proj = _open_project(project, component, "main")
        proj.build_primitive(
            "annulus",
            r_in=r_in,
            r_out=r_out,
            z0=z0,
            z1=z1,
            ni=ni,
            nj=nj,
            nk=nk,
            theta_start=theta_start,
            theta_end=theta_end,
            axis=axis,
            center=tuple(center),
            i_law=i_law,
            j_law=j_law,
            k_law=k_law,
            i_beta=i_beta,
            j_beta=j_beta,
            k_beta=k_beta,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    _emit_primitive_writeback(
        proj, journal, output, default_name="annulus", adjacency_json=adjacency_json
    )
    _save_project(proj, project)


# ----- shore primitive flat-caps -----------------------------------------


@primitive_app.command("flat-caps")
def primitive_flat_caps(
    ctx: typer.Context,
    config: Annotated[
        Path | None,
        typer.Option(
            "-c",
            "--config",
            help='TOML primitive config (kind = "primitive", type = "flat_caps").',
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "-o",
            "--output",
            help="Output stem; writes 5 files: <stem>_center.geo + <stem>_sub_*.geo.",
        ),
    ] = None,
    r_out: Annotated[float, typer.Option("--r-out", help="Outer radius.")] = 2.0,
    z0: Annotated[float, typer.Option("--z0", help="Axial start (along --axis).")] = 0.0,
    z1: Annotated[float, typer.Option("--z1", help="Axial end (along --axis).")] = 1.0,
    axis: Annotated[str, typer.Option("--axis", help=_AXIS_HELP)] = "z",
    center: Annotated[
        tuple[float, float, float],
        typer.Option("--center"),
    ] = (0.0, 0.0, 0.0),
    inner_frac: Annotated[
        float,
        typer.Option(
            "--inner-frac",
            help=(
                "Ratio r_inner/r_out - the central square's corner radius "
                "as a fraction of the outer disk radius (must lie in (0, 1))."
            ),
        ),
    ] = 0.5,
    ni: Annotated[int, typer.Option("--ni", help="Radial node count for each sub-block.")] = 20,
    nc: Annotated[
        int, typer.Option("--nc", help="Side node count for the central square (shared).")
    ] = 20,
    nk: Annotated[int, typer.Option(help="Axial node count (k-axis).")] = 30,
    i_law: Annotated[str, typer.Option("--i-spacing")] = "uniform",
    c_law: Annotated[
        str,
        typer.Option("--c-spacing", help="Distribution along the central-square / arc direction."),
    ] = "uniform",
    k_law: Annotated[str, typer.Option("--k-spacing")] = "uniform",
    i_beta: Annotated[float, typer.Option("--i-beta")] = 3.0,
    c_beta: Annotated[float, typer.Option("--c-beta")] = 3.0,
    k_beta: Annotated[float, typer.Option("--k-beta")] = 3.0,
    adjacency_json: Annotated[
        Path | None,
        typer.Option(
            "--adjacency-json",
            help="Also write a .adjacency.json sidecar (for shore cc-par).",
        ),
    ] = None,
    journal: _JOURNAL = None,
    project: _PROJECT = None,
    component: _COMPONENT = None,
) -> None:
    """Build a 5-block flat-caps cylinder (no central hole - Chimera-free background)."""
    cfg = _load_prim_config(config)
    if cfg is not None and cfg.geometry.type != "flat_caps":
        _err.print(
            f"[red]Error:[/red] config geometry.type='{cfg.geometry.type}' "
            "does not match 'flat_caps'."
        )
        raise typer.Exit(1) from None

    pick = _pick_prim(cfg, ctx)
    g = cfg.geometry if cfg else None
    grid = cfg.grid if cfg else None
    sp = cfg.spacing if cfg else None

    # In the flat_caps config the central-square direction lives in
    # ``[grid].nj`` and ``[spacing].j_law`` for visual consistency
    # with box / annulus.  We map them here onto the kwargs that
    # ``FlatCapsCylinder.build`` actually expects (``nc``, ``c_law``).
    r_out = pick("r_out", r_out, g.r_out if g else None)
    z0 = pick("z0", z0, g.z0 if g else None)
    z1 = pick("z1", z1, g.z1 if g else None)
    axis = pick("axis", axis, g.axis if g else None)
    center = pick("center", center, g.center if g else None)
    inner_frac = pick("inner_frac", inner_frac, g.inner_frac if g else None)
    ni = pick("ni", ni, grid.ni if grid else None)
    nc = pick("nc", nc, grid.nj if grid else None)
    nk = pick("nk", nk, grid.nk if grid else None)
    i_law = pick("i_law", i_law, sp.i_law if sp else None)
    c_law = pick("c_law", c_law, sp.j_law if sp else None)
    k_law = pick("k_law", k_law, sp.k_law if sp else None)
    i_beta = pick("i_beta", i_beta, sp.i_beta if sp else None)
    c_beta = pick("c_beta", c_beta, sp.j_beta if sp else None)
    k_beta = pick("k_beta", k_beta, sp.k_beta if sp else None)
    output = pick("output", output, cfg.output if cfg else None)

    try:
        proj = _open_project(project, component, "main")
        proj.build_primitive(
            "flat_caps",
            r_out=r_out,
            z0=z0,
            z1=z1,
            ni=ni,
            nc=nc,
            nk=nk,
            inner_frac=inner_frac,
            axis=axis,
            center=tuple(center),
            i_law=i_law,
            c_law=c_law,
            k_law=k_law,
            i_beta=i_beta,
            c_beta=c_beta,
            k_beta=k_beta,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    _emit_primitive_writeback(
        proj, journal, output, default_name="flat_caps", adjacency_json=adjacency_json
    )
    _save_project(proj, project)


# ----- shore primitive h-grid-rect ---------------------------------------


@primitive_app.command("h-grid-rect")
def primitive_h_grid_rect(
    output: Annotated[
        Path | None,
        typer.Option("-o", "--output", help="Output .geo file."),
    ] = None,
    length: Annotated[
        float, typer.Option("--length", help="Duct length along the streamwise axis.")
    ] = 1.0,
    width: Annotated[
        float, typer.Option("--width", help="Cross-section width (along block j).")
    ] = 0.1,
    height: Annotated[
        float, typer.Option("--height", help="Cross-section height (along block k).")
    ] = 0.1,
    ni: Annotated[int, typer.Option(help="Streamwise node count (i-axis).")] = 40,
    nj: Annotated[int, typer.Option(help="Cross-section width node count (j-axis).")] = 20,
    nk: Annotated[int, typer.Option(help="Cross-section height node count (k-axis).")] = 20,
    streamwise_axis: Annotated[
        str,
        typer.Option(
            "--streamwise-axis",
            help="World axis the duct is aligned with: 'x' (default), 'y', or 'z'.",
        ),
    ] = "x",
    i_law: Annotated[
        str, typer.Option("--i-spacing", help="Streamwise spacing: uniform | tanh | tanh2.")
    ] = "uniform",
    j_law: Annotated[str, typer.Option("--j-spacing", help="Cross-section j spacing.")] = "uniform",
    k_law: Annotated[str, typer.Option("--k-spacing", help="Cross-section k spacing.")] = "uniform",
    i_beta: Annotated[float, typer.Option("--i-beta")] = 3.0,
    j_beta: Annotated[float, typer.Option("--j-beta")] = 3.0,
    k_beta: Annotated[float, typer.Option("--k-beta")] = 3.0,
    adjacency_json: Annotated[
        Path | None,
        typer.Option(
            "--adjacency-json",
            help="Also write a .adjacency.json sidecar (for shore cc-par).",
        ),
    ] = None,
    journal: _JOURNAL = None,
    project: _PROJECT = None,
    component: _COMPONENT = None,
) -> None:
    """Build a rectangular-duct H-grid (single block, internal flow).

    Default BCs: i-faces FREE (streamwise inflow / outflow), the four
    cross-section faces (j_lo, j_hi, k_lo, k_hi) WALL.  Use case:
    prismatic ducts; the analytic baseline for internal flow.
    """
    try:
        proj = _open_project(project, component, "main")
        proj.build_primitive(
            "h_grid_rect",
            length=length,
            width=width,
            height=height,
            ni=ni,
            nj=nj,
            nk=nk,
            streamwise_axis=streamwise_axis,
            i_law=i_law,
            j_law=j_law,
            k_law=k_law,
            i_beta=i_beta,
            j_beta=j_beta,
            k_beta=k_beta,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    _emit_primitive_writeback(
        proj, journal, output, default_name="h_grid_rect", adjacency_json=adjacency_json
    )
    _save_project(proj, project)


# ----- shore primitive h-grid-tube ---------------------------------------


@primitive_app.command("h-grid-tube")
def primitive_h_grid_tube(
    output: Annotated[
        Path | None,
        typer.Option("-o", "--output", help="Output .geo file."),
    ] = None,
    length: Annotated[
        float, typer.Option("--length", help="Tube length along the streamwise axis.")
    ] = 1.0,
    radius: Annotated[float, typer.Option("--radius", help="Outer-wall radius.")] = 0.05,
    n_axial: Annotated[int, typer.Option("--n-axial", help="Streamwise node count.")] = 40,
    n_circ: Annotated[
        int, typer.Option("--n-circ", help="Circumferential node count (j-axis, periodic).")
    ] = 32,
    n_radial: Annotated[
        int,
        typer.Option(
            "--n-radial",
            help="Radial node count (k-axis: 0 = centerline, last = wall).",
        ),
    ] = 8,
    streamwise_axis: Annotated[
        str,
        typer.Option(
            "--streamwise-axis",
            help="World axis the tube is aligned with: 'x' (default), 'y', or 'z'.",
        ),
    ] = "x",
    i_law: Annotated[
        str, typer.Option("--i-spacing", help="Streamwise spacing: uniform | tanh | tanh2.")
    ] = "uniform",
    k_law: Annotated[
        str,
        typer.Option(
            "--k-spacing",
            help="Radial spacing (centerline → wall).  Use tanh / tanh2 to cluster near the wall.",
        ),
    ] = "uniform",
    i_beta: Annotated[float, typer.Option("--i-beta")] = 3.0,
    k_beta: Annotated[float, typer.Option("--k-beta")] = 3.0,
    adjacency_json: Annotated[
        Path | None,
        typer.Option(
            "--adjacency-json",
            help="Also write a .adjacency.json sidecar (for shore cc-par).",
        ),
    ] = None,
    journal: _JOURNAL = None,
    project: _PROJECT = None,
    component: _COMPONENT = None,
) -> None:
    """Build a circular-bore H-grid (single block, j-periodic).

    The cross-section is a 2D polar grid: ``j`` circumferential (periodic),
    ``k`` radial from the centerline (k=0, polar singular axis) to the
    outer wall (k=n_radial-1).  Default BCs: i-faces FREE, k_lo and k_hi
    WALL.  Use case: straight cylindrical ducts, simple SRM chambers.

    The k=0 polar axis collapses every angular node to a single point —
    accept the singularity (Xall handles it) or use ``h-grid-rect`` if
    you want a non-singular cross-section.
    """
    try:
        proj = _open_project(project, component, "main")
        proj.build_primitive(
            "h_grid_tube",
            length=length,
            radius=radius,
            n_axial=n_axial,
            n_circ=n_circ,
            n_radial=n_radial,
            streamwise_axis=streamwise_axis,
            i_law=i_law,
            k_law=k_law,
            i_beta=i_beta,
            k_beta=k_beta,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    _emit_primitive_writeback(
        proj, journal, output, default_name="h_grid_tube", adjacency_json=adjacency_json
    )
    _save_project(proj, project)


# ---------------------------------------------------------------------------
# shore cap-tip
# ---------------------------------------------------------------------------


@app.command(name="cap-tip")
def cap_tip(
    body_geo: Annotated[
        Path,
        typer.Argument(
            help=(
                "Path to the C-grid body block's .geo file (typically "
                "'<stem>_body.geo' from `shore mesh --topology cgrid`)."
            ),
        ),
    ],
    tip: Annotated[
        str,
        typer.Option(
            "--tip",
            help="Which spanwise tip to cap: 'jlo' (j=0 station) or 'jhi' (j=n_stations-1).",
        ),
    ] = "jlo",
    le_frac: Annotated[
        float,
        typer.Option(
            "--le-frac",
            help="Chord fraction for the LE-cut (0 < le_frac < te_frac; default 0.05).",
        ),
    ] = 0.05,
    te_frac: Annotated[
        float,
        typer.Option(
            "--te-frac",
            help=(
                "Chord fraction for the TE-cut (le_frac < te_frac < 1; default 0.95).  "
                "Splits the chord band into a 'mid' region (cap_upper / cap_lower) "
                "and a 'TE' region (cap_te_upper / cap_te_lower).  Lower te_frac "
                "values give the TE region more chordwise resolution."
            ),
        ),
    ] = 0.95,
    te_thickness: Annotated[
        float | None,
        typer.Option(
            "--te-thickness",
            help=(
                "Auto-round synthetic TE thickness (absolute, in body units).  "
                "Mutually exclusive with --te-thickness-frac.  Only used when "
                "the input cross-section has a sharp TE."
            ),
        ),
    ] = None,
    te_thickness_frac: Annotated[
        float,
        typer.Option(
            "--te-thickness-frac",
            help=(
                "Auto-round synthetic TE thickness, as a fraction of chord "
                "(default 0.002 = 0.2%).  Only used when the input cross-section "
                "has a sharp TE.  The body wrap is left untouched; only the cap "
                "is rounded for its own structural soundness."
            ),
        ),
    ] = 0.002,
    no_auto_round_te: Annotated[
        bool,
        typer.Option(
            "--no-auto-round-te",
            help=(
                "Disable auto-rounding.  Sharp-TE inputs are then rejected with "
                "an error rather than rounded.  The default (auto-round on) is "
                "almost always what you want."
            ),
        ),
    ] = False,
    tip_inset: Annotated[
        float | None,
        typer.Option(
            "--tip-inset",
            help=(
                "Spanwise wall thickness of the cap.  Default = "
                "1e-3 * chord (auto-detected from the cross-section)."
            ),
        ),
    ] = None,
    n_cap_j: Annotated[
        int,
        typer.Option("--n-cap-j", help="Cap j-axis (airfoil-thickness) node count."),
    ] = 12,
    nk: Annotated[
        int | None,
        typer.Option(
            "--nk",
            help=(
                "Cap k-axis (spanwise) node count.  Default: match the "
                "body block's k-axis node count."
            ),
        ),
    ] = None,
    cap_j_spacing: Annotated[
        str,
        typer.Option(
            "--cap-j-spacing",
            help="Cap j-axis spacing law: uniform | tanh | tanh2.  Default 'tanh'.",
        ),
    ] = "tanh",
    cap_j_beta: Annotated[
        float, typer.Option("--cap-j-beta", help="Clustering strength for tanh / tanh2.")
    ] = 4.0,
    span_axis: Annotated[
        str,
        typer.Option(
            "--span-axis",
            help="World axis the wing spans: 'x', 'y', or 'z' (default 'z').",
        ),
    ] = "z",
    output: Annotated[
        Path | None,
        typer.Option("-o", "--output", help="Output stem (default: 'cap_<tip>')."),
    ] = None,
    adjacency_json: Annotated[
        Path | None,
        typer.Option(
            "--adjacency-json",
            help="Also write the .adjacency.json sidecar (for shore cc-par).",
        ),
    ] = None,
    journal: _JOURNAL = None,
) -> None:
    """Build a wing-tip cap from a C-grid body block's .geo file.

    The cap is an independent 3-block (sharp TE) or 4-block (blunt TE)
    H-style structured component that closes the wing's airfoil cross-
    section at one spanwise extreme.  It composes against the C-grid
    body wrap as a chimera component via ``shore cc-par`` — the cap's
    k_lo face is WALL (the actual wing-tip closure), k_hi is FREE
    (chimera fringe receiving from the body wrap).

    Sharp / blunt is **inferred** from the body block's TE shape: if
    the body's i=0 and i=ni-1 nodes coincide at every spanwise station,
    sharp.  Otherwise blunt — which currently requires the full
    wake_base block; emit a clear error pointing the user at the
    Python API in that case (CLI blunt support pending).
    """
    if tip not in ("jlo", "jhi"):
        _err.print(f"[red]Error:[/red] --tip must be 'jlo' or 'jhi', got {tip!r}.")
        raise typer.Exit(1) from None
    if cap_j_spacing not in ("uniform", "tanh", "tanh2"):
        _err.print(
            f"[red]Error:[/red] --cap-j-spacing must be 'uniform' / 'tanh' / 'tanh2', "
            f"got {cap_j_spacing!r}."
        )
        raise typer.Exit(1) from None
    if span_axis not in ("x", "y", "z"):
        _err.print("[red]Error:[/red] --span-axis must be 'x', 'y', or 'z'.")
        raise typer.Exit(1) from None
    if not (0.0 < le_frac < 1.0):
        _err.print(f"[red]Error:[/red] --le-frac must be in (0, 1), got {le_frac}.")
        raise typer.Exit(1) from None
    if not (0.0 < te_frac < 1.0):
        _err.print(f"[red]Error:[/red] --te-frac must be in (0, 1), got {te_frac}.")
        raise typer.Exit(1) from None
    if not (le_frac < te_frac):
        _err.print(f"[red]Error:[/red] --le-frac ({le_frac}) must be < --te-frac ({te_frac}).")
        raise typer.Exit(1) from None
    if te_thickness is not None and te_thickness <= 0:
        _err.print(f"[red]Error:[/red] --te-thickness must be > 0, got {te_thickness}.")
        raise typer.Exit(1) from None
    if not body_geo.exists():
        _err.print(f"[red]Error:[/red] body geo file not found: {body_geo}")
        raise typer.Exit(1) from None

    # ------------------------------------------------------------
    # Load the body block, extract the spanwise-extreme cross-section
    # at the wall (k=0).
    # ------------------------------------------------------------
    proj = Project()
    try:
        proj.build_tip_cap(
            body_geo,
            tip=tip,
            le_frac=le_frac,
            te_frac=te_frac,
            span_axis=span_axis,
            tip_inset=tip_inset,
            nk=nk,
            n_cap_j=n_cap_j,
            cap_j_spacing=cap_j_spacing,
            cap_j_beta=cap_j_beta,
            auto_round_te=not no_auto_round_te,
            te_thickness=te_thickness,
            te_thickness_frac=te_thickness_frac,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    resolved = proj.last_result or {}
    if "tip_inset" in resolved:
        _console.print(f"[dim]auto tip_inset = {resolved['tip_inset']:.6g}[/dim]")
    if "nk" in resolved:
        _console.print(f"[dim]auto nk = {resolved['nk']} (matches body)[/dim]")

    out_stem = output if output is not None else Path(f"cap_{tip}")
    proj.write_geo(out_stem, per_block=True)
    assert proj.mesh is not None
    for blk in proj.mesh.blocks:
        _console.print(f"[green]wrote[/green] {out_stem.parent / out_stem.name}_{blk.label}.geo")

    if adjacency_json is not None:
        proj.write_adjacency(adjacency_json)
        _console.print(f"[green]wrote[/green] {adjacency_json}")
    if journal is not None:
        proj.save_journal(journal)
        _console.print(f"[green]wrote[/green] {journal}")


# ---------------------------------------------------------------------------
# shore split
# ---------------------------------------------------------------------------


def _blocks_from_adjacency_shapes(adjacency_path):
    """Build lightweight wired HexBlocks from an adjacency sidecar.

    Nodes are zero-filled at the sidecar's shape (only extents matter for
    cut propagation / floor checking — coordinates are irrelevant here), and
    faces are wired from the sidecar.  Used by `shore balance` to floor-check
    seam-propagated cuts without reading the (large) .geo node data.
    """
    import numpy as np

    from shore.io.adjacency import read_adjacency_json, wire_blocks_from_adjacency
    from shore.mesh.block import HexBlock

    payload = read_adjacency_json(adjacency_path)
    blocks = []
    for b in payload["blocks"]:
        ni, nj, nk = b["shape"]
        blocks.append(HexBlock(label=b["label"], nodes=np.zeros((ni, nj, nk, 3))))
    wire_blocks_from_adjacency(blocks, payload, source=str(adjacency_path))
    return blocks


@app.command()
def split(
    inputs: Annotated[
        list[Path],
        typer.Argument(
            help=(
                "Input .geo files to split.  Block labels are derived from each "
                "file's stem (e.g. 'sub0' from 'sphere_sub0.geo')."
            ),
        ),
    ],
    config: Annotated[
        Path,
        typer.Option(
            "-c",
            "--config",
            help="TOML split config (kind = 'split') describing one or more BlockSplit directives.",
        ),
    ],
    output: Annotated[
        Path | None,
        typer.Option(
            "-o",
            "--output",
            help=(
                "Output stem.  Writes one '<stem>_<label>.geo' per chunk plus "
                "'<stem>.adjacency.json' describing the SHARED / DIRICHLET graph."
            ),
        ),
    ] = None,
    adjacency: Annotated[
        list[Path] | None,
        typer.Option(
            "--adjacency",
            help=(
                "Input .adjacency.json sidecar(s) describing the pre-existing "
                "inter-block seams of the input .geo files.  Repeatable (one per "
                "component).  Without this, .geo files carry only geometry, so "
                "every face defaults to FREE and the topology's SHARED/DIRICHLET "
                "seams are LOST through the split (they become FREE -> chimera in "
                "cc.par -> overset hunts donors across a 1-to-1 seam).  Pass each "
                "component's sidecar to preserve its seams onto the right chunks."
            ),
        ),
    ] = None,
    verbose: _VERBOSE = False,
    journal: _JOURNAL = None,
) -> None:
    """Split structured blocks along ijk vertex indices, preserving adjacency.

    Reads existing ``.geo`` files, applies the splits described in the
    TOML config (one ``[[splits]]`` entry per block), writes new
    ``.geo`` files and a JSON sidecar with the adjacency graph.

    Pass ``--adjacency`` with each input component's sidecar so the
    pre-existing inter-block seams (SHARED / DIRICHLET) are wired onto
    the input blocks before splitting and thus preserved on the output
    chunks; otherwise they are lost (the .geo files are geometry only).
    """
    if not inputs:
        _err.print("[red]Error:[/red] no input .geo files given.")
        raise typer.Exit(1) from None
    for p in inputs:
        if not p.exists():
            _err.print(f"[red]Error:[/red] file not found: {p}")
            raise typer.Exit(1) from None
        if p.suffix.lower() != ".geo":
            _err.print(f"[red]Error:[/red] input is not a .geo file: {p}")
            raise typer.Exit(1) from None

    try:
        cfg = load_split_config(config)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    # Blocks from the input .geo files (labels: file stems minus their common
    # ``<base>_`` prefix, e.g. ``sphere_sub0.geo`` -> ``sub0``), with the
    # pre-existing inter-block seams wired from the input sidecar(s) so the
    # split preserves them (the geometry-only .geo blocks are all FREE).
    from shore.pipeline import common_base_stem

    base_stem = common_base_stem([p.stem for p in inputs])
    proj = Project()
    try:
        proj.load_blocks(inputs, adjacency or ())
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    assert proj.mesh is not None
    if verbose:
        _console.print(f"Read {len(proj.mesh.blocks)} input block(s):")
        for blk in proj.mesh.blocks:
            _console.print(f"  {blk.label:24s} {blk.nodes.shape}")
        if adjacency:
            n_seam = sum(
                1
                for blk in proj.mesh.blocks
                for f in blk.faces()
                if f.bc.name in ("SHARED", "DIRICHLET", "PERIODIC")
            )
            _console.print(f"Wired {n_seam} connection face(s) from {len(adjacency)} sidecar(s).")

    # split_blocks handles several axes per block and propagates cuts across
    # SHARED/DIRICHLET seams, so the directives go straight through.
    directives = [(d.label, d.axis, d.at) for d in cfg.splits]
    try:
        proj.split(directives)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    # Resolve output stem (CLI > config > default).
    if output is None:
        output = (
            cfg.output if cfg.output is not None else (inputs[0].parent / (base_stem or "split"))
        )

    # One .geo per output block, plus the adjacency JSON.
    proj.write_geo(output, per_block=True)
    assert proj.mesh is not None
    if verbose:
        for blk in proj.mesh.blocks:
            _console.print(f"  wrote {output.parent / f'{output.name}_{blk.label}.geo'}")
    adj_path = output.parent / f"{output.name}.adjacency.json"
    proj.write_adjacency(adj_path)
    if journal is not None:
        proj.save_journal(journal)
    _console.print(
        f"[green]✓[/green] Wrote {len(proj.mesh.blocks)} blocks "
        f"+ {adj_path.name} ({len(directives)} split directive(s))"
    )


# ---------------------------------------------------------------------------
# shore grd
# ---------------------------------------------------------------------------


@app.command(name="grd")
def grd(
    inputs: Annotated[
        list[Path],
        typer.Argument(
            help=(
                "Input .geo files to pack, one per block.  Block ordering in "
                "the .grd matches the order given here."
            ),
        ),
    ],
    output: Annotated[
        Path,
        typer.Option("-o", "--output", help="Output .grd file path."),
    ],
    journal: _JOURNAL = None,
) -> None:
    """Pack one or more .geo files into a single Xall binary .grd file.

    The .grd format is Fortran unformatted sequential binary: a block-count
    record, per-block cell-count triples (ni-1, nj-1, nk-1), then three
    coordinate records (X, Y, Z) per block in Fortran column-major order.

    Block ordering in the .grd matches the order of the input files, which
    must agree with the ordering used when generating the companion cc.par.
    """
    if not inputs:
        _err.print("[red]Error:[/red] no input .geo files given.")
        raise typer.Exit(1) from None
    for p in inputs:
        if not p.exists():
            _err.print(f"[red]Error:[/red] file not found: {p}")
            raise typer.Exit(1) from None
        if p.suffix.lower() != ".geo":
            _err.print(f"[red]Error:[/red] input is not a .geo file: {p}")
            raise typer.Exit(1) from None

    proj = Project()
    try:
        proj.write_grd(output, inputs)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    if journal is not None:
        proj.save_journal(journal)

    _console.print(f"[green]✓[/green] Wrote {proj.last_result} ({len(inputs)} blocks)")


# ---------------------------------------------------------------------------
# shore cc-par
# ---------------------------------------------------------------------------


@app.command(name="cc-par")
def cc_par(
    adjacency: Annotated[
        Path,
        typer.Argument(
            help="Input .adjacency.json sidecar (produced by shore split or shore mesh).",
        ),
    ],
    output: Annotated[
        Path,
        typer.Argument(
            help="Output cc.par path.",
        ),
    ],
    grd_filename: Annotated[
        str,
        typer.Option(
            "--grd",
            help="Name of the .grd file Xall should pair with this cc.par "
            "(written into header line 1).",
        ),
    ],
    output_basename: Annotated[
        str,
        typer.Option(
            "--output-basename",
            help="Output basename Xall uses for derived files (header line 2).",
        ),
    ] = "cc.",
    free_face_bc: Annotated[
        int,
        typer.Option(
            "--free-bc",
            help="BC code emitted for FREE faces on every block that has no "
            "entry in --meta.  Default 40 (chimera cell-center, body-fitted/"
            "buffer meshes).  Use a negative natural BC for far-field background "
            "blocks: -3 inflow, -4 inoutflow, -9 extrapolation.",
        ),
    ] = 40,
    block_meta_json: Annotated[
        str | None,
        typer.Option(
            "--meta",
            help="Per-block BlockMeta overrides as a JSON object or a path to a "
            "JSON file.  Keys are block labels; values are objects with any "
            "subset of {level, group, priority, comment, free_face_bc}.  "
            "Missing keys inherit the block's defaults (level=1, group=0, "
            "priority=1, free_face_bc=--free-bc).  "
            'Example: \'{"cap_north": {"level": 2, "free_face_bc": -9}}\'',
        ),
    ] = None,
    boxes_json: Annotated[
        str | None,
        typer.Option(
            "--boxes",
            help="Hexahedral marker boxes for the cc.par boxes section, as a "
            "JSON list or a path to a JSON file.  Each entry is an object "
            "with keys 'associated_block' (block label) and 'vertices' "
            "(list of 8 [x, y, z] triples in canonical i-fastest order).  "
            "See shore.io.boxes for the vertex-order convention.",
        ),
    ] = None,
    multigrid_levels: Annotated[
        tuple[int, int],
        typer.Option(
            "--mgl",
            help="Multigrid levels: 'N_levels finest_level_active'.  Header line 7.",
        ),
    ] = (4, 1),
    boundary_layer_thickness: Annotated[
        float,
        typer.Option(
            "--boundary-layer-thickness",
            help="Header boundary-layer thickness (line 9).",
        ),
    ] = -0.1,
    numerical_beach: Annotated[
        float,
        typer.Option(
            "--numerical-beach",
            help="Header numerical-beach extension (line 10).",
        ),
    ] = -1.0,
    journal: _JOURNAL = None,
) -> None:
    """Convert a SHORE .adjacency.json sidecar to an Xall cc.par file.

    Maps SHORE's per-face boundary conditions to Xall cc.par codes:
    WALL -> 1 (parete), FREE -> --free-bc (default 40), SHARED / DIRICHLET
    -> connection orientation codes with paired patch indices.  Cubed-sphere
    cap-equator cross-axis seams use the hard-coded orientation table from
    shore.io.cc_par.

    Use --meta to override per-block metadata (level, group, priority,
    free_face_bc) for selected blocks.  Any block not listed inherits the
    global --free-bc default.
    """
    import json as _json

    from shore.exceptions import ParameterError

    if not adjacency.exists():
        _err.print(f"[red]Error:[/red] {adjacency} not found.")
        raise typer.Exit(1)

    # Parse --meta: try inline JSON first; if that fails treat the value as
    # a file path.  This order avoids the ambiguity where a valid JSON string
    # that also happens to match a filename in the working directory would
    # silently be misinterpreted as a path.
    meta_overrides: dict[str, dict] = {}
    if block_meta_json is not None:
        try:
            meta_overrides = _json.loads(block_meta_json)
        except _json.JSONDecodeError:
            meta_path = Path(block_meta_json)
            if not meta_path.exists():
                _err.print(
                    f"[red]Error:[/red] --meta is not valid JSON and not an existing file: "
                    f"{block_meta_json!r}"
                )
                raise typer.Exit(1) from None
            try:
                meta_overrides = _json.loads(meta_path.read_text())
            except _json.JSONDecodeError as exc:
                _err.print(f"[red]Error:[/red] --meta file {meta_path} is not valid JSON: {exc}")
                raise typer.Exit(1) from None
        if not isinstance(meta_overrides, dict):
            _err.print(
                "[red]Error:[/red] --meta must be a JSON object mapping block labels to "
                f"BlockMeta fields, got {type(meta_overrides).__name__}"
            )
            raise typer.Exit(1) from None

    # Parse --boxes the same way as --meta: inline JSON first, then file path.
    boxes_payload: list[dict] | None = None
    if boxes_json is not None:
        try:
            boxes_payload = _json.loads(boxes_json)
        except _json.JSONDecodeError:
            boxes_path = Path(boxes_json)
            if not boxes_path.exists():
                _err.print(
                    f"[red]Error:[/red] --boxes is not valid JSON and not an existing "
                    f"file: {boxes_json!r}"
                )
                raise typer.Exit(1) from None
            try:
                boxes_payload = _json.loads(boxes_path.read_text())
            except _json.JSONDecodeError as exc:
                _err.print(f"[red]Error:[/red] --boxes file {boxes_path} is not valid JSON: {exc}")
                raise typer.Exit(1) from None
        if not isinstance(boxes_payload, list):
            _err.print(
                "[red]Error:[/red] --boxes must be a JSON list of "
                "{associated_block, vertices} entries, got "
                f"{type(boxes_payload).__name__}"
            )
            raise typer.Exit(1) from None

    if boxes_payload is not None:
        from shore.io.boxes import box_from_json

        try:
            for entry in boxes_payload:
                box_from_json(entry)  # validate up front, for a --boxes error message
        except (ParameterError, TypeError, ValueError) as exc:
            _err.print(f"[red]Error:[/red] --boxes: {_esc(str(exc))}")
            raise typer.Exit(1) from None

    proj = Project()
    try:
        proj.write_cc_par(
            adjacency,
            output,
            grd_filename=grd_filename,
            output_basename=output_basename,
            free_face_bc=free_face_bc,
            multigrid_levels=multigrid_levels,
            boundary_layer_thickness=boundary_layer_thickness,
            numerical_beach=numerical_beach,
            meta=meta_overrides,
            boxes=boxes_payload or (),
        )
    except (ParameterError, ValueError) as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    if journal is not None:
        proj.save_journal(journal)
    _console.print(f"[green]✓[/green] Wrote {proj.last_result}")


# ---------------------------------------------------------------------------
# shore grd-merge
# ---------------------------------------------------------------------------


@app.command(name="grd-merge")
def grd_merge(
    inputs: Annotated[
        list[Path],
        typer.Argument(
            help=(
                "Two-or-more input .grd files, in the desired output order.  "
                "Use the same ordering when running shore cc-par-merge."
            ),
        ),
    ],
    output: Annotated[
        Path,
        typer.Option("-o", "--output", help="Output merged .grd file path."),
    ],
) -> None:
    """Concatenate multiple .grd files into one for an Xall multi-mesh run.

    Used when an Xall job combines several independent Chimera components
    (e.g. wall + background): the merged .grd lists every block from
    inputs[0] first, then inputs[1], etc., with the header block count
    rewritten to the sum.  Block payloads are passed through byte-for-byte.
    """
    from shore.io.grd import merge_grd

    if len(inputs) < 2:
        _err.print("[red]Error:[/red] need at least 2 input .grd files.")
        raise typer.Exit(1) from None
    for p in inputs:
        if not p.exists():
            _err.print(f"[red]Error:[/red] file not found: {p}")
            raise typer.Exit(1) from None

    try:
        out = merge_grd(inputs, output)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    _console.print(f"[green]✓[/green] Wrote {out} ({len(inputs)} inputs merged)")


# ---------------------------------------------------------------------------
# shore cc-par-merge
# ---------------------------------------------------------------------------


@app.command(name="cc-par-merge")
def cc_par_merge(
    inputs: Annotated[
        list[Path],
        typer.Argument(
            help=(
                "Two-or-more input cc.par files, in the desired output order.  "
                "Must match the ordering used for shore grd-merge."
            ),
        ),
    ],
    output: Annotated[
        Path,
        typer.Option("-o", "--output", help="Output merged cc.par path."),
    ],
    grd_filename: Annotated[
        str,
        typer.Option(
            "--grd",
            help="Name of the merged .grd file (written into header line 1, "
            "overriding inputs[0]'s value).",
        ),
    ],
) -> None:
    """Concatenate multiple cc.par files into one for an Xall multi-mesh run.

    Block indices and patch family pointers in inputs[1:] are shifted by
    the cumulative block / patch counts of the preceding inputs so the
    merged file is internally consistent.  Header fields other than the
    .grd filename are taken from inputs[0]; if any input disagrees, a
    warning is printed but the merge proceeds.

    Block-section comments are auto-prefixed with each input's stem (e.g.
    'wall.sub0', 'background.center') so blocks from different topologies
    are distinguishable in the merged file.
    """
    from shore.io.cc_par import merge_cc_par

    if len(inputs) < 2:
        _err.print("[red]Error:[/red] need at least 2 input cc.par files.")
        raise typer.Exit(1) from None
    for p in inputs:
        if not p.exists():
            _err.print(f"[red]Error:[/red] file not found: {p}")
            raise typer.Exit(1) from None

    try:
        out = merge_cc_par(inputs, output, grd_filename=grd_filename)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    _console.print(f"[green]✓[/green] Wrote {out} ({len(inputs)} inputs merged)")


# ---------------------------------------------------------------------------
# shore proc-input
# ---------------------------------------------------------------------------


@app.command(name="proc-input")
def proc_input(
    grd: Annotated[
        Path,
        typer.Argument(help="Path to the .grd file (block ordering and weights)."),
    ],
    output: Annotated[
        Path,
        typer.Argument(help="Output proc.input path."),
    ],
    np: Annotated[
        int,
        typer.Option(
            "--np",
            help=(
                "Number of MPI ranks (>= 1).  Default 1 (serial).  When > 1, "
                "blocks are distributed by greedy weight-balanced assignment "
                "across range(np) unless overridden via --meta."
            ),
        ),
    ] = 1,
    block_meta_json: Annotated[
        str | None,
        typer.Option(
            "--meta",
            help=(
                "Per-block overrides for proc / group / body.  Inline JSON or "
                "path to a JSON file.  Keys are 1-based block indices "
                "(strings); values are objects with any subset of "
                "{proc, group, body}.  Missing fields use defaults: proc from "
                "greedy balance, group=0, body=0."
            ),
        ),
    ] = None,
    header_comment: Annotated[
        str,
        typer.Option(
            "--header-comment",
            help="Free-form text written on header line 1 of proc.input.",
        ),
    ] = "generated by SHORE",
    journal: _JOURNAL = None,
) -> None:
    """Write a proc.input file from a .grd, with greedy weight-balanced ranks.

    proc.input is the runtime work-distribution descriptor Xall reads at MPI
    startup: it maps each block to an MPI rank and carries the per-block
    group / body indices the solver's force/moment integration needs.

    Block ordering and weights come from the .grd (cell counts ni*nj*nk per
    block).  Without --meta, proc=0 for --np 1 and a greedy weight-balanced
    distribution otherwise; group=0 and body=0 for every block.

    Use --meta to override per-block values:

        --meta '{"1": {"proc": 0, "body": 1}, "5": {"group": 2, "body": 1}}'

    A balance summary is printed on stderr for --np > 1.
    """
    if np < 1:
        _err.print(f"[red]Error:[/red] --np must be >= 1, got {np}")
        raise typer.Exit(1) from None
    if not grd.exists():
        _err.print(f"[red]Error:[/red] .grd file not found: {grd}")
        raise typer.Exit(1) from None

    import json as _json

    # Parse --meta: inline JSON, then path fallback.  Same convention as
    # `shore cc-par --meta`.
    overrides: dict[str, dict] = {}
    if block_meta_json is not None:
        try:
            overrides = _json.loads(block_meta_json)
        except _json.JSONDecodeError:
            meta_path = Path(block_meta_json)
            if not meta_path.exists():
                _err.print(
                    f"[red]Error:[/red] --meta is not valid JSON and not an existing file: "
                    f"{block_meta_json!r}"
                )
                raise typer.Exit(1) from None
            try:
                overrides = _json.loads(meta_path.read_text())
            except _json.JSONDecodeError as exc:
                _err.print(f"[red]Error:[/red] --meta file {meta_path} is not valid JSON: {exc}")
                raise typer.Exit(1) from None
        if not isinstance(overrides, dict):
            _err.print(
                "[red]Error:[/red] --meta must be a JSON object mapping block-index "
                f"strings to {{proc, group, body}} dicts, got {type(overrides).__name__}"
            )
            raise typer.Exit(1) from None

    proc_assignment: dict[int, int] = {}
    group_assignment: dict[int, int] = {}
    body_assignment: dict[int, int] = {}
    for k_str, v in overrides.items():
        try:
            k = int(k_str)
        except (TypeError, ValueError):
            _err.print(f"[red]Error:[/red] --meta key {k_str!r} is not an integer block index.")
            raise typer.Exit(1) from None
        if not isinstance(v, dict):
            _err.print(
                f"[red]Error:[/red] --meta value for block {k} must be a JSON object, "
                f"got {type(v).__name__}"
            )
            raise typer.Exit(1) from None
        if "proc" in v:
            proc_assignment[k] = int(v["proc"])
        if "group" in v:
            group_assignment[k] = int(v["group"])
        if "body" in v:
            body_assignment[k] = int(v["body"])

    proj = Project()
    try:
        proj.write_proc_input(grd, output, np=np, meta=overrides, header_comment=header_comment)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    out, summary = proj.last_result
    if journal is not None:
        proj.save_journal(journal)

    _console.print(f"[green]✓[/green] Wrote {out}")
    if np > 1:
        _err.print(summary.format())
        if summary.imbalance > 0.05:
            _err.print(
                f"[yellow]Warning:[/yellow] load imbalance is "
                f"{summary.imbalance * 100:.1f}% (> 5%).  "
                "Run `shore balance` to plan a topology-safe set of splits "
                "that brings the imbalance within tolerance."
            )


# ---------------------------------------------------------------------------
# shore proc-input-merge
# ---------------------------------------------------------------------------


@app.command(name="proc-input-merge")
def proc_input_merge(
    inputs: Annotated[
        list[Path],
        typer.Argument(
            help=(
                "Two-or-more input proc.input files, in the desired output order.  "
                "Must match the ordering used for shore grd-merge."
            ),
        ),
    ],
    output: Annotated[
        Path,
        typer.Option("-o", "--output", help="Output merged proc.input path."),
    ],
    grd_path: Annotated[
        Path | None,
        typer.Option(
            "--grd",
            help=(
                "Optional path to the merged .grd.  When given, the printed "
                "balance summary uses real cell counts and shows imbalance.  "
                "Without it, only per-rank block counts are reported "
                "(proc.input does not store per-block weights)."
            ),
        ),
    ] = None,
    header_comment: Annotated[
        str,
        typer.Option(
            "--header-comment",
            help="Free-form text written on header line 1 of the merged file.",
        ),
    ] = "merged by SHORE",
) -> None:
    """Concatenate multiple proc.input files into one.

    Block indices in inputs[1:] are shifted by the cumulative block count
    of the preceding inputs.  group, body, and proc columns pass through
    verbatim — the merger does NOT re-run load balancing or renumber bodies
    across components.

    Invoke shore proc-input-merge alongside shore grd-merge with the same
    input ordering so the merged proc.input matches the merged .grd.
    Pass --grd to read the merged .grd's per-block cell counts so the
    balance summary printed at the end shows real cell-count imbalance.
    """
    if len(inputs) < 2:
        _err.print("[red]Error:[/red] need at least 2 input proc.input files.")
        raise typer.Exit(1) from None
    for p in inputs:
        if not p.exists():
            _err.print(f"[red]Error:[/red] file not found: {p}")
            raise typer.Exit(1) from None
    if grd_path is not None and not grd_path.exists():
        _err.print(f"[red]Error:[/red] .grd file not found: {grd_path}")
        raise typer.Exit(1) from None

    from shore.io.proc_input import merge_proc_input

    try:
        out, summary = merge_proc_input(
            inputs,
            output,
            grd_path=grd_path,
            header_comment=header_comment,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    _console.print(f"[green]✓[/green] Wrote {out} ({len(inputs)} inputs merged)")
    if summary.np > 1:
        _err.print(summary.format())


# ---------------------------------------------------------------------------
# shore balance
# ---------------------------------------------------------------------------


@app.command(name="balance")
def balance(
    grd: Annotated[
        Path,
        typer.Argument(help="Path to the .grd file (block ordering and weights)."),
    ],
    output: Annotated[
        Path,
        typer.Option("-o", "--output", help="Output split.toml path."),
    ],
    np: Annotated[
        int,
        typer.Option(
            "--np",
            help="Number of MPI ranks (>= 1).  The planner balances against this rank count.",
        ),
    ],
    blocks_per_rank: Annotated[
        float,
        typer.Option(
            "--blocks-per-rank",
            help=(
                "Target average blocks per rank (default 1.0).  Higher values give "
                "finer granularity and better balance, at more blocks (and more "
                "overset donor-search work)."
            ),
        ),
    ] = 1.0,
    ngr: Annotated[
        int,
        typer.Option(
            "--ngr",
            help=(
                "Multigrid level count the mesh must support (default 1).  Every "
                "chunk axis is kept divisible by 2^(ngr-1)."
            ),
        ),
    ] = 1,
    min_coarse_cells: Annotated[
        int,
        typer.Option(
            "--min-coarse-cells",
            help=(
                "Minimum cells per axis on the coarsest (level-ngr) grid (default 8).  "
                "Every chunk axis is kept >= min_coarse_cells * 2^(ngr-1) fine cells."
            ),
        ),
    ] = 8,
    tolerance: Annotated[
        float,
        typer.Option(
            "--tolerance",
            help=(
                "Target fractional imbalance (max-min)/mean for the 'converged' "
                "report (default 0.05).  The planner reports achieved imbalance "
                "regardless; it never fails on an unreachable tolerance."
            ),
        ),
    ] = 0.05,
    adjacency: Annotated[
        Path | None,
        typer.Option(
            "--adjacency",
            help=(
                "Path to a .adjacency.json sidecar.  STRONGLY recommended: it lets "
                "the planner enforce the multigrid floor against seam-PROPAGATED "
                "cuts (shore split mirrors cuts across SHARED/DIRICHLET seams).  "
                "Without it the floor is only checked per block in isolation."
            ),
        ),
    ] = None,
    plan_format: Annotated[
        str,
        typer.Option(
            "--format",
            help=(
                "Output format: 'toml' (default, for the PRE-overset `shore split`) "
                "or 'bsplit-par' (for the POST-overset external `bsplit`, which "
                "splits framed overset outputs cc.0N/cc.0N.grd).  'bsplit-par' "
                "requires a grd-only plan (do not pass --adjacency)."
            ),
        ),
    ] = "toml",
) -> None:
    """Plan a multi-way equal-partition split for MPI load balancing.

    Partitions every block into ~equal, ~cubic chunks sized so the total
    lands near np*blocks-per-rank, keeping every chunk axis >=
    min-coarse-cells * 2^(ngr-1) cells (divisible by 2^(ngr-1)).  Cut counts
    are chosen per-axis so seam-linked blocks (equal extent) get identical
    cuts — conformal by construction; `shore split` then propagates and
    executes them.  Writes a split.toml in the `shore split` schema.

    Typical use:

        shore balance wall.grd -o splits.toml --np 16 --ngr 3 \\
            --blocks-per-rank 8 --adjacency wall.adjacency.json
        shore split -c splits.toml wall_*.geo --adjacency wall.adjacency.json
        shore grd wall_split_*.geo -o wall.grd
        shore proc-input wall.grd wall.proc.input --np 16

    Best-effort: when the floor cannot be met (e.g. an input block already
    below the floor), it warns and emits the best plan.
    """
    if plan_format not in ("toml", "bsplit-par"):
        _err.print(
            f"[red]Error:[/red] --format must be 'toml' or 'bsplit-par', got {plan_format!r}"
        )
        raise typer.Exit(1) from None
    if plan_format == "bsplit-par" and adjacency is not None:
        _err.print(
            "[red]Error:[/red] --format bsplit-par needs a grd-only plan; "
            "do not pass --adjacency (bsplit needs grd block numbers)."
        )
        raise typer.Exit(1) from None
    if not grd.exists():
        _err.print(f"[red]Error:[/red] .grd file not found: {grd}")
        raise typer.Exit(1) from None
    if adjacency is not None and not adjacency.exists():
        _err.print(f"[red]Error:[/red] adjacency file not found: {adjacency}")
        raise typer.Exit(1) from None

    from shore.balance import plan_balance

    # Build lightweight wired blocks from the sidecar so the planner can
    # floor-check seam-propagated cuts (no .geo node data needed — only
    # shapes + the seam graph).
    sim_blocks = None
    if adjacency is not None:
        try:
            sim_blocks = _blocks_from_adjacency_shapes(adjacency)
        except Exception as exc:
            _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
            raise typer.Exit(1) from None

    try:
        plan = plan_balance(
            grd,
            np=np,
            blocks_per_rank=blocks_per_rank,
            ngr=ngr,
            min_coarse_cells=min_coarse_cells,
            tolerance=tolerance,
            adjacency_path=adjacency,
            blocks=sim_blocks,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    if plan_format == "bsplit-par":
        try:
            text = plan.to_bsplit_par(ngr=ngr)
        except Exception as exc:
            _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
            raise typer.Exit(1) from None
    else:
        text = plan.to_toml()
    output.write_text(text, encoding="utf-8")
    _console.print(f"[green]✓[/green] Wrote {output} ({plan_format})")

    blocks_touched = len({s[0] for s in plan.splits})
    _err.print(
        f"shore balance: {plan.initial_summary.n_blocks} blocks "
        f"(imbalance {plan.initial_summary.imbalance * 100:.1f}%) → "
        f"{plan.final_summary.n_blocks} blocks "
        f"(imbalance {plan.final_summary.imbalance * 100:.2f}%) "
        f"via {len(plan.splits)} cut(s) across {blocks_touched} block(s)"
    )
    if not plan.converged:
        _err.print(
            f"[yellow]Note:[/yellow] achieved imbalance "
            f"{plan.final_summary.imbalance * 100:.2f}% > tolerance "
            f"{tolerance * 100:.1f}% (floor-limited or finer --blocks-per-rank needed)."
        )


# ---------------------------------------------------------------------------
# shore boxes-from-stl
# ---------------------------------------------------------------------------


@app.command(name="boxes-from-stl")
def boxes_from_stl_cmd(
    stl: Annotated[
        Path,
        typer.Argument(
            help="Body STL file to enclose in marker boxes.",
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "-o",
            "--output",
            help="Output JSON file in the shore cc-par --boxes schema.",
        ),
    ],
    associated_block: Annotated[
        str,
        typer.Option(
            "--block",
            help="SHORE block label every generated box references.  Resolved "
            "to a 1-based index when shore cc-par writes the box section.",
        ),
    ],
    max_outward_gap: Annotated[
        float,
        typer.Option(
            "--max-outward-gap",
            help="Max distance (STL units) from any box-surface point to the "
            "nearest STL surface, measured outward only.  Strictly positive; "
            "should be smaller than the body's near-mesh wall-normal thickness.  "
            "voxel-fill chooses voxel size = 2 * gap / sqrt(3) so every "
            "voxel-corner stays within the gap of the body surface.",
        ),
    ],
    method: Annotated[
        str,
        typer.Option(
            "--method",
            help="Box-generation algorithm: 'voxel-fill' (default — voxelise "
            "the body's AABB, emit each inside-body voxel as a Box; full "
            "body-volume coverage), 'aabb' (single world AABB, best effort), "
            "'obb' (single oriented bounding box, best effort).",
        ),
    ] = "voxel-fill",
    face_samples_per_axis: Annotated[
        int,
        typer.Option(
            "--face-samples-per-axis",
            help="aabb / obb only.  Box-surface sampling density for the "
            "outward-gap post-check: each face sampled on a k-by-k grid "
            "(total = 8 corners + 6*k*k samples per box).  Default 4 is "
            "fast and accurate; raise to 6-8 for sharply concave geometries.  "
            "Ignored by voxel-fill.",
        ),
    ] = 4,
) -> None:
    """Generate marker boxes from a body STL, ready for shore cc-par --boxes.

    Output is a JSON list of {associated_block, vertices} entries —
    the same schema cc-par --boxes consumes.  Pipe the two together
    in a typical workflow:

        shore boxes-from-stl sphere.stl --block sub0 --max-outward-gap 0.10 \\
                             -o sphere_boxes.json
        shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \\
                     --boxes sphere_boxes.json
    """
    import json as _json

    from shore.boxes import boxes_from_stl
    from shore.exceptions import MeshInputError, ParameterError
    from shore.io.boxes import box_to_json

    if not stl.exists():
        _err.print(f"[red]Error:[/red] STL file not found: {stl}")
        raise typer.Exit(1) from None

    if method not in ("aabb", "obb", "voxel-fill"):
        _err.print(
            f"[red]Error:[/red] --method must be 'aabb', 'obb', or 'voxel-fill'; got {method!r}"
        )
        raise typer.Exit(1) from None

    try:
        boxes = boxes_from_stl(
            stl,
            associated_block,
            max_outward_gap=max_outward_gap,
            method=method,  # type: ignore[arg-type]
            face_samples_per_axis=face_samples_per_axis,
        )
    except (ParameterError, MeshInputError) as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    payload = [box_to_json(box) for box in boxes]
    output.write_text(_json.dumps(payload, indent=2) + "\n")
    _console.print(f"[green]✓[/green] Wrote {output} ({len(boxes)} box(es), method={method!r})")


# ---------------------------------------------------------------------------
# shore boxes-export
# ---------------------------------------------------------------------------


@app.command(name="boxes-export")
def boxes_export(
    boxes_json: Annotated[
        Path,
        typer.Argument(
            help=(
                "Input boxes JSON file in the shore cc-par --boxes / shore boxes-from-stl schema."
            ),
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "-o",
            "--output",
            help="Output .vtm path for ParaView.",
        ),
    ],
) -> None:
    """Export marker boxes to a PyVista MultiBlock (.vtm) for ParaView.

    Reads the same JSON format ``shore cc-par --boxes`` consumes and
    ``shore boxes-from-stl`` produces; writes one ``vtkUnstructuredGrid``
    cell (``VTK_HEXAHEDRON``) per box wrapped in a MultiBlock so each
    box can be selected individually in ParaView.  Each cell carries
    ``associated_block`` (string) and ``box_index`` (int) cell data
    for colouring or filtering.

    Closes the loop alongside ``shore boxes-from-stl``::

        shore boxes-from-stl sphere.stl --block sub0 --max-outward-gap 0.05 \\
                             -o sphere_boxes.json
        shore boxes-export sphere_boxes.json -o sphere_boxes.vtm
    """
    import json as _json

    from shore.exceptions import ParameterError
    from shore.io.boxes import box_from_json

    if not boxes_json.exists():
        _err.print(f"[red]Error:[/red] {boxes_json} not found.")
        raise typer.Exit(1) from None

    try:
        payload = _json.loads(boxes_json.read_text())
    except _json.JSONDecodeError as exc:
        _err.print(f"[red]Error:[/red] {boxes_json} is not valid JSON: {exc}")
        raise typer.Exit(1) from None
    if not isinstance(payload, list):
        _err.print(
            "[red]Error:[/red] boxes JSON must be a list of "
            "{associated_block, vertices} entries, got "
            f"{type(payload).__name__}"
        )
        raise typer.Exit(1) from None
    if not payload:
        _err.print(f"[red]Error:[/red] {boxes_json} contains zero boxes.")
        raise typer.Exit(1) from None

    try:
        boxes = [box_from_json(entry) for entry in payload]
    except (ParameterError, TypeError, ValueError) as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    try:
        from shore.io.vtk import write_boxes_vtm
    except ImportError as exc:
        _err.print(
            f"[red]Error:[/red] VTK export requires pyvista — install "
            f"with `pip install shore[vis]`.  ({_esc(str(exc))})"
        )
        raise typer.Exit(1) from None

    out = write_boxes_vtm(boxes, output)
    _console.print(f"[green]✓[/green] Wrote {out} ({len(boxes)} box(es))")


# ---------------------------------------------------------------------------
# shore profile  (2D airfoil section coordinates)
# ---------------------------------------------------------------------------


profile_app = typer.Typer(
    name="profile",
    help="2D airfoil section coordinate generators (NACA 4-digit, 5-digit).",
    no_args_is_help=True,
)
app.add_typer(profile_app)


@profile_app.command("naca")
def profile_naca(
    designation: Annotated[
        str,
        typer.Argument(help="NACA digit string (e.g. '0012', '2412', '23012')."),
    ],
    output: Annotated[
        Path,
        typer.Option("-o", "--output", help="Output .dat file (Selig format)."),
    ],
    n_pts: Annotated[
        int,
        typer.Option("--n-pts", help="Chord samples per side (LE..TE)."),
    ] = 80,
    closed_te: Annotated[
        bool,
        typer.Option(
            "--closed-te/--open-te",
            help="Use the closed-TE thickness polynomial.",
        ),
    ] = True,
    name: Annotated[
        str | None,
        typer.Option("--name", help="Optional name comment on the .dat header."),
    ] = None,
) -> None:
    """Generate a NACA airfoil section and write it as a Selig .dat file.

    The output is a ``(2 * n_pts - 1, 2)`` polyline ordered upper-TE
    -> LE -> lower-TE, written as whitespace-separated ``x y`` rows.
    Compatible with the UIUC airfoil database format.

    Examples
    --------

        shore profile naca 0012 -o naca0012.dat
        shore profile naca 2412 --n-pts 120 -o naca2412.dat
        shore profile naca 23012 -o naca23012.dat --name "NACA 23012"
    """
    from shore.exceptions import ParameterError
    from shore.profiles import naca, write_section

    try:
        coords = naca(designation, n_pts=n_pts, closed_te=closed_te)
    except ParameterError as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    name_label = name if name is not None else f"NACA {designation}"
    write_section(output, coords, name=name_label)
    _console.print(
        f"[green]✓[/green] Wrote {output}  ({coords.shape[0]} points, "
        f"{'closed' if closed_te else 'open'} TE)"
    )


# ---------------------------------------------------------------------------
# shore body  (3D STL bodies)
# ---------------------------------------------------------------------------


body_app = typer.Typer(
    name="body",
    help="3D STL body generators (analytic primitives + section extrusion).",
    no_args_is_help=True,
)
app.add_typer(body_app)


def _write_stl(mesh_obj, output: Path, *, label: str) -> None:
    """Export a trimesh.Trimesh to *output* and print a summary."""
    output.parent.mkdir(parents=True, exist_ok=True)
    mesh_obj.export(str(output))
    _console.print(
        f"[green]✓[/green] Wrote {output}  "
        f"({len(mesh_obj.vertices)} vertices, {len(mesh_obj.faces)} faces, "
        f"{label})"
    )


@body_app.command("sphere")
def body_sphere(
    output: Annotated[Path, typer.Option("-o", "--output", help="Output .stl file.")],
    radius: Annotated[float, typer.Option("--radius", help="Sphere radius.")] = 1.0,
    subdivisions: Annotated[
        int,
        typer.Option(
            "--subdivisions",
            help="Icosphere subdivision level (0..6+, each level x4 triangles).",
        ),
    ] = 3,
) -> None:
    """Build an icosahedral sphere STL centred at the origin."""
    from shore.exceptions import ParameterError
    from shore.shapes import sphere

    try:
        mesh_obj = sphere(radius=radius, subdivisions=subdivisions)
    except ParameterError as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    _write_stl(mesh_obj, output, label="sphere")


@body_app.command("cylinder")
def body_cylinder(
    output: Annotated[Path, typer.Option("-o", "--output", help="Output .stl file.")],
    radius: Annotated[float, typer.Option("--radius", help="Cylinder radius.")] = 1.0,
    height: Annotated[float, typer.Option("--height", help="Cylinder height (z-extent).")] = 1.0,
    sections: Annotated[
        int, typer.Option("--sections", help="Segments around the circumference.")
    ] = 32,
) -> None:
    """Build a right circular cylinder STL, axis along z."""
    from shore.exceptions import ParameterError
    from shore.shapes import cylinder

    try:
        mesh_obj = cylinder(radius=radius, height=height, sections=sections)
    except ParameterError as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    _write_stl(mesh_obj, output, label="cylinder")


@body_app.command("torus")
def body_torus(
    output: Annotated[Path, typer.Option("-o", "--output", help="Output .stl file.")],
    r_major: Annotated[
        float, typer.Option("--r-major", help="Major radius (centre of tube).")
    ] = 1.0,
    r_minor: Annotated[float, typer.Option("--r-minor", help="Minor radius (tube radius).")] = 0.25,
    major_sections: Annotated[
        int, typer.Option("--major-sections", help="Segments around the major circle.")
    ] = 32,
    minor_sections: Annotated[
        int, typer.Option("--minor-sections", help="Segments around the tube.")
    ] = 32,
) -> None:
    """Build a torus STL (axis of revolution = z)."""
    from shore.exceptions import ParameterError
    from shore.shapes import torus

    try:
        mesh_obj = torus(
            r_major=r_major,
            r_minor=r_minor,
            major_sections=major_sections,
            minor_sections=minor_sections,
        )
    except ParameterError as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    _write_stl(mesh_obj, output, label="torus")


@body_app.command("box")
def body_box(
    output: Annotated[Path, typer.Option("-o", "--output", help="Output .stl file.")],
    extents: Annotated[
        tuple[float, float, float],
        typer.Option("--extents", help="Side lengths (dx dy dz)."),
    ] = (1.0, 1.0, 1.0),
) -> None:
    """Build a rectangular box STL centred at the origin."""
    from shore.exceptions import ParameterError
    from shore.shapes import box as _box

    try:
        mesh_obj = _box(extents=extents)
    except ParameterError as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    _write_stl(mesh_obj, output, label="box")


@body_app.command("capsule")
def body_capsule(
    output: Annotated[Path, typer.Option("-o", "--output", help="Output .stl file.")],
    radius: Annotated[float, typer.Option("--radius", help="Cap and cylinder radius.")] = 0.5,
    height: Annotated[float, typer.Option("--height", help="Cylindrical middle length.")] = 1.0,
    sections: Annotated[
        int, typer.Option("--sections", help="Segments per circumferential / meridional axis.")
    ] = 16,
) -> None:
    """Build a capsule STL (cylinder + hemispherical caps), axis along z."""
    from shore.exceptions import ParameterError
    from shore.shapes import capsule

    try:
        mesh_obj = capsule(radius=radius, height=height, sections=sections)
    except ParameterError as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None
    _write_stl(mesh_obj, output, label="capsule")


@body_app.command("extrude")
def body_extrude(
    section: Annotated[
        Path,
        typer.Option(
            "--section",
            help="Selig .dat section file (use 'shore profile naca' to generate one).",
        ),
    ],
    output: Annotated[Path, typer.Option("-o", "--output", help="Output .stl file.")],
    span: Annotated[float, typer.Option("--span", help="Total wing span.")] = 1.0,
    n_stations: Annotated[
        int, typer.Option("--n-stations", help="Span stations including root and tip.")
    ] = 8,
    taper: Annotated[float, typer.Option("--taper", help="Tip-to-root chord ratio.")] = 1.0,
    sweep_deg: Annotated[
        float, typer.Option("--sweep-deg", help="Sweep angle of the LE in degrees.")
    ] = 0.0,
    twist_deg: Annotated[
        float, typer.Option("--twist-deg", help="Geometric twist at tip in degrees.")
    ] = 0.0,
    span_axis: Annotated[
        str, typer.Option("--span-axis", help="Span axis: 'x', 'y', or 'z'.")
    ] = "z",
    end_caps: Annotated[
        bool,
        typer.Option(
            "--end-caps/--open-ends",
            help="Add centroid-fan caps at the span ends (watertight).",
        ),
    ] = False,
) -> None:
    """Extrude a 2D Selig section into a 3D wing STL.

    Reads ``--section`` as a Selig .dat (see `shore profile naca`),
    builds a wing by stacking copies along ``--span-axis`` with
    optional taper, sweep, and twist, and writes an STL.

    Examples
    --------

        # Constant-section straight wing:
        shore profile naca 0012 -o naca0012.dat
        shore body extrude --section naca0012.dat --span 1.0 --n-stations 8 -o wing.stl

        # Tapered swept twisted wing with watertight ends:
        shore body extrude --section naca0012.dat --span 2.0 --n-stations 12 \\
            --taper 0.6 --sweep-deg 25 --twist-deg -3 --end-caps -o wing.stl
    """
    from shore.exceptions import MeshInputError, ParameterError
    from shore.profiles import read_section
    from shore.shapes import extrude_section

    try:
        coords = read_section(section)
    except (MeshInputError, ParameterError) as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    try:
        mesh_obj = extrude_section(
            coords,
            span=span,
            n_stations=n_stations,
            taper=taper,
            sweep_deg=sweep_deg,
            twist_deg=twist_deg,
            span_axis=span_axis,
            end_caps=end_caps,
        )
    except ParameterError as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    label = f"wing, {n_stations} stations"
    if taper != 1.0:
        label += f", taper={taper}"
    if sweep_deg != 0.0:
        label += f", sweep={sweep_deg}°"
    if twist_deg != 0.0:
        label += f", twist={twist_deg}°"
    if end_caps:
        label += ", capped"
    _write_stl(mesh_obj, output, label=label)
