"""SHORE command-line interface."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import click
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
    CubedSphereTopology,
    Mesh,
    OGridTopology,
    Spacing1D,
    check_geo,
    export_geo,
    export_surface,
    info_geo,
    info_surface,
    view_geo_file,
    view_surface_file,
)
from shore.primitive import (
    AnnularCylinderMesh,
    BoxMesh,
    FlatCapsCylinder,
)

app = typer.Typer(
    name="shore",
    help="SHORE - Structured Hyperbolic Overlapping Remeshing Engine.",
    add_completion=False,
    no_args_is_help=True,
)

_console = Console()
_err = Console(stderr=True)


# ---------------------------------------------------------------------------
# Shared option aliases
# ---------------------------------------------------------------------------

_NI = Annotated[int, typer.Option(help="Latitude resolution of the O-grid surface.")]
_NJ = Annotated[int, typer.Option(help="Longitude resolution of the O-grid surface.")]
_NK = Annotated[int, typer.Option(help="Number of wall-normal layers (includes k=0).")]
_DS = Annotated[float, typer.Option(help="First-layer thickness (same units as the STL).")]
_GROWTH = Annotated[float, typer.Option(help="Geometric growth ratio between successive layers.")]
_SMOOTH = Annotated[float, typer.Option(help="Laplacian smooth strength per layer [0-1].")]
_SITERS = Annotated[int, typer.Option(help="Laplacian sweeps per layer.")]
_BLEND_K = Annotated[
    int,
    typer.Option(
        "--blend-normals-k",
        help="Layers over which to blend STL surface normals into grid normals (0 = off).",
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


def _param_source(name: str) -> click.core.ParameterSource | None:
    """Return how *name* was supplied for the current click invocation, or None."""
    try:
        ctx = click.get_current_context()
    except RuntimeError:
        return None
    return ctx.get_parameter_source(name)


def _make_picker(cfg):  # type: ignore[no-untyped-def]
    """Build a (name, cli_value, cfg_value) -> resolved-value function.

    Precedence: explicit CLI flag wins; otherwise config (when present and
    not None); otherwise the CLI default.
    """

    def pick(name: str, cli_value, cfg_value):  # type: ignore[no-untyped-def]
        if cfg is None:
            return cli_value
        src = _param_source(name)
        if src == click.core.ParameterSource.COMMANDLINE:
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
    input: Annotated[
        Path | None,
        typer.Argument(help="Input STL file (optional when --config supplies it)."),
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
            help="Topology to build: 'cubed_sphere' (6 blocks, default) or 'ogrid' (1 block).",
        ),
    ] = "cubed_sphere",
    verbose: _VERBOSE = False,
) -> None:
    """Full pipeline: STL -> structured volume -> geo files.

    Two topologies are supported:

    - ``cubed_sphere`` (default): 6-block decomposition (4 equatorial
      sub-blocks + 2 polar caps); writes 6 ``<stem>_<label>.geo`` files.
    - ``ogrid``: single-block lat-lon grid, j-periodic; writes one
      ``<stem>.geo`` file.

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
    _pick = _make_picker(cfg)
    output = _pick("output", output, cfg.output if cfg else None)
    ni = _pick("ni", ni, cfg.grid.ni if cfg else None)
    nj = _pick("nj", nj, cfg.grid.nj if cfg else None)
    nk = _pick("nk", nk, cfg.grid.nk if cfg else None)
    ds = _pick("ds", ds, cfg.march.ds if cfg else None)
    growth = _pick("growth", growth, cfg.march.growth if cfg else None)
    smooth = _pick("smooth", smooth, cfg.march.smooth if cfg else None)
    smooth_iters = _pick("smooth_iters", smooth_iters, cfg.march.smooth_iters if cfg else None)
    blend_normals_k = _pick(
        "blend_normals_k", blend_normals_k, cfg.march.blend_normals_k if cfg else None
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
    if topology not in ("cubed_sphere", "ogrid"):
        _err.print(
            f"[red]Error:[/red] unknown --topology {topology!r}; "
            "expected 'cubed_sphere' or 'ogrid'."
        )
        raise typer.Exit(1) from None

    # theta_cap is special: CLI exposes two flags (--theta-cap and
    # --theta-cap-deg) but the config has a single resolved radian
    # value.  Apply the same picker logic to whichever CLI flag was set.
    if cfg is not None and cfg.grid.theta_cap is not None:
        cli_set_rad = _param_source("theta_cap") == click.core.ParameterSource.COMMANDLINE
        cli_set_deg = _param_source("theta_cap_deg") == click.core.ParameterSource.COMMANDLINE
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
    spacing_k = _build_spacing_k(
        volume_k_spacing, volume_k_beta, volume_k_thickness, nk, ds, growth
    )
    out_stem = output or _default_output(input, "")

    # --- load body ---
    if verbose:
        with _progress() as prog:
            t = prog.add_task("Loading STL   ", total=1, info="")
            body = Body.from_stl(input)
            prog.update(t, advance=1, info=f"{len(body.mesh.faces)} triangles")
    else:
        body = Body.from_stl(input)

    # --- build topology + march ---
    _theta_topo = theta if theta is not None else 0.05 * np.pi
    if topology == "cubed_sphere":
        topo: CubedSphereTopology | OGridTopology = CubedSphereTopology()
    else:
        topo = OGridTopology()
    try:
        blocks = topo.build(
            surface=None,
            nk=nk,
            ni=ni,
            nj=nj,
            spacing_k=spacing_k,
            projector=body.projector,
            theta_cap=_theta_topo,
            i_spacing=surface_i_spacing,
            i_beta=surface_i_beta,
            j_spacing=surface_j_spacing,
            j_beta=surface_j_beta,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    m = Mesh(body=body, blocks=blocks, topology=topo)

    steps = spacing_k.build_steps(nk)

    if verbose:
        with _progress() as prog:
            t = prog.add_task("Marching      ", total=nk - 1, info="")

            def _cb(k: int, jmin: float) -> None:
                prog.update(t, advance=1, info=f"jmin={jmin:.2e}")

            try:
                m.march(
                    ds=ds,
                    growth=growth,
                    smooth_strength=smooth,
                    smooth_iters=smooth_iters,
                    steps=steps,
                    layer_callback=_cb,
                    blend_normals_k=blend_normals_k,
                )
            except Exception as exc:
                _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
                raise typer.Exit(1) from None
    else:
        try:
            m.march(
                ds=ds,
                growth=growth,
                smooth_strength=smooth,
                smooth_iters=smooth_iters,
                steps=steps,
                blend_normals_k=blend_normals_k,
            )
        except Exception as exc:
            _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
            raise typer.Exit(1) from None

    # --- write ---
    if topology == "cubed_sphere":
        m.write_geo(out_stem)
        if verbose:
            nc = nj // 4 + 1
            _console.print(
                f"[green]✓[/green] Wrote 6 .geo files  (sub: {ni}x{nc}x{nk}, caps: {nc}x{nc}x{nk})"
            )
    else:
        # Single-block O-grid: write one ``<stem>.geo`` (no _<label>).
        from shore.io.geo import write_geo

        single_path = out_stem if out_stem.suffix == ".geo" else out_stem.with_suffix(".geo")
        write_geo(single_path, m.blocks[0].nodes)
        if verbose:
            _console.print(f"[green]✓[/green] Wrote {single_path}  ({ni}x{nj}x{nk})")


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
    input: Annotated[Path, typer.Argument(help="Input STL file.")],
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
            t = prog.add_task("Loading STL   ", total=1, info="")
            body = Body.from_stl(input)
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
        body = Body.from_stl(input)
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
                "One or more input files. A single .geo writes a .vts; a "
                "single .npz writes a .vts; multiple .geo files are bundled "
                "into a single .vtm MultiBlock."
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
    if len(inputs) == 1:
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
    primitive_obj: BoxMesh | AnnularCylinderMesh | FlatCapsCylinder,
    output: Path | None,
    default_name: str,
) -> None:
    """Resolve output path and call write_geo on the primitive."""
    out = output if output is not None else Path(f"{default_name}.geo")
    primitive_obj.write_geo(out)
    if isinstance(primitive_obj, FlatCapsCylinder):
        stem = out.with_suffix("") if out.suffix == ".geo" else out
        for blk in primitive_obj.blocks:
            _console.print(f"[green]wrote[/green] {stem}_{blk.label}.geo")
    else:
        target = out if out.suffix == ".geo" else out.with_suffix(".geo")
        _console.print(f"[green]wrote[/green] {target}")


def _load_prim_config(path: Path | None) -> PrimitiveConfig | None:
    if path is None:
        return None
    try:
        return load_primitive_config(path)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None


def _pick_prim(cfg: PrimitiveConfig | None):  # type: ignore[no-untyped-def]
    """Same precedence picker as ``shore mesh``: explicit CLI > config > default."""

    def pick(name: str, cli_value, cfg_value):  # type: ignore[no-untyped-def]
        if cfg is None:
            return cli_value
        src = _param_source(name)
        if src == click.core.ParameterSource.COMMANDLINE:
            return cli_value
        if cfg_value is not None:
            return cfg_value
        return cli_value

    return pick


# ----- shore primitive box ------------------------------------------------


@primitive_app.command("box")
def primitive_box(
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
    ni: _NI = 40,
    nj: _NJ = 40,
    nk: _NK = 30,
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
) -> None:
    """Build a Cartesian box primitive (single block)."""
    cfg = _load_prim_config(config)
    if cfg is not None and cfg.geometry.type != "box":
        _err.print(
            f"[red]Error:[/red] config geometry.type='{cfg.geometry.type}' does not match 'box'."
        )
        raise typer.Exit(1) from None

    pick = _pick_prim(cfg)
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
        prim = BoxMesh.build(
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

    _emit_primitive_writeback(prim, output, default_name="box")


# ----- shore primitive annulus -------------------------------------------


@primitive_app.command("annulus")
def primitive_annulus(
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
    ni: _NI = 20,
    nj: _NJ = 64,
    nk: _NK = 30,
    i_law: Annotated[str, typer.Option("--i-spacing")] = "uniform",
    j_law: Annotated[str, typer.Option("--j-spacing")] = "uniform",
    k_law: Annotated[str, typer.Option("--k-spacing")] = "uniform",
    i_beta: Annotated[float, typer.Option("--i-beta")] = 3.0,
    j_beta: Annotated[float, typer.Option("--j-beta")] = 3.0,
    k_beta: Annotated[float, typer.Option("--k-beta")] = 3.0,
) -> None:
    """Build an annular-cylinder primitive (single block, optionally periodic-j)."""
    cfg = _load_prim_config(config)
    if cfg is not None and cfg.geometry.type != "annulus":
        _err.print(
            f"[red]Error:[/red] config geometry.type='{cfg.geometry.type}' "
            "does not match 'annulus'."
        )
        raise typer.Exit(1) from None

    pick = _pick_prim(cfg)
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
        prim = AnnularCylinderMesh.build(
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

    _emit_primitive_writeback(prim, output, default_name="annulus")


# ----- shore primitive flat-caps -----------------------------------------


@primitive_app.command("flat-caps")
def primitive_flat_caps(
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
    nk: _NK = 30,
    i_law: Annotated[str, typer.Option("--i-spacing")] = "uniform",
    c_law: Annotated[
        str,
        typer.Option("--c-spacing", help="Distribution along the central-square / arc direction."),
    ] = "uniform",
    k_law: Annotated[str, typer.Option("--k-spacing")] = "uniform",
    i_beta: Annotated[float, typer.Option("--i-beta")] = 3.0,
    c_beta: Annotated[float, typer.Option("--c-beta")] = 3.0,
    k_beta: Annotated[float, typer.Option("--k-beta")] = 3.0,
) -> None:
    """Build a 5-block flat-caps cylinder (no central hole - Chimera-free background)."""
    cfg = _load_prim_config(config)
    if cfg is not None and cfg.geometry.type != "flat_caps":
        _err.print(
            f"[red]Error:[/red] config geometry.type='{cfg.geometry.type}' "
            "does not match 'flat_caps'."
        )
        raise typer.Exit(1) from None

    pick = _pick_prim(cfg)
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
        prim = FlatCapsCylinder.build(
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

    _emit_primitive_writeback(prim, output, default_name="flat_caps")


# ---------------------------------------------------------------------------
# shore split
# ---------------------------------------------------------------------------


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
    verbose: _VERBOSE = False,
) -> None:
    """Split structured blocks along ijk vertex indices, preserving adjacency.

    Reads existing ``.geo`` files, applies the splits described in the
    TOML config (one ``[[splits]]`` entry per block), writes new
    ``.geo`` files and a JSON sidecar with the adjacency graph.
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

    # Build HexBlocks from the input .geo files.  Block label = filename
    # stem; if a stem looks like ``<base>_<label>``, the suffix is used.
    from shore.io.adjacency import write_adjacency_json
    from shore.io.geo import read_geo, write_geo
    from shore.mesh.block import HexBlock
    from shore.split import BlockSplit, split_blocks

    # Determine block labels.  When all inputs share a common
    # ``<base>_`` prefix, strip it so the labels match the topology
    # builder's labels (e.g. ``sphere_sub0.geo`` → label ``sub0``).
    # Otherwise fall back to the full stem of each input.
    stems = [p.stem for p in inputs]
    if len(stems) >= 2:
        # Longest common ``<base>_`` prefix among all stems, ending
        # at the last underscore boundary.
        common = stems[0]
        for s in stems[1:]:
            i = 0
            limit = min(len(common), len(s))
            while i < limit and common[i] == s[i]:
                i += 1
            common = common[:i]
        # Trim back to the last underscore (so the prefix is full segments).
        if "_" in common:
            base_stem = common.rsplit("_", 1)[0] + "_"
        else:
            base_stem = ""
    else:
        # Single input: strip up to and including the *last* underscore
        # if present, so e.g. ``background_box.geo`` → label ``box``.
        base_stem = stems[0].rsplit("_", 1)[0] + "_" if "_" in stems[0] else ""

    blocks: list[HexBlock] = []
    for p, stem in zip(inputs, stems, strict=True):
        nodes = read_geo(p)
        label = stem[len(base_stem) :] if base_stem and stem.startswith(base_stem) else stem
        blocks.append(HexBlock(label=label, nodes=nodes))

    if verbose:
        _console.print(f"Read {len(blocks)} input block(s):")
        for blk in blocks:
            _console.print(f"  {blk.label:24s} {blk.nodes.shape}")

    directives = [BlockSplit(label=d.label, axis=d.axis, at=d.at) for d in cfg.splits]
    try:
        new_blocks = split_blocks(blocks, directives)
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {_esc(str(exc))}")
        raise typer.Exit(1) from None

    # Resolve output stem (CLI > config > default).
    if output is None:
        output = (
            cfg.output if cfg.output is not None else (inputs[0].parent / (base_stem or "split"))
        )

    out_dir = output.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write one .geo per output block, plus the adjacency JSON.
    for blk in new_blocks:
        target = out_dir / f"{output.name}_{blk.label}.geo"
        write_geo(target, blk.nodes)
        if verbose:
            _console.print(f"  wrote {target}")

    adj_path = out_dir / f"{output.name}.adjacency.json"
    write_adjacency_json(new_blocks, adj_path)
    _console.print(
        f"[green]✓[/green] Wrote {len(new_blocks)} blocks "
        f"+ {adj_path.name} ({len(directives)} split directive(s))"
    )


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
            help="Name of the .grd file Xnavis should pair with this cc.par "
            "(written into header line 1).",
        ),
    ],
    output_basename: Annotated[
        str,
        typer.Option(
            "--output-basename",
            help="Output basename Xnavis uses for derived files (header line 2).",
        ),
    ] = "cc.",
    free_face_bc: Annotated[
        int,
        typer.Option(
            "--free-bc",
            help="BC code emitted for FREE faces.  Default 40 (chimera "
            "cell-center, body-fitted/buffer meshes).  Use a negative natural "
            "BC for far-field background blocks: -3 inflow, -4 inoutflow, -9 "
            "extrapolation.",
        ),
    ] = 40,
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
) -> None:
    """Convert a SHORE .adjacency.json sidecar to an Xnavis cc.par file.

    Maps SHORE's per-face boundary conditions to Xnavis cc.par codes:
    WALL -> 1 (parete), FREE -> --free-bc (default 40), SHARED / DIRICHLET
    -> connection orientation codes with paired patch indices.  Cubed-sphere
    cap-equator cross-axis seams use the hard-coded orientation table from
    shore.io.cc_par.
    """
    from shore.io.adjacency import read_adjacency_json
    from shore.io.cc_par import BlockMeta, write_cc_par

    if not adjacency.exists():
        _err.print(f"[red]Error:[/red] {adjacency} not found.")
        raise typer.Exit(1)

    payload = read_adjacency_json(adjacency)
    meta = {
        blk["label"]: BlockMeta(comment=blk["label"], free_face_bc=free_face_bc)
        for blk in payload["blocks"]
    }
    out = write_cc_par(
        payload,
        output,
        grd_filename=grd_filename,
        output_basename=output_basename,
        multigrid_levels=multigrid_levels,
        boundary_layer_thickness=boundary_layer_thickness,
        numerical_beach=numerical_beach,
        block_metadata=meta,
    )
    _console.print(f"[green]✓[/green] Wrote {out}")
