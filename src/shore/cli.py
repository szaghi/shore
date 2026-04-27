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

from shore.config import MeshConfig, dump_template, load_config
from shore.mesh import (
    Body,
    CapsTopology,
    Mesh,
    Spacing1D,
    check_geo,
    export_geo,
    export_surface,
    info_geo,
    info_surface,
    view_geo_file,
    view_surface_file,
)

app = typer.Typer(
    name="shore",
    help="SHORE — Structured Hyperbolic Overlapping Remeshing Engine.",
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
        help="[DEPRECATED, Phase 4] Ignored — cap interior is now built by "
        "Shirley + stereographic projection with C1 seam continuity.",
        hidden=True,
    ),
]
_CAP_SQUARENESS = Annotated[
    float | None,
    typer.Option(
        "--cap-squareness",
        help="[DEPRECATED] Ignored — cap k=0 is now built by gnomonic "
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
    verbose: _VERBOSE = False,
) -> None:
    """Full pipeline: STL → 6-block sphere O-grid volume → geo files.

    Parameter precedence (highest first): explicit CLI flag → config file → built-in default.
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
    # Gnomonic algorithm projects the surface inside topology.build through
    # the cached BVH; no separate body.project_caps call is needed.
    _theta_topo = theta if theta is not None else 0.05 * np.pi
    topo = CapsTopology()
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
            t = prog.add_task("Marching caps ", total=nk - 1, info="")

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
    m.write_geo(out_stem)

    if verbose:
        nc = nj // 4 + 1
        _console.print(
            f"[green]✓[/green] Wrote 6 .geo files  (sub: {ni}x{nc}x{nk}, caps: {nc}x{nc}x{nk})"
        )


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
        _console.print(f"  grid          ni={cfg.grid.ni} nj={cfg.grid.nj} nk={cfg.grid.nk}")
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

    table = Table(title=f"Jacobian check — {input}", show_header=True, header_style="bold cyan")
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
    input: Annotated[Path, typer.Argument(help="Input .geo or .npz file.")],
    output: Annotated[Path | None, typer.Option("-o", "--output", help="Output .vts file.")] = None,
    no_jacobian: Annotated[
        bool,
        typer.Option("--no-jacobian", help="Skip Jacobian cell data (geo only)."),
    ] = False,
) -> None:
    """Export a .geo volume or .npz surface to VTK XML (.vts) for ParaView."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

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
