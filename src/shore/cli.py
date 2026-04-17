"""SHORE command-line interface."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import numpy as np
import typer
from rich.console import Console
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

from shore.io.geo import read_geo, write_geo
from shore.surface.io import load_stl
from shore.surface.remesh import project_sphere_to_surface
from shore.volume.hyperbolic import march
from shore.volume.jacobian import jacobian_per_layer
from shore.volume.spacing import build_k_steps

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

_NI       = Annotated[int,   typer.Option(help="Latitude resolution of the O-grid surface.")]
_NJ       = Annotated[int,   typer.Option(help="Longitude resolution of the O-grid surface.")]
_NK       = Annotated[int,   typer.Option(help="Number of wall-normal layers (includes k=0).")]
_DS       = Annotated[float, typer.Option(help="First-layer thickness (same units as the STL).")]
_GROWTH   = Annotated[float, typer.Option(help="Geometric growth ratio between successive layers.")]
_SMOOTH   = Annotated[float, typer.Option(help="Laplacian smooth strength per layer [0–1].")]
_SITERS   = Annotated[int,   typer.Option(help="Laplacian sweeps per layer.")]
_VERBOSE  = Annotated[bool,  typer.Option("--verbose", "-v", help="Show rich progress bar.")]
_TOPOLOGY = Annotated[str,   typer.Option("--topology", help="'single' (default), 'double' (Yin-Yang), or 'caps' (6-block sphere polar cap).")]
_THETA_CAP = Annotated[float | None, typer.Option("--theta-cap", help="Polar exclusion angle in radians: the equatorial band covers [theta_cap, pi-theta_cap]. Cap boundary parallel radius = sin(theta_cap). Default 0.05*pi (~9 deg). Applies to 'single' and 'caps' topologies.")]
_THETA_CAP_DEG = Annotated[float | None, typer.Option("--theta-cap-deg", help="Convenience alias for --theta-cap in degrees. Mutually exclusive with --theta-cap.")]
_THETA_MARGIN = Annotated[float | None, typer.Option("--theta-margin", help="Deprecated alias for --theta-cap (radians). Use --theta-cap instead.")]
_CAP_PULL       = Annotated[float, typer.Option("--cap-pull", help="[caps] Blend factor [0,1) doming the cap interior toward the pole. 0=off (default). Values 0.3-0.5 improve interior cell quality.")]
_CAP_SQUARENESS = Annotated[float, typer.Option("--cap-squareness", help="[caps] Blend factor [0,1] straightening each cap boundary edge from circular arc toward geodesic chord (inscribed square). 0=circular (default), 1=full chord (~90 deg corner angle).")]

# --- spacing options ---
_SURFACE_I_SPACING = Annotated[str,         typer.Option("--surface-i-spacing", help="Latitude spacing law: 'uniform' (default), 'tanh' (cluster near theta_cap), 'tanh2' (both ends).")]
_SURFACE_I_BETA    = Annotated[float,        typer.Option("--surface-i-beta",    help="Clustering strength for --surface-i-spacing tanh/tanh2 (default 3.0).")]
_SURFACE_J_SPACING = Annotated[str,         typer.Option("--surface-j-spacing", help="Longitude spacing law: 'uniform' (default), 'tanh2' (cluster toward --j-cluster-lon).")]
_SURFACE_J_BETA    = Annotated[float,        typer.Option("--surface-j-beta",    help="Clustering strength for --surface-j-spacing tanh2 (default 3.0).")]
_J_CLUSTER_LON     = Annotated[float | None, typer.Option("--j-cluster-lon",     help="Longitude in radians to cluster toward with --surface-j-spacing tanh2.")]
_VOLUME_K_SPACING   = Annotated[str,         typer.Option("--volume-k-spacing",   help="Wall-normal spacing law: 'geometric' (default), 'tanh' (cluster near wall), 'tanh2' (both ends).")]
_VOLUME_K_BETA      = Annotated[float,        typer.Option("--volume-k-beta",      help="Clustering strength for --volume-k-spacing tanh/tanh2 (default 3.0).")]
_VOLUME_K_THICKNESS = Annotated[float | None, typer.Option("--volume-k-thickness", help="Total wall-normal thickness. Required for --volume-k-spacing tanh/tanh2.")]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _progress() -> Progress:
    """Build the standard SHORE progress bar."""
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


def _make_march_callback(progress: Progress, task_id: int):
    """Return a layer_callback that updates *progress*."""
    def callback(k: int, jmin: float, step: float) -> None:
        progress.update(
            task_id,
            advance=1,
            info=f"jmin={jmin:.2e}  Δs={step:.2e}",
        )
    return callback


def _default_output(input_path: Path, suffix: str) -> Path:
    return input_path.parent / (input_path.stem + suffix)


def _resolve_spacing(
    surface_i_spacing: str,
    surface_i_beta: float,
    surface_j_spacing: str,
    surface_j_beta: float,
    j_cluster_lon: float | None,
    volume_k_spacing: str,
    volume_k_beta: float,
    volume_k_thickness: float | None,
    nk: int,
    ds: float,
    growth: float,
) -> tuple[dict, np.ndarray]:
    """Validate spacing params, build surface kwargs and k-steps array."""
    if volume_k_spacing in {"tanh", "tanh2"} and volume_k_thickness is None:
        _err.print(
            "[red]Error:[/red] --volume-k-thickness is required when "
            f"--volume-k-spacing={volume_k_spacing}."
        )
        raise typer.Exit(1) from None

    surface_kwargs: dict = {}
    if surface_i_spacing != "uniform":
        surface_kwargs["i_spacing"] = surface_i_spacing
        surface_kwargs["i_beta"]    = surface_i_beta
    if surface_j_spacing != "uniform":
        surface_kwargs["j_spacing"]     = surface_j_spacing
        surface_kwargs["j_beta"]        = surface_j_beta
        surface_kwargs["j_cluster_lon"] = j_cluster_lon

    try:
        k_steps = build_k_steps(
            nk, ds, growth,
            spacing=volume_k_spacing,
            beta=volume_k_beta,
            total_thickness=volume_k_thickness,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(1) from None

    return surface_kwargs, k_steps


# ---------------------------------------------------------------------------
# shore mesh
# ---------------------------------------------------------------------------

@app.command()
def mesh(
    input: Annotated[Path, typer.Argument(help="Input STL file.")],
    output:       Annotated[Path | None, typer.Option("-o", "--output", help="Output .geo file.")] = None,
    ni:           _NI      = 40,
    nj:           _NJ      = 60,
    nk:           _NK      = 30,
    ds:           _DS      = 1e-3,
    growth:       _GROWTH  = 1.15,
    smooth:       _SMOOTH  = 0.2,
    smooth_iters: _SITERS  = 2,
    topology:     _TOPOLOGY = "single",
    theta_cap:     _THETA_CAP     = None,
    theta_cap_deg: _THETA_CAP_DEG = None,
    theta_margin:  _THETA_MARGIN  = None,
    cap_pull:        _CAP_PULL        = 0.0,
    cap_squareness:  _CAP_SQUARENESS  = 0.0,
    surface_i_spacing: _SURFACE_I_SPACING = "uniform",
    surface_i_beta:    _SURFACE_I_BETA    = 3.0,
    surface_j_spacing: _SURFACE_J_SPACING = "uniform",
    surface_j_beta:    _SURFACE_J_BETA    = 3.0,
    j_cluster_lon:     _J_CLUSTER_LON     = None,
    volume_k_spacing:   _VOLUME_K_SPACING   = "geometric",
    volume_k_beta:      _VOLUME_K_BETA      = 3.0,
    volume_k_thickness: _VOLUME_K_THICKNESS = None,
    verbose:      _VERBOSE = False,
) -> None:
    """Full pipeline: STL → structured O-grid volume → geo."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    if topology not in {"single", "double", "caps"}:
        _err.print(f"[red]Error:[/red] unknown topology '{topology}'. Choose 'single', 'double', or 'caps'.")
        raise typer.Exit(1) from None

    n_theta = sum(x is not None for x in (theta_cap, theta_cap_deg, theta_margin))
    if n_theta > 1:
        _err.print("[red]Error:[/red] supply at most one of --theta-cap, --theta-cap-deg, --theta-margin.")
        raise typer.Exit(1) from None
    if theta_cap_deg is not None:
        theta_cap = np.deg2rad(theta_cap_deg)
    elif theta_margin is not None:
        theta_cap = theta_margin

    # Yin-Yang default: pi/6 if user gave nothing
    if topology == "double" and theta_cap is None:
        theta_cap = np.pi / 6

    out = output or _default_output(input, ".geo")

    surface_kwargs, k_steps = _resolve_spacing(
        surface_i_spacing, surface_i_beta,
        surface_j_spacing, surface_j_beta, j_cluster_lon,
        volume_k_spacing, volume_k_beta, volume_k_thickness,
        nk, ds, growth,
    )

    # --- stage 1: load ---
    if verbose:
        with _progress() as prog:
            t = prog.add_task("Loading STL   ", total=1, info="")
            mesh_obj = load_stl(input)
            prog.update(t, advance=1, info=f"{len(mesh_obj.faces)} triangles")
    else:
        mesh_obj = load_stl(input)

    if topology == "single":
        _run_single_topology(mesh_obj, out, ni, nj, nk, ds, growth, smooth, smooth_iters, theta_cap, verbose, surface_kwargs, k_steps)
    elif topology == "caps":
        _run_caps_topology(mesh_obj, out, ni, nj, nk, ds, growth, smooth, smooth_iters, theta_cap, cap_pull, cap_squareness, verbose, surface_kwargs, k_steps)
    else:
        _run_double_topology(mesh_obj, out, ni, nj, nk, ds, growth, smooth, smooth_iters, theta_cap, verbose, surface_kwargs, k_steps)


def _run_single_topology(
    mesh_obj,
    out: Path,
    ni: int,
    nj: int,
    nk: int,
    ds: float,
    growth: float,
    smooth: float,
    smooth_iters: int,
    theta_cap: float | None,
    verbose: bool,
    surface_kwargs: dict,
    k_steps: np.ndarray,
) -> None:
    """Run the single O-grid topology pipeline."""
    # --- stage 1: project ---
    if verbose:
        with _progress() as prog:
            t = prog.add_task("Projecting    ", total=1, info="")
            k0 = project_sphere_to_surface(mesh_obj, ni, nj, theta_cap=theta_cap, **surface_kwargs)
            prog.update(t, advance=1, info=f"({ni}×{nj}) surface layer")
    else:
        k0 = project_sphere_to_surface(mesh_obj, ni, nj, theta_cap=theta_cap, **surface_kwargs)

    # --- stage 2: march ---
    if verbose:
        with _progress() as prog:
            t = prog.add_task("Marching      ", total=nk - 1, info="")
            grid = march(
                k0, nk, ds, growth, smooth, smooth_iters,
                steps=k_steps,
                layer_callback=_make_march_callback(prog, t),
            )
    else:
        grid = march(k0, nk, ds, growth, smooth, smooth_iters, steps=k_steps)

    # --- stage 3: write ---
    write_geo(out, grid)

    if verbose:
        _console.print(f"[green]✓[/green] Wrote {out}  ({ni}×{nj}×{nk})")


def _run_caps_topology(
    mesh_obj,
    out: Path,
    ni: int,
    nj: int,
    nk: int,
    ds: float,
    growth: float,
    smooth: float,
    smooth_iters: int,
    theta_cap: float | None,
    cap_pull: float,
    cap_squareness: float,
    verbose: bool,
    surface_kwargs: dict,
    k_steps: np.ndarray,
) -> None:
    """Run the 6-block sphere polar cap topology pipeline."""
    from shore.topology.sphere_cap import split_ogrid_with_polar_caps

    # --- stage 1: project ---
    if verbose:
        with _progress() as prog:
            t = prog.add_task("Projecting    ", total=1, info="")
            k0 = project_sphere_to_surface(mesh_obj, ni, nj, theta_cap=theta_cap, **surface_kwargs)
            prog.update(t, advance=1, info=f"({ni}×{nj}) surface layer")
    else:
        k0 = project_sphere_to_surface(mesh_obj, ni, nj, theta_cap=theta_cap, **surface_kwargs)

    # --- stage 2: march all 6 blocks ---
    nc = nj // 4 + 1
    if verbose:
        with _progress() as prog:
            t = prog.add_task("Marching caps ", total=nk - 1, info="")

            def _cb(k: int, jmin: float) -> None:
                prog.update(t, advance=1, info=f"jmin={jmin:.2e}")

            blocks = split_ogrid_with_polar_caps(
                k0, nk, ds, growth, smooth, smooth_iters,
                cap_pull=cap_pull, cap_squareness=cap_squareness,
                steps=k_steps, layer_callback=_cb,
            )
    else:
        blocks = split_ogrid_with_polar_caps(
            k0, nk, ds, growth, smooth, smooth_iters,
            cap_pull=cap_pull, cap_squareness=cap_squareness,
            steps=k_steps,
        )

    # --- stage 3: write ---
    stem = out.stem
    labels = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]
    for label, block in zip(labels, blocks, strict=True):
        geo_path = out.parent / f"{stem}_{label}.geo"
        write_geo(geo_path, block)

    if verbose:
        _console.print(
            f"[green]✓[/green] Wrote 6 blocks  "
            f"(sub: {ni}×{nc}×{nk}, caps: {nc}×{nc}×{nk})"
        )


def _run_double_topology(
    mesh_obj,
    out: Path,
    ni: int,
    nj: int,
    nk: int,
    ds: float,
    growth: float,
    smooth: float,
    smooth_iters: int,
    theta_cap: float,
    verbose: bool,
    surface_kwargs: dict,
    k_steps: np.ndarray,
) -> None:
    """Run the double O-grid (Yin-Yang) topology pipeline."""
    from shore.surface.yinyang import project_yang_to_surface

    # --- stage 1: project both ---
    if verbose:
        with _progress() as prog:
            t = prog.add_task("Projecting Yin ", total=1, info="")
            k0_yin = project_sphere_to_surface(mesh_obj, ni, nj, theta_cap=theta_cap, **surface_kwargs)
            prog.update(t, advance=1, info=f"({ni}×{nj}) surface layer")

        with _progress() as prog:
            t = prog.add_task("Projecting Yang", total=1, info="")
            k0_yang = project_yang_to_surface(mesh_obj, ni, nj, theta_margin=theta_cap, **surface_kwargs)
            prog.update(t, advance=1, info=f"({ni}×{nj}) surface layer")
    else:
        k0_yin = project_sphere_to_surface(mesh_obj, ni, nj, theta_cap=theta_cap, **surface_kwargs)
        k0_yang = project_yang_to_surface(mesh_obj, ni, nj, theta_margin=theta_cap, **surface_kwargs)

    # --- stage 2: march both ---
    stem = out.stem.removesuffix("_yin").removesuffix("_yang")
    yin_out = out.parent / f"{stem}_yin.geo"
    yang_out = out.parent / f"{stem}_yang.geo"

    if verbose:
        with _progress() as prog:
            t = prog.add_task("Marching Yin  ", total=nk - 1, info="")
            yin_grid = march(
                k0_yin, nk, ds, growth, smooth, smooth_iters,
                steps=k_steps, layer_callback=_make_march_callback(prog, t),
            )
        with _progress() as prog:
            t = prog.add_task("Marching Yang ", total=nk - 1, info="")
            yang_grid = march(
                k0_yang, nk, ds, growth, smooth, smooth_iters,
                steps=k_steps, layer_callback=_make_march_callback(prog, t),
            )
    else:
        yin_grid = march(k0_yin, nk, ds, growth, smooth, smooth_iters, steps=k_steps)
        yang_grid = march(k0_yang, nk, ds, growth, smooth, smooth_iters, steps=k_steps)

    # --- stage 3: write ---
    write_geo(yin_out, yin_grid)
    write_geo(yang_out, yang_grid)

    if verbose:
        _console.print(f"[green]✓[/green] Wrote {yin_out} and {yang_out}  ({ni}×{nj}×{nk})")


# ---------------------------------------------------------------------------
# shore project
# ---------------------------------------------------------------------------

@app.command()
def project(
    input: Annotated[Path, typer.Argument(help="Input STL file.")],
    output:        Annotated[Path | None, typer.Option("-o", "--output", help="Output .npz file.")] = None,
    ni:            _NI            = 40,
    nj:            _NJ            = 60,
    theta_cap:     _THETA_CAP     = None,
    theta_cap_deg: _THETA_CAP_DEG = None,
    theta_margin:  _THETA_MARGIN  = None,
    surface_i_spacing: _SURFACE_I_SPACING = "uniform",
    surface_i_beta:    _SURFACE_I_BETA    = 3.0,
    surface_j_spacing: _SURFACE_J_SPACING = "uniform",
    surface_j_beta:    _SURFACE_J_BETA    = 3.0,
    j_cluster_lon:     _J_CLUSTER_LON     = None,
    verbose:       _VERBOSE       = False,
) -> None:
    """Project a structured sphere grid onto an STL surface (k=0 layer only)."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    n_theta = sum(x is not None for x in (theta_cap, theta_cap_deg, theta_margin))
    if n_theta > 1:
        _err.print("[red]Error:[/red] supply at most one of --theta-cap, --theta-cap-deg, --theta-margin.")
        raise typer.Exit(1) from None
    if theta_cap_deg is not None:
        theta_cap = np.deg2rad(theta_cap_deg)
    elif theta_margin is not None:
        theta_cap = theta_margin

    surface_kwargs: dict = {}
    if surface_i_spacing != "uniform":
        surface_kwargs["i_spacing"] = surface_i_spacing
        surface_kwargs["i_beta"]    = surface_i_beta
    if surface_j_spacing != "uniform":
        surface_kwargs["j_spacing"]     = surface_j_spacing
        surface_kwargs["j_beta"]        = surface_j_beta
        surface_kwargs["j_cluster_lon"] = j_cluster_lon

    out = output or _default_output(input, "_k0.npz")

    if verbose:
        with _progress() as prog:
            t = prog.add_task("Loading STL   ", total=1, info="")
            mesh_obj = load_stl(input)
            prog.update(t, advance=1, info=f"{len(mesh_obj.faces)} triangles")

        with _progress() as prog:
            t = prog.add_task("Projecting    ", total=1, info="")
            k0 = project_sphere_to_surface(mesh_obj, ni, nj, theta_cap=theta_cap, **surface_kwargs)
            prog.update(t, advance=1, info=f"({ni}×{nj}) surface layer")
    else:
        mesh_obj = load_stl(input)
        k0 = project_sphere_to_surface(mesh_obj, ni, nj, theta_cap=theta_cap, **surface_kwargs)

    np.savez_compressed(out, k0=k0)

    if verbose:
        _console.print(f"[green]✓[/green] Wrote {out}  ({ni}×{nj}×3)")


# ---------------------------------------------------------------------------
# shore march
# ---------------------------------------------------------------------------

@app.command(name="march")
def march_cmd(
    input: Annotated[Path, typer.Argument(help="Input .npz file (k=0 layer from `shore project`).")],
    output:       Annotated[Path | None, typer.Option("-o", "--output", help="Output .geo file.")] = None,
    nk:           _NK      = 30,
    ds:           _DS      = 1e-3,
    growth:       _GROWTH  = 1.15,
    smooth:       _SMOOTH  = 0.2,
    smooth_iters: _SITERS  = 2,
    volume_k_spacing:   _VOLUME_K_SPACING   = "geometric",
    volume_k_beta:      _VOLUME_K_BETA      = 3.0,
    volume_k_thickness: _VOLUME_K_THICKNESS = None,
    verbose:      _VERBOSE = False,
) -> None:
    """Extrude a structured surface layer (.npz) into a volume (.geo)."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    if volume_k_spacing in {"tanh", "tanh2"} and volume_k_thickness is None:
        _err.print(
            "[red]Error:[/red] --volume-k-thickness is required when "
            f"--volume-k-spacing={volume_k_spacing}."
        )
        raise typer.Exit(1) from None

    out = output or _default_output(input, ".geo")

    data = np.load(input)
    if "k0" not in data:
        _err.print("[red]Error:[/red] .npz file has no 'k0' array. Was it created by `shore project`?")
        raise typer.Exit(1) from None

    k0 = data["k0"]

    try:
        k_steps = build_k_steps(
            nk, ds, growth,
            spacing=volume_k_spacing,
            beta=volume_k_beta,
            total_thickness=volume_k_thickness,
        )
    except Exception as exc:
        _err.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(1) from None

    if verbose:
        with _progress() as prog:
            t = prog.add_task("Marching      ", total=nk - 1, info="")
            grid = march(
                k0, nk, ds, growth, smooth, smooth_iters,
                steps=k_steps, layer_callback=_make_march_callback(prog, t),
            )
    else:
        grid = march(k0, nk, ds, growth, smooth, smooth_iters, steps=k_steps)

    ni, nj, *_ = grid.shape
    write_geo(out, grid)

    if verbose:
        _console.print(f"[green]✓[/green] Wrote {out}  ({ni}×{nj}×{nk})")



# ---------------------------------------------------------------------------
# shore info
# ---------------------------------------------------------------------------

@app.command()
def info(
    input: Annotated[Path, typer.Argument(help="Input .geo file.")],
) -> None:
    """Print metadata of a structured .geo file."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    grid = read_geo(input)
    ni, nj, nk, _ = grid.shape

    xyz_min = grid.reshape(-1, 3).min(axis=0)
    xyz_max = grid.reshape(-1, 3).max(axis=0)

    table = Table(title=str(input), show_header=True, header_style="bold cyan")
    table.add_column("Property", style="dim")
    table.add_column("Value")

    table.add_row("Dimensions", f"ni={ni}  nj={nj}  nk={nk}")
    table.add_row("Points", f"{ni * nj * nk:,}")
    table.add_row("Cells", f"{(ni-1) * (nj-1) * (nk-1):,}")
    table.add_row("x bounds", f"[{xyz_min[0]:.6g},  {xyz_max[0]:.6g}]")
    table.add_row("y bounds", f"[{xyz_min[1]:.6g},  {xyz_max[1]:.6g}]")
    table.add_row("z bounds", f"[{xyz_min[2]:.6g},  {xyz_max[2]:.6g}]")

    _console.print(table)


# ---------------------------------------------------------------------------
# shore check
# ---------------------------------------------------------------------------

@app.command()
def check(
    input: Annotated[Path, typer.Argument(help="Input .geo file.")],
    warn_threshold: Annotated[float, typer.Option(help="Warn if jmin drops below this value.")] = 1e-6,
) -> None:
    """Validate cell Jacobians layer-by-layer in a .geo file.

    Exits with code 1 if any layer has a negative Jacobian.
    """
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    grid = read_geo(input)
    ni, nj, nk, _ = grid.shape

    table = Table(title=f"Jacobian check — {input}", show_header=True, header_style="bold cyan")
    table.add_column("Layer", justify="right", style="dim")
    table.add_column("jmin", justify="right")
    table.add_column("Status")

    any_bad = False
    jmins = jacobian_per_layer(grid)   # (nk-1,)

    for k, jmin in enumerate(jmins):
        if jmin <= 0.0:
            status = "[bold red]INVERTED[/bold red]"
            any_bad = True
        elif jmin < warn_threshold:
            status = "[yellow]warn[/yellow]"
        else:
            status = "[green]ok[/green]"
        table.add_row(str(k + 1), f"{jmin:.4e}", status)

    _console.print(table)

    if any_bad:
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
    no_jacobian: Annotated[bool, typer.Option("--no-jacobian", help="Skip Jacobian cell data (geo only).")] = False,
) -> None:
    """Export a .geo volume or .npz surface layer to VTK XML (.vts) for ParaView."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    try:
        from shore.io.vtk import write_surface_vts, write_vts  # noqa: PLC0415
    except ImportError as exc:
        _err.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(1) from None

    suffix = input.suffix.lower()

    if suffix == ".geo":
        out = output or _default_output(input, ".vts")
        grid = read_geo(input)
        write_vts(out, grid, with_jacobian=not no_jacobian)
        _console.print(f"[green]✓[/green] Wrote {out}")

    elif suffix == ".npz":
        out = output or _default_output(input, ".vts")
        data = np.load(input)
        if "k0" not in data:
            _err.print("[red]Error:[/red] .npz has no 'k0' array.")
            raise typer.Exit(1) from None
        write_surface_vts(out, data["k0"])
        _console.print(f"[green]✓[/green] Wrote {out}")

    else:
        _err.print(f"[red]Error:[/red] unsupported extension '{suffix}'. Expected .geo or .npz.")
        raise typer.Exit(1) from None


# ---------------------------------------------------------------------------
# shore view
# ---------------------------------------------------------------------------

@app.command()
def view(
    input: Annotated[Path, typer.Argument(help="Input .geo or .npz file.")],
    colorby: Annotated[str, typer.Option(help="Colour field: 'Jacobian' (geo) or 'radius' (npz) or 'none'.")] = "",
    no_edges: Annotated[bool, typer.Option("--no-edges", help="Hide cell edges.")] = False,
    opacity: Annotated[float, typer.Option(help="Grid opacity [0–1] (geo only).")] = 0.6,
) -> None:
    """Open an interactive PyVista window to inspect a .geo or .npz file."""
    if not input.exists():
        _err.print(f"[red]Error:[/red] file not found: {input}")
        raise typer.Exit(1) from None

    try:
        from shore.vis import view_geo, view_surface  # noqa: PLC0415
    except ImportError as exc:
        _err.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(1) from None

    suffix = input.suffix.lower()

    if suffix == ".geo":
        grid = read_geo(input)
        cb = colorby if colorby else "Jacobian"
        view_geo(grid, colorby=cb, show_edges=not no_edges, opacity=opacity)

    elif suffix == ".npz":
        data = np.load(input)
        if "k0" not in data:
            _err.print("[red]Error:[/red] .npz has no 'k0' array.")
            raise typer.Exit(1) from None
        cb = colorby if colorby else "radius"
        view_surface(data["k0"], colorby=cb, show_edges=not no_edges)

    else:
        _err.print(f"[red]Error:[/red] unsupported extension '{suffix}'. Expected .geo or .npz.")
        raise typer.Exit(1) from None
