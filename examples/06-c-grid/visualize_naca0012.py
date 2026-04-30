"""Standalone example: C-grid mesh around a NACA 0012 airfoil.

Builds a NACA 0012 STL from the analytic 4-digit profile, runs the
SHORE C-grid pipeline (TE detection -> span slice -> C-curve build ->
hyperbolic march), and writes a ``.geo`` plus a ParaView ``.vtm``.

By default everything lands under ``examples_out/`` (cwd-relative):

    examples_out/
        naca0012.stl
        naca0012.geo
        naca0012.vtm
        naca0012/                 # block payloads referenced by the .vtm

Run from the project root after activating .venv:

    # Generate a small C-grid and write naca0012.{stl,geo,vtm}:
    python examples/06-c-grid/visualize_naca0012.py --generate

    # Larger C-grid with finer wake clustering:
    python examples/06-c-grid/visualize_naca0012.py --generate \\
        --ni-body 160 --ni-wake 40 --n-stations 8 \\
        --wake-length 25 --wake-growth 1.20 \\
        --nk 30 --ds 5e-4 --growth 1.15

    # Visualise an existing .geo (any topology):
    python examples/06-c-grid/visualize_naca0012.py naca0012.geo

C-grid layout
-------------
A C-grid wraps a sharp-TE body (airfoil, hydrofoil) with a single
hexahedral block and adds two wake branches downstream of the TE.
The i-axis sweeps as one continuous parametric line:

    upper-wake-far-tip  ->  upper-TE  ->  suction  ->  LE
                                                       |
                                                       v
    lower-wake-far-tip  <-  lower-TE  <-  pressure   <-

The j-axis is the span direction (one node per station); the k-axis
is the wall-normal extrusion.  See :mod:`shore.surface.cgrid` for the
detection / slicing / wake-construction details.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# NACA 0012 STL builder
# ---------------------------------------------------------------------------


def _make_naca0012_stl(
    path: Path,
    *,
    n_section_pts: int = 80,
    n_stations: int = 8,
    chord: float = 1.0,
    span: float = 0.5,
) -> None:
    """Write a NACA 0012 wing STL with sharp TE and centroid-fan end caps.

    Uses :func:`shore.profiles.naca` for the analytic 2D section and
    :func:`shore.shapes.extrude_section` for the 3D extrusion.
    """
    from shore.profiles import naca
    from shore.shapes import extrude_section

    section = naca("0012", n_pts=n_section_pts // 2) * chord
    mesh = extrude_section(
        section,
        span=span,
        n_stations=n_stations,
        end_caps=True,
    )
    mesh.export(str(path))


# ---------------------------------------------------------------------------
# Mesh generation
# ---------------------------------------------------------------------------


def generate_cgrid(
    stl_path: Path,
    geo_path: Path,
    *,
    ni_body: int,
    ni_wake: int,
    n_stations: int,
    wake_length: float,
    wake_growth: float,
    span_axis: str,
    nk: int,
    ds: float,
    growth: float,
    smooth: float,
    smooth_iters: int,
) -> Path:
    """Build the STL (if missing), march the C-grid, write the .geo."""
    from shore.io.geo import write_geo
    from shore.mesh import Body, CGridTopology, Mesh, Spacing1D

    if not stl_path.exists():
        print(f"Generating {stl_path} ...")
        _make_naca0012_stl(stl_path, n_stations=max(8, n_stations + 2))
        print("  STL written.")

    if geo_path.exists():
        print(f"  {geo_path} already exists, skipping marching.")
        return geo_path

    print(
        f"Building C-grid: ni_body={ni_body} ni_wake={ni_wake} "
        f"n_stations={n_stations} wake_length={wake_length}"
    )
    body = Body.from_stl(stl_path, skip_hygiene=True)
    topo = CGridTopology()
    spacing_k = Spacing1D(law="geometric", ds=ds, growth=growth)
    blocks = topo.build(
        surface=None,
        nk=nk,
        mesh=body.mesh,
        spacing_k=spacing_k,
        ni_body=ni_body,
        ni_wake=ni_wake,
        n_stations=n_stations,
        wake_length=wake_length,
        wake_growth=wake_growth,
        span_axis=span_axis,
    )
    m = Mesh(body=body, blocks=blocks, topology=topo)

    def _cb(k: int, jmin: float) -> None:
        print(f"  layer {k:3d}  jmin={jmin:.3e}", end="\r")

    m.march(
        ds=ds,
        growth=growth,
        smooth_strength=smooth,
        smooth_iters=smooth_iters,
        layer_callback=_cb,
    )
    print()

    write_geo(geo_path, m.blocks[0].nodes)
    print(f"  {geo_path} written  shape={m.blocks[0].nodes.shape[:3]}")
    return geo_path


def export_to_vtk(geo_path: Path) -> Path:
    """Write a single-block .vtm next to the .geo."""
    from shore.io.geo import read_geo
    from shore.io.vtk import write_caps_vtm

    vtm_path = geo_path.with_suffix(".vtm")
    grid = read_geo(geo_path)
    write_caps_vtm(vtm_path, [grid], labels=["cgrid"], with_jacobian=True)
    print(f"  {vtm_path} written.")
    return vtm_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build / visualise a SHORE C-grid around NACA 0012.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--generate",
        action="store_true",
        help="Generate naca0012.stl + naca0012.geo",
    )
    parser.add_argument("--ni-body", type=int, default=120)
    parser.add_argument("--ni-wake", type=int, default=30)
    parser.add_argument("--n-stations", type=int, default=5)
    parser.add_argument("--wake-length", type=float, default=15.0)
    parser.add_argument("--wake-growth", type=float, default=1.15)
    parser.add_argument("--span-axis", default="z", choices=["x", "y", "z"])
    parser.add_argument("--nk", type=int, default=20)
    parser.add_argument("--ds", type=float, default=5e-4)
    parser.add_argument("--growth", type=float, default=1.15)
    parser.add_argument("--smooth", type=float, default=0.1)
    parser.add_argument("--smooth-iters", type=int, default=1)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("examples_out"),
        help="Output directory; all files land here.",
    )
    parser.add_argument(
        "input",
        type=Path,
        nargs="?",
        default=None,
        help=".geo file to visualise",
    )
    parser.add_argument("--no-view", action="store_true", help="Generate only, skip the viewer")

    args = parser.parse_args()
    out_dir: Path = args.output
    out_dir.mkdir(parents=True, exist_ok=True)
    stl_path = out_dir / "naca0012.stl"
    geo_path = out_dir / "naca0012.geo"

    if args.generate:
        try:
            generate_cgrid(
                stl_path,
                geo_path,
                ni_body=args.ni_body,
                ni_wake=args.ni_wake,
                n_stations=args.n_stations,
                wake_length=args.wake_length,
                wake_growth=args.wake_growth,
                span_axis=args.span_axis,
                nk=args.nk,
                ds=args.ds,
                growth=args.growth,
                smooth=args.smooth,
                smooth_iters=args.smooth_iters,
            )
        except Exception as exc:  # pragma: no cover - example script
            print(f"Error: {exc}", file=sys.stderr)
            sys.exit(1)

        export_to_vtk(geo_path)

        if args.no_view:
            return
        target = geo_path
    else:
        target = args.input
        if target is None:
            parser.error("Pass --generate or a .geo file to visualise.")

    try:
        from shore.mesh import view_geo_file
    except ImportError:
        print(
            "PyVista not installed; install with `pip install shore[vis]`.",
            file=sys.stderr,
        )
        sys.exit(0)
    view_geo_file(target)


if __name__ == "__main__":
    main()
