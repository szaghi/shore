"""Standalone example: build a flat-caps cylinder primitive and export it.

Run from the project root after activating .venv:

    # Generate flat_caps_*.geo (5 files) + flat_caps.vtm; open the viewer.
    python examples/02-primitives/visualize_flat_caps.py

    # Generate only — skip the viewer.
    python examples/02-primitives/visualize_flat_caps.py --no-view

    # Save a screenshot instead of opening a window.
    python examples/02-primitives/visualize_flat_caps.py --screenshot flat_caps.png

    # Drive the parameters from the TOML config:
    python examples/02-primitives/visualize_flat_caps.py --config examples/02-primitives/config_primitive_flat_caps.toml

The flat-caps primitive is a 5-block conforming cylindrical mesh: a
central square block plus 4 trapezoidal sub-blocks.  It eliminates
the Chimera fringe a hollow annular cylinder would otherwise need at
its inner boundary — useful for high-order schemes, conservation, or
GPU-friendly memory access.

Use case: cylindrical background mesh when you want to avoid Chimera
fringes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def generate_flat_caps(
    output_stem: Path,
    r_out: float,
    z0: float,
    z1: float,
    ni: int,
    nc: int,
    nk: int,
    inner_frac: float,
    axis: str,
    center: tuple[float, float, float],
    i_law: str,
    c_law: str,
    k_law: str,
    i_beta: float,
    c_beta: float,
    k_beta: float,
    write_vtk: bool,
) -> tuple[list[Path], Path | None]:
    """Build the flat-caps cylinder primitive and write its .geo files.

    Returns
    -------
    (geo_paths, vtm_path or None)
    """
    from shore.primitive import FlatCapsCylinder

    print(f"Building flat-caps cylinder primitive (ni={ni}, nc={nc}, nk={nk}) ...")
    print(f"  r_out = {r_out},  z = [{z0}, {z1}],  inner_frac = {inner_frac}")
    print(f"  axis = {axis}")
    print(f"  spacing: i={i_law}, c={c_law}, k={k_law}")

    fc = FlatCapsCylinder.build(
        r_out=r_out,
        z0=z0,
        z1=z1,
        ni=ni,
        nc=nc,
        nk=nk,
        inner_frac=inner_frac,
        axis=axis,
        center=center,
        i_law=i_law,
        c_law=c_law,
        k_law=k_law,
        i_beta=i_beta,
        c_beta=c_beta,
        k_beta=k_beta,
    )

    # Multi-block primitive: write_geo treats *output_stem* as a stem.
    fc.write_geo(output_stem)
    geo_paths = [output_stem.parent / f"{output_stem.name}_{blk.label}.geo" for blk in fc.blocks]
    for p in geo_paths:
        print(f"  wrote {p}")

    vtm_path: Path | None = None
    if write_vtk:
        try:
            vtm_path = output_stem.with_suffix(".vtm")
            fc.write_vtk(vtm_path)
            print(f"  wrote {vtm_path}")
        except ImportError as exc:
            print(f"Warning: VTK export skipped ({exc}).", file=sys.stderr)
            vtm_path = None

    return geo_paths, vtm_path


def _load_from_config(path: Path) -> dict:
    from shore.config import load_primitive_config

    cfg = load_primitive_config(path)
    if cfg.geometry.type != "flat_caps":
        raise SystemExit(
            f"Error: config geometry.type='{cfg.geometry.type}' (expected 'flat_caps')."
        )
    g = cfg.geometry
    # Per the CLI mapping in src/shore/cli.py, the central-square
    # direction lives in [grid].nj and [spacing].j_law.
    return {
        "r_out": g.r_out if g.r_out is not None else 2.0,
        "z0": g.z0 if g.z0 is not None else 0.0,
        "z1": g.z1 if g.z1 is not None else 4.0,
        "inner_frac": g.inner_frac if g.inner_frac is not None else 0.5,
        "axis": g.axis,
        "center": tuple(g.center),
        "ni": cfg.grid.ni,
        "nc": cfg.grid.nj,
        "nk": cfg.grid.nk,
        "i_law": cfg.spacing.i_law,
        "c_law": cfg.spacing.j_law,
        "k_law": cfg.spacing.k_law,
        "i_beta": cfg.spacing.i_beta,
        "c_beta": cfg.spacing.j_beta,
        "k_beta": cfg.spacing.k_beta,
        "_output_from_cfg": cfg.output,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a 5-block flat-caps cylinder primitive and export to .geo / .vtm.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="TOML primitive config (kind='primitive', type='flat_caps').",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("flat_caps"),
        help="Output stem; writes <stem>_<label>.geo (5 files) and <stem>.vtm.",
    )
    parser.add_argument("--r-out", type=float, default=2.0)
    parser.add_argument("--z0", type=float, default=0.0)
    parser.add_argument("--z1", type=float, default=4.0)
    parser.add_argument(
        "--inner-frac",
        type=float,
        default=0.5,
        help="Central square corner radius / r_out (in (0, 1)).",
    )
    parser.add_argument(
        "--axis", default="z", help="'x' / 'y' / 'z' or '<vx> <vy> <vz>' (3-vector)."
    )
    parser.add_argument(
        "--center", type=float, nargs=3, default=(0.0, 0.0, 0.0), metavar=("X", "Y", "Z")
    )
    parser.add_argument("--ni", type=int, default=16, help="Radial node count for each sub-block.")
    parser.add_argument(
        "--nc", type=int, default=24, help="Side node count for the central square (shared)."
    )
    parser.add_argument("--nk", type=int, default=30)
    parser.add_argument("--i-spacing", default="uniform", choices=["uniform", "tanh", "tanh2"])
    parser.add_argument(
        "--c-spacing",
        default="uniform",
        choices=["uniform", "tanh", "tanh2"],
        help="Distribution along the central-square / arc direction.",
    )
    parser.add_argument("--k-spacing", default="uniform", choices=["uniform", "tanh", "tanh2"])
    parser.add_argument("--i-beta", type=float, default=3.0)
    parser.add_argument("--c-beta", type=float, default=3.0)
    parser.add_argument("--k-beta", type=float, default=3.0)
    parser.add_argument(
        "--no-vtk",
        action="store_true",
        help="Skip the .vtm MultiBlock export (only the 5 .geo files are written).",
    )
    parser.add_argument("--no-view", action="store_true")
    parser.add_argument("--screenshot", type=Path, default=None)
    parser.add_argument("--no-edges", action="store_true")
    parser.add_argument("--opacity", type=float, default=0.6)

    args = parser.parse_args()

    kwargs = {
        "r_out": args.r_out,
        "z0": args.z0,
        "z1": args.z1,
        "inner_frac": args.inner_frac,
        "axis": args.axis,
        "center": tuple(args.center),
        "ni": args.ni,
        "nc": args.nc,
        "nk": args.nk,
        "i_law": args.i_spacing,
        "c_law": args.c_spacing,
        "k_law": args.k_spacing,
        "i_beta": args.i_beta,
        "c_beta": args.c_beta,
        "k_beta": args.k_beta,
    }
    output = args.output
    if args.config is not None:
        cfg = _load_from_config(args.config)
        cfg_output = cfg.pop("_output_from_cfg")
        kwargs.update(cfg)
        if cfg_output is not None and args.output == Path("flat_caps"):
            output = cfg_output

    geo_paths, _vtm_path = generate_flat_caps(
        output_stem=output,
        write_vtk=not args.no_vtk,
        **kwargs,
    )

    if args.no_view and args.screenshot is None:
        return

    # Visualise the *center* block by default (representative of the assembly).
    center_path = next((p for p in geo_paths if p.stem.endswith("_center")), geo_paths[0])
    try:
        from shore.io.geo import read_geo
        from shore.vis import save_screenshot, view_geo
    except ImportError as exc:
        print(f"Visualisation unavailable: {exc}.  Install with `pip install shore[vis]`.")
        return

    grid = read_geo(center_path)
    if args.screenshot is not None:
        save_screenshot(grid, str(args.screenshot))
        print(f"  screenshot saved to {args.screenshot}")
    else:
        print(f"  showing {center_path.name}; the .vtm bundles all 5 blocks for ParaView.")
        view_geo(grid, colorby="Jacobian", show_edges=not args.no_edges, opacity=args.opacity)


if __name__ == "__main__":
    main()
