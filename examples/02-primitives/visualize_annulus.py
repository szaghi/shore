"""Standalone example: build an annular cylinder primitive and export it.

Run from the project root after activating .venv:

    # Generate annulus.geo + annulus.vtm with the defaults; open the viewer.
    python examples/02-primitives/visualize_annulus.py

    # Partial sector (half-pipe).
    python examples/02-primitives/visualize_annulus.py --theta-end 3.14159

    # Generate only — skip the viewer.
    python examples/02-primitives/visualize_annulus.py --no-view

    # Drive the parameters from the TOML config:
    python examples/02-primitives/visualize_annulus.py --config examples/02-primitives/config_primitive_annulus.toml

The annulus is a single block with optional theta-periodicity.  A full
annulus is periodic-j (no duplicate wrap node); a partial sector
(``--theta-end < 2*pi``) is non-periodic.

Use case: buffer between a near-body sphere/cylinder mesh and an
outer Cartesian background.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def generate_annulus(
    output: Path,
    r_in: float,
    r_out: float,
    z0: float,
    z1: float,
    ni: int,
    nj: int,
    nk: int,
    axis: str,
    center: tuple[float, float, float],
    theta_start: float,
    theta_end: float | None,
    i_law: str,
    j_law: str,
    k_law: str,
    i_beta: float,
    j_beta: float,
    k_beta: float,
    write_vtk: bool,
) -> tuple[Path, Path | None]:
    """Build the annular-cylinder primitive and write its .geo (and .vtm).

    Returns
    -------
    (geo_path, vtm_path or None)
    """
    from shore.primitive import AnnularCylinderMesh

    print(f"Building annular cylinder primitive ({ni}x{nj}x{nk}) ...")
    print(f"  r_in = {r_in},  r_out = {r_out},  z = [{z0}, {z1}]")
    if theta_end is None:
        print("  full annulus (j periodic)")
    else:
        print(f"  partial sector: theta in [{theta_start}, {theta_end}] (j non-periodic)")
    print(f"  axis = {axis}")
    print(f"  spacing: i={i_law}, j={j_law}, k={k_law}")

    ann = AnnularCylinderMesh.build(
        r_in=r_in,
        r_out=r_out,
        z0=z0,
        z1=z1,
        ni=ni,
        nj=nj,
        nk=nk,
        axis=axis,
        center=center,
        theta_start=theta_start,
        theta_end=theta_end,
        i_law=i_law,
        j_law=j_law,
        k_law=k_law,
        i_beta=i_beta,
        j_beta=j_beta,
        k_beta=k_beta,
    )

    geo_path = output.with_suffix(".geo")
    ann.write_geo(geo_path)
    print(f"  wrote {geo_path}")

    vtm_path: Path | None = None
    if write_vtk:
        try:
            vtm_path = output.with_suffix(".vtm")
            ann.write_vtk(vtm_path)
            print(f"  wrote {vtm_path}")
        except ImportError as exc:
            print(f"Warning: VTK export skipped ({exc}).", file=sys.stderr)
            vtm_path = None

    return geo_path, vtm_path


def _load_from_config(path: Path) -> dict:
    from shore.config import load_primitive_config

    cfg = load_primitive_config(path)
    if cfg.geometry.type != "annulus":
        raise SystemExit(f"Error: config geometry.type='{cfg.geometry.type}' (expected 'annulus').")
    g = cfg.geometry
    return {
        "r_in": g.r_in if g.r_in is not None else 0.5,
        "r_out": g.r_out if g.r_out is not None else 2.0,
        "z0": g.z0 if g.z0 is not None else 0.0,
        "z1": g.z1 if g.z1 is not None else 4.0,
        "axis": g.axis,
        "center": tuple(g.center),
        "theta_start": g.theta_start,
        "theta_end": g.theta_end,
        "ni": cfg.grid.ni,
        "nj": cfg.grid.nj,
        "nk": cfg.grid.nk,
        "i_law": cfg.spacing.i_law,
        "j_law": cfg.spacing.j_law,
        "k_law": cfg.spacing.k_law,
        "i_beta": cfg.spacing.i_beta,
        "j_beta": cfg.spacing.j_beta,
        "k_beta": cfg.spacing.k_beta,
        "_output_from_cfg": cfg.output,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate an annular-cylinder primitive and export to .geo / .vtm.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="TOML primitive config (kind='primitive', type='annulus').",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("annulus"),
        help="Output stem; writes <stem>.geo and <stem>.vtm.",
    )
    parser.add_argument("--r-in", type=float, default=0.5)
    parser.add_argument("--r-out", type=float, default=2.0)
    parser.add_argument("--z0", type=float, default=0.0)
    parser.add_argument("--z1", type=float, default=4.0)
    parser.add_argument(
        "--axis", default="z", help="'x' / 'y' / 'z' or '<vx> <vy> <vz>' (3-vector)."
    )
    parser.add_argument(
        "--center", type=float, nargs=3, default=(0.0, 0.0, 0.0), metavar=("X", "Y", "Z")
    )
    parser.add_argument("--theta-start", type=float, default=0.0, help="Sector start in radians.")
    parser.add_argument(
        "--theta-end",
        type=float,
        default=None,
        help="Sector end in radians; omit for a full annulus.",
    )
    parser.add_argument("--ni", type=int, default=20)
    parser.add_argument("--nj", type=int, default=64)
    parser.add_argument("--nk", type=int, default=30)
    parser.add_argument("--i-spacing", default="uniform", choices=["uniform", "tanh", "tanh2"])
    parser.add_argument("--j-spacing", default="uniform", choices=["uniform", "tanh", "tanh2"])
    parser.add_argument("--k-spacing", default="uniform", choices=["uniform", "tanh", "tanh2"])
    parser.add_argument("--i-beta", type=float, default=3.0)
    parser.add_argument("--j-beta", type=float, default=3.0)
    parser.add_argument("--k-beta", type=float, default=3.0)
    parser.add_argument(
        "--no-vtk",
        action="store_true",
        help="Skip the .vtm MultiBlock export (only .geo is written).",
    )
    parser.add_argument("--no-view", action="store_true")
    parser.add_argument("--screenshot", type=Path, default=None)
    parser.add_argument("--no-edges", action="store_true")
    parser.add_argument("--opacity", type=float, default=0.6)

    args = parser.parse_args()

    kwargs = {
        "r_in": args.r_in,
        "r_out": args.r_out,
        "z0": args.z0,
        "z1": args.z1,
        "axis": args.axis,
        "center": tuple(args.center),
        "theta_start": args.theta_start,
        "theta_end": args.theta_end,
        "ni": args.ni,
        "nj": args.nj,
        "nk": args.nk,
        "i_law": args.i_spacing,
        "j_law": args.j_spacing,
        "k_law": args.k_spacing,
        "i_beta": args.i_beta,
        "j_beta": args.j_beta,
        "k_beta": args.k_beta,
    }
    output = args.output
    if args.config is not None:
        cfg = _load_from_config(args.config)
        cfg_output = cfg.pop("_output_from_cfg")
        kwargs.update(cfg)
        if cfg_output is not None and args.output == Path("annulus"):
            output = cfg_output

    geo_path, _vtm_path = generate_annulus(
        output=output,
        write_vtk=not args.no_vtk,
        **kwargs,
    )

    if args.no_view and args.screenshot is None:
        return

    try:
        from shore.io.geo import read_geo
        from shore.vis import save_screenshot, view_geo
    except ImportError as exc:
        print(f"Visualisation unavailable: {exc}.  Install with `pip install shore[vis]`.")
        return

    grid = read_geo(geo_path)
    if args.screenshot is not None:
        save_screenshot(grid, str(args.screenshot))
        print(f"  screenshot saved to {args.screenshot}")
    else:
        view_geo(grid, colorby="Jacobian", show_edges=not args.no_edges, opacity=args.opacity)


if __name__ == "__main__":
    main()
