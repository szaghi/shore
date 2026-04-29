"""Standalone example: build a Cartesian box primitive and export it.

Run from the project root after activating .venv:

    # Generate box.geo + box.vtm with the defaults; open the viewer.
    python examples/02-primitives/visualize_box.py

    # Generate only — skip the interactive viewer (CI-friendly).
    python examples/02-primitives/visualize_box.py --no-view

    # Save a screenshot instead of opening a window.
    python examples/02-primitives/visualize_box.py --screenshot box.png

    # Drive the parameters from the TOML config:
    python examples/02-primitives/visualize_box.py --config examples/02-primitives/config_primitive_box.toml

The box is a single axis-aligned hexahedral block with independent
spacing laws per axis.  Use case: outermost Cartesian background of a
Chimera assembly, or a simple buffer when geometry permits axis-aligned
cells.

For visualisation/Jupyter use, see examples/01-body-fitted/tour.ipynb.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def generate_box(
    output: Path,
    vmin: tuple[float, float, float],
    vmax: tuple[float, float, float],
    ni: int,
    nj: int,
    nk: int,
    i_law: str,
    j_law: str,
    k_law: str,
    i_beta: float,
    j_beta: float,
    k_beta: float,
    write_vtk: bool,
) -> tuple[Path, Path | None]:
    """Build the box primitive and write its .geo (and optionally .vtm).

    Returns
    -------
    (geo_path, vtm_path or None)
    """
    from shore.primitive import BoxMesh

    print(f"Building box primitive ({ni}x{nj}x{nk}) ...")
    print(f"  vmin = {vmin},  vmax = {vmax}")
    print(f"  spacing: i={i_law}, j={j_law}, k={k_law}")

    box = BoxMesh.build(
        vmin=vmin,
        vmax=vmax,
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

    geo_path = output.with_suffix(".geo")
    box.write_geo(geo_path)
    print(f"  wrote {geo_path}")

    vtm_path: Path | None = None
    if write_vtk:
        try:
            vtm_path = output.with_suffix(".vtm")
            box.write_vtk(vtm_path)
            print(f"  wrote {vtm_path}")
        except ImportError as exc:
            print(f"Warning: VTK export skipped ({exc}).", file=sys.stderr)
            vtm_path = None

    return geo_path, vtm_path


def _load_from_config(path: Path) -> dict:
    """Read a primitive TOML config and return the kwargs for generate_box."""
    from shore.config import load_primitive_config

    cfg = load_primitive_config(path)
    if cfg.geometry.type != "box":
        raise SystemExit(f"Error: config geometry.type='{cfg.geometry.type}' (expected 'box').")
    return {
        "vmin": cfg.geometry.vmin or (0.0, 0.0, 0.0),
        "vmax": cfg.geometry.vmax or (1.0, 1.0, 1.0),
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
        description="Generate a Cartesian box primitive and export to .geo / .vtm.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="TOML primitive config (kind = 'primitive', type = 'box').",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("box"),
        help="Output stem; writes <stem>.geo and <stem>.vtm.",
    )
    parser.add_argument(
        "--min",
        dest="vmin",
        type=float,
        nargs=3,
        default=(0.0, 0.0, 0.0),
        metavar=("X", "Y", "Z"),
        help="Minimum corner of the box.",
    )
    parser.add_argument(
        "--max",
        dest="vmax",
        type=float,
        nargs=3,
        default=(10.0, 5.0, 3.0),
        metavar=("X", "Y", "Z"),
        help="Maximum corner of the box.",
    )
    parser.add_argument("--ni", type=int, default=40)
    parser.add_argument("--nj", type=int, default=30)
    parser.add_argument("--nk", type=int, default=20)
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
    parser.add_argument(
        "--no-view", action="store_true", help="Generate only; skip the interactive viewer."
    )
    parser.add_argument(
        "--screenshot",
        type=Path,
        default=None,
        help="Save a PNG screenshot instead of opening a window.",
    )
    parser.add_argument("--no-edges", action="store_true", help="Hide cell edges in viewer.")
    parser.add_argument(
        "--opacity", type=float, default=0.6, help="Grid opacity [0-1] in the viewer."
    )

    args = parser.parse_args()

    # Resolve parameters: CLI flags only, or override CLI defaults from a config.
    kwargs = {
        "vmin": tuple(args.vmin),
        "vmax": tuple(args.vmax),
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
        if cfg_output is not None and args.output == Path("box"):
            output = cfg_output

    geo_path, _vtm_path = generate_box(
        output=output,
        write_vtk=not args.no_vtk,
        **kwargs,
    )

    if args.no_view and args.screenshot is None:
        return

    # --- visualise ---
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
