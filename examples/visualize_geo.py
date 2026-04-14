"""Standalone visualisation script for SHORE outputs.

Run from the project root after activating .venv:

    # Generate a unit-sphere test case and open the viewer in one step:
    python examples/visualize_geo.py --generate

    # Visualise an existing file:
    python examples/visualize_geo.py sphere.geo
    python examples/visualize_geo.py sphere_k0.npz
    python examples/visualize_geo.py sphere.geo --screenshot out.png

For Jupyter use, see examples/visualize.ipynb.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Sphere generation helper
# ---------------------------------------------------------------------------

def generate_sphere(
    stl_path: Path,
    geo_path: Path,
    npz_path: Path,
    ni: int,
    nj: int,
    nk: int,
    ds: float,
    growth: float,
) -> None:
    """Create sphere.stl, sphere_k0.npz, and sphere.geo using the SHORE API."""
    import trimesh

    from shore.io.geo import write_geo
    from shore.surface.io import load_stl
    from shore.surface.remesh import project_sphere_to_surface
    from shore.volume.hyperbolic import march

    # 1 — STL
    if not stl_path.exists():
        print(f"Generating {stl_path} …")
        mesh = trimesh.creation.icosphere(subdivisions=4, radius=1.0)
        mesh.export(str(stl_path))
        print(f"  {len(mesh.faces)} triangles written.")
    else:
        print(f"  {stl_path} already exists, skipping STL generation.")

    # 2 — surface projection (k=0)
    if not npz_path.exists():
        print(f"Projecting {ni}×{nj} surface layer → {npz_path} …")
        mesh = load_stl(stl_path)
        k0 = project_sphere_to_surface(mesh, ni=ni, nj=nj)
        np.savez_compressed(str(npz_path), k0=k0)
        print(f"  Surface layer shape: {k0.shape}.")
    else:
        print(f"  {npz_path} already exists, skipping projection.")
        k0 = np.load(str(npz_path))["k0"]

    # 3 — hyperbolic marching
    if not geo_path.exists():
        print(f"Marching {nk} layers → {geo_path} …")

        def _cb(k: int, jmin: float, step: float) -> None:
            print(f"  layer {k:3d}  jmin={jmin:.3e}  step={step:.3e}", end="\r")

        grid = march(k0, nk=nk, ds=ds, growth=growth, layer_callback=_cb)
        print()  # newline after \r
        write_geo(str(geo_path), grid)
        print(f"  {geo_path} written — shape {grid.shape[:3]}.")
    else:
        print(f"  {geo_path} already exists, skipping marching.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualise a SHORE .geo or .npz file with PyVista.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # -- generation options --------------------------------------------------
    gen = parser.add_argument_group("sphere generation (--generate)")
    gen.add_argument(
        "--generate", action="store_true",
        help="Generate sphere.stl + sphere_k0.npz + sphere.geo before visualising",
    )
    gen.add_argument("--ni",     type=int,   default=40,   help="Latitude points")
    gen.add_argument("--nj",     type=int,   default=60,   help="Longitude points")
    gen.add_argument("--nk",     type=int,   default=20,   help="Wall-normal layers")
    gen.add_argument("--ds",     type=float, default=0.02, help="First-layer thickness")
    gen.add_argument("--growth", type=float, default=1.10, help="Layer growth ratio")

    # -- visualisation options -----------------------------------------------
    parser.add_argument(
        "input", type=Path, nargs="?", default=None,
        help=".geo or .npz file to visualise (default: sphere.geo when --generate is set)",
    )
    parser.add_argument(
        "--colorby", default="",
        help="Colour field: 'Jacobian' (geo default), 'radius' (npz default), or 'none'",
    )
    parser.add_argument(
        "--screenshot", type=Path, default=None,
        help="Save a PNG screenshot instead of opening an interactive window",
    )
    parser.add_argument("--no-edges",  action="store_true", help="Hide cell edges")
    parser.add_argument("--no-view",   action="store_true", help="Generate only, skip viewer")
    parser.add_argument("--opacity",   type=float, default=0.6, help="Grid opacity [0-1]")

    args = parser.parse_args()

    # -- generation ----------------------------------------------------------
    if args.generate:
        stl_path = Path("sphere.stl")
        npz_path = Path("sphere_k0.npz")
        geo_path = Path("sphere.geo")
        generate_sphere(
            stl_path, geo_path, npz_path,
            ni=args.ni, nj=args.nj, nk=args.nk,
            ds=args.ds, growth=args.growth,
        )
        if args.no_view:
            return
        if args.input is None:
            args.input = geo_path

    # -- resolve input -------------------------------------------------------
    if args.input is None:
        parser.error("Provide an input file or use --generate.")

    if not args.input.exists():
        print(f"Error: file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    suffix = args.input.suffix.lower()

    # -- visualise -----------------------------------------------------------
    if suffix == ".geo":
        from shore.io.geo import read_geo
        from shore.vis import save_screenshot, view_geo

        grid = read_geo(args.input)
        colorby = args.colorby or "Jacobian"

        if args.screenshot:
            save_screenshot(grid, str(args.screenshot))
            print(f"Screenshot saved to {args.screenshot}")
        else:
            view_geo(
                grid,
                colorby=colorby,
                show_edges=not args.no_edges,
                opacity=args.opacity,
            )

    elif suffix == ".npz":
        from shore.vis import view_surface

        data = np.load(args.input)
        if "k0" not in data:
            print("Error: .npz has no 'k0' array.", file=sys.stderr)
            sys.exit(1)

        colorby = args.colorby or "radius"
        view_surface(data["k0"], colorby=colorby, show_edges=not args.no_edges)

    else:
        print(
            f"Error: unsupported extension '{suffix}'. Expected .geo or .npz.",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
