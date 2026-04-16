"""Standalone visualisation script for SHORE outputs.

Run from the project root after activating .venv:

    # Generate a unit-sphere test case and open the viewer in one step:
    python examples/visualize_geo.py --generate

    # Generate and export to VTK:
    python examples/visualize_geo.py --generate --vtk

    # Generate Yin-Yang double O-grid and export as MultiBlock:
    python examples/visualize_geo.py --generate --topology double --vtk

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
    npz_path: Path | None,
    ni: int,
    nj: int,
    nk: int,
    ds: float,
    growth: float,
    topology: str = "single",
) -> tuple[Path, ...]:
    """Create sphere.stl, sphere_k0.npz, and sphere.geo using the SHORE API.

    Returns tuple of generated .geo paths.
    """
    import trimesh
    import numpy as np

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

    mesh = load_stl(stl_path)

    if topology == "single":
        # 2a — surface projection (k=0) for single O-grid
        if npz_path and not npz_path.exists():
            print(f"Projecting {ni}×{nj} surface layer → {npz_path} …")
            k0 = project_sphere_to_surface(mesh, ni=ni, nj=nj)
            np.savez_compressed(str(npz_path), k0=k0)
            print(f"  Surface layer shape: {k0.shape}.")
        elif npz_path:
            print(f"  {npz_path} already exists, skipping projection.")
            k0 = np.load(str(npz_path))["k0"]
        else:
            k0 = project_sphere_to_surface(mesh, ni=ni, nj=nj)

        # 3a — hyperbolic marching for single O-grid
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

        return (geo_path,)

    else:  # double topology (Yin-Yang)
        from shore.surface.yinyang import project_yang_to_surface

        theta_margin = np.pi / 6

        # 2b — project both Yin and Yang
        k0_yin = project_sphere_to_surface(mesh, ni=ni, nj=nj, theta_margin=theta_margin)
        k0_yang = project_yang_to_surface(mesh, ni=ni, nj=nj, theta_margin=theta_margin)

        stem = geo_path.stem.removesuffix("_yin").removesuffix("_yang")
        yin_geo_path = geo_path.parent / f"{stem}_yin.geo"
        yang_geo_path = geo_path.parent / f"{stem}_yang.geo"

        # 3b — march both Yin and Yang
        if not yin_geo_path.exists():
            print(f"Marching {nk} layers Yin → {yin_geo_path} …")
            yin_grid = march(k0_yin, nk=nk, ds=ds, growth=growth)
            write_geo(str(yin_geo_path), yin_grid)
            print(f"  {yin_geo_path} written — shape {yin_grid.shape[:3]}.")
        else:
            print(f"  {yin_geo_path} already exists, skipping marching.")

        if not yang_geo_path.exists():
            print(f"Marching {nk} layers Yang → {yang_geo_path} …")
            yang_grid = march(k0_yang, nk=nk, ds=ds, growth=growth)
            write_geo(str(yang_geo_path), yang_grid)
            print(f"  {yang_geo_path} written — shape {yang_grid.shape[:3]}.")
        else:
            print(f"  {yang_geo_path} already exists, skipping marching.")

        return (yin_geo_path, yang_geo_path)


def export_to_vtk(
    geo_paths: tuple[Path, ...],
    topology: str = "single",
) -> list[Path]:
    """Export .geo files to VTK format.

    For single topology: writes a .vts file for each .geo path.
    For double topology (Yin-Yang): writes individual .vts files for each
    block plus a combined .vtm MultiBlock for ParaView.

    Returns list of generated VTK file paths.
    """
    from shore.io.geo import read_geo
    from shore.io.vtk import write_vts, write_double_vtm

    vtk_paths: list[Path] = []

    if topology == "single":
        for geo_path in geo_paths:
            vts_path = geo_path.with_suffix(".vts")
            print(f"Exporting {geo_path} → {vts_path} …")
            grid = read_geo(geo_path)
            write_vts(vts_path, grid, with_jacobian=True)
            print(f"  {vts_path} written.")
            vtk_paths.append(vts_path)

    else:  # double topology (Yin-Yang)
        # Write individual .vts files for each block
        for geo_path in geo_paths:
            vts_path = geo_path.with_suffix(".vts")
            print(f"Exporting {geo_path} → {vts_path} …")
            grid = read_geo(geo_path)
            write_vts(vts_path, grid, with_jacobian=True)
            print(f"  {vts_path} written.")
            vtk_paths.append(vts_path)

        # Write combined MultiBlock (.vtm)
        yin_path = geo_paths[0]
        yang_path = geo_paths[1]
        stem = yin_path.stem.removesuffix("_yin")
        vtm_path = yin_path.parent / f"{stem}.vtm"

        print(f"Exporting combined Yin-Yang MultiBlock → {vtm_path} …")
        yin_grid = read_geo(yin_path)
        yang_grid = read_geo(yang_path)
        write_double_vtm(vtm_path, yin_grid, yang_grid, with_jacobian=True)
        print(f"  {vtm_path} written.")
        vtk_paths.append(vtm_path)

    return vtk_paths


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
    gen.add_argument(
        "--topology", default="single",
        choices=["single", "double"],
        help="'single' (default) for single O-grid, 'double' for Yin-Yang",
    )
    gen.add_argument(
        "--vtk", action="store_true",
        help="Export generated .geo file(s) to VTK format after generation. "
             "For single: writes .vts. For double: writes .vts per block plus .vtm MultiBlock.",
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
        npz_path = Path("sphere_k0.npz") if args.topology == "single" else None
        geo_path = Path("sphere.geo")
        geo_paths = generate_sphere(
            stl_path, geo_path, npz_path,
            ni=args.ni, nj=args.nj, nk=args.nk,
            ds=args.ds, growth=args.growth,
            topology=args.topology,
        )

        # VTK export after generation
        if args.vtk:
            try:
                export_to_vtk(geo_paths, topology=args.topology)
            except Exception as e:
                print(f"Warning: VTK export failed: {e}", file=sys.stderr)

        if args.no_view:
            return
        if args.input is None:
            args.input = geo_paths[0]  # default to first geo file

    # -- VTK export for existing files ---------------------------------------
    if args.vtk and args.input is not None:
        if args.input.suffix == ".geo":
            geo_paths = (args.input,)
            try:
                # Determine topology from filename
                topo = "double" if "_yin" in args.input.stem or "_yang" in args.input.stem else "single"
                if topo == "double":
                    # For Yin or Yang block, find the counterpart
                    stem = args.input.stem.removesuffix("_yin").removesuffix("_yang")
                    if "_yin" in args.input.stem:
                        yin_path = args.input
                        yang_path = args.input.parent / f"{stem}_yang.geo"
                        geo_paths = (yin_path, yang_path) if yang_path.exists() else geo_paths
                    elif "_yang" in args.input.stem:
                        yang_path = args.input
                        yin_path = args.input.parent / f"{stem}_yin.geo"
                        geo_paths = (yin_path, yang_path) if yin_path.exists() else geo_paths
                export_to_vtk(geo_paths, topology=topo)
            except Exception as e:
                print(f"Warning: VTK export failed: {e}", file=sys.stderr)
        elif args.input.suffix == ".npz":
            try:
                from shore.io.vtk import write_surface_vts
                vts_path = args.input.with_suffix(".vts")
                print(f"Exporting {args.input} → {vts_path} …")
                data = np.load(args.input)
                if "k0" in data:
                    write_surface_vts(vts_path, data["k0"])
                    print(f"  {vts_path} written.")
                else:
                    print(f"Warning: .npz has no 'k0' array.", file=sys.stderr)
            except Exception as e:
                print(f"Warning: VTK export failed: {e}", file=sys.stderr)
        elif args.generate:
            pass  # VTK already done during generation
        else:
            print(f"Warning: --vtk requires a .geo or .npz input.", file=sys.stderr)

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