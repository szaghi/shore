"""Standalone visualisation script for SHORE outputs.

Run from the project root after activating .venv:

    # Generate a unit-sphere test case and open the viewer in one step:
    python examples/visualize_geo.py --generate

    # Generate and export to VTK:
    python examples/visualize_geo.py --generate --vtk

    # Generate Yin-Yang double O-grid and export as MultiBlock:
    python examples/visualize_geo.py --generate --topology double --vtk

    # Tanh-spaced boundary layer with fixed total thickness:
    python examples/visualize_geo.py --generate \\
        --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 2.0

    # Latitude clustering (toward cap boundary) + longitude clustering:
    python examples/visualize_geo.py --generate \\
        --surface-i-spacing tanh --surface-i-beta 4.0 \\
        --surface-j-spacing tanh2 --surface-j-beta 5.0 --j-cluster-lon 0.0


    # Caps topology with interface spacing correction:
    python examples/visualize_geo.py --generate --topology caps \
        --correct-spacing --cap-pull 0.4 --cap-squareness 0.5
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
    theta_cap: float | None = None,
    cap_pull: float = 0.0,
    cap_squareness: float = 0.0,
    correct_spacing: bool = False,
    smooth: float = 0.2,
    smooth_iters: int = 2,
    surface_i_spacing: str = "uniform",
    surface_i_beta: float = 3.0,
    surface_j_spacing: str = "uniform",
    surface_j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
    volume_k_spacing: str = "geometric",
    volume_k_beta: float = 3.0,
    volume_k_thickness: float | None = None,
) -> tuple[Path, ...]:
    """Create sphere.stl, sphere_k0.npz, and sphere.geo using the SHORE API.

    Returns tuple of generated .geo paths.
    """
    import trimesh

    from shore.io.geo import write_geo
    from shore.surface.io import load_stl
    from shore.surface.remesh import project_sphere_to_surface
    from shore.volume.hyperbolic import march
    from shore.volume.spacing import build_k_steps

    # --- build surface kwargs (only pass non-defaults) ---
    surface_kwargs: dict = {}
    if surface_i_spacing != "uniform":
        surface_kwargs["i_spacing"] = surface_i_spacing
        surface_kwargs["i_beta"] = surface_i_beta
    if surface_j_spacing != "uniform":
        surface_kwargs["j_spacing"] = surface_j_spacing
        surface_kwargs["j_beta"] = surface_j_beta
        surface_kwargs["j_cluster_lon"] = j_cluster_lon

    # --- build k-steps ---
    if volume_k_spacing in {"tanh", "tanh2"} and volume_k_thickness is None:
        print("Error: --volume-k-thickness is required for tanh/tanh2 volume spacing.",
              file=sys.stderr)
        sys.exit(1)
    k_steps = build_k_steps(
        nk, ds, growth,
        spacing=volume_k_spacing,
        beta=volume_k_beta,
        total_thickness=volume_k_thickness,
    )

    # 1 — STL
    if not stl_path.exists():
        print(f"Generating {stl_path} ...")
        mesh = trimesh.creation.icosphere(subdivisions=4, radius=1.0)
        mesh.export(str(stl_path))
        print(f"  {len(mesh.faces)} triangles written.")
    else:
        print(f"  {stl_path} already exists, skipping STL generation.")

    mesh = load_stl(stl_path)

    if topology == "single":
        # 2a — surface projection (k=0) for single O-grid
        if npz_path and not npz_path.exists():
            print(f"Projecting {ni}x{nj} surface layer -> {npz_path} ...")
            k0 = project_sphere_to_surface(mesh, ni=ni, nj=nj, theta_cap=theta_cap,
                                           **surface_kwargs)
            np.savez_compressed(str(npz_path), k0=k0)
            print(f"  Surface layer shape: {k0.shape}.")
        elif npz_path:
            print(f"  {npz_path} already exists, skipping projection.")
            k0 = np.load(str(npz_path))["k0"]
        else:
            k0 = project_sphere_to_surface(mesh, ni=ni, nj=nj, theta_cap=theta_cap,
                                           **surface_kwargs)

        # 3a — hyperbolic marching for single O-grid
        if not geo_path.exists():
            print(f"Marching {nk} layers -> {geo_path} ...")

            def _cb(k: int, jmin: float, step: float) -> None:
                print(f"  layer {k:3d}  jmin={jmin:.3e}  step={step:.3e}", end="\r")

            grid = march(k0, nk=nk, ds=ds, growth=growth,
                         smooth_strength=smooth, smooth_iters=smooth_iters,
                         steps=k_steps, layer_callback=_cb)
            print()  # newline after \r
            write_geo(str(geo_path), grid)
            print(f"  {geo_path} written -- shape {grid.shape[:3]}.")
        else:
            print(f"  {geo_path} already exists, skipping marching.")

        return (geo_path,)

    elif topology == "caps":
        from shore.topology.sphere_cap import split_ogrid_with_polar_caps

        # 2c — surface projection (k=0) for caps topology
        k0 = project_sphere_to_surface(mesh, ni=ni, nj=nj, theta_cap=theta_cap,
                                       **surface_kwargs)

        # 3c — split O-grid with polar caps
        stem = geo_path.stem
        labels = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]
        cap_geo_paths = [geo_path.parent / f"{stem}_{lbl}.geo" for lbl in labels]

        if not cap_geo_paths[0].exists():
            print(f"Generating 6-block sphere-cap topology ({nk} layers) ...")

            def _cb(k: int, jmin: float) -> None:
                print(f"  layer {k:3d}  jmin={jmin:.3e}", end="\r")

            blocks = split_ogrid_with_polar_caps(
                k0, nk=nk, ds=ds, growth=growth,
                smooth_strength=smooth, smooth_iters=smooth_iters,
                cap_pull=cap_pull, cap_squareness=cap_squareness,
                steps=k_steps, layer_callback=_cb,
            )
            print()
            if correct_spacing:
                from shore.topology.sphere_cap import correct_interface_spacing
                print("Correcting interface spacing ...")
                correct_interface_spacing(blocks)
                print("  Done.")
            for lbl, p, blk in zip(labels, cap_geo_paths, blocks, strict=True):
                write_geo(str(p), blk)
                print(f"  {p} written -- shape {blk.shape[:3]}.")
        else:
            print(f"  {cap_geo_paths[0]} already exists, skipping generation.")

        return tuple(cap_geo_paths)

    else:  # double topology (Yin-Yang)
        from shore.surface.yinyang import project_yang_to_surface

        theta_margin = np.pi / 6

        # 2b — project both Yin and Yang
        k0_yin = project_sphere_to_surface(mesh, ni=ni, nj=nj, theta_cap=theta_cap,
                                           **surface_kwargs)
        k0_yang = project_yang_to_surface(mesh, ni=ni, nj=nj, theta_margin=theta_cap,
                                          **surface_kwargs)

        stem = geo_path.stem.removesuffix("_yin").removesuffix("_yang")
        yin_geo_path = geo_path.parent / f"{stem}_yin.geo"
        yang_geo_path = geo_path.parent / f"{stem}_yang.geo"

        # 3b — march both Yin and Yang
        if not yin_geo_path.exists():
            print(f"Marching {nk} layers Yin -> {yin_geo_path} ...")
            yin_grid = march(k0_yin, nk=nk, ds=ds, growth=growth,
                             smooth_strength=smooth, smooth_iters=smooth_iters,
                             steps=k_steps)
            write_geo(str(yin_geo_path), yin_grid)
            print(f"  {yin_geo_path} written -- shape {yin_grid.shape[:3]}.")
        else:
            print(f"  {yin_geo_path} already exists, skipping marching.")

        if not yang_geo_path.exists():
            print(f"Marching {nk} layers Yang -> {yang_geo_path} ...")
            yang_grid = march(k0_yang, nk=nk, ds=ds, growth=growth,
                              smooth_strength=smooth, smooth_iters=smooth_iters,
                              steps=k_steps)
            write_geo(str(yang_geo_path), yang_grid)
            print(f"  {yang_geo_path} written -- shape {yang_grid.shape[:3]}.")
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
    from shore.io.vtk import write_double_vtm, write_vts

    vtk_paths: list[Path] = []

    if topology == "single":
        for geo_path in geo_paths:
            vts_path = geo_path.with_suffix(".vts")
            print(f"Exporting {geo_path} → {vts_path} …")
            grid = read_geo(geo_path)
            write_vts(vts_path, grid, with_jacobian=True)
            print(f"  {vts_path} written.")
            vtk_paths.append(vts_path)

    elif topology == "caps":
        from shore.io.vtk import write_caps_vtm

        # Individual .vts files for each of the 6 blocks
        blocks = []
        for geo_path in geo_paths:
            grid = read_geo(geo_path)
            blocks.append(grid)
            vts_path = geo_path.with_suffix(".vts")
            print(f"Exporting {geo_path} → {vts_path} …")
            from shore.io.vtk import open_to_pyvista
            open_to_pyvista(grid, with_jacobian=True).save(str(vts_path))
            print(f"  {vts_path} written.")
            vtk_paths.append(vts_path)

        # Combined MultiBlock
        stem = geo_paths[0].stem.rsplit("_sub0", 1)[0]
        vtm_path = geo_paths[0].parent / f"{stem}_caps.vtm"
        print(f"Exporting combined sphere-cap MultiBlock → {vtm_path} …")
        write_caps_vtm(vtm_path, blocks, with_jacobian=True)
        print(f"  {vtm_path} written.")
        vtk_paths.append(vtm_path)

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
        choices=["single", "double", "caps"],
        help="'single' (default) for single O-grid, 'double' for Yin-Yang, 'caps' for 6-block sphere-cap",
    )
    gen.add_argument(
        "--vtk", action="store_true",
        help="Export generated .geo file(s) to VTK format after generation. "
             "For single: writes .vts. For double: writes .vts per block plus .vtm MultiBlock.",
    )
    gen.add_argument("--ni",     type=int,   default=40,   help="Latitude points")
    gen.add_argument("--nj",     type=int,   default=60,   help="Longitude points")
    gen.add_argument("--nk",     type=int,   default=20,   help="Wall-normal layers")
    gen.add_argument("--ds",     type=float, default=0.02, help="First-layer thickness (geometric)")
    gen.add_argument("--growth", type=float, default=1.10, help="Layer growth ratio (geometric)")
    gen.add_argument("--smooth", type=float, default=0.2,  help="Laplacian smooth strength [0-1]")
    gen.add_argument("--smooth-iters", type=int, default=2, help="Laplacian sweeps per layer")

    # surface spacing
    gen.add_argument(
        "--surface-i-spacing", default="uniform",
        choices=["uniform", "tanh", "tanh2"],
        help="Latitude spacing law",
    )
    gen.add_argument("--surface-i-beta", type=float, default=3.0,
                     help="Clustering strength for latitude tanh/tanh2")
    gen.add_argument(
        "--surface-j-spacing", default="uniform",
        choices=["uniform", "tanh2"],
        help="Longitude spacing law",
    )
    gen.add_argument("--surface-j-beta", type=float, default=3.0,
                     help="Clustering strength for longitude tanh2")
    gen.add_argument("--j-cluster-lon", type=float, default=None,
                     help="Longitude to cluster toward (radians, tanh2 only)")

    # volume spacing
    gen.add_argument(
        "--volume-k-spacing", default="geometric",
        choices=["geometric", "tanh", "tanh2"],
        help="Wall-normal spacing law",
    )
    gen.add_argument("--volume-k-beta", type=float, default=3.0,
                     help="Clustering strength for tanh/tanh2")
    gen.add_argument("--volume-k-thickness", type=float, default=None,
                     help="Total wall-normal thickness (required for tanh/tanh2)")
    theta_grp = gen.add_mutually_exclusive_group()
    theta_grp.add_argument(
        "--theta-cap", type=float, default=None, metavar="RAD",
        help="Polar exclusion angle in radians (default 0.05*pi ≈ 9°). "
             "Cap boundary parallel radius = sin(theta_cap) on a unit sphere.",
    )
    theta_grp.add_argument(
        "--theta-cap-deg", type=float, default=None, metavar="DEG",
        help="Polar exclusion angle in degrees (alias for --theta-cap).",
    )
    gen.add_argument(
        "--cap-pull", type=float, default=0.0, metavar="PULL",
        help="[caps topology] Blend factor [0,1) doming the cap interior toward the pole. "
             "0=off (default). Values 0.3-0.5 improve interior cell quality.",
    )
    gen.add_argument(
        "--cap-squareness", type=float, default=0.0, metavar="SQ",
        help="[caps topology] Blend factor [0,1] straightening each cap boundary edge "
             "from circular arc toward geodesic chord (inscribed square). "
             "0=circular (default), 1=full chord (~90 deg corner angle).",
    )
    gen.add_argument(
        "--correct-spacing", action="store_true",
        help="[caps topology] Equalise cell widths at all block interfaces after marching.",
    )

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
        # Resolve theta_cap from mutually-exclusive group
        theta_cap: float | None = None
        if args.theta_cap_deg is not None:
            theta_cap = np.deg2rad(args.theta_cap_deg)
        elif args.theta_cap is not None:
            theta_cap = args.theta_cap

        stl_path = Path("sphere.stl")
        npz_path = Path("sphere_k0.npz") if args.topology == "single" else None
        geo_path = Path("sphere.geo")
        geo_paths = generate_sphere(
            stl_path, geo_path, npz_path,
            ni=args.ni, nj=args.nj, nk=args.nk,
            ds=args.ds, growth=args.growth,
            topology=args.topology,
            theta_cap=theta_cap,
            cap_pull=args.cap_pull,
            cap_squareness=args.cap_squareness,
            correct_spacing=args.correct_spacing,
            smooth=args.smooth,
            smooth_iters=args.smooth_iters,
            surface_i_spacing=args.surface_i_spacing,
            surface_i_beta=args.surface_i_beta,
            surface_j_spacing=args.surface_j_spacing,
            surface_j_beta=args.surface_j_beta,
            j_cluster_lon=args.j_cluster_lon,
            volume_k_spacing=args.volume_k_spacing,
            volume_k_beta=args.volume_k_beta,
            volume_k_thickness=args.volume_k_thickness,
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
                    print("Warning: .npz has no 'k0' array.", file=sys.stderr)
            except Exception as e:
                print(f"Warning: VTK export failed: {e}", file=sys.stderr)
        elif args.generate:
            pass  # VTK already done during generation
        else:
            print("Warning: --vtk requires a .geo or .npz input.", file=sys.stderr)

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