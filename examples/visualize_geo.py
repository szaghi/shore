"""Standalone visualisation script for SHORE outputs.

Run from the project root after activating .venv:

    # Generate a unit-sphere caps mesh and open the viewer:
    python examples/visualize_geo.py --generate

    # Generate and export all 6 blocks to VTK MultiBlock:
    python examples/visualize_geo.py --generate --vtk

    # Tanh-spaced boundary layer with fixed total thickness:
    python examples/visualize_geo.py --generate \\
        --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 2.0

    # Latitude clustering (toward cap boundary) + longitude clustering:
    python examples/visualize_geo.py --generate \\
        --surface-i-spacing tanh --surface-i-beta 4.0 \\
        --surface-j-spacing tanh2 --surface-j-beta 5.0

    # Cap quality controls:
    python examples/visualize_geo.py --generate \\
        --cap-pull 0.4 --cap-squareness 0.5

    # Visualise an existing file:
    python examples/visualize_geo.py sphere_sub0.geo
    python examples/visualize_geo.py sphere_k0.npz
    python examples/visualize_geo.py sphere_sub0.geo --screenshot out.png

For Jupyter use, see examples/visualize.ipynb.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def generate_sphere(
    stl_path: Path,
    npz_path: Path,
    geo_stem: Path,
    ni: int,
    nj: int,
    nk: int,
    ds: float,
    growth: float,
    theta_cap: float | None = None,
    cap_pull: float = 0.0,
    cap_squareness: float = 0.0,
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
    blend_normals_k: int = 0,
) -> list[Path]:
    """Create sphere.stl, sphere_k0.npz, and 6 sphere_<label>.geo files.

    Returns list of generated .geo paths.
    """
    import trimesh

    from shore.mesh import Body, CapsTopology, Mesh, Spacing1D

    if volume_k_spacing in {"tanh", "tanh2"} and volume_k_thickness is None:
        print(
            "Error: --volume-k-thickness is required for tanh/tanh2 volume spacing.",
            file=sys.stderr,
        )
        sys.exit(1)

    spacing_k = Spacing1D(
        law=volume_k_spacing,
        ds=ds,
        growth=growth,
        beta=volume_k_beta,
        total_thickness=volume_k_thickness,
    )

    # --- STL ---
    if not stl_path.exists():
        print(f"Generating {stl_path} ...")
        sphere = trimesh.creation.icosphere(subdivisions=4, radius=1.0)
        sphere.export(str(stl_path))
        print(f"  {len(sphere.faces)} triangles written.")
    else:
        print(f"  {stl_path} already exists, skipping STL generation.")

    body = Body.from_stl(stl_path)

    # --- surface projection ---
    first_geo = geo_stem.parent / f"{geo_stem.name}_sub0.geo"
    if not npz_path.exists():
        print(f"Projecting {ni}×{nj} unified surface -> {npz_path} ...")
        k0, nc_cap = body.project_caps(
            ni, nj,
            theta_cap=theta_cap,
            i_spacing=surface_i_spacing, i_beta=surface_i_beta,
            j_spacing=surface_j_spacing, j_beta=surface_j_beta,
            j_cluster_lon=j_cluster_lon,
        )
        Body.save_surface(k0, npz_path)
        print(f"  Surface layer shape: {k0.shape}, nc_cap={nc_cap}.")
    else:
        print(f"  {npz_path} already exists, skipping projection.")
        nc_cap = None  # will be recomputed from shape when loading

    # --- march ---
    labels = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]
    geo_paths = [geo_stem.parent / f"{geo_stem.name}_{lbl}.geo" for lbl in labels]

    if not first_geo.exists():
        print(f"Generating 6-block caps topology ({nk} layers) ...")
        k0 = Body.load_surface(npz_path)

        if nc_cap is None:
            # Derive nc_cap from the unified surface shape: ni_total = ni + 2*(nc_cap-1)
            ni_total = k0.shape[0]
            nc_cap = (ni_total - ni) // 2 + 1

        topo = CapsTopology(cap_pull=cap_pull, cap_squareness=cap_squareness)
        blocks = topo.build(k0, nk, nc_cap=nc_cap, spacing_k=spacing_k, mesh=body.mesh)
        m = Mesh(body=body, blocks=blocks, topology=topo)

        steps = spacing_k.build_steps(nk)

        def _cb(k: int, jmin: float) -> None:
            print(f"  layer {k:3d}  jmin={jmin:.3e}", end="\r")

        m.march(
            ds=ds, growth=growth,
            smooth_strength=smooth, smooth_iters=smooth_iters,
            steps=steps, layer_callback=_cb,
            blend_normals_k=blend_normals_k,
        )
        print()
        m.write_geo(geo_stem)
        for p, blk in zip(geo_paths, m.blocks, strict=True):
            print(f"  {p} written -- shape {blk.nodes.shape[:3]}.")
    else:
        print(f"  {first_geo} already exists, skipping marching.")

    return geo_paths


def export_to_vtk(geo_paths: list[Path]) -> Path:
    """Export 6 .geo block files to a single VTK MultiBlock .vtm file.

    Returns the written .vtm path.
    """
    from shore.io.geo import read_geo
    from shore.io.vtk import write_caps_vtm

    stem = geo_paths[0].stem.rsplit("_sub0", 1)[0]
    vtm_path = geo_paths[0].parent / f"{stem}.vtm"

    print(f"Exporting 6 blocks -> {vtm_path} ...")
    blocks = [read_geo(p) for p in geo_paths]
    write_caps_vtm(vtm_path, blocks, with_jacobian=True)
    print(f"  {vtm_path} written.")
    return vtm_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualise a SHORE .geo or .npz file with PyVista.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    gen = parser.add_argument_group("sphere generation (--generate)")
    gen.add_argument("--generate", action="store_true",
                     help="Generate sphere.stl + sphere_k0.npz + 6 sphere_<label>.geo files")
    gen.add_argument("--vtk", action="store_true",
                     help="Export generated .geo files to a VTK MultiBlock .vtm")
    gen.add_argument("--ni",     type=int,   default=40,   help="Latitude points")
    gen.add_argument("--nj",     type=int,   default=60,   help="Longitude points (multiple of 4)")
    gen.add_argument("--nk",     type=int,   default=20,   help="Wall-normal layers")
    gen.add_argument("--ds",     type=float, default=0.02, help="First-layer thickness")
    gen.add_argument("--growth", type=float, default=1.10, help="Layer growth ratio")
    gen.add_argument("--smooth", type=float, default=0.2,  help="Laplacian smooth strength [0-1]")
    gen.add_argument("--smooth-iters", type=int, default=2, help="Laplacian sweeps per layer")

    gen.add_argument("--surface-i-spacing", default="uniform",
                     choices=["uniform", "tanh", "tanh2"], help="Latitude spacing law")
    gen.add_argument("--surface-i-beta", type=float, default=3.0)
    gen.add_argument("--surface-j-spacing", default="uniform",
                     choices=["uniform", "tanh2"], help="Longitude spacing law")
    gen.add_argument("--surface-j-beta", type=float, default=3.0)
    gen.add_argument("--j-cluster-lon", type=float, default=None,
                     help="Longitude to cluster toward (radians, tanh2 only)")

    gen.add_argument("--volume-k-spacing", default="geometric",
                     choices=["geometric", "tanh", "tanh2"], help="Wall-normal spacing law")
    gen.add_argument("--volume-k-beta",      type=float, default=3.0)
    gen.add_argument("--volume-k-thickness", type=float, default=None,
                     help="Total wall-normal thickness (required for tanh/tanh2)")

    theta_grp = gen.add_mutually_exclusive_group()
    theta_grp.add_argument("--theta-cap",     type=float, default=None, metavar="RAD",
                            help="Polar exclusion angle in radians (default 0.05*pi ≈ 9°)")
    theta_grp.add_argument("--theta-cap-deg", type=float, default=None, metavar="DEG",
                            help="Polar exclusion angle in degrees")

    gen.add_argument("--cap-pull",       type=float, default=0.0,
                     help="Dome factor [0,1) for cap interior toward pole")
    gen.add_argument("--cap-squareness", type=float, default=0.0,
                     help="Straighten cap boundary edges toward geodesic chord [0,1]")
    gen.add_argument("--blend-normals-k", type=int, default=0,
                     help="Layers over which to blend STL surface normals (0 = off)")

    parser.add_argument("input", type=Path, nargs="?", default=None,
                        help=".geo or .npz file to visualise")
    parser.add_argument("--colorby",    default="",
                        help="Colour field: 'Jacobian' (geo), 'radius' (npz), or 'none'")
    parser.add_argument("--screenshot", type=Path, default=None,
                        help="Save a PNG screenshot instead of opening interactive window")
    parser.add_argument("--no-edges",   action="store_true", help="Hide cell edges")
    parser.add_argument("--no-view",    action="store_true", help="Generate only, skip viewer")
    parser.add_argument("--opacity",    type=float, default=0.6, help="Grid opacity [0-1]")

    args = parser.parse_args()

    import numpy as np

    # --- generation ---
    geo_paths: list[Path] = []
    if args.generate:
        theta_cap: float | None = None
        if args.theta_cap_deg is not None:
            theta_cap = np.deg2rad(args.theta_cap_deg)
        elif args.theta_cap is not None:
            theta_cap = args.theta_cap

        geo_paths = generate_sphere(
            stl_path=Path("sphere.stl"),
            npz_path=Path("sphere_k0.npz"),
            geo_stem=Path("sphere"),
            ni=args.ni, nj=args.nj, nk=args.nk,
            ds=args.ds, growth=args.growth,
            theta_cap=theta_cap,
            cap_pull=args.cap_pull,
            cap_squareness=args.cap_squareness,
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
            blend_normals_k=args.blend_normals_k,
        )

        if args.vtk:
            try:
                export_to_vtk(geo_paths)
            except Exception as exc:
                print(f"Warning: VTK export failed: {exc}", file=sys.stderr)

        if args.no_view:
            return
        if args.input is None:
            args.input = geo_paths[0]   # default: first sub-block

    # --- VTK export for an existing .geo or .npz (no --generate) ---
    if args.vtk and not args.generate and args.input is not None:
        try:
            from shore.mesh import export_geo, export_surface
            if args.input.suffix == ".geo":
                out = export_geo(args.input)
                print(f"Exported {out}")
            elif args.input.suffix == ".npz":
                out = export_surface(args.input)
                print(f"Exported {out}")
            else:
                print("Warning: --vtk requires a .geo or .npz input.", file=sys.stderr)
        except Exception as exc:
            print(f"Warning: VTK export failed: {exc}", file=sys.stderr)

    # --- resolve input ---
    if args.input is None:
        parser.error("Provide an input file or use --generate.")

    if not args.input.exists():
        print(f"Error: file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    suffix = args.input.suffix.lower()

    # --- visualise ---
    if suffix == ".geo":
        from shore.io.geo import read_geo
        from shore.vis import save_screenshot, view_geo

        grid    = read_geo(args.input)
        colorby = args.colorby or "Jacobian"

        if args.screenshot:
            save_screenshot(grid, str(args.screenshot))
            print(f"Screenshot saved to {args.screenshot}")
        else:
            view_geo(grid, colorby=colorby, show_edges=not args.no_edges, opacity=args.opacity)

    elif suffix == ".npz":
        from shore.mesh import Body
        from shore.vis import view_surface

        k0      = Body.load_surface(args.input)
        colorby = args.colorby or "radius"
        view_surface(k0, colorby=colorby, show_edges=not args.no_edges)

    else:
        print(f"Error: unsupported extension '{suffix}'. Expected .geo or .npz.",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
