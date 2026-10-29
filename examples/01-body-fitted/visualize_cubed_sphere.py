"""Standalone example: 6-block cubed-sphere mesh around a unit sphere.

When ``--generate`` is given, this script runs the full SHORE
body-fitted pipeline against an internally-built unit icosphere
through the cubed-sphere topology, writing all 6 ``.geo`` blocks
(``sphere_sub0..3.geo`` + ``sphere_cap_north/south.geo``) and
optionally a ``.vtm`` MultiBlock for ParaView.  The cap and equator
spacing flags below are specific to this topology — see
``examples/02-primitives/visualize_box.py`` / ``visualize_annulus.py`` /
``visualize_flat_caps.py`` for the parametric primitives.

When called *without* ``--generate``, the script is a generic
``.geo`` / ``.npz`` PyVista viewer — it will read and display any
structured grid produced by SHORE, regardless of topology.

Run from the project root after activating .venv:

    # Generate a cubed-sphere mesh on a unit icosphere and open the viewer:
    python examples/01-body-fitted/visualize_cubed_sphere.py --generate

    # Generate and export all 6 blocks to VTK MultiBlock:
    python examples/01-body-fitted/visualize_cubed_sphere.py --generate --vtk

    # Tanh-spaced boundary layer with fixed total thickness:
    python examples/01-body-fitted/visualize_cubed_sphere.py --generate \\
        --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 2.0

    # Latitude clustering (toward cap boundary) + longitude clustering:
    python examples/01-body-fitted/visualize_cubed_sphere.py --generate \\
        --surface-i-spacing tanh --surface-i-beta 4.0 \\
        --surface-j-spacing tanh2 --surface-j-beta 5.0

    # Visualise an existing file (any topology):
    python examples/01-body-fitted/visualize_cubed_sphere.py sphere_sub0.geo
    python examples/01-body-fitted/visualize_cubed_sphere.py sphere_k0.npz
    python examples/01-body-fitted/visualize_cubed_sphere.py sphere_sub0.geo --screenshot out.png

Cap/equator C1 seam matching
----------------------------
The cap k=0 layer is built by gnomonic flat-square projection (a tangent-plane
square radially projected onto the body).  The equator pins its first/last
meridional step to the cap's first inward step at every seam vertex, giving
C1 spacing continuity at the seam.  This pinning is only stable when the cap
and equator radial steps are within a factor of ~2 of each other:

    theta_cap / (nc - 1)   ~   (pi - 2*theta_cap) / (ni - 1)

with ``nc = nj // 4 + 1``.  If they differ by more than 2x the equator falls
back to its unpinned i_spacing law (a UserWarning is emitted).  The default
``--ni 40 --nj 60 --theta-cap (pi/6)`` lands inside the C1 window.

For Jupyter use, see examples/01-body-fitted/tour.ipynb.
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
    smooth: float = 0.2,
    smooth_iters: int = 2,
    surface_i_spacing: str = "uniform",
    surface_i_beta: float = 3.0,
    surface_j_spacing: str = "uniform",
    surface_j_beta: float = 3.0,
    volume_k_spacing: str = "geometric",
    volume_k_beta: float = 3.0,
    volume_k_thickness: float | None = None,
) -> list[Path]:
    """Create sphere.stl + 6 sphere_<label>.geo files via the
    cubed-sphere topology.

    Returns list of generated .geo paths.
    """
    import trimesh

    from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

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

    # --- march ---
    labels = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]
    geo_paths = [geo_stem.parent / f"{geo_stem.name}_{lbl}.geo" for lbl in labels]
    first_geo = geo_paths[0]

    if not first_geo.exists():
        print(f"Generating 6-block cubed-sphere topology ({ni}x{nj}, {nk} layers) ...")
        m = Mesh.from_stl(
            stl_path,
            topology=CubedSphereTopology(),
            ni=ni, nj=nj, nk=nk,
            theta_cap=theta_cap,
            i_spacing=surface_i_spacing, i_beta=surface_i_beta,
            j_spacing=surface_j_spacing, j_beta=surface_j_beta,
            spacing_k=spacing_k,
        )

        steps = spacing_k.build_steps(nk)

        def _cb(k: int, jmin: float) -> None:
            print(f"  layer {k:3d}  jmin={jmin:.3e}", end="\r")

        m.march(
            ds=ds, growth=growth,
            smooth_strength=smooth, smooth_iters=smooth_iters,
            steps=steps, layer_callback=_cb,
        )
        print()

        # Save the k=0 cap surface as a quick reference (single block).
        if not npz_path.exists():
            from shore.mesh import Body
            Body.save_surface(m.blocks[4].nodes[:, :, 0, :], npz_path)
            print(f"  k=0 cap_north surface saved → {npz_path}")

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
    # Deprecated: gnomonic cap is symmetric in 4 quadrants, no cluster_lon.
    gen.add_argument("--j-cluster-lon", type=float, default=None,
                     help=argparse.SUPPRESS)

    gen.add_argument("--volume-k-spacing", default="geometric",
                     choices=["geometric", "tanh", "tanh2"], help="Wall-normal spacing law")
    gen.add_argument("--volume-k-beta",      type=float, default=3.0)
    gen.add_argument("--volume-k-thickness", type=float, default=None,
                     help="Total wall-normal thickness (required for tanh/tanh2)")

    theta_grp = gen.add_mutually_exclusive_group()
    theta_grp.add_argument("--theta-cap",     type=float, default=None, metavar="RAD",
                            help="Polar exclusion angle in radians "
                                 "(default pi/6 = 30 deg). For C1 seam "
                                 "matching, choose ni so that "
                                 "theta_cap/(nc-1) is within a factor of 2 "
                                 "of (pi-2*theta_cap)/(ni-1).")
    theta_grp.add_argument("--theta-cap-deg", type=float, default=None, metavar="DEG",
                            help="Polar exclusion angle in degrees")

    # Deprecated cap-quality knobs: accepted for backwards-compat, ignored.
    gen.add_argument("--cap-squareness", type=float, default=None,
                     help=argparse.SUPPRESS)
    gen.add_argument("--cap-pull", type=float, default=None,
                     help=argparse.SUPPRESS)

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

    if args.cap_pull is not None:
        import warnings as _warn
        _warn.warn(
            "--cap-pull is deprecated and ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
    if args.cap_squareness is not None:
        import warnings as _warn
        _warn.warn(
            "--cap-squareness is deprecated and ignored: cap k=0 is now built "
            "by gnomonic flat-square projection.",
            DeprecationWarning,
            stacklevel=2,
        )
    if args.j_cluster_lon is not None:
        import warnings as _warn
        _warn.warn(
            "--j-cluster-lon is deprecated and ignored.",
            DeprecationWarning,
            stacklevel=2,
        )

    # --- generation ---
    geo_paths: list[Path] = []
    if args.generate:
        # Default theta_cap = pi/6 (30 deg).  Combined with the default
        # ni=40/nj=60 this lands inside the cap/equator C1 pinning window
        # (cap radial step ~2 deg, equator step ~3 deg, ratio 1.5x).
        if args.theta_cap_deg is not None:
            theta_cap: float = float(np.deg2rad(args.theta_cap_deg))
        elif args.theta_cap is not None:
            theta_cap = float(args.theta_cap)
        else:
            theta_cap = float(np.pi / 6.0)

        geo_paths = generate_sphere(
            stl_path=Path("sphere.stl"),
            npz_path=Path("sphere_k0.npz"),
            geo_stem=Path("sphere"),
            ni=args.ni, nj=args.nj, nk=args.nk,
            ds=args.ds, growth=args.growth,
            theta_cap=theta_cap,
            smooth=args.smooth,
            smooth_iters=args.smooth_iters,
            surface_i_spacing=args.surface_i_spacing,
            surface_i_beta=args.surface_i_beta,
            surface_j_spacing=args.surface_j_spacing,
            surface_j_beta=args.surface_j_beta,
            volume_k_spacing=args.volume_k_spacing,
            volume_k_beta=args.volume_k_beta,
            volume_k_thickness=args.volume_k_thickness,
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
