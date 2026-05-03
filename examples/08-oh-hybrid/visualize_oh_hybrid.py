"""Standalone example: 12-block OH hybrid (cubed-sphere wrap + H-channel).

The OH topology is a Style-2 conforming hybrid: 6 cubed-sphere wrap
blocks (marched outward from the STL surface) plus 6 channel-frustum
blocks (ruled-surface from the wrap envelope to the user-supplied
axis-aligned bounding box).  Wrap and channel are SHARED-by-view at
the seam — they share one contiguous numpy buffer per pair, so
``wrap.k_hi`` and ``channel.k_lo`` are bit-identical with no copy.

Default block layout:

    sub0..sub3        4 equatorial wrap sub-blocks (marched)
    cap_north         polar cap wrap (marched, faces -z)
    cap_south         polar cap wrap (marched, faces +z)
    chan_sub0..3      4 equatorial channel frustums (filled post-march)
    chan_north        north-cap channel (faces -z box face)
    chan_south        south-cap channel (faces +z box face)

The wrap is constructed with ``theta_cap = pi/4`` and ``phi_offset = pi/4``
so each of the 6 wrap blocks faces exactly one of the bounding box's
6 faces — the key geometric trick that makes a 1:1 wrap-to-channel
mapping possible.  The channel-channel lateral seams sit on the cube
*edges*, with bit-exact memory aliasing.

When ``--generate`` is given, this script runs the full SHORE pipeline
through the OH topology against an internally-built unit icosphere,
writing all 12 ``.geo`` blocks and a ``.vtm`` MultiBlock for ParaView.

All outputs land in ``examples_out/`` by default; pass ``-o /tmp/run`` to
override.  Run from the project root after activating .venv:

    # Generate a default OH hybrid mesh and open the viewer:
    python examples/08-oh-hybrid/visualize_oh_hybrid.py --generate

    # Generate and export to VTK MultiBlock (ParaView-friendly):
    python examples/08-oh-hybrid/visualize_oh_hybrid.py --generate --vtk

    # Larger / denser channels with tanh clustering toward the wrap seam:
    python examples/08-oh-hybrid/visualize_oh_hybrid.py --generate \\
        --box-half 5.0 --nk-channel 20 \\
        --channel-k-law tanh --channel-k-beta 4.0

    # Visualise an already-generated block:
    python examples/08-oh-hybrid/visualize_oh_hybrid.py examples_out/oh_sub0.geo
    python examples/08-oh-hybrid/visualize_oh_hybrid.py examples_out/oh_chan_north.geo

Note
----
``OHTopology`` is a Python-only API at the moment.  ``shore mesh``
(the CLI) currently dispatches only ``cubed_sphere`` and ``ogrid``;
the OH topology has no ``--topology oh-hybrid`` flag yet.  Drive it
from Python via :class:`shore.mesh.OHTopology`, as this script does.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


# Block layout returned by OHTopology.build — kept in sync with
# shore.mesh.topology.OHTopology so writers know each block's role.
_OH_BLOCK_LABELS = [
    "sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south",
    "chan_sub0", "chan_sub1", "chan_sub2", "chan_sub3",
    "chan_north", "chan_south",
]


def generate_oh_hybrid(
    stl_path: Path,
    geo_stem: Path,
    *,
    ni: int,
    nj: int,
    nk_wrap: int,
    nk_channel: int,
    ds: float,
    growth: float,
    box_vmin: tuple[float, float, float],
    box_vmax: tuple[float, float, float],
    channel_k_law: str,
    channel_k_beta: float,
    smooth: float = 0.2,
    smooth_iters: int = 2,
    blend_normals_k: int = 0,
) -> list[Path]:
    """Create sphere.stl + 12 oh_<label>.geo files via the OH hybrid topology.

    Returns the list of generated .geo paths in OHTopology block order.
    """
    import trimesh

    from shore.mesh import Mesh, OHTopology, Spacing1D

    spacing_k = Spacing1D(law="geometric", ds=ds, growth=growth)

    # --- STL ---
    if not stl_path.exists():
        print(f"Generating {stl_path} ...")
        sphere = trimesh.creation.icosphere(subdivisions=4, radius=1.0)
        sphere.export(str(stl_path))
        print(f"  {len(sphere.faces)} triangles written.")
    else:
        print(f"  {stl_path} already exists, skipping STL generation.")

    geo_paths = [geo_stem.parent / f"{geo_stem.name}_{lbl}.geo" for lbl in _OH_BLOCK_LABELS]
    first_geo = geo_paths[0]

    if first_geo.exists():
        print(f"  {first_geo} already exists, skipping marching.")
        return geo_paths

    # --- OH hybrid build + march ---
    print(
        f"Generating 12-block OH hybrid topology "
        f"({ni}x{nj}, wrap={nk_wrap} channel={nk_channel}) "
        f"box=[{box_vmin}, {box_vmax}] ..."
    )
    m = Mesh.from_stl(
        stl_path,
        topology=OHTopology(
            box_vmin=box_vmin,
            box_vmax=box_vmax,
            nk_channel=nk_channel,
            k_law_channel=channel_k_law,
            k_beta_channel=channel_k_beta,
        ),
        ni=ni,
        nj=nj,
        nk=nk_wrap,
        spacing_k=spacing_k,
    )

    def _cb(k: int, jmin: float) -> None:
        print(f"  wrap layer {k:3d}  jmin={jmin:.3e}", end="\r")

    m.march(
        ds=ds, growth=growth,
        smooth_strength=smooth, smooth_iters=smooth_iters,
        layer_callback=_cb,
        blend_normals_k=blend_normals_k,
    )
    print()
    print("  channel blocks filled by ruled-surface from wrap envelope to box.")

    m.write_geo(geo_stem)
    for p, blk in zip(geo_paths, m.blocks, strict=True):
        print(f"  {p} written -- shape {blk.nodes.shape[:3]}.")

    return geo_paths


def export_to_vtk(geo_stem: Path, geo_paths: list[Path]) -> Path:
    """Export the 12 OH blocks to a single PyVista MultiBlock .vtm."""
    from shore.io.geo import read_geo
    from shore.io.vtk import write_caps_vtm

    vtm_path = geo_stem.parent / f"{geo_stem.name}.vtm"
    print(f"Exporting 12 blocks -> {vtm_path} ...")
    blocks = [read_geo(p) for p in geo_paths]
    write_caps_vtm(vtm_path, blocks, with_jacobian=True, labels=_OH_BLOCK_LABELS)
    print(f"  {vtm_path} written.")
    return vtm_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build / visualise a SHORE OH hybrid mesh (cubed-sphere wrap + H-channel).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    gen = parser.add_argument_group("OH hybrid generation (--generate)")
    gen.add_argument("--generate", action="store_true",
                     help="Generate sphere.stl + 12 oh_<label>.geo files.")
    gen.add_argument("--vtk", action="store_true",
                     help="Export the 12 .geo blocks to a single .vtm MultiBlock.")
    gen.add_argument("--ni",         type=int,   default=24,   help="Latitude points per wrap block.")
    gen.add_argument("--nj",         type=int,   default=32,   help="Longitude points (multiple of 4).")
    gen.add_argument("--nk-wrap",    type=int,   default=20,   help="Wrap wall-normal layers.")
    gen.add_argument("--nk-channel", type=int,   default=10,   help="Channel wall-normal layers.")
    gen.add_argument("--ds",         type=float, default=0.02, help="Wrap first-layer thickness.")
    gen.add_argument("--growth",     type=float, default=1.10, help="Wrap geometric growth ratio.")
    gen.add_argument("--smooth",     type=float, default=0.2,  help="Laplacian smooth strength [0-1].")
    gen.add_argument("--smooth-iters", type=int, default=2,    help="Laplacian sweeps per layer.")

    gen.add_argument("--box-half", type=float, default=3.0,
                     help="Half-extent of the cube bounding box (centred on origin).")
    gen.add_argument("--channel-k-law", default="uniform",
                     choices=["uniform", "tanh", "tanh2"],
                     help="Channel wall-normal spacing law (k=0 is the wrap seam, "
                          "k=last is the box face).")
    gen.add_argument("--channel-k-beta", type=float, default=3.0,
                     help="Clustering strength for tanh / tanh2 channel spacing.")
    gen.add_argument("--blend-normals-k", type=int, default=0,
                     help="Layers over which to blend STL surface normals (0 = off).")

    parser.add_argument("input", type=Path, nargs="?", default=None,
                        help=".geo file to visualise (any single OH block).")
    parser.add_argument("--colorby", default="",
                        help="Colour field: 'Jacobian' (default for .geo) or 'none'.")
    parser.add_argument("--screenshot", type=Path, default=None,
                        help="Save a PNG screenshot instead of opening the interactive window.")
    parser.add_argument("--no-edges", action="store_true", help="Hide cell edges.")
    parser.add_argument("--no-view", action="store_true", help="Generate only, skip viewer.")
    parser.add_argument("--opacity", type=float, default=0.6, help="Grid opacity [0-1].")
    parser.add_argument("-o", "--output", type=Path, default=Path("examples_out"),
                        help="Output directory; .stl / .geo / .vtm files land here.")

    args = parser.parse_args()
    out_dir: Path = args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    geo_paths: list[Path] = []
    if args.generate:
        half = float(args.box_half)
        if half <= 1.0:
            print(
                "Error: --box-half must be > 1.0 (the unit-sphere body radius).",
                file=sys.stderr,
            )
            sys.exit(1)
        box_vmin = (-half, -half, -half)
        box_vmax = (+half, +half, +half)

        geo_paths = generate_oh_hybrid(
            stl_path=out_dir / "sphere.stl",
            geo_stem=out_dir / "oh",
            ni=args.ni, nj=args.nj,
            nk_wrap=args.nk_wrap, nk_channel=args.nk_channel,
            ds=args.ds, growth=args.growth,
            box_vmin=box_vmin, box_vmax=box_vmax,
            channel_k_law=args.channel_k_law,
            channel_k_beta=args.channel_k_beta,
            smooth=args.smooth, smooth_iters=args.smooth_iters,
            blend_normals_k=args.blend_normals_k,
        )

        if args.vtk:
            try:
                export_to_vtk(out_dir / "oh", geo_paths)
            except Exception as exc:
                print(f"Warning: VTK export failed: {exc}", file=sys.stderr)

        if args.no_view:
            return
        if args.input is None:
            args.input = geo_paths[0]  # default: first wrap sub-block

    # --- viewer / standalone export ---
    if args.input is None:
        parser.error("Provide a .geo file or use --generate.")

    if not args.input.exists():
        print(f"Error: file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if args.input.suffix.lower() != ".geo":
        print(
            f"Error: unsupported extension '{args.input.suffix}'. Expected .geo.",
            file=sys.stderr,
        )
        sys.exit(1)

    from shore.io.geo import read_geo
    from shore.vis import save_screenshot, view_geo

    grid = read_geo(args.input)
    colorby = args.colorby or "Jacobian"

    if args.screenshot:
        save_screenshot(grid, str(args.screenshot))
        print(f"Screenshot saved to {args.screenshot}")
    else:
        view_geo(grid, colorby=colorby, show_edges=not args.no_edges, opacity=args.opacity)


if __name__ == "__main__":
    main()
