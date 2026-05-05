"""Standalone example: 12-/14-block CH hybrid (C-grid wrap + H-channel).

The CH topology is a sharp-/blunt-TE airfoil analogue of OH: a C-grid
wrap (3 or 4 blocks: ``upper`` / ``body`` / ``lower`` / optional
``wake_base``) plus 9 (sharp) or 10 (blunt) H-channel blocks bridging
the wrap envelope to a chord/normal/span-aligned far-field box.

The 5 ``k_hi``-side channel blocks are SHARED-by-view with their wrap
partners (one contiguous numpy buffer of depth ``nk_wrap+nk_channel-1``
per wrap segment); the 4 (or 5) outflow sub-blocks at the rear frame
are conformal SHARED to the wrap.upper.i_lo / lower.i_hi faces (wake
sub-blocks, nk = nk_wrap) and to chan_upper.i_lo / chan_lower.i_hi
(top sub-blocks, nk = nk_channel) -- proper structured-grid coupling
across the entire wake-far-tip seam.

Default block layout (sharp TE):

    upper, body, lower                                # wrap (marched)
    chan_upper, chan_body_top, chan_body_front,
        chan_body_bot, chan_lower                     # k_hi-side channels
    chan_outflow_top_upper, chan_outflow_wake_upper,
    chan_outflow_wake_lower, chan_outflow_top_lower   # outflow sub-blocks

For blunt-TE bodies a ``wake_base`` wrap block + a
``chan_outflow_wake_base`` outflow sub-block are added.

The CH topology is exposed via the ``shore mesh --topology ch`` CLI as
well; this script drives it from Python so the channel construction
is visible.

All outputs land in ``examples_out/`` by default; pass ``-o /tmp/run`` to
override.  Run from the project root after activating .venv:

    # Generate a CH mesh on a NACA 0012 + open the viewer:
    python examples/09-ch-hybrid/visualize_ch_hybrid.py --generate

    # Generate and export to VTK MultiBlock (ParaView-friendly):
    python examples/09-ch-hybrid/visualize_ch_hybrid.py --generate --vtk

    # Larger / denser channels with tanh clustering toward the wrap seam:
    python examples/09-ch-hybrid/visualize_ch_hybrid.py --generate \\
        --chord-min -3.0 --chord-max 18.0 \\
        --normal-min -3.5 --normal-max 3.5 \\
        --nk-channel 16 --channel-k-law tanh --channel-k-beta 4.0

    # Visualise an already-generated block:
    python examples/09-ch-hybrid/visualize_ch_hybrid.py examples_out/ch_body.geo
    python examples/09-ch-hybrid/visualize_ch_hybrid.py examples_out/ch_chan_outflow_upper.geo
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def generate_ch_hybrid(
    stl_path: Path,
    geo_stem: Path,
    *,
    ni_body: int,
    ni_wake: int,
    n_stations: int,
    wake_length: float,
    nk_wrap: int,
    nk_channel: int,
    ds: float,
    growth: float,
    chord_min: float,
    chord_max: float,
    normal_min: float,
    normal_max: float,
    channel_k_law: str,
    channel_k_beta: float,
    smooth: float = 0.1,
    smooth_iters: int = 1,
) -> list[Path]:
    """Create naca0012.stl + 10 ch_<label>.geo files via the CH hybrid topology.

    Returns the generated .geo paths in CHTopology block order.
    """
    from shore.mesh import Body, CHTopology, Mesh, Spacing1D
    from shore.profiles import naca
    from shore.shapes import extrude_section

    # --- STL ---
    if not stl_path.exists():
        print(f"Generating {stl_path} (NACA 0012, sharp TE) ...")
        section = naca("0012", n_pts=40, closed_te=True)
        body_mesh = extrude_section(
            section, span=0.5, n_stations=n_stations, span_axis="z", end_caps=True
        )
        body_mesh.export(str(stl_path))
        print(f"  {len(body_mesh.faces)} triangles written.")
    else:
        print(f"  {stl_path} already exists, skipping STL generation.")

    body = Body.from_stl(stl_path, skip_hygiene=True)

    # --- CH hybrid build + march ---
    print(
        f"Generating CH hybrid topology "
        f"(ni_body={ni_body}, ni_wake={ni_wake}, n_stations={n_stations}, "
        f"wrap={nk_wrap} channel={nk_channel}) "
        f"frame chord=[{chord_min},{chord_max}] normal=[{normal_min},{normal_max}] ..."
    )
    spacing_k = Spacing1D(law="geometric", ds=ds, growth=growth)
    topo = CHTopology(
        chord_min=chord_min,
        chord_max=chord_max,
        normal_min=normal_min,
        normal_max=normal_max,
        nk_channel=nk_channel,
        k_law_channel=channel_k_law,
        k_beta_channel=channel_k_beta,
    )
    blocks = topo.build(
        surface=None,
        nk=nk_wrap,
        mesh=body.mesh,
        spacing_k=spacing_k,
        ni_body=ni_body,
        ni_wake=ni_wake,
        n_stations=n_stations,
        wake_length=wake_length,
    )
    m = Mesh(body=body, blocks=blocks, topology=topo)

    def _cb(k: int, jmin: float) -> None:
        print(f"  wrap layer {k:3d}  jmin={jmin:.3e}", end="\r")

    m.march(
        ds=ds,
        growth=growth,
        smooth_strength=smooth,
        smooth_iters=smooth_iters,
        layer_callback=_cb,
    )
    print()
    print("  channel blocks filled by ruled-surface from wrap envelope to frame.")

    geo_stem.parent.mkdir(parents=True, exist_ok=True)
    m.write_geo(geo_stem)
    geo_paths = [geo_stem.parent / f"{geo_stem.name}_{b.label}.geo" for b in m.blocks]
    for p, blk in zip(geo_paths, m.blocks, strict=True):
        print(f"  {p} written -- shape {blk.nodes.shape[:3]}.")

    return geo_paths


def export_to_vtk(geo_stem: Path, geo_paths: list[Path]) -> Path:
    """Export the CH blocks to a single PyVista MultiBlock .vtm."""
    from shore.io.geo import read_geo
    from shore.io.vtk import write_caps_vtm

    vtm_path = geo_stem.parent / f"{geo_stem.name}.vtm"
    print(f"Exporting {len(geo_paths)} blocks -> {vtm_path} ...")
    arrays = [read_geo(p) for p in geo_paths]
    labels = [p.stem.split("_", 1)[1] for p in geo_paths]
    write_caps_vtm(vtm_path, arrays, with_jacobian=True, labels=labels)
    print(f"  {vtm_path} written.")
    return vtm_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build / visualise a SHORE CH hybrid mesh (C-grid wrap + H-channel far-field).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    gen = parser.add_argument_group("CH hybrid generation (--generate)")
    gen.add_argument(
        "--generate", action="store_true", help="Generate naca0012.stl + ch_<label>.geo files."
    )
    gen.add_argument(
        "--vtk", action="store_true", help="Export the .geo blocks to a single .vtm MultiBlock."
    )
    gen.add_argument(
        "--ni-body", type=int, default=80, help="Body wrap node count (upper-TE -> LE -> lower-TE)."
    )
    gen.add_argument("--ni-wake", type=int, default=20, help="Per-branch wake node count.")
    gen.add_argument("--n-stations", type=int, default=5, help="Span stations.")
    gen.add_argument("--wake-length", type=float, default=15.0, help="Wake length.")
    gen.add_argument("--nk-wrap", type=int, default=12, help="Wrap layers.")
    gen.add_argument("--nk-channel", type=int, default=8, help="Channel layers.")
    gen.add_argument("--ds", type=float, default=5e-4, help="First-layer thickness.")
    gen.add_argument("--growth", type=float, default=1.08, help="Wrap growth ratio.")
    gen.add_argument("--smooth", type=float, default=0.0, help="Laplacian strength.")
    gen.add_argument("--smooth-iters", type=int, default=0, help="Laplacian sweeps.")

    gen.add_argument("--chord-min", type=float, default=-2.0)
    gen.add_argument("--chord-max", type=float, default=20.0)
    gen.add_argument("--normal-min", type=float, default=-2.5)
    gen.add_argument("--normal-max", type=float, default=2.5)
    gen.add_argument("--channel-k-law", default="uniform", choices=["uniform", "tanh", "tanh2"])
    gen.add_argument("--channel-k-beta", type=float, default=3.0)

    parser.add_argument(
        "input",
        type=Path,
        nargs="?",
        default=None,
        help=".geo file to visualise (any single CH block).",
    )
    parser.add_argument(
        "--colorby", default="", help="Colour field: 'Jacobian' (default) or 'none'."
    )
    parser.add_argument("--screenshot", type=Path, default=None)
    parser.add_argument("--no-edges", action="store_true")
    parser.add_argument("--no-view", action="store_true")
    parser.add_argument("--opacity", type=float, default=0.6)
    parser.add_argument("-o", "--output", type=Path, default=Path("examples_out"))

    args = parser.parse_args()
    out_dir: Path = args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    geo_paths: list[Path] = []
    if args.generate:
        geo_paths = generate_ch_hybrid(
            stl_path=out_dir / "naca0012.stl",
            geo_stem=out_dir / "ch",
            ni_body=args.ni_body,
            ni_wake=args.ni_wake,
            n_stations=args.n_stations,
            wake_length=args.wake_length,
            nk_wrap=args.nk_wrap,
            nk_channel=args.nk_channel,
            ds=args.ds,
            growth=args.growth,
            chord_min=args.chord_min,
            chord_max=args.chord_max,
            normal_min=args.normal_min,
            normal_max=args.normal_max,
            channel_k_law=args.channel_k_law,
            channel_k_beta=args.channel_k_beta,
            smooth=args.smooth,
            smooth_iters=args.smooth_iters,
        )

        if args.vtk:
            try:
                export_to_vtk(out_dir / "ch", geo_paths)
            except Exception as exc:
                print(f"Warning: VTK export failed: {exc}", file=sys.stderr)

        if args.no_view:
            return
        if args.input is None:
            args.input = geo_paths[1]  # default: body wrap

    # --- viewer ---
    if args.input is None:
        parser.error("Provide a .geo file or use --generate.")

    if not args.input.exists():
        print(f"Error: file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if args.input.suffix.lower() != ".geo":
        print(
            f"Error: unsupported extension '{args.input.suffix}'. Expected .geo.", file=sys.stderr
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
