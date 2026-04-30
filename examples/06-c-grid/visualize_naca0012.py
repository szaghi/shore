"""Standalone example: C-grid mesh around a NACA 0012 airfoil.

Builds a NACA 0012 STL from the analytic 4-digit profile, runs the
SHORE C-grid pipeline (TE detection -> span slice -> C-curve build ->
hyperbolic march), and writes ``.geo`` files plus a ParaView ``.vtm``.

By default everything lands under ``examples_out/`` (cwd-relative).
With the default 3-block layout (Xall-correct):

    examples_out/
        naca0012.stl              # input wing STL
        naca0012_upper.geo        # upper-wake block
        naca0012_body.geo         # body wrap block (TE -> suction -> LE -> pressure -> TE)
        naca0012_lower.geo        # lower-wake block
        naca0012.vtm              # ParaView MultiBlock (3 blocks)
        naca0012/                 # block payloads referenced by the .vtm

With ``--cgrid-blocks 1`` (legacy single-block layout):

    examples_out/
        naca0012.stl
        naca0012.geo              # one concatenated block (wake cut as i-overlap)
        naca0012.vtm

Run from the project root after activating .venv:

    # Generate the 3-block C-grid (default) and write the .geo / .vtm files:
    python examples/06-c-grid/visualize_naca0012.py --generate

    # Larger C-grid with finer wake clustering:
    python examples/06-c-grid/visualize_naca0012.py --generate \\
        --ni-body 160 --ni-wake 40 --n-stations 8 \\
        --wake-length 25 --wake-growth 1.20 \\
        --nk 30 --ds 5e-4 --growth 1.15

    # Legacy single-block layout (Xall has trouble with the wake-cut overlap):
    python examples/06-c-grid/visualize_naca0012.py --generate --cgrid-blocks 1

    # Visualise an existing .geo (any topology):
    python examples/06-c-grid/visualize_naca0012.py naca0012_body.geo

C-grid layout (3-block default)
-------------------------------
Three hex blocks with SHARED seams at the upper-TE column, the
lower-TE column, and the wake cut.  The wake cut SHARED face has an
i-reversal AxisMap on the partner side, the manifestation of the
"j-counter flips direction" pathology that the 1-block layout's
self-overlap encoded silently.  Xall sees clean face-to-face
adjacency on all three seams and applies BCs correctly.

i-axes per block:

    upper:  upper-wake-far-tip -> upper-TE
    body:   upper-TE -> suction -> LE -> pressure -> lower-TE
    lower:  lower-TE -> lower-wake-far-tip

The j-axis is the span direction (one node per station); the k-axis
is the wall-normal extrusion.  See :mod:`shore.surface.cgrid` and
:class:`shore.mesh.CGridTopology` for the detection / slicing /
wake-construction / SHARED-seam details.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# NACA 0012 STL builder
# ---------------------------------------------------------------------------


def _make_naca0012_stl(
    path: Path,
    *,
    n_section_pts: int = 80,
    n_stations: int = 8,
    chord: float = 1.0,
    span: float = 0.5,
) -> None:
    """Write a NACA 0012 wing STL with sharp TE and centroid-fan end caps.

    Uses :func:`shore.profiles.naca` for the analytic 2D section and
    :func:`shore.shapes.extrude_section` for the 3D extrusion.
    """
    from shore.profiles import naca
    from shore.shapes import extrude_section

    section = naca("0012", n_pts=n_section_pts // 2) * chord
    mesh = extrude_section(
        section,
        span=span,
        n_stations=n_stations,
        end_caps=True,
    )
    mesh.export(str(path))


# ---------------------------------------------------------------------------
# Mesh generation
# ---------------------------------------------------------------------------


def generate_cgrid(
    stl_path: Path,
    geo_stem: Path,
    *,
    ni_body: int,
    ni_wake: int,
    n_stations: int,
    wake_length: float,
    wake_growth: float,
    span_axis: str,
    nk: int,
    ds: float,
    growth: float,
    smooth: float,
    smooth_iters: int,
    cgrid_blocks: int = 3,
) -> list[Path]:
    """Build the STL (if missing), march the C-grid, write the .geo files.

    Returns the list of ``.geo`` paths written (one per block).
    """
    from shore.mesh import Body, CGridTopology, Mesh, Spacing1D

    if not stl_path.exists():
        print(f"Generating {stl_path} ...")
        _make_naca0012_stl(stl_path, n_stations=max(8, n_stations + 2))
        print("  STL written.")

    print(
        f"Building C-grid (blocks={cgrid_blocks}): ni_body={ni_body} "
        f"ni_wake={ni_wake} n_stations={n_stations} wake_length={wake_length}"
    )
    body = Body.from_stl(stl_path, skip_hygiene=True)
    topo = CGridTopology()
    spacing_k = Spacing1D(law="geometric", ds=ds, growth=growth)
    blocks = topo.build(
        surface=None,
        nk=nk,
        mesh=body.mesh,
        spacing_k=spacing_k,
        ni_body=ni_body,
        ni_wake=ni_wake,
        n_stations=n_stations,
        wake_length=wake_length,
        wake_growth=wake_growth,
        span_axis=span_axis,
        blocks=cgrid_blocks,
    )
    m = Mesh(body=body, blocks=blocks, topology=topo)

    def _cb(k: int, jmin: float) -> None:
        print(f"  layer {k:3d}  jmin={jmin:.3e}", end="\r")

    m.march(
        ds=ds,
        growth=growth,
        smooth_strength=smooth,
        smooth_iters=smooth_iters,
        layer_callback=_cb,
    )
    print()

    geo_stem.parent.mkdir(parents=True, exist_ok=True)
    if cgrid_blocks == 1:
        from shore.io.geo import write_geo

        single = geo_stem.with_suffix(".geo")
        write_geo(single, m.blocks[0].nodes)
        print(f"  {single} written  shape={m.blocks[0].nodes.shape[:3]}")
        return [single]
    # Multi-block: <stem>_<label>.geo via Mesh.write_geo.
    m.write_geo(geo_stem)
    paths: list[Path] = []
    for blk in m.blocks:
        p = geo_stem.parent / f"{geo_stem.name}_{blk.label}.geo"
        print(f"  {p} written  shape={blk.nodes.shape[:3]}")
        paths.append(p)
    return paths


def export_to_vtk(geo_paths: list[Path], vtm_path: Path) -> Path:
    """Bundle one or more .geo files into a single ParaView .vtm."""
    from shore.io.geo import read_geo
    from shore.io.vtk import write_caps_vtm

    grids = [read_geo(p) for p in geo_paths]
    labels = [p.stem.split("_")[-1] if "_" in p.stem else "cgrid" for p in geo_paths]
    write_caps_vtm(vtm_path, grids, labels=labels, with_jacobian=True)
    print(f"  {vtm_path} written  ({len(grids)} block(s))")
    return vtm_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build / visualise a SHORE C-grid around NACA 0012.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--generate",
        action="store_true",
        help="Generate naca0012.stl + naca0012.geo",
    )
    parser.add_argument("--ni-body", type=int, default=120)
    parser.add_argument("--ni-wake", type=int, default=30)
    parser.add_argument("--n-stations", type=int, default=5)
    parser.add_argument("--wake-length", type=float, default=15.0)
    parser.add_argument("--wake-growth", type=float, default=1.15)
    parser.add_argument("--span-axis", default="z", choices=["x", "y", "z"])
    parser.add_argument("--nk", type=int, default=20)
    parser.add_argument("--ds", type=float, default=5e-4)
    parser.add_argument("--growth", type=float, default=1.15)
    parser.add_argument("--smooth", type=float, default=0.1)
    parser.add_argument("--smooth-iters", type=int, default=1)
    parser.add_argument(
        "--cgrid-blocks",
        type=int,
        default=3,
        choices=[1, 3],
        help="Block layout: 3 (default, Xall-correct) or 1 (legacy single-block).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("examples_out"),
        help="Output directory; all files land here.",
    )
    parser.add_argument(
        "input",
        type=Path,
        nargs="?",
        default=None,
        help=".geo file to visualise",
    )
    parser.add_argument("--no-view", action="store_true", help="Generate only, skip the viewer")

    args = parser.parse_args()
    out_dir: Path = args.output
    out_dir.mkdir(parents=True, exist_ok=True)
    stl_path = out_dir / "naca0012.stl"
    geo_stem = out_dir / "naca0012"
    vtm_path = out_dir / "naca0012.vtm"

    if args.generate:
        try:
            geo_paths = generate_cgrid(
                stl_path,
                geo_stem,
                ni_body=args.ni_body,
                ni_wake=args.ni_wake,
                n_stations=args.n_stations,
                wake_length=args.wake_length,
                wake_growth=args.wake_growth,
                span_axis=args.span_axis,
                nk=args.nk,
                ds=args.ds,
                growth=args.growth,
                smooth=args.smooth,
                smooth_iters=args.smooth_iters,
                cgrid_blocks=args.cgrid_blocks,
            )
        except Exception as exc:  # pragma: no cover - example script
            print(f"Error: {exc}", file=sys.stderr)
            sys.exit(1)

        export_to_vtk(geo_paths, vtm_path)

        if args.no_view:
            return
        # Default viewer target: first emitted .geo (body block in the
        # 3-block layout, or the sole .geo in the 1-block layout).
        target = geo_paths[0] if args.cgrid_blocks == 1 else geo_paths[1]
    else:
        target = args.input
        if target is None:
            parser.error("Pass --generate or a .geo file to visualise.")

    try:
        from shore.mesh import view_geo_file
    except ImportError:
        print(
            "PyVista not installed; install with `pip install shore[vis]`.",
            file=sys.stderr,
        )
        sys.exit(0)
    view_geo_file(target)


if __name__ == "__main__":
    main()
