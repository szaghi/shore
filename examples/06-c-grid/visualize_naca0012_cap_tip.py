"""Standalone example: C-grid + tip caps around a NACA 0012 wing (ADR 0059).

Closes both spanwise tips of a C-grid wing with overset caps on the tip
faces. When ``cap_tips != "none"``, the C-grid auto-promotes its sharp-TE
3-block layout to the 4-block layout (with a ``wake_base`` block bridging
the TE thickness gap), the body's TE is rounded over the full span by
``cap_te_thickness_frac`` (default 0.2% chord), and the first and last
stations move onto the tips. The cap blocks come from the same single
``CGridTopology.build(...)`` call.

Outputs (default ``examples_out/``):

    naca0012.dat                 # 2D Selig section (intermediate)
    naca0012.stl                 # 3D wing STL (sharp TE, with end caps)
    wing_upper.geo               # C-grid upper-wake block
    wing_body.geo                # C-grid body wrap (rounded TE)
    wing_lower.geo               # C-grid lower-wake block
    wing_wake_base.geo           # C-grid wake-base block (TE bridge)
    wing_cap_jlo_core.geo        # jlo cap blocks (5 of them)
    wing_cap_jlo_upper.geo
    wing_cap_jlo_le.geo
    wing_cap_jlo_lower.geo
    wing_cap_jlo_te.geo
    wing_cap_jhi_*.geo           # jhi cap blocks (5)
    wing_capped.vtm              # ParaView MultiBlock (14 blocks)

Run from the project root after activating .venv:

    python examples/06-c-grid/visualize_naca0012_cap_tip.py --generate
    python examples/06-c-grid/visualize_naca0012_cap_tip.py --generate -o /tmp/run

The bash sibling (``visualize_naca0012_cap_tip.sh``) does the same via
``shore mesh --cap-tips both`` and is the canonical CLI reference.

Cap-block layout per tip
------------------------
Five hex blocks per tip cover the tip face, in an O-H layout:

- ``cap_<tip>_core``  — inside the section, inset from both sides.
- ``cap_<tip>_upper`` — from the core to the i = 0 side's wall.
- ``cap_<tip>_le``    — from the core round the leading edge.
- ``cap_<tip>_lower`` — from the core to the other side's wall.
- ``cap_<tip>_te``    — from the core across the TE base.

All cap blocks have ``k_lo = WALL`` (the tip face) and grow outwards; the
faces between them are SHARED, the others FREE (overset fringes). The
cap's edge is the body's tip-station wall nodes, bit for bit.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# NACA 0012 STL + conformal-cap C-grid in one call
# ---------------------------------------------------------------------------


def _make_naca0012_stl(
    path: Path,
    *,
    n_section_pts: int = 80,
    n_stations: int = 8,
    chord: float = 1.0,
    span: float = 0.5,
) -> None:
    from shore.profiles import naca
    from shore.shapes import extrude_section

    section = naca("0012", n_pts=n_section_pts // 2) * chord
    mesh = extrude_section(section, span=span, n_stations=n_stations, end_caps=True)
    mesh.export(str(path))


def _build_capped_cgrid(
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
    cap_tips: str,
    cap_le_frac: float,
    cap_flare_deg: float,
    cap_te_thickness_frac: float,
    cap_tip_inset: float,
    cap_n_cap_j: int,
) -> list[Path]:
    from shore.mesh import Body, CGridTopology, Mesh, Spacing1D

    if not stl_path.exists():
        print(f"Generating {stl_path} ...")
        _make_naca0012_stl(stl_path, n_stations=max(8, n_stations + 2))

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
        blocks=3,
        cap_tips=cap_tips,
        cap_le_frac=cap_le_frac,
        cap_flare_deg=cap_flare_deg,
        cap_te_thickness_frac=cap_te_thickness_frac,
        cap_tip_inset=cap_tip_inset,
        cap_n_cap_j=cap_n_cap_j,
    )
    m = Mesh(body=body, blocks=blocks, topology=topo)
    m.march(ds=ds, growth=growth, smooth_strength=smooth, smooth_iters=smooth_iters)
    m.write_geo(geo_stem)
    return [geo_stem.parent / f"{geo_stem.name}_{blk.label}.geo" for blk in blocks]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a NACA 0012 C-grid + tip caps and write outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    gen = parser.add_argument_group("generation")
    gen.add_argument("--generate", action="store_true", help="Build the mesh + caps.")
    gen.add_argument("-o", "--output", type=Path, default=Path("examples_out"),
                     help="Output directory (default: examples_out/).")

    gen.add_argument("--ni-body", type=int, default=120)
    gen.add_argument("--ni-wake", type=int, default=30)
    gen.add_argument("--n-stations", type=int, default=5)
    gen.add_argument("--wake-length", type=float, default=15.0)
    gen.add_argument("--wake-growth", type=float, default=1.15)
    gen.add_argument("--span-axis", default="z", choices=["x", "y", "z"])
    gen.add_argument("--nk", type=int, default=20)
    gen.add_argument("--ds", type=float, default=5e-4)
    gen.add_argument("--growth", type=float, default=1.15)
    gen.add_argument("--smooth", type=float, default=0.1)
    gen.add_argument("--smooth-iters", type=int, default=1)

    gen.add_argument("--cap-tips", default="both", choices=["none", "jlo", "jhi", "both"])
    gen.add_argument("--cap-le-frac", type=float, default=0.05)
    gen.add_argument("--cap-flare-deg", type=float, default=30.0)
    gen.add_argument("--cap-te-thickness-frac", type=float, default=0.002)
    gen.add_argument("--cap-tip-inset", type=float, default=0.0)
    gen.add_argument("--cap-n-cap-j", type=int, default=12)

    parser.add_argument("input", type=Path, nargs="?", default=None,
                        help="A .geo file to visualise (skip --generate).")
    parser.add_argument("--no-view", action="store_true", help="Generate only, skip viewer.")

    args = parser.parse_args()
    out_dir: Path = args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    geo_paths: list[Path] = []
    if args.generate:
        try:
            from shore.io.geo import read_geo
        except Exception as exc:
            print(f"Error: {exc}", file=sys.stderr)
            sys.exit(1)

        geo_paths = _build_capped_cgrid(
            stl_path=out_dir / "naca0012.stl",
            geo_stem=out_dir / "wing",
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
            cap_tips=args.cap_tips,
            cap_le_frac=args.cap_le_frac,
            cap_flare_deg=args.cap_flare_deg,
            cap_te_thickness_frac=args.cap_te_thickness_frac,
            cap_tip_inset=args.cap_tip_inset,
            cap_n_cap_j=args.cap_n_cap_j,
        )
        for p in geo_paths:
            print(f"  wrote {p}")

        # ParaView .vtm with all blocks.
        try:
            from shore.io.vtk import write_caps_vtm

            arrays = [read_geo(p) for p in geo_paths]
            labels = [p.stem.replace("wing_", "") for p in geo_paths]
            write_caps_vtm(
                out_dir / "wing_capped.vtm",
                arrays,
                labels=labels,
                with_jacobian=True,
            )
            print(f"  wrote {out_dir / 'wing_capped.vtm'}")
        except Exception as exc:
            print(f"  (skipped .vtm: {exc})", file=sys.stderr)

        if args.no_view:
            return

    target = args.input if args.input is not None else (
        geo_paths[0] if geo_paths else None
    )
    if target is None:
        parser.error("Pass --generate or a .geo file to visualise.")

    try:
        from shore.mesh import view_geo_file
    except ImportError:
        print("PyVista not installed; install with `pip install shore-mesh[vis]`.", file=sys.stderr)
        sys.exit(0)
    view_geo_file(target)


if __name__ == "__main__":
    main()
