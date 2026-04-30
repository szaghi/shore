"""Standalone example: H-grid primitives (rectangular duct + circular tube).

Builds two H-grid variants for internal-flow / SRM-style cases:

    examples_out/
        h_grid_rect.geo      # 1m-long rectangular duct, 0.1m x 0.2m cross-section
        h_grid_rect.vtm      # ParaView MultiBlock (1 block)
        h_grid_tube.geo      # 0.5m-long circular tube, radius 0.05m
        h_grid_tube.vtm      # ParaView MultiBlock (1 block)

Run from the project root after activating .venv:

    python examples/07-h-grid/visualize_h_grid.py --generate
    python examples/07-h-grid/visualize_h_grid.py --generate -o /tmp/out
    python examples/07-h-grid/visualize_h_grid.py --generate --no-view
    python examples/07-h-grid/visualize_h_grid.py examples_out/h_grid_rect.geo

Both variants are single-block primitives with default BCs that match
the typical internal-flow configuration: streamwise faces FREE
(inflow / outflow), cross-section walls WALL.  See
:mod:`shore.primitive.h_grid_rect` and
:mod:`shore.primitive.h_grid_tube` for the API; see
[`shore primitive h-grid-rect / h-grid-tube`](/reference/primitive)
for the CLI surface.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def generate_rect(out_path: Path) -> Path:
    """Build the rectangular-duct H-grid and write .geo + .vtm."""
    from shore.primitive import HGridRectMesh

    print(f"Building rectangular H-grid → {out_path}.geo")
    hg = HGridRectMesh.build(
        length=1.0,
        width=0.1,
        height=0.2,
        ni=80,
        nj=20,
        nk=24,
        streamwise_axis="x",
        # Cluster cross-section nodes near the duct walls (boundary layer).
        j_law="tanh2",
        j_beta=4.0,
        k_law="tanh2",
        k_beta=4.0,
    )
    hg.write_geo(out_path)
    hg.write_vtk(out_path)
    print(f"  shape={hg.blocks[0].nodes.shape}")
    return out_path.with_suffix(".geo")


def generate_tube(out_path: Path) -> Path:
    """Build the circular-tube H-grid and write .geo + .vtm."""
    from shore.primitive import HGridTubeMesh

    print(f"Building circular-bore H-grid → {out_path}.geo")
    ht = HGridTubeMesh.build(
        length=0.5,
        radius=0.05,
        n_axial=60,
        n_circ=32,
        n_radial=10,
        streamwise_axis="x",
        # Cluster radial nodes near the outer wall via a one-sided
        # tanh — the polar singular axis at k=0 is naturally clustered
        # by the geometry (all angular nodes collapse there).
        k_law="tanh",
        k_beta=3.0,
    )
    ht.write_geo(out_path)
    ht.write_vtk(out_path)
    print(f"  shape={ht.blocks[0].nodes.shape}")
    return out_path.with_suffix(".geo")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build / visualise SHORE H-grid primitives (rect + tube).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--generate",
        action="store_true",
        help="Generate h_grid_rect.geo + h_grid_tube.geo (and matching .vtm files).",
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
        help=".geo file to visualise.",
    )
    parser.add_argument(
        "--no-view",
        action="store_true",
        help="Generate only, skip the viewer.",
    )

    args = parser.parse_args()
    out_dir: Path = args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.generate:
        try:
            rect_path = generate_rect(out_dir / "h_grid_rect")
            generate_tube(out_dir / "h_grid_tube")
        except Exception as exc:  # pragma: no cover - example script
            print(f"Error: {exc}", file=sys.stderr)
            sys.exit(1)

        if args.no_view:
            return
        target = rect_path
    else:
        target = args.input
        if target is None:
            parser.error("Pass --generate or a .geo file to visualise.")

    try:
        from shore.mesh import view_geo_file
    except ImportError:
        print(
            "PyVista not installed; install with `pip install shore-mesh[vis]`.",
            file=sys.stderr,
        )
        sys.exit(0)
    view_geo_file(target)


if __name__ == "__main__":
    main()
