"""End-to-end example: per-edge i/j spacing on a BoxMesh primitive.

Box, Annulus, and FlatCaps are SHORE's *parametric* primitives — they
build their nodes directly from analytic geometry, no STL projection
and no hyperbolic marching.  Of the three, **Box** is the natural place
to demonstrate per-edge i/j spacing because every axis is straight and
free of geometric coupling: the user can request a different spacing
law on any of the 4 parallel edges along any axis, and the result is a
transfinite-interpolation (TFI) blend between the four edge
distributions.

What this example demonstrates
------------------------------
Two boxes side-by-side:

1. **Tensor-product baseline** — uniform along all three axes.  Every
   parallel edge of every axis carries the same ``Spacing1D``;
   construction collapses to the legacy ``np.meshgrid`` tensor
   product.

2. **Per-edge i-pin** — the i-edge at ``(j=0, k=0)`` is overridden to
   use one-sided tanh clustering toward the wall (small first cell);
   the other 3 i-edges remain uniform.  TFI blends the four i-edge
   distributions across the (j, k) face, producing a smooth
   non-tensor-product mesh in which the bottom-front i-row is densely
   packed near ``vmin[0]`` while the top-back i-row stays uniform.

The print-out compares the first 3 i-coordinates at three corners of
the (j, k) face — the pinned corner, the opposite corner, and the
midpoint — so you can verify TFI is doing what it claims.

Why this matters
----------------
For an overset (Chimera) assembly, the *background* mesh is often
parametric.  A typical pattern is "uniform background everywhere
except near the body, where I want clustering on one side."
Pre-existing meshers force you to either accept uniform-everywhere or
restructure the block.  Per-edge spacing lets you express the
boundary-layer-style request as a single edge override, and the TFI
takes care of the rest.

The j and k directions accept the same kind of per-edge override; this
example exercises i for clarity but the mechanism is identical for the
other axes.

Note: per-edge i/j is **not** wired on the cubed-sphere or
FlatCapsCylinder primitives — see those primitives' docstrings for the
geometric constraints that prevent per-edge control.  See
``per_edge_k_spacing_cubed_sphere.py`` for the cubed-sphere k-only
case.

Run from the project root after activating .venv:

    # Generate, print the comparison, and open the viewer.
    python examples/03-per-edge-spacing/per_edge_ij_box.py

    # Generate only — skip the viewer (CI-friendly).
    python examples/03-per-edge-spacing/per_edge_ij_box.py --no-view

    # Custom output directory.
    python examples/03-per-edge-spacing/per_edge_ij_box.py -o my_run --no-view
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Geometry / discretisation defaults
# ---------------------------------------------------------------------------

_VMIN = (0.0, 0.0, 0.0)
_VMAX = (1.0, 1.0, 1.0)

_NI = 21
_NJ = 11
_NK = 11

_TANH_BETA = 4.0


# ---------------------------------------------------------------------------
# Build helpers
# ---------------------------------------------------------------------------


def build_baseline():
    """Build a uniform-everywhere BoxMesh — the tensor-product baseline."""
    from shore.primitive import BoxMesh

    return BoxMesh.build(vmin=_VMIN, vmax=_VMAX, ni=_NI, nj=_NJ, nk=_NK, label="box_baseline")


def build_per_edge_pinned():
    """Build a BoxMesh with one i-edge overridden to tanh clustering."""
    from shore.mesh.spacing import Spacing1D
    from shore.primitive import BoxMesh

    return BoxMesh.build(
        vmin=_VMIN,
        vmax=_VMAX,
        ni=_NI,
        nj=_NJ,
        nk=_NK,
        # Pin the i-axis edge at (j=0, k=0) to one-sided tanh, dense at
        # vmin[0]=0.  The other 3 parallel i-edges remain at the default
        # uniform i_law.  TFI blends the four edge distributions over the
        # (j, k) face.
        edge_pins={"j_lo_k_lo": Spacing1D(law="tanh", beta=_TANH_BETA)},
        label="box_pinned",
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def i_row(box, j: int, k: int) -> np.ndarray:
    """Return the (ni,) array of x-coordinates along the i-row at (j, k)."""
    return box.blocks[0].nodes[:, j, k, 0]


def print_comparison(baseline, pinned) -> None:
    print()
    print("First 5 x-coordinates along three i-rows of the box")
    print("-----------------------------------------------------")
    print(f"{'corner':>22}  {'baseline (uniform)':>22}  {'pinned (TFI blend)':>22}")
    rows = [
        ("j=0, k=0  (pinned)", 0, 0),
        ("j=mid, k=mid  (blend)", _NJ // 2, _NK // 2),
        ("j=last, k=last  (uni)", _NJ - 1, _NK - 1),
    ]
    for label, j, k in rows:
        x_base = i_row(baseline, j, k)[:5]
        x_pin = i_row(pinned, j, k)[:5]
        print(f"  {label:>20s}")
        for idx, (a, b) in enumerate(zip(x_base, x_pin, strict=True)):
            print(f"    i={idx}: {a:>20.6f}  {b:>20.6f}")
    print()
    print(
        "At (j=0, k=0) the pinned distribution is exactly tanh — first "
        "interval is far below 1/(ni-1)."
    )
    print(
        "At (j=last, k=last) all 4 weights collapse onto the uniform-i edges "
        "→ result is identical to the baseline."
    )
    print(
        "At (j=mid, k=mid) TFI blends 25% of the tanh distribution into "
        "75% uniform.  The first cell is intermediate."
    )
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build two BoxMeshes — a uniform tensor-product baseline and a "
            "TFI mesh with one i-edge overridden to tanh clustering — and "
            "print the first 5 x-coords at three corners of the (j, k) "
            "face to demonstrate the blend."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("examples_out"),
        help="Output directory; .geo + .vtm files land here.",
    )
    parser.add_argument(
        "--no-view",
        action="store_true",
        help="Skip the interactive viewer.",
    )
    parser.add_argument(
        "--screenshot",
        type=Path,
        default=None,
        help="Save a PNG screenshot instead of opening a window.",
    )
    args = parser.parse_args()
    out_dir: Path = args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=== per_edge_ij_box ===")
    print(f"output directory: {out_dir.resolve()}")
    print(f"box vmin={_VMIN}, vmax={_VMAX}, shape=({_NI}, {_NJ}, {_NK})")
    print()

    print("Run 1: tensor-product baseline (uniform everywhere).")
    baseline = build_baseline()
    print(
        "Run 2: per-edge override on i-edge 'j_lo_k_lo' "
        f"(tanh, beta={_TANH_BETA}); other i-edges stay uniform."
    )
    pinned = build_per_edge_pinned()

    print_comparison(baseline, pinned)

    print("Writing outputs:")
    baseline.write_geo(out_dir / "box_baseline")
    pinned.write_geo(out_dir / "box_pinned")
    print(f"  wrote {out_dir / 'box_baseline.geo'}")
    print(f"  wrote {out_dir / 'box_pinned.geo'}")
    try:
        baseline.write_vtk(out_dir / "box_baseline.vtm")
        pinned.write_vtk(out_dir / "box_pinned.vtm")
        print(f"  wrote {out_dir / 'box_baseline.vtm'}")
        print(f"  wrote {out_dir / 'box_pinned.vtm'}")
    except ImportError as exc:
        print(f"  (skipped .vtm: {exc})")

    if args.no_view and args.screenshot is None:
        return

    try:
        from shore.io.geo import read_geo
        from shore.vis import save_screenshot, view_geo
    except ImportError as exc:
        print(f"\nVisualisation unavailable: {exc}.  Install with `pip install shore[vis]`.")
        return

    view_target = out_dir / "box_pinned.geo"
    grid = read_geo(view_target)
    if args.screenshot is not None:
        save_screenshot(grid, str(args.screenshot))
        print(f"  screenshot saved to {args.screenshot}")
    else:
        print(f"\nShowing {view_target.name}; the .vtm bundles the same data for ParaView.")
        view_geo(grid, colorby="Jacobian", show_edges=True, opacity=0.6)


if __name__ == "__main__":
    main()
