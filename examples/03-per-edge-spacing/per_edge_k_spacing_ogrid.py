"""End-to-end example: per-edge k-spacing on a single-block O-grid.

The O-grid is the simpler peer of the 6-block cubed-sphere topology:
one ``(ni, nj, nk, 3)`` hexahedral block, j-periodic, with polar
singularities at ``i = 0`` and ``i = ni - 1``.  Per-edge spacing on
this topology is *only* meaningful on the wall-normal (k) direction:

- **i-direction**: all four parallel i-edges (``j_lo_k_lo`` etc.) lie
  at adjacent azimuths under the j-wrap, not at distinct transverse
  corners, so per-edge i has no real degree of freedom over the
  per-axis ``i_spacing`` argument.

- **j-direction**: j is periodic — the "first" and "last" cells of
  any j-edge are the same cell, so per-edge j is geometrically
  meaningless.

- **k-direction**: 4 parallel k-edges that pass through marching at
  the same k-index of a single block.  Replacing their ``Spacing1D``
  field with a shared user instance switches ``Mesh.march()`` from the
  call-site ``(ds, growth)`` to the user-supplied schedule — same
  identity-check mechanism as on the cubed sphere, but with only one
  block to coordinate.

This example mirrors ``per_edge_k_spacing_cubed_sphere.py`` but on the
O-grid: build the same sphere mesh twice (legacy ``(ds, growth)`` vs
per-edge override) and print the first three layer thicknesses
side-by-side.  The override is the value that lands in the actual
mesh; the call-site ``ds`` / ``growth`` are silently ignored once the
4 k-edges of the single O-grid block are user-overridden.

For per-edge i/j on parametric **primitives** (Box / Annulus), see
``per_edge_ij_box.py``.

Run from the project root after activating .venv:

    # Generate, print the comparison, and open the viewer.
    python examples/03-per-edge-spacing/per_edge_k_spacing_ogrid.py

    # Generate only — skip the viewer (CI-friendly).
    python examples/03-per-edge-spacing/per_edge_k_spacing_ogrid.py --no-view

    # Custom output directory.
    python examples/03-per-edge-spacing/per_edge_k_spacing_ogrid.py -o my_run --no-view
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Defaults — small mesh for a fast smoke test, well inside the regime where
# both schedules produce all-positive Jacobians.
# ---------------------------------------------------------------------------

_BODY_RADIUS = 0.5
_ICOSPHERE_SUBDIVISIONS = 3

_NI = 12
_NJ = 24
_NK = 8
_THETA_CAP_DEG = 15.0

_LEGACY_DS = 0.05
_LEGACY_GROWTH = 1.10

_OVERRIDE_DS = 0.01
_OVERRIDE_GROWTH = 1.20


def make_sphere_stl(out_dir: Path) -> Path:
    import trimesh

    out_dir.mkdir(parents=True, exist_ok=True)
    stl_path = out_dir / "sphere.stl"
    sphere = trimesh.creation.icosphere(subdivisions=_ICOSPHERE_SUBDIVISIONS, radius=_BODY_RADIUS)
    sphere.export(str(stl_path))
    return stl_path


def build_ogrid(stl_path: Path):
    """Build the unsplit single-block O-grid (no marching yet)."""
    from shore.mesh import Mesh, OGridTopology

    return Mesh.from_stl(
        stl_path,
        topology=OGridTopology(),
        ni=_NI,
        nj=_NJ,
        nk=_NK,
        theta_cap=float(np.deg2rad(_THETA_CAP_DEG)),
    )


def apply_per_edge_override(mesh, ds: float, growth: float) -> None:
    """Stamp one shared ``Spacing1D`` onto all 4 k-edges of the O-grid block."""
    from shore.mesh.spacing import Spacing1D

    override = Spacing1D(law="geometric", ds=ds, growth=growth)
    blk = mesh.blocks[0]
    for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
        blk.edge(name).spacing = override


def first_layer_thicknesses(mesh, n: int = 3) -> list[float]:
    """Radial step lengths near the south pole row (i=0, j=0)."""
    blk = mesh.blocks[0]
    out: list[float] = []
    for k in range(n):
        d = blk.nodes[0, 0, k + 1, :] - blk.nodes[0, 0, k, :]
        out.append(float(np.linalg.norm(d)))
    return out


def print_comparison(legacy: list[float], override: list[float]) -> None:
    print()
    print("First-layer thicknesses (block 0, i=0, j=0)")
    print("---------------------------------------------------")
    print(f"{'k':>3}  {'legacy (ds, growth)':>22}  {'per-edge override':>22}")
    for k, (a, b) in enumerate(zip(legacy, override, strict=True)):
        print(f"{k + 1:>3}  {a:>22.6f}  {b:>22.6f}")
    print()
    print(
        f"legacy schedule:    ds={_LEGACY_DS}, growth={_LEGACY_GROWTH} "
        f"  →  expected first ≈ {_LEGACY_DS}"
    )
    print(
        f"override schedule:  ds={_OVERRIDE_DS}, growth={_OVERRIDE_GROWTH}"
        f"  →  expected first ≈ {_OVERRIDE_DS}"
    )
    print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a single-block O-grid sphere mesh twice — once via "
            "Mesh.march(ds, growth) and once via per-edge k-spacing — and "
            "print the first three layer thicknesses side-by-side."
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

    print("=== per_edge_k_spacing_ogrid ===")
    print(f"output directory: {out_dir.resolve()}")
    print(f"O-grid: ni={_NI}, nj={_NJ}, nk={_NK}, theta_cap={_THETA_CAP_DEG} deg")
    print()

    stl_path = make_sphere_stl(out_dir)
    print(f"  wrote {stl_path}")

    print()
    print("Run 1: legacy Mesh.march(ds, growth) — no per-edge override.")
    mesh_legacy = build_ogrid(stl_path)
    mesh_legacy.march(ds=_LEGACY_DS, growth=_LEGACY_GROWTH, smooth_strength=0.2, smooth_iters=2)
    legacy_thicknesses = first_layer_thicknesses(mesh_legacy)

    print("Run 2: per-edge k-spacing override (Mesh.march ds/growth ignored).")
    mesh_override = build_ogrid(stl_path)
    apply_per_edge_override(mesh_override, ds=_OVERRIDE_DS, growth=_OVERRIDE_GROWTH)
    mesh_override.march(ds=_LEGACY_DS, growth=_LEGACY_GROWTH, smooth_strength=0.2, smooth_iters=2)
    override_thicknesses = first_layer_thicknesses(mesh_override)

    print_comparison(legacy_thicknesses, override_thicknesses)

    print("Writing outputs:")
    legacy_stem = out_dir / "ogrid_legacy"
    override_stem = out_dir / "ogrid_override"
    mesh_legacy.write_geo(legacy_stem)
    mesh_override.write_geo(override_stem)
    print(f"  wrote {legacy_stem}_*.geo  (1 file)")
    print(f"  wrote {override_stem}_*.geo  (1 file)")

    try:
        mesh_override.write_vtk(out_dir / "ogrid_override.vtm")
        print(f"  wrote {out_dir / 'ogrid_override.vtm'}")
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

    view_target = out_dir / "ogrid_override_ogrid.geo"
    grid = read_geo(view_target)
    if args.screenshot is not None:
        save_screenshot(grid, str(args.screenshot))
        print(f"  screenshot saved to {args.screenshot}")
    else:
        print(f"\nShowing {view_target.name}.")
        view_geo(grid, colorby="Jacobian", show_edges=True, opacity=0.6)


if __name__ == "__main__":
    main()
