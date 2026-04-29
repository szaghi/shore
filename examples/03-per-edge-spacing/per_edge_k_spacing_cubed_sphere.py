"""End-to-end example: per-edge k-spacing on a cubed-sphere mesh.

Of the three SHORE topologies that take a body-fitted hyperbolic march
— cubed sphere and O-grid — the **wall-normal (k) direction** is the
only axis where per-edge spacing has a real degree of freedom on a
cubed sphere.  i and j are constrained by the cap / equator gnomonic
construction (the C1 seam pinning is what gives the topology its
quality, and it's geometric, not user-configurable).  k, by contrast,
is the user's call: how thick is the first layer, how fast does the
mesh grow into the far field?

This example shows the per-edge data-model path for that one knob.
Instead of passing ``ds`` / ``growth`` globally to ``Mesh.march()``,
you replace the four k-edges' ``spacing`` field on every block with a
shared ``Spacing1D`` instance, and the marcher picks it up.

How the override is detected
----------------------------
``HexBlock`` constructs every k-edge with the same default
``Spacing1D`` instance.  An edge whose ``spacing`` is no longer that
instance has been replaced by the caller.  When all 4 k-edges of every
block carry the *same* user instance, the marcher uses that schedule;
mixed states fall back to the call-site ``(ds, growth)`` and emit a
``UserWarning``.

Why all 6 blocks must agree
---------------------------
The cubed sphere has *k-coupled* seams: cap / equator DIRICHLET pins
at the i-boundary and sub-sub SHARED views at the j-boundary all live
at the same k-index, so the layer-by-layer node positions across all 6
blocks must remain coincident.  A different k-schedule on one block
would push its layer position out of registration.  ``Mesh._resolve_k_steps``
enforces this: it reads block 0's k-edges as canonical and warns if any
other block disagrees.

What the script does
--------------------
Builds the *same* sphere mesh twice — once via the legacy
``(ds, growth)`` path and once via the per-edge override — and prints
the first three layer thicknesses side-by-side so you can verify the
override landed.

For per-edge i/j on **primitives** (Box / Annulus), see
``examples/03-per-edge-spacing/per_edge_ij_box.py``.  For per-edge k on **O-grid**, see
``examples/03-per-edge-spacing/per_edge_k_spacing_ogrid.py``.

Run from the project root after activating .venv:

    # Generate, print the comparison, and open the viewer.
    python examples/03-per-edge-spacing/per_edge_k_spacing_cubed_sphere.py

    # Generate only — skip the viewer (CI-friendly).
    python examples/03-per-edge-spacing/per_edge_k_spacing_cubed_sphere.py --no-view

    # Custom output directory.
    python examples/03-per-edge-spacing/per_edge_k_spacing_cubed_sphere.py -o my_run --no-view
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Defaults — small mesh that runs fast and stays inside the C1 cap-equator
# pinning window.  See ``docs/topologies/cubed-sphere.md`` for what those
# parameters mean.
# ---------------------------------------------------------------------------

_BODY_RADIUS = 0.5
_ICOSPHERE_SUBDIVISIONS = 3

_NI = 30
_NJ = 40
_NK = 12
_THETA_CAP_DEG = 30.0

# Two distinct schedules so the printout makes the override unambiguous.
_LEGACY_DS = 0.05
_LEGACY_GROWTH = 1.10

_OVERRIDE_DS = 0.01
_OVERRIDE_GROWTH = 1.20


# ---------------------------------------------------------------------------
# Mesh construction helpers
# ---------------------------------------------------------------------------


def make_sphere_stl(out_dir: Path) -> Path:
    """Write an in-memory icosphere to ``<out_dir>/sphere.stl``."""
    import trimesh

    out_dir.mkdir(parents=True, exist_ok=True)
    stl_path = out_dir / "sphere.stl"
    sphere = trimesh.creation.icosphere(subdivisions=_ICOSPHERE_SUBDIVISIONS, radius=_BODY_RADIUS)
    sphere.export(str(stl_path))
    return stl_path


def build_mesh(stl_path: Path):
    """Build the unsplit 6-block cubed-sphere mesh (no marching yet)."""
    from shore.mesh import CubedSphereTopology, Mesh

    return Mesh.from_stl(
        stl_path,
        topology=CubedSphereTopology(),
        ni=_NI,
        nj=_NJ,
        nk=_NK,
        theta_cap=float(np.deg2rad(_THETA_CAP_DEG)),
    )


def apply_per_edge_override(mesh, ds: float, growth: float) -> None:
    """Replace every block's 4 k-edges with one shared ``Spacing1D``.

    The resolver looks for *identity*: a single ``Spacing1D`` instance
    stamped onto all 4 k-edges of every block tells the marcher
    "this is the user's per-block k-schedule".  Different instances on
    different edges (even with identical values) trigger the warning
    path and a fallback to ``(ds, growth)``.
    """
    from shore.mesh.spacing import Spacing1D

    override = Spacing1D(law="geometric", ds=ds, growth=growth)

    k_edges = ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi")
    for blk in mesh.blocks:
        for name in k_edges:
            blk.edge(name).spacing = override


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def first_layer_thicknesses(mesh, n: int = 3) -> list[float]:
    """Return the first *n* radial step lengths at block 0's (i=0, j=0) corner."""
    blk = mesh.blocks[0]
    out: list[float] = []
    for k in range(n):
        d = blk.nodes[0, 0, k + 1, :] - blk.nodes[0, 0, k, :]
        out.append(float(np.linalg.norm(d)))
    return out


def print_comparison(legacy: list[float], override: list[float]) -> None:
    print()
    print("First-layer thicknesses (block 0, corner i=0, j=0)")
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


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a cubed-sphere mesh twice — once via Mesh.march(ds, growth) "
            "and once via the per-edge k-spacing override — and print the "
            "first three layer thicknesses side-by-side."
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

    print("=== per_edge_k_spacing ===")
    print(f"output directory: {out_dir.resolve()}")
    print(f"mesh: ni={_NI}, nj={_NJ}, nk={_NK}, theta_cap={_THETA_CAP_DEG} deg")
    print()

    stl_path = make_sphere_stl(out_dir)
    print(f"  wrote {stl_path}")

    # --- Run 1: legacy (ds, growth) path ---
    print()
    print("Run 1: legacy Mesh.march(ds, growth) — no per-edge override.")
    mesh_legacy = build_mesh(stl_path)
    mesh_legacy.march(ds=_LEGACY_DS, growth=_LEGACY_GROWTH, smooth_strength=0.2, smooth_iters=2)
    legacy_thicknesses = first_layer_thicknesses(mesh_legacy)

    # --- Run 2: per-edge override path ---
    print("Run 2: per-edge k-spacing override (Mesh.march ds/growth ignored).")
    mesh_override = build_mesh(stl_path)
    apply_per_edge_override(mesh_override, ds=_OVERRIDE_DS, growth=_OVERRIDE_GROWTH)
    # Pass the *legacy* ds/growth here on purpose.  If the resolver works,
    # they should be ignored in favour of the per-edge schedule.
    mesh_override.march(ds=_LEGACY_DS, growth=_LEGACY_GROWTH, smooth_strength=0.2, smooth_iters=2)
    override_thicknesses = first_layer_thicknesses(mesh_override)

    print_comparison(legacy_thicknesses, override_thicknesses)

    # --- Outputs ---
    print("Writing outputs:")
    legacy_stem = out_dir / "sphere_legacy"
    override_stem = out_dir / "sphere_override"
    mesh_legacy.write_geo(legacy_stem)
    mesh_override.write_geo(override_stem)
    print(f"  wrote {legacy_stem}_*.geo  (6 files)")
    print(f"  wrote {override_stem}_*.geo  (6 files)")

    try:
        mesh_override.write_vtk(out_dir / "sphere_override.vtm")
        print(f"  wrote {out_dir / 'sphere_override.vtm'}")
    except ImportError as exc:
        print(f"  (skipped .vtm: {exc})")

    if args.no_view and args.screenshot is None:
        return

    # --- Visualise the overridden mesh (the one with the small ds) ---
    try:
        from shore.io.geo import read_geo
        from shore.vis import save_screenshot, view_geo
    except ImportError as exc:
        print(f"\nVisualisation unavailable: {exc}.  Install with `pip install shore[vis]`.")
        return

    view_target = out_dir / "sphere_override_sub0.geo"
    grid = read_geo(view_target)
    if args.screenshot is not None:
        save_screenshot(grid, str(args.screenshot))
        print(f"  screenshot saved to {args.screenshot}")
    else:
        print(f"\nShowing {view_target.name}; the .vtm bundles all 6 blocks for ParaView.")
        view_geo(grid, colorby="Jacobian", show_edges=True, opacity=0.6)


if __name__ == "__main__":
    main()
