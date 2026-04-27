"""End-to-end example: split a cubed-sphere mesh for parallel decomposition.

The cubed-sphere topology produces 6 blocks at the natural topological
scale (4 equatorial sub-blocks + 2 polar caps).  CFD solvers and
parallel decomposers often want *more, smaller* blocks — for load
balancing across MPI ranks, GPU streaming-multiprocessor occupancy,
halo-exchange granularity, or cache fit.  ``shore.split`` is the
post-processing pass that subdivides any list of HexBlocks while
recording every internal and external SHARED adjacency in a JSON
sidecar.

This example demonstrates the canonical "split for parallel" workflow:

1. Build a small cubed-sphere mesh on an in-memory icosphere.
2. Split each of the 4 equatorial sub-blocks along the wall-normal
   (k) axis at the same vertex index.  Splitting all 4 together is
   *required* — the j-periodic ring shares each sub_q's j_lo with
   sub_(q-1)'s j_hi, so splitting only one along k would create an
   irreconcilable mismatch on the j-seams.
3. Write 10 ``.geo`` files (8 sub-block chunks + 2 unchanged caps),
   one ``assembly.vtm`` MultiBlock for ParaView, and one
   ``.adjacency.json`` describing every face's BC and SHARED partner.
4. Print a tabular adjacency summary so you can verify the equatorial
   ring is intact across the splits.

Splitting a single primitive (box, annulus, flat-caps) is the simpler
case — each primitive's blocks have at most a few SHARED faces and
the cross-block constraint never kicks in.  See ``shore split --help``
or ``docs/reference/split.md`` for the per-primitive recipe.

Run from the project root after activating .venv:

    # Generate everything; open the viewer at the end.
    python examples/split_cubed_sphere.py

    # Generate only — skip the viewer (CI-friendly).
    python examples/split_cubed_sphere.py --no-view

    # Custom output directory.
    python examples/split_cubed_sphere.py -o my_split --no-view

    # Drive the splits from a TOML config (CLI flags still override).
    python examples/split_cubed_sphere.py \\
        --config examples/config_split_cubed_sphere.toml --no-view
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Defaults — tuned for a fast end-to-end smoke test that still produces a
# visually interesting split.  See ``docs/guide/assembly.md`` for the
# rationale behind these ni / nj / nk choices.
# ---------------------------------------------------------------------------

_BODY_RADIUS = 0.5
_ICOSPHERE_SUBDIVISIONS = 3

_NI = 30
_NJ = 40
_NK = 12
_THETA_CAP_DEG = 30.0
_DS = 0.01
_GROWTH = 1.1

# Split each sub-block along k at this single vertex index.  Mid-block
# is the most common request for parallel decomposition.
_K_CUT = _NK // 2  # = 6


# ---------------------------------------------------------------------------
# Mesh construction
# ---------------------------------------------------------------------------


def make_sphere_stl(out_dir: Path) -> Path:
    """Write an in-memory icosphere to ``<out_dir>/sphere.stl``."""
    import trimesh

    out_dir.mkdir(parents=True, exist_ok=True)
    stl_path = out_dir / "sphere.stl"
    sphere = trimesh.creation.icosphere(subdivisions=_ICOSPHERE_SUBDIVISIONS, radius=_BODY_RADIUS)
    sphere.export(str(stl_path))
    return stl_path


def build_cubed_sphere(stl_path: Path):
    """Build the unsplit 6-block cubed-sphere mesh."""
    import numpy as np

    from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

    spacing_k = Spacing1D(law="geometric", ds=_DS, growth=_GROWTH)
    print(f"Building cubed-sphere mesh ({_NI}x{_NJ}, nk={_NK}) ...")
    print(f"  body radius {_BODY_RADIUS} (diameter {2 * _BODY_RADIUS})")
    print(f"  theta_cap = {_THETA_CAP_DEG} deg")
    print(f"  wall-normal: geometric ds={_DS} growth={_GROWTH}")
    m = Mesh.from_stl(
        stl_path,
        topology=CubedSphereTopology(),
        ni=_NI,
        nj=_NJ,
        nk=_NK,
        theta_cap=float(np.deg2rad(_THETA_CAP_DEG)),
        spacing_k=spacing_k,
    )
    m.march(ds=_DS, growth=_GROWTH, smooth_strength=0.2, smooth_iters=2)
    return m


# ---------------------------------------------------------------------------
# Split application
# ---------------------------------------------------------------------------


def apply_splits(blocks: list, splits: list | None = None) -> list:
    """Apply paired k-splits across the 4 equatorial sub-blocks.

    When *splits* is None, builds the canonical "split each sub at
    k=_K_CUT" directive list.  When passed (typically loaded from a
    TOML config), uses it verbatim.
    """
    from shore.split import BlockSplit, split_blocks

    if splits is None:
        splits = [BlockSplit(label=f"sub{q}", axis="k", at=[_K_CUT]) for q in range(4)]

    print(f"Applying {len(splits)} split directive(s):")
    for sp in splits:
        print(f"  {sp.label} along {sp.axis} at {list(sp.at)}")

    new_blocks = split_blocks(blocks, splits)
    print(f"  → {len(blocks)} blocks before, {len(new_blocks)} blocks after")
    return new_blocks


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


def write_geo_files(blocks: list, out_stem: Path) -> list[Path]:
    """Write each block to ``<stem>_<label>.geo`` and return the paths."""
    from shore.io.geo import write_geo

    paths: list[Path] = []
    for blk in blocks:
        p = out_stem.parent / f"{out_stem.name}_{blk.label}.geo"
        write_geo(p, blk.nodes)
        paths.append(p)
    print(f"  wrote {len(paths)} .geo files under {out_stem.parent}/")
    return paths


def write_adjacency(blocks: list, out_path: Path) -> Path:
    """Write the adjacency JSON sidecar."""
    from shore.io.adjacency import write_adjacency_json

    p = write_adjacency_json(blocks, out_path)
    print(f"  wrote {p}")
    return p


def write_assembly_vtm(geo_paths: list[Path], vtm_path: Path) -> Path | None:
    """Bundle every .geo into one ParaView MultiBlock .vtm."""
    try:
        from shore.io.geo import read_geo
        from shore.io.vtk import write_caps_vtm
    except ImportError as exc:
        print(f"Warning: VTK export skipped ({exc}).", file=sys.stderr)
        return None

    blocks = [read_geo(p) for p in geo_paths]
    labels = [p.stem.split("_", 1)[1] if "_" in p.stem else p.stem for p in geo_paths]
    write_caps_vtm(vtm_path, blocks, with_jacobian=True, labels=labels)
    print(f"  wrote {vtm_path}  ({len(blocks)} blocks)")
    return vtm_path


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def print_adjacency_summary(adj_path: Path) -> None:
    """Print the per-block adjacency graph from the JSON sidecar.

    Uses :func:`shore.io.adjacency.summarise_adjacency` so the output
    format stays aligned with whatever the library produces.
    """
    from shore.io.adjacency import read_adjacency_json, summarise_adjacency

    payload = read_adjacency_json(adj_path)
    print()
    print(f"Adjacency graph ({len(payload['blocks'])} blocks):")
    print(summarise_adjacency(payload))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _load_split_directives_from_toml(path: Path) -> list:
    """Read a ``shore split`` TOML config into BlockSplit directives."""
    from shore.config import load_split_config
    from shore.split import BlockSplit

    cfg = load_split_config(path)
    return [BlockSplit(label=d.label, axis=d.axis, at=d.at) for d in cfg.splits]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a cubed-sphere mesh, split each sub-block along k for "
            "parallel decomposition, and write .geo + .vtm + .adjacency.json."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("examples_out"),
        help="Output directory; all .geo, .vtm, .adjacency.json files land here.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional TOML split config (kind = 'split').  When omitted, the "
        "script applies its built-in defaults: split each of sub0..sub3 "
        "along k at vertex index nk//2.",
    )
    parser.add_argument(
        "--no-vtm",
        action="store_true",
        help="Skip the bundled .vtm assembly export.",
    )
    parser.add_argument(
        "--no-view",
        action="store_true",
        help="Generate only; skip the interactive viewer.",
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

    print("=== split_cubed_sphere ===")
    print(f"output directory: {out_dir.resolve()}")

    # 1. STL + cubed-sphere mesh
    stl_path = make_sphere_stl(out_dir)
    print(f"  wrote {stl_path}")
    m = build_cubed_sphere(stl_path)

    # 2. Splits
    if args.config is not None:
        splits = _load_split_directives_from_toml(args.config)
    else:
        splits = None
    new_blocks = apply_splits(m.blocks, splits)

    # 3. Outputs
    print("Writing outputs:")
    out_stem = out_dir / "sphere"
    geo_paths = write_geo_files(new_blocks, out_stem)
    adj_path = write_adjacency(new_blocks, out_dir / "sphere.adjacency.json")
    if not args.no_vtm:
        write_assembly_vtm(geo_paths, out_dir / "sphere.vtm")

    # 4. Tabular summary
    print_adjacency_summary(adj_path)

    if args.no_view and args.screenshot is None:
        return

    # 5. Visualise (PyVista on one of the chunks; the .vtm has them all)
    try:
        from shore.io.geo import read_geo
        from shore.vis import save_screenshot, view_geo
    except ImportError as exc:
        print(f"\nVisualisation unavailable: {exc}.  Install with `pip install shore[vis]`.")
        return

    view_target = out_dir / "sphere_sub0__k0.geo"
    grid = read_geo(view_target)
    if args.screenshot is not None:
        save_screenshot(grid, str(args.screenshot))
        print(f"  screenshot saved to {args.screenshot}")
    else:
        print(f"\nShowing {view_target.name}; the .vtm bundles all 10 chunks for ParaView.")
        view_geo(grid, colorby="Jacobian", show_edges=True, opacity=0.6)


if __name__ == "__main__":
    main()
