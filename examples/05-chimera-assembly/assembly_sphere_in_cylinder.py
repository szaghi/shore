"""First-tentative Chimera assembly: sphere body inside a cylindrical domain.

This example exercises the topology + primitive layer end-to-end by
producing two structured meshes that **share a domain in world
coordinates** but are generated and written independently:

- **Near-body**: a 6-block cubed-sphere O-grid wrapping a unit-radius
  sphere (diameter 1.0, centred at the origin), marched outward to a
  total wall-normal extent of ~0.2 (so the outer cubed-sphere shell
  sits at radius ~0.7).

- **Background**: a 5-block flat-caps cylinder spanning the full
  domain ``[-10, 10]^2`` in the (x, y) plane and ``z ∈ [-10, 50]``,
  axis along z.  No central hole — the sub-blocks plus the central
  square block fill the disk conformally, so the background does
  not introduce a Chimera fringe at its inner boundary.

The two meshes overlap in the cylindrical region ``r < ~0.7`` near
the origin; resolving that overlap (hole-cutting, fringe
identification, donor search) is the **solver's** responsibility, not
this example's.  SHORE just generates the two meshes side by side.

Output
------

By default everything lands under ``examples_out/`` (cwd-relative):

::

    examples_out/
      fuselage_sub0.geo  ..  fuselage_sub3.geo
      fuselage_cap_north.geo
      fuselage_cap_south.geo                 # 6 cubed-sphere blocks
      background_center.geo
      background_sub_e/n/w/s.geo             # 5 flat-caps blocks
      assembly.vtm                           # 11-block ParaView MultiBlock

The ``.vtm`` bundles all 11 blocks (6 near-body + 5 background) into
one file you can open in ParaView to verify the assembly geometry.

Run from the project root after activating .venv:

    # Generate everything; open the viewer on the assembly.
    python examples/05-chimera-assembly/assembly_sphere_in_cylinder.py

    # Generate only — skip the interactive viewer.
    python examples/05-chimera-assembly/assembly_sphere_in_cylinder.py --no-view

    # Custom output directory.
    python examples/05-chimera-assembly/assembly_sphere_in_cylinder.py -o my_run

    # Drive the two pieces from TOML configs (CLI flags still override):
    python examples/05-chimera-assembly/assembly_sphere_in_cylinder.py \\
        --near-body-config examples/05-chimera-assembly/config_assembly_near_body.toml \\
        --background-config examples/05-chimera-assembly/config_assembly_background.toml
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Defaults — match the scenario described in the docstring above.
# ---------------------------------------------------------------------------

# Body
_BODY_RADIUS = 0.5
_ICOSPHERE_SUBDIVISIONS = 4

# Near-body cubed-sphere mesh.  ``ni`` / ``nj`` chosen so the cap radial
# step (theta_cap/(nc-1)) and the equator step ((pi-2*theta_cap)/(ni-1))
# land inside the per-j C1 pinning window (ratio in [0.5, 2.0]).
#
# Wall-normal: geometric law (uniform aspect ratio across layers, no
# extreme first-cell clustering).  ``ds=0.01`` gives a first cell ~0.3x
# the surface cell scale, ``growth=1.1`` over ``nk=13`` layers reaches
# a total thickness of ~0.21 — outer cubed-sphere shell at r ≈ 0.7.
_NB_NI = 30
_NB_NJ = 40
_NB_NK = 13
_NB_THETA_CAP_DEG = 30.0
_NB_DS = 0.01
_NB_GROWTH = 1.1

# Background flat-caps cylinder
_BG_R_OUT = 10.0
_BG_Z0 = -10.0
_BG_Z1 = 50.0
_BG_INNER_FRAC = 0.4  # central square corners at r = 4.0
_BG_NI = 16  # radial sub-block resolution
_BG_NC = 24  # central square / arc resolution (shared)
_BG_NK = 60  # axial — long extent needs more layers


# ---------------------------------------------------------------------------
# Body / near-body mesh
# ---------------------------------------------------------------------------


def make_sphere_stl(out_dir: Path, radius: float = _BODY_RADIUS) -> Path:
    """Write an in-memory icosphere to ``<out_dir>/sphere.stl``.

    Parameters
    ----------
    out_dir :
        Directory to write into; created if missing.
    radius :
        Sphere radius (default 0.5 → diameter 1.0).
    """
    import trimesh

    out_dir.mkdir(parents=True, exist_ok=True)
    stl_path = out_dir / "sphere.stl"
    sphere = trimesh.creation.icosphere(subdivisions=_ICOSPHERE_SUBDIVISIONS, radius=radius)
    sphere.export(str(stl_path))
    return stl_path


def build_near_body(
    stl_path: Path,
    out_stem: Path,
    cfg: dict | None = None,
) -> list[Path]:
    """Build the cubed-sphere near-body mesh.

    Returns the list of 6 ``.geo`` paths.  When *cfg* is supplied
    (loaded from a TOML primitive-style mesh config) its values
    override the hardcoded defaults; CLI overrides still beat both via
    the calling script.
    """
    import numpy as np

    from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

    if cfg is None:
        cfg = {}

    ni = cfg.get("ni", _NB_NI)
    nj = cfg.get("nj", _NB_NJ)
    nk = cfg.get("nk", _NB_NK)
    theta_cap_deg = cfg.get("theta_cap_deg", _NB_THETA_CAP_DEG)
    ds = cfg.get("ds", _NB_DS)
    growth = cfg.get("growth", _NB_GROWTH)

    spacing_k = Spacing1D(law="geometric", ds=ds, growth=growth)
    total_thick = float(spacing_k.build_steps(nk).sum())
    print(f"Building near-body cubed-sphere ({ni}x{nj}, nk={nk}) ...")
    print(f"  theta_cap = {theta_cap_deg} deg")
    print(
        f"  wall-normal: geometric ds={ds} growth={growth} "
        f"(total thickness ~{total_thick:.3f}, outer at r ~ {_BODY_RADIUS + total_thick:.3f})"
    )

    m = Mesh.from_stl(
        stl_path,
        topology=CubedSphereTopology(),
        ni=ni,
        nj=nj,
        nk=nk,
        theta_cap=float(np.deg2rad(theta_cap_deg)),
        spacing_k=spacing_k,
    )
    m.march(ds=ds, growth=growth, smooth_strength=0.2, smooth_iters=2)

    m.write_geo(out_stem)
    labels = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]
    geo_paths = [out_stem.parent / f"{out_stem.name}_{lbl}.geo" for lbl in labels]
    for p in geo_paths:
        print(f"  wrote {p}")
    return geo_paths


# ---------------------------------------------------------------------------
# Background mesh
# ---------------------------------------------------------------------------


def build_background(
    out_stem: Path,
    cfg: dict | None = None,
) -> list[Path]:
    """Build the flat-caps cylinder background mesh.

    Returns the list of 5 ``.geo`` paths.
    """
    from shore.primitive import FlatCapsCylinder

    if cfg is None:
        cfg = {}

    r_out = cfg.get("r_out", _BG_R_OUT)
    z0 = cfg.get("z0", _BG_Z0)
    z1 = cfg.get("z1", _BG_Z1)
    inner_frac = cfg.get("inner_frac", _BG_INNER_FRAC)
    ni = cfg.get("ni", _BG_NI)
    nc = cfg.get("nc", _BG_NC)
    nk = cfg.get("nk", _BG_NK)

    print(f"Building background flat-caps cylinder (ni={ni}, nc={nc}, nk={nk}) ...")
    print(f"  r_out = {r_out},  z = [{z0}, {z1}],  inner_frac = {inner_frac}")

    fc = FlatCapsCylinder.build(
        r_out=r_out,
        z0=z0,
        z1=z1,
        ni=ni,
        nc=nc,
        nk=nk,
        inner_frac=inner_frac,
        axis="z",
    )
    fc.write_geo(out_stem)
    geo_paths = [out_stem.parent / f"{out_stem.name}_{blk.label}.geo" for blk in fc.blocks]
    for p in geo_paths:
        print(f"  wrote {p}")
    return geo_paths


# ---------------------------------------------------------------------------
# Assembly bundle
# ---------------------------------------------------------------------------


def write_assembly_vtm(
    near_body_paths: list[Path],
    background_paths: list[Path],
    vtm_path: Path,
) -> Path:
    """Bundle all 11 blocks into a single PyVista MultiBlock ``.vtm``.

    Block keys are derived from each ``.geo`` filename (stem) and
    prefixed with ``nb_`` / ``bg_`` so the two meshes are visually
    grouped in ParaView.
    """
    from shore.io.geo import read_geo
    from shore.io.vtk import write_caps_vtm

    blocks = []
    labels = []
    for p in near_body_paths:
        blocks.append(read_geo(p))
        labels.append(f"nb_{p.stem.split('_', 1)[1]}")  # drop the file stem prefix
    for p in background_paths:
        blocks.append(read_geo(p))
        labels.append(f"bg_{p.stem.split('_', 1)[1]}")

    write_caps_vtm(vtm_path, blocks, with_jacobian=True, labels=labels)
    print(f"  wrote {vtm_path}  ({len(blocks)} blocks: {', '.join(labels)})")
    return vtm_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _load_optional_config(path: Path | None, expected_kind: str) -> dict | None:
    """Read a TOML config and return a flat dict of relevant keys.

    *expected_kind* is ``"mesh"`` for the near-body cubed-sphere config
    or ``"primitive"`` for the background flat-caps config.  Returns
    None when *path* is None.
    """
    if path is None:
        return None
    if expected_kind == "mesh":
        from shore.config import load_config

        cfg = load_config(path)
        return {
            "ni": cfg.grid.ni,
            "nj": cfg.grid.nj,
            "nk": cfg.grid.nk,
            "theta_cap_deg": (
                None
                if cfg.grid.theta_cap is None
                else float(cfg.grid.theta_cap) * 180.0 / 3.141592653589793
            ),
            "k_thickness": cfg.volume.k_thickness,
            "k_beta": cfg.volume.k_beta,
        }
    if expected_kind == "primitive":
        from shore.config import load_primitive_config

        cfg = load_primitive_config(path)
        if cfg.geometry.type != "flat_caps":
            raise SystemExit(
                f"Error: background config geometry.type='{cfg.geometry.type}' "
                "(expected 'flat_caps')."
            )
        g = cfg.geometry
        return {
            "r_out": g.r_out if g.r_out is not None else _BG_R_OUT,
            "z0": g.z0 if g.z0 is not None else _BG_Z0,
            "z1": g.z1 if g.z1 is not None else _BG_Z1,
            "inner_frac": g.inner_frac if g.inner_frac is not None else _BG_INNER_FRAC,
            "ni": cfg.grid.ni,
            "nc": cfg.grid.nj,
            "nk": cfg.grid.nk,
        }
    raise ValueError(f"unknown expected_kind: {expected_kind}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=("First-tentative assembly: sphere body in a flat-caps cylinder background."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("examples_out"),
        help="Output directory; all .geo and .vtm files land here.",
    )
    parser.add_argument(
        "--near-body-config",
        type=Path,
        default=None,
        help="Optional TOML config for the near-body cubed-sphere mesh (``shore mesh`` schema).",
    )
    parser.add_argument(
        "--background-config",
        type=Path,
        default=None,
        help="Optional TOML config for the background flat-caps cylinder "
        "(``shore primitive`` schema with type='flat_caps').",
    )
    parser.add_argument(
        "--no-view",
        action="store_true",
        help="Generate only; skip the interactive viewer.",
    )
    parser.add_argument(
        "--no-vtm",
        action="store_true",
        help="Skip the bundled .vtm assembly export.",
    )
    parser.add_argument(
        "--screenshot",
        type=Path,
        default=None,
        help="Save a PNG screenshot instead of opening a window (implies --no-view).",
    )

    args = parser.parse_args()

    out_dir: Path = args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. STL ---
    print("=== assembly_sphere_in_cylinder ===")
    print(f"output directory: {out_dir.resolve()}")
    stl_path = make_sphere_stl(out_dir)
    print(f"  wrote {stl_path} (icosphere, radius {_BODY_RADIUS})")

    # --- 2. Near-body cubed-sphere ---
    nb_cfg = _load_optional_config(args.near_body_config, "mesh")
    nb_paths = build_near_body(stl_path, out_dir / "fuselage", cfg=nb_cfg)

    # --- 3. Background flat-caps cylinder ---
    bg_cfg = _load_optional_config(args.background_config, "primitive")
    bg_paths = build_background(out_dir / "background", cfg=bg_cfg)

    # --- 4. Bundled assembly .vtm ---
    if not args.no_vtm:
        try:
            write_assembly_vtm(nb_paths, bg_paths, out_dir / "assembly.vtm")
        except ImportError as exc:
            print(f"Warning: VTK export skipped ({exc}).", file=sys.stderr)

    if args.no_view and args.screenshot is None:
        return

    # --- 5. Visualise ---
    try:
        from shore.io.geo import read_geo
        from shore.vis import save_screenshot, view_geo
    except ImportError as exc:
        print(
            f"Visualisation unavailable: {exc}.  Install with `pip install shore[vis]`.",
        )
        return

    # Default view: a cap block from the near-body (it sits inside the
    # background and its smaller scale is the most visually informative
    # single block).  For the full assembly use the .vtm in ParaView.
    view_target = out_dir / "fuselage_cap_north.geo"
    grid = read_geo(view_target)
    if args.screenshot is not None:
        save_screenshot(grid, str(args.screenshot))
        print(f"  screenshot saved to {args.screenshot}")
    else:
        print(f"  showing {view_target.name}; the .vtm bundles all 11 blocks for ParaView.")
        view_geo(grid, colorby="Jacobian", show_edges=True, opacity=0.6)


if __name__ == "__main__":
    main()
