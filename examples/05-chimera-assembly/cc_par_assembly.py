"""Generate the full Xall solver bundle from the sphere-in-cylinder assembly.

This example is a direct continuation of ``assembly_sphere_in_cylinder.py``.
It re-builds the same two meshes, writes their adjacency sidecars, exports
VTK MultiBlock files for ParaView inspection, then produces every file the
Xall pipeline needs:

- ``cc.par`` — chimera descriptor consumed by the overset preprocessor.
- ``proc.input`` — runtime work-distribution descriptor consumed directly
  by Xall at MPI startup.

Each is written per-component (wall + background) and then **merged** into
a single ``assembly.{grd, cc.par, proc.input}`` triple for an Xall run
that consumes both Chimera components as one 11-block job.

The two meshes are independent Chimera components — they share a domain
in world-space but have no SHARED faces between them.  Each gets its own
files first; the merge step concatenates them with the right block /
patch / proc index shifts so a single Xall job can consume the assembly.

The near-body cubed-sphere uses per-block ``BlockMeta`` to assign different
``free_face_bc`` values to wall blocks vs. cap blocks, demonstrating that
``write_cc_par``'s ``block_metadata`` argument (and the CLI ``--meta``
flag) give full per-block control:

- ``sub{0..3}`` — equatorial sub-blocks, chimera handoff → BC ``40`` (offgen).
- ``cap_north / cap_south`` — polar cap blocks, same → BC ``40``.
- Background blocks — far-field extrapolation → BC ``-9``.

The ``proc.input`` step shows the load-balancing summary on stderr and
sets ``body=1`` on every wall block (so Xall's force / moment integration
treats them as one physical body); background blocks stay at ``body=0``.

Output
------

::

    examples_out/
      wall_sub{0..3}.geo
      wall_cap_north.geo / wall_cap_south.geo    # near-body .geo files
      wall.vtm                                   # near-body ParaView MultiBlock
      wall.grd                                   # near-body Xall binary grid
      background_{center,sub_e,sub_n,sub_w,sub_s}.geo
      background.vtm                             # background ParaView MultiBlock
      background.grd                             # background Xall binary grid
      wall.adjacency.json
      background.adjacency.json
      wall.cc.par                                # includes marker boxes from sphere.stl
      background.cc.par
      wall.proc.input                            # MPI rank assignments + body=1
      background.proc.input                      # MPI rank assignments
      sphere_boxes.vtm                           # standalone box visualisation
      assembly.grd                               # merged: background + wall
      assembly.cc.par                            # merged: background + wall
      assembly.proc.input                        # merged: background + wall
      # When --balance is passed the wall.* and background.* files
      # above are *regenerated* on a larger split-balanced block list,
      # and the planner's TOML is preserved alongside:
      wall.splits.toml                           # shore balance plan (wall)
      background.splits.toml                     # shore balance plan (background)

Run from the project root after activating .venv::

    python examples/05-chimera-assembly/cc_par_assembly.py
    python examples/05-chimera-assembly/cc_par_assembly.py -o my_run
    python examples/05-chimera-assembly/cc_par_assembly.py --np 16
    python examples/05-chimera-assembly/cc_par_assembly.py --np 16 --balance
    python examples/05-chimera-assembly/cc_par_assembly.py --view      # open interactive viewer
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Mesh parameters — match assembly_sphere_in_cylinder.py defaults exactly.
# ---------------------------------------------------------------------------

_BODY_RADIUS = 0.5
_ICOSPHERE_SUBDIVISIONS = 4

_NB_NI = 30
_NB_NJ = 40
_NB_NK = 13
_NB_THETA_CAP_DEG = 30.0
_NB_DS = 0.01
_NB_GROWTH = 1.1

_BG_R_OUT = 10.0
_BG_Z0 = -10.0
_BG_Z1 = 50.0
_BG_INNER_FRAC = 0.4
_BG_NI = 16
_BG_NC = 24
_BG_NK = 60


# ---------------------------------------------------------------------------
# Mesh builders  (identical to assembly_sphere_in_cylinder.py)
# ---------------------------------------------------------------------------


def make_sphere_stl(out_dir: Path) -> Path:
    """Write an icosphere to ``<out_dir>/sphere.stl``."""
    import trimesh

    out_dir.mkdir(parents=True, exist_ok=True)
    stl_path = out_dir / "sphere.stl"
    sphere = trimesh.creation.icosphere(subdivisions=_ICOSPHERE_SUBDIVISIONS, radius=_BODY_RADIUS)
    sphere.export(str(stl_path))
    return stl_path


def build_near_body(stl_path: Path, out_stem: Path) -> list:
    """Build the 6-block cubed-sphere near-body mesh and write .geo files."""
    import numpy as np

    from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

    spacing_k = Spacing1D(law="geometric", ds=_NB_DS, growth=_NB_GROWTH)
    print(f"Building near-body cubed-sphere ({_NB_NI}x{_NB_NJ}, nk={_NB_NK}) ...")
    m = Mesh.from_stl(
        stl_path,
        topology=CubedSphereTopology(),
        ni=_NB_NI,
        nj=_NB_NJ,
        nk=_NB_NK,
        theta_cap=float(np.deg2rad(_NB_THETA_CAP_DEG)),
        spacing_k=spacing_k,
    )
    m.march(ds=_NB_DS, growth=_NB_GROWTH, smooth_strength=0.2, smooth_iters=2)
    m.write_geo(out_stem)
    labels = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]
    for lbl in labels:
        print(f"  wrote {out_stem.parent / (out_stem.name + '_' + lbl + '.geo')}")
    return m.blocks


def build_background(out_stem: Path) -> list:
    """Build the 5-block flat-caps cylinder background mesh and write .geo files."""
    from shore.primitive import FlatCapsCylinder

    print(f"Building background flat-caps cylinder (ni={_BG_NI}, nc={_BG_NC}, nk={_BG_NK}) ...")
    fc = FlatCapsCylinder.build(
        r_out=_BG_R_OUT,
        z0=_BG_Z0,
        z1=_BG_Z1,
        ni=_BG_NI,
        nc=_BG_NC,
        nk=_BG_NK,
        inner_frac=_BG_INNER_FRAC,
        axis="z",
    )
    fc.write_geo(out_stem)
    for blk in fc.blocks:
        print(f"  wrote {out_stem.parent / (out_stem.name + '_' + blk.label + '.geo')}")
    return fc.blocks


# ---------------------------------------------------------------------------
# VTK export
# ---------------------------------------------------------------------------


def write_vtm(
    blocks: list,
    vtm_path: Path,
    prefix: str,
    boxes: list | None = None,
) -> Path | None:
    """Bundle *blocks* (optionally + *boxes*) into a ParaView MultiBlock .vtm file."""
    try:
        from shore.io.geo import read_geo
        from shore.io.vtk import write_caps_vtm
    except ImportError as exc:
        print(f"  Warning: VTK export skipped ({exc}).", file=sys.stderr)
        return None

    geo_dir = vtm_path.parent
    stem = vtm_path.stem  # e.g. "wall"
    grids = []
    labels = []
    for blk in blocks:
        p = geo_dir / f"{stem}_{blk.label}.geo"
        grids.append(read_geo(p))
        labels.append(blk.label)
    write_caps_vtm(vtm_path, grids, with_jacobian=True, labels=labels, boxes=boxes)
    preview = ", ".join(labels[:3]) + ("..." if len(labels) > 3 else "")
    n_boxes = len(boxes) if boxes else 0
    box_suffix = f" + {n_boxes} box(es)" if n_boxes else ""
    print(f"  wrote {vtm_path}  ({len(grids)} blocks: {preview}{box_suffix})")
    return vtm_path


# ---------------------------------------------------------------------------
# Binary .grd writer
# ---------------------------------------------------------------------------


def write_grd_file(blocks: list, out_dir: Path, stem: str) -> Path:
    """Pack all blocks into one Xall binary ``.grd`` file, reading from disk.

    Resolves each block's ``.geo`` path as ``<out_dir>/<stem>_<label>.geo``
    and delegates to :func:`shore.io.grd.write_grd`.  Standalone: usable
    after ``shore split`` writes ``.geo`` files without in-memory blocks.

    Parameters
    ----------
    blocks :
        List of HexBlock objects (labels used to locate the .geo files).
    out_dir :
        Directory containing the ``.geo`` files.
    stem :
        Shared stem prefix (e.g. ``"wall"``).
    """
    from shore.io.geo import read_geo
    from shore.io.grd import write_grd

    geo_paths = [out_dir / f"{stem}_{blk.label}.geo" for blk in blocks]
    out = write_grd(out_dir / f"{stem}.grd", geo_paths)
    total_nodes = sum(
        g.shape[0] * g.shape[1] * g.shape[2] for g in (read_geo(p) for p in geo_paths)
    )
    print(f"  wrote {out}  ({len(geo_paths)} blocks, {total_nodes:,} nodes total)")
    return out


# ---------------------------------------------------------------------------
# Adjacency + cc.par writers
# ---------------------------------------------------------------------------


def write_adjacency(blocks: list, path: Path) -> Path:
    """Serialise the block adjacency graph to a JSON sidecar."""
    from shore.io.adjacency import write_adjacency_json

    p = write_adjacency_json(blocks, path)
    print(f"  wrote {p}")
    return p


def write_cc_par_file(
    adj_path: Path,
    ccpar_path: Path,
    *,
    grd_filename: str,
    block_metadata: dict,
    boxes: list | None = None,
) -> Path:
    """Convert an adjacency JSON sidecar to a ``cc.par`` file.

    Parameters
    ----------
    adj_path :
        Path to the ``.adjacency.json`` sidecar.
    ccpar_path :
        Destination ``cc.par`` path.
    grd_filename :
        Name of the ``.grd`` file that Xall pairs with this ``cc.par``
        (written into header line 1 verbatim).
    block_metadata :
        ``{label: BlockMeta}`` dict giving per-block level / group / priority /
        free_face_bc.  Any block not present gets ``BlockMeta()`` defaults.
    boxes :
        Optional list of :class:`shore.io.boxes.Box` instances written
        as the cc.par boxes section.  Default: no boxes.
    """
    from shore.io.adjacency import read_adjacency_json
    from shore.io.cc_par import write_cc_par

    payload = read_adjacency_json(adj_path)
    out = write_cc_par(
        payload,
        ccpar_path,
        grd_filename=grd_filename,
        block_metadata=block_metadata,
        boxes=boxes,
    )
    n_blocks = len(payload["blocks"])
    n_boxes = len(boxes) if boxes else 0
    print(f"  wrote {out}  ({n_blocks} blocks, {n_boxes} boxes)")
    return out


# ---------------------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------------------


def print_cc_par_summary(ccpar_path: Path) -> None:
    """Print the header and patch-count line from a cc.par file."""
    lines = ccpar_path.read_text().splitlines()
    print(f"\n--- {ccpar_path.name} (first 25 lines) ---")
    for ln in lines[:25]:
        print(f"  {ln}")
    if len(lines) > 25:
        print(f"  ... ({len(lines) - 25} more lines)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build sphere-in-cylinder assembly meshes and convert each "
            "mesh to an Xall cc.par chimera pre-processor input file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("examples_out"),
        help="Output directory; all files land here.",
    )
    parser.add_argument(
        "--view",
        action="store_true",
        default=False,
        help="Open an interactive PyVista viewer at the end.",
    )
    parser.add_argument(
        "--screenshot",
        type=Path,
        default=None,
        help="Save a PNG screenshot instead of opening a window (implies --view).",
    )
    parser.add_argument(
        "--np",
        type=int,
        default=4,
        help=(
            "Number of MPI ranks for proc.input load balancing.  "
            "Default 4 (the wall and background each have enough blocks "
            "for a non-trivial 4-way distribution).  Use 1 for serial."
        ),
    )
    parser.add_argument(
        "--balance",
        action="store_true",
        default=False,
        help=(
            "After the initial proc.input pass, run shore balance "
            "to plan splits that bring the imbalance below tolerance, "
            "apply them with shore split, then regenerate .grd, "
            "adjacency, cc.par, and proc.input on the larger block "
            "list.  Demonstrates the full plan-then-execute cycle."
        ),
    )
    parser.add_argument(
        "--balance-tolerance",
        type=float,
        default=0.05,
        help=(
            "Target imbalance for shore balance, as a fraction "
            "((max - min) / mean).  Default 0.05 = 5%%.  Only used "
            "when --balance is set."
        ),
    )

    args = parser.parse_args()
    out_dir: Path = args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=== cc_par_assembly ===")
    print(f"output directory: {out_dir.resolve()}")

    # 1. STL
    stl_path = make_sphere_stl(out_dir)
    print(f"  wrote {stl_path}")

    # 2. Near-body cubed-sphere mesh
    nb_blocks = build_near_body(stl_path, out_dir / "wall")

    # 3. Background flat-caps cylinder mesh
    bg_blocks = build_background(out_dir / "background")

    # 4. Generate marker boxes from the body STL.
    #
    # The boxes flag cells inside the sphere body as "internal wall"
    # for the chimera hole-cutter — when the background mesh's cells
    # land inside one of these boxes they get marked as donors-only
    # rather than active flow cells.  The boxes therefore must:
    #
    #   - contain the STL surface (so every body-interior background
    #     cell gets caught), and
    #   - fit inside the wall mesh's near-mesh outer extent (so they
    #     don't spuriously hole-cut legitimate near-body cells).
    #
    # ``max_outward_gap`` is the crucial knob: max distance from any
    # box-surface point to the STL surface (outward only).  The wall
    # is ~0.214 thick here (ds=0.01, growth=1.1, nk=13), so 0.10
    # leaves a comfortable safety margin against the wall outer
    # extent at radius 0.5+0.214≈0.71.
    #
    # We use method='voxel-fill' (the default): voxelise the body's
    # AABB at side length 2*gap/sqrt(3), test each voxel-centre with
    # ``trimesh.contains``, emit each strictly-inside-body voxel as
    # an axis-aligned Box.  This guarantees body-volume coverage
    # (every body-interior cell-centre lies in some box) while
    # keeping every box vertex within ``max_outward_gap`` of the STL
    # surface.  Single-OBB methods cannot satisfy both constraints
    # for non-cubic bodies (a sphere's bounding cube corners sit at
    # r*(sqrt(3)-1) past the body surface).
    from shore.boxes import boxes_from_stl

    print("\nGenerating marker boxes from the body STL:")
    body_boxes = boxes_from_stl(
        stl_path,
        associated_block="sub0",
        max_outward_gap=0.10,
        method="voxel-fill",
    )
    print(f"  {len(body_boxes)} box(es) from method='voxel-fill' (max_outward_gap=0.10)")

    # 5. VTK MultiBlock exports.
    #
    # ``wall.vtm`` layers the marker boxes onto the cubed-sphere blocks
    # so ParaView shows the box's outward gap to the body and inward
    # gap to the near-mesh side-by-side.  ``sphere_boxes.vtm`` is a
    # standalone export of just the boxes, useful when comparing
    # against the merged ``assembly.vtm`` in a separate ParaView pane.
    print("\nWriting VTK MultiBlock files:")
    write_vtm(nb_blocks, out_dir / "wall.vtm", prefix="wall", boxes=body_boxes)
    write_vtm(bg_blocks, out_dir / "background.vtm", prefix="background")

    try:
        from shore.io.vtk import write_boxes_vtm

        sphere_boxes_vtm = write_boxes_vtm(body_boxes, out_dir / "sphere_boxes.vtm")
        print(f"  wrote {sphere_boxes_vtm}  (standalone, {len(body_boxes)} box(es))")
    except ImportError:
        # pyvista not installed — write_vtm already warned for this
        # invocation, so we silently skip the standalone too.
        pass

    # 6. Binary .grd files for Xall
    print("\nWriting binary .grd files:")
    write_grd_file(nb_blocks, out_dir, "wall")
    write_grd_file(bg_blocks, out_dir, "background")

    # 7. Adjacency JSON sidecars
    print("\nWriting adjacency sidecars:")
    nb_adj = write_adjacency(nb_blocks, out_dir / "wall.adjacency.json")
    bg_adj = write_adjacency(bg_blocks, out_dir / "background.adjacency.json")

    # 8. cc.par files — one per independent Chimera component.
    #
    # Near-body: per-block BlockMeta assigns chimera cell-centre BC (40) to
    # every sub-block and cap; this is also the write_cc_par default, but we
    # spell it out explicitly to show the block_metadata interface.
    #
    # Background: all blocks get far-field extrapolation BC (-9).
    from shore.io.cc_par import BlockMeta

    nb_meta = {
        "sub0": BlockMeta(comment="sub0", free_face_bc=40),
        "sub1": BlockMeta(comment="sub1", free_face_bc=40),
        "sub2": BlockMeta(comment="sub2", free_face_bc=40),
        "sub3": BlockMeta(comment="sub3", free_face_bc=40),
        "cap_north": BlockMeta(level=1, comment="cap_north", free_face_bc=40),
        "cap_south": BlockMeta(level=1, comment="cap_south", free_face_bc=40),
    }
    bg_meta = {blk.label: BlockMeta(comment=blk.label, free_face_bc=-9) for blk in bg_blocks}

    print("\nWriting cc.par files:")
    nb_ccpar = write_cc_par_file(
        nb_adj,
        out_dir / "wall.cc.par",
        grd_filename="wall.grd",
        block_metadata=nb_meta,
        boxes=body_boxes,
    )
    bg_ccpar = write_cc_par_file(
        bg_adj,
        out_dir / "background.cc.par",
        grd_filename="background.grd",
        block_metadata=bg_meta,
    )

    # 9. Per-component proc.input — Xall's runtime work-distribution descriptor.
    #
    # cc.par is consumed by the overset preprocessor (which produces the runtime
    # artefacts the solver reads); proc.input is consumed directly by Xall at MPI
    # startup to learn block→rank assignments and per-block group / body indices.
    #
    # Block ordering and weights come from each component's .grd; without
    # explicit assignments, ranks are distributed by greedy weight-balanced
    # descent across range(np).  body=1 marks every wall block as belonging
    # to the same physical body for force / moment integration; the
    # background blocks stay at body=0 (no rigid-body motion).
    from shore.io.proc_input import write_proc_input

    print(f"\nWriting proc.input files (np={args.np}):")
    nb_proc_path, nb_proc_summary = write_proc_input(
        out_dir / "wall.grd",
        out_dir / "wall.proc.input",
        np=args.np,
        body_assignment=dict.fromkeys(range(1, len(nb_blocks) + 1), 1),
    )
    print(f"  wrote {nb_proc_path}")
    if args.np > 1:
        for line in nb_proc_summary.format().splitlines():
            print(f"    {line}")

    bg_proc_path, bg_proc_summary = write_proc_input(
        out_dir / "background.grd",
        out_dir / "background.proc.input",
        np=args.np,
    )
    print(f"  wrote {bg_proc_path}")
    if args.np > 1:
        for line in bg_proc_summary.format().splitlines():
            print(f"    {line}")

    # 10. Optional load-balance pass.
    #
    # The initial proc.input above does its best with the natural
    # 6-block / 5-block layout.  At --np=4 the wall imbalance is ~30%
    # and the background ~70% — both well above the 5% tolerance shore
    # uses by default — because the cubed sphere has 6 nearly-equal
    # blocks (poor 4-way fit) and the flat-caps cylinder has one large
    # central block that dwarfs its 4 lateral sub-blocks.
    #
    # `shore balance` plans splits that bring imbalance under tolerance;
    # `shore split` applies the plan; we then regenerate everything
    # downstream of the block list (.grd, adjacency, cc.par, proc.input).
    # Plan-then-execute: every step is reviewable on its own.
    if args.balance:
        from shore.balance import plan_balance

        print(f"\n=== Load-balance plan (tolerance={args.balance_tolerance:.1%}) ===")
        print("  Plan-only mode: this writes <stem>.splits.toml for each")
        print("  component but does *not* apply the splits.  Apply with:")
        print("    shore split -c <stem>.splits.toml <stem>_*.geo")
        print("  then regenerate .grd / cc.par / proc.input on the larger")
        print("  block list.  Plan-then-execute keeps every step reviewable")
        print("  and lets you tune the tolerance before mutating files.")
        print()
        print("  Note: applying the plan to the wall + background of *this*")
        print("  particular assembly hits two known shore.split v1 limitations")
        print("  (cap-equator DIRICHLET partner-pointer rewiring; cross-axis")
        print("  SHARED seams in the flat-caps cylinder).  Both are tracked")
        print("  on the roadmap.  The TOMLs themselves are correct and the")
        print("  planner has all the information needed to act once the")
        print("  splitter catches up.")

        for stem in ("wall", "background"):
            print(f"\n--- {stem} ---")
            plan = plan_balance(
                out_dir / f"{stem}.grd",
                np=args.np,
                tolerance=args.balance_tolerance,
                max_iterations=200,
                adjacency_path=out_dir / f"{stem}.adjacency.json",
            )
            print(
                f"  plan: {plan.initial_summary.imbalance * 100:.1f}% "
                f"\u2192 {plan.final_summary.imbalance * 100:.1f}% "
                f"after {plan.iterations} iteration(s) "
                f"({len(plan.splits)} cut(s) across "
                f"{len({s[0] for s in plan.splits})} block(s)) "
                f"[{'converged' if plan.converged else 'capped'}]"
            )
            if not plan.splits:
                print(f"  {stem}: already balanced; no plan written")
                continue
            plan_path = out_dir / f"{stem}.splits.toml"
            plan_path.write_text(plan.to_toml())
            print(f"  wrote {plan_path}")

    # 11. Merge into a single Xall multi-mesh job.
    #
    # The wall (6 blocks) and background (5 blocks) are independent
    # Chimera components — they overlap in world space but share no
    # SHARED faces.  Xall consumes them as one .grd + one cc.par with
    # 11 blocks; the merge utilities concatenate the two outputs and
    # shift cc.par block / family indices to keep the merged file
    # internally consistent.
    #
    # Order: background first, wall second.  Either is valid; the
    # convention matters because Xall's overset solver uses block
    # ordering to break level / priority ties.
    from shore.io.cc_par import merge_cc_par
    from shore.io.grd import merge_grd
    from shore.io.proc_input import merge_proc_input

    print("\nMerging into one assembly grid + cc.par + proc.input:")
    assembly_grd = merge_grd(
        [out_dir / "background.grd", out_dir / "wall.grd"],
        out_dir / "assembly.grd",
    )
    print(f"  wrote {assembly_grd}  (5 + 6 = 11 blocks)")
    assembly_ccpar = merge_cc_par(
        [out_dir / "background.cc.par", out_dir / "wall.cc.par"],
        out_dir / "assembly.cc.par",
        grd_filename="assembly.grd",
    )
    print(f"  wrote {assembly_ccpar}  (block comments auto-prefixed background.* / wall.*)")
    # Pass grd_path=assembly_grd so the merged summary carries real
    # cell-count weights (the merger reads them from the merged .grd).
    # Without it the summary would show only per-rank block counts.
    assembly_proc, assembly_proc_summary = merge_proc_input(
        [out_dir / "background.proc.input", out_dir / "wall.proc.input"],
        out_dir / "assembly.proc.input",
        grd_path=assembly_grd,
    )
    print(f"  wrote {assembly_proc}")
    if assembly_proc_summary.np > 1:
        for line in assembly_proc_summary.format().splitlines():
            print(f"    {line}")

    # 12. Spot-check printout
    print_cc_par_summary(nb_ccpar)
    print_cc_par_summary(bg_ccpar)
    print_cc_par_summary(assembly_ccpar)

    # 13. Equivalent CLI commands (for reference)
    import json

    nb_meta_json = json.dumps({k: {"free_face_bc": v.free_face_bc} for k, v in nb_meta.items()})
    bg_meta_json = json.dumps({k: {"free_face_bc": v.free_face_bc} for k, v in bg_meta.items()})
    print("\nEquivalent CLI commands:")
    print(
        f"  shore boxes-from-stl {stl_path} --block sub0 --max-outward-gap 0.10 "
        f"-o {out_dir / 'sphere_boxes.json'}"
    )
    print(f"  shore boxes-export {out_dir / 'sphere_boxes.json'} -o {out_dir / 'sphere_boxes.vtm'}")
    print(
        f"  shore cc-par {nb_adj} {nb_ccpar} --grd wall.grd --meta '{nb_meta_json}' "
        f"--boxes {out_dir / 'sphere_boxes.json'}"
    )
    print(f"  shore cc-par {bg_adj} {bg_ccpar} --grd background.grd --meta '{bg_meta_json}'")
    nb_proc_meta_json = json.dumps({str(i): {"body": 1} for i in range(1, len(nb_blocks) + 1)})
    print(
        f"  shore proc-input {out_dir / 'wall.grd'} {nb_proc_path} "
        f"--np {args.np} --meta '{nb_proc_meta_json}'"
    )
    print(f"  shore proc-input {out_dir / 'background.grd'} {bg_proc_path} --np {args.np}")
    if args.balance:
        print("  # Re-balance: plan + apply + repack + re-emit cc.par/proc.input")
        for stem in ("wall", "background"):
            print(
                f"  shore balance {out_dir / f'{stem}.grd'} "
                f"-o {out_dir / f'{stem}.splits.toml'} --np {args.np} "
                f"--tolerance {args.balance_tolerance} "
                f"--adjacency {out_dir / f'{stem}.adjacency.json'}"
            )
            print(
                f"  shore split -c {out_dir / f'{stem}.splits.toml'} "
                f"{out_dir / f'{stem}_*.geo'}  # then re-run grd / cc-par / proc-input"
            )
    print(
        f"  shore grd-merge {out_dir / 'background.grd'} {out_dir / 'wall.grd'} -o {assembly_grd}"
    )
    print(
        f"  shore cc-par-merge {out_dir / 'background.cc.par'} {out_dir / 'wall.cc.par'} "
        f"-o {assembly_ccpar} --grd assembly.grd"
    )
    print(
        f"  shore proc-input-merge {out_dir / 'background.proc.input'} "
        f"{out_dir / 'wall.proc.input'} -o {assembly_proc} --grd {assembly_grd}"
    )

    if not args.view and args.screenshot is None:
        return

    # 14. Optional visualisation (one near-body cap block)
    try:
        from shore.io.geo import read_geo
        from shore.vis import save_screenshot, view_geo
    except ImportError as exc:
        print(
            f"\nVisualisation unavailable: {exc}.  Install with `pip install shore[vis]`.",
            file=sys.stderr,
        )
        return

    view_target = out_dir / "wall_cap_north.geo"
    grid = read_geo(view_target)
    if args.screenshot is not None:
        save_screenshot(grid, str(args.screenshot))
        print(f"  screenshot saved to {args.screenshot}")
    else:
        print(f"\nShowing {view_target.name}.")
        view_geo(grid, colorby="Jacobian", show_edges=True, opacity=0.6)


if __name__ == "__main__":
    main()
