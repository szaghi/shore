"""Headless spike checks (ADR 0006, items 2 and 3) on an offscreen canvas.

1. Upload + first frame and orbit frame time for an N-node assembly. Offscreen
   ``draw()`` includes the GPU->CPU readback, so the frame times are an upper
   bound; ``app.py`` measures the on-screen rate.
2. The float32 trap: pick a wall cell of height ``ds_wall`` at a world offset
   of 1e3, in each origin mode, and check the returned (block, i, j, k).

Results go to ``results/<host>-<adapter>.json``.

Usage::

    python bench_headless.py --nodes 10000000
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import socket
import time
from pathlib import Path

import numpy as np
import pygfx as gfx
from interaction import final_distance, wall_cell, zoom_to
from rendercanvas.offscreen import RenderCanvas
from synthetic import make_assembly
from viewport import OriginMode, PygfxViewport

SIZE = (1280, 720)


def _adapter_summary() -> str:
    return gfx.renderers.wgpu.get_shared().adapter.summary


def bench_frames(n_nodes: int, n_frames: int) -> dict:
    """Upload, first-frame and orbit timings for one assembly."""
    t0 = time.perf_counter()
    blocks = make_assembly(n_nodes=n_nodes)
    t_build = time.perf_counter() - t0

    canvas = RenderCanvas(size=SIZE)
    vp = PygfxViewport(canvas, origin_mode="floating")
    t0 = time.perf_counter()
    for b in blocks:
        vp.add_block(b)
    vp.show_all()
    canvas.draw()
    t_first = time.perf_counter() - t0

    target = vp.target.copy()
    eye0 = np.asarray(vp.camera.local.position, dtype=np.float64) - target
    times = []
    for f in range(n_frames):
        ang = 2.0 * np.pi * f / n_frames
        c, s = np.cos(ang), np.sin(ang)
        rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        vp.camera.local.position = target + rot @ eye0
        vp.camera.look_at(target)
        t = time.perf_counter()
        canvas.draw()
        times.append(1e3 * (time.perf_counter() - t))
    a = np.asarray(times[2:])  # drop pipeline-compilation warm-up
    return {
        "n_nodes": sum(b.n_nodes for b in blocks),
        "n_blocks": len(blocks),
        "surface_nodes": int(sum(sum(g.size // 3 for g in _faces(b.nodes)) for b in blocks)),
        "build_s": round(t_build, 3),
        "upload_first_frame_s": round(t_first, 3),
        "frame_ms_mean": round(float(a.mean()), 2),
        "frame_ms_p95": round(float(np.percentile(a, 95)), 2),
        "frame_ms_max": round(float(a.max()), 2),
        "note": "offscreen draw() includes GPU->CPU readback: upper bound",
    }


def _faces(nodes: np.ndarray) -> list[np.ndarray]:
    return [
        nodes[0],
        nodes[-1],
        nodes[:, 0],
        nodes[:, -1],
        nodes[:, :, 0],
        nodes[:, :, -1],
    ]


def precision_trap(ds_wall: float, offset: float, k_cells: tuple[int, ...]) -> dict:
    """Pick i_lo-face cells of height ``ds_wall`` near the wall in each mode."""
    blocks = make_assembly(n_nodes=200_000, ds_wall=ds_wall, offset=(offset, 0.0, 0.0))
    block = blocks[len(blocks) // 2]
    nodes = block.nodes
    nj = nodes.shape[1]
    b = nj // 2
    out: dict[str, dict] = {}
    modes: tuple[OriginMode, ...] = ("absolute", "block", "floating")
    for mode in modes:
        canvas = RenderCanvas(size=SIZE)
        vp = PygfxViewport(canvas, origin_mode=mode)
        for blk in blocks:
            vp.add_block(blk)
        hits = []
        for k in k_cells:
            quad = np.array(
                [nodes[0, b, k], nodes[0, b + 1, k], nodes[0, b + 1, k + 1], nodes[0, b, k + 1]]
            )
            centre = quad.mean(axis=0)
            # look at the i_lo face head-on; 1e-6 cell ~ 20 px tall at this distance
            normal = np.cross(quad[1] - quad[0], quad[3] - quad[0])
            distance = ds_wall * SIZE[1] / (20.0 * 2.0 * np.tan(np.radians(22.5)))
            vp.look_at_point(centre, distance, view_dir=normal)
            canvas.draw()
            x, y = vp.project(centre)
            res = vp.pick(x, y)
            # coincident seam faces z-fight, so either block may answer: judge by
            # the float64 centre of the picked quad, not by the block name
            err = None if res is None else float(np.linalg.norm(res.xyz - centre))
            hits.append(
                {
                    "k": k,
                    "expected": [block.name, "i_lo", [0, b, k]],
                    "got": None if res is None else [res.block, res.face, list(res.cell)],
                    "centre_err": err,
                    "ok": err is not None and err < 0.5 * ds_wall,
                    "pick_ms": None if res is None else round(res.elapsed_ms, 2),
                }
            )
        out[mode] = {"all_ok": all(h["ok"] for h in hits), "picks": hits}
    return {"ds_wall": ds_wall, "offset": offset, "modes": out}


EPS32 = float(np.finfo(np.float32).eps)
SWEEP_MODES: tuple[OriginMode, ...] = ("floating", "tiled", "rebase")


def extent_sweep(
    radii: tuple[float, ...],
    ds_over_extent: tuple[float, ...],
    offset: float,
    aspect_ratio: float,
    modes: tuple[OriginMode, ...] = SWEEP_MODES,
) -> list[dict]:
    """Pick trap over block extent x wall spacing, per origin mode.

    The wall cells picked have a fixed tangential/normal aspect ratio
    (``aspect_ratio``, realistic boundary layers are 1e3-1e5): polar angle is
    clustered about each block's centre node. Without that, a coarse test
    assembly yields aspect ratios of ~1e9, and no float32 vertex scheme can
    represent a quad thinner than ~eps32 times its own width.

    ``res32`` is ``eps32 * max|p - block centroid|``, the vertex resolution of
    the single-origin-per-block (``floating``) scheme. The assembly uses
    ``n_phi=5`` so ``i_lo`` faces are not coordinate planes. 18 picks per row:
    3 blocks x the 2 cells adjacent to the cluster centre x k = 0..2.
    """
    rows = []
    for radius in radii:
        for ratio in ds_over_extent:
            # block extent ~ 2.85 * radius for this assembly (n_phi=5, thickness 2r)
            ds_wall = ratio * 2.85 * radius
            blocks = make_assembly(
                n_nodes=200_000,
                n_phi=5,
                n_theta=4,
                radius=radius,
                thickness=2.0 * radius,
                ds_wall=ds_wall,
                offset=(offset, 0.0, 0.0),
                aspect_ratio=aspect_ratio,
            )
            sample = (blocks[1], blocks[8], blocks[17])
            extent = max(float(np.ptp(b.nodes.reshape(-1, 3), axis=0).max()) for b in sample)
            targets = []
            res32 = 0.0
            for blk in sample:
                nodes = blk.nodes
                centroid = nodes.reshape(-1, 3).mean(axis=0)
                c = (nodes.shape[1] - 1) // 2
                for b in (c - 1, c):
                    for k in (0, 1, 2):
                        quad = np.array(
                            [
                                nodes[0, b, k],
                                nodes[0, b + 1, k],
                                nodes[0, b + 1, k + 1],
                                nodes[0, b, k + 1],
                            ]
                        )
                        res32 = max(res32, EPS32 * float(np.abs(quad - centroid).max()))
                        targets.append(quad)
            distance = ds_wall * SIZE[1] / (20.0 * 2.0 * np.tan(np.radians(22.5)))
            for mode in modes:
                canvas = RenderCanvas(size=SIZE)
                vp = PygfxViewport(canvas, origin_mode=mode)
                t0 = time.perf_counter()
                for blk in blocks:
                    vp.add_block(blk)
                setup_s = time.perf_counter() - t0
                n_ok = 0
                draw_ms, rebase_ms = [], []
                for quad in targets:
                    centre = quad.mean(axis=0)
                    normal = np.cross(quad[1] - quad[0], quad[3] - quad[0])
                    vp.look_at_point(centre, distance, view_dir=normal)
                    rebase_ms.append(vp.last_rebase_ms)
                    t = time.perf_counter()
                    canvas.draw()
                    draw_ms.append(1e3 * (time.perf_counter() - t))
                    res = vp.pick(*vp.project(centre))
                    err = None if res is None else float(np.linalg.norm(res.xyz - centre))
                    n_ok += int(err is not None and err < 0.5 * ds_wall)
                row = {
                    "mode": mode,
                    "radius": radius,
                    "block_extent": round(extent, 6),
                    "ds_wall": ds_wall,
                    "ds_over_extent": ds_wall / extent,
                    "aspect_ratio": aspect_ratio,
                    "ds_over_res32": ds_wall / res32,
                    "picks_ok": n_ok,
                    "picks": len(targets),
                    "meshes": vp.n_pieces,
                    "tiles_capped": vp.tiles_capped,
                    "setup_s": round(setup_s, 3),
                    "draw_ms_median": round(float(np.median(draw_ms[1:])), 1),
                    "rebase_ms_median": round(float(np.median(rebase_ms)), 2),
                }
                rows.append(row)
                cols = [
                    f"{mode:8s}",
                    f"extent {extent:7.3g}",
                    f"ds/extent {row['ds_over_extent']:7.1e}",
                    f"ds/res32 {row['ds_over_res32']:7.2g}",
                    f"picks {n_ok:2d}/{len(targets)}",
                    f"meshes {vp.n_pieces:6d}{'*' if vp.tiles_capped else ' '}",
                    f"setup {setup_s:6.2f} s",
                    f"draw {row['draw_ms_median']:6.1f} ms",
                    f"rebase {row['rebase_ms_median']:6.2f} ms",
                ]
                print("  ".join(cols), flush=True)
    return rows


def rebase_policy_check(offset: float, aspect_ratio: float) -> list[dict]:
    """Zoom like a user onto a wall cell; the trigger alone must keep it exact.

    Extent ~285, ``ds/extent = 1e-9``: the regime where one origin per block
    fails (10-12/18 in the sweep). Factor 100 (default) vs inf (trigger off).
    A zoom does not probe the factor bound: it keeps zooming after the
    trigger fires, so it always ends with the origin at the target. See
    ``drift_tolerance_check`` for the bound itself.
    """
    radius = 100.0
    ds_wall = 1e-9 * 2.85 * radius
    blocks = make_assembly(
        n_nodes=200_000,
        n_phi=5,
        n_theta=4,
        radius=radius,
        thickness=2.0 * radius,
        ds_wall=ds_wall,
        offset=(offset, 0.0, 0.0),
        aspect_ratio=aspect_ratio,
    )
    rows = []
    for factor in (100.0, float("inf")):
        canvas = RenderCanvas(size=SIZE)
        vp = PygfxViewport(canvas, origin_mode="rebase", rebase_factor=factor)
        for blk in blocks:
            vp.add_block(blk)
        vp.show_all()
        canvas.draw()
        n_ok = 0
        for blk in (blocks[1], blocks[8], blocks[17]):
            cell = wall_cell(blk)
            for _ in zoom_to(vp, cell, d_start=3.0 * radius, n_steps=60):
                canvas.draw()
            res = vp.pick(*vp.project(cell.centre))
            err = None if res is None else float(np.linalg.norm(res.xyz - cell.centre))
            n_ok += int(err is not None and err < 0.5 * cell.height)
        rows.append(
            {
                "rebase_factor": factor,
                "picks_ok": n_ok,
                "picks": 3,
                "rebases": vp.rebase_count,
                "last_rebase_ms": round(vp.last_rebase_ms, 2),
            }
        )
        print(
            f"rebase trigger factor {factor:>5}: picks {n_ok}/3, rebases {vp.rebase_count}",
            flush=True,
        )
    return rows


def drift_tolerance_check(offset: float, aspect_ratio: float) -> list[dict]:
    """Pick a 20 px wall cell with the origin ``D`` view widths away (no trigger).

    Predicted vertex error is ``eps32 * D * V`` world units, i.e.
    ``eps32 * D * 720`` px at 720 px view height. Picks tolerate ~half a cell
    (10 px); the visual bound used for the trigger factor is 0.1 px.
    """
    radius = 100.0
    blocks = make_assembly(
        n_nodes=200_000,
        n_phi=5,
        n_theta=4,
        radius=radius,
        thickness=2.0 * radius,
        ds_wall=1e-9 * 2.85 * radius,
        offset=(offset, 0.0, 0.0),
        aspect_ratio=aspect_ratio,
    )
    sample = [wall_cell(b) for b in (blocks[1], blocks[8], blocks[17])]
    rows = []
    for views in (1e2, 1e3, 1e4, 1e5, 1e6):
        canvas = RenderCanvas(size=SIZE)
        vp = PygfxViewport(canvas, origin_mode="rebase", rebase_factor=float("inf"))
        for blk in blocks:
            vp.add_block(blk)
        n_ok = 0
        for cell in sample:
            vp.look_at_point(cell.centre, final_distance(cell), view_dir=cell.normal)
            tangent = np.cross(cell.normal, [0.0, 0.0, 1.0])
            tangent /= np.linalg.norm(tangent)
            vp.set_origin(cell.centre + views * vp.view_width() * tangent)
            canvas.draw()
            res = vp.pick(*vp.project(cell.centre))
            err = None if res is None else float(np.linalg.norm(res.xyz - cell.centre))
            n_ok += int(err is not None and err < 0.5 * cell.height)
        predicted_px = EPS32 * views * SIZE[1]
        rows.append({"drift_views": views, "predicted_px": predicted_px, "picks_ok": n_ok})
        print(
            f"origin drift {views:7.0e} views: predicted error {predicted_px:8.3f} px,"
            f" picks {n_ok}/3",
            flush=True,
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--nodes", type=int, default=10_000_000)
    ap.add_argument("--frames", type=int, default=60)
    ap.add_argument("--ds-wall", type=float, default=1e-6)
    ap.add_argument("--offset", type=float, default=1e3)
    ap.add_argument("--skip-frames", action="store_true", help="only run the precision trap")
    ap.add_argument("--skip-sweep", action="store_true", help="skip the extent/ds_wall sweep")
    ap.add_argument("--aspect-ratio", type=float, default=1e4, help="wall-cell AR in the sweep")
    args = ap.parse_args()

    adapter = _adapter_summary()
    print(f"adapter: {adapter}")
    result = {
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "adapter": adapter,
        "versions": {"pygfx": gfx.__version__},
    }
    if not args.skip_frames:
        result["frames"] = bench_frames(args.nodes, args.frames)
        print(json.dumps(result["frames"], indent=2))
    result["precision"] = precision_trap(args.ds_wall, args.offset, k_cells=(0, 1, 2, 5))
    for mode, r in result["precision"]["modes"].items():
        print(f"precision trap [{mode:8s}]: {'PASS' if r['all_ok'] else 'FAIL'}")
    result["rebase_policy"] = rebase_policy_check(args.offset, args.aspect_ratio)
    result["drift_tolerance"] = drift_tolerance_check(args.offset, args.aspect_ratio)
    if not args.skip_sweep:
        result["extent_sweep"] = extent_sweep(
            radii=(1.0, 100.0),
            ds_over_extent=(1e-7, 1e-8, 1e-9, 1e-10),
            offset=args.offset,
            aspect_ratio=args.aspect_ratio,
        )

    out = Path(__file__).parent / "results"
    out.mkdir(exist_ok=True)
    tag = re.sub(r"[^A-Za-z0-9]+", "-", f"{result['host']}-{adapter}")[:80].strip("-")
    path = out / f"{tag}.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(f"written {path}")


if __name__ == "__main__":
    main()
