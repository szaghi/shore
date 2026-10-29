"""Measure ADR 0006's GUI budgets on the real GUI (ADR 0037).

Run on a machine with a GPU (not WSL2, whose wgpu adapter is the CPU)::

    python scripts/gui_bench.py --nodes 10000000 --out bench

1. a project with one box primitive of about ``--nodes`` nodes is built and
   saved (its grids); the timed part is ``MainWindow.open_project(...,
   replay=False)`` up to the first frame (budget 2 s);
2. orbit: the camera is turned in 1-degree steps, each frame drawn
   synchronously (offscreen, including the GPU->CPU readback: an upper
   bound; budget 16 ms);
3. pick at the image centre (budget 50 ms);
4. a hyperbolic march on the worker thread while a 10 ms QTimer measures
   the UI's lateness (budget 100 ms).

Prints one JSON object. The frame times are meaningless on a CPU adapter
(``(CPU)`` in the adapter line).
"""

from __future__ import annotations

import argparse
import json
import os
import time
import warnings
from pathlib import Path

import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=int, default=10_000_000)
    parser.add_argument("--out", type=Path, default=Path("gui-bench"))
    parser.add_argument("--frames", type=int, default=120)
    args = parser.parse_args()
    warnings.simplefilter("ignore")

    import trimesh
    import wgpu
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication

    from shore.gui.window import MainWindow
    from shore.project import Project

    app = QApplication.instance() or QApplication([])
    adapter = wgpu.gpu.request_adapter_sync(power_preference="high-performance").summary
    result: dict[str, object] = {"adapter": adapter}

    project_dir = args.out / f"box{args.nodes}"
    n = round(args.nodes ** (1.0 / 3.0))
    if not (project_dir / "project.toml").is_file():
        p = Project()
        p.build_primitive("box", vmin=(0.0, 0.0, 0.0), vmax=(1.0, 1.0, 1.0), ni=n, nj=n, nk=n)
        p.save(project_dir)
    result["nodes"] = n**3

    # 1. load to first frame
    w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    t0 = time.perf_counter()
    w.open_project(str(project_dir), replay=False)
    w.viewport.draw_now()
    result["load_first_frame_s"] = round(time.perf_counter() - t0, 3)

    # 2. orbit
    vp = w.viewport
    times = []
    for step in range(args.frames):
        angle = np.radians(1.0)
        pos = np.asarray(vp.camera.local.position, np.float64)
        c, s = np.cos(angle), np.sin(angle)
        vp.camera.local.position = (c * pos[0] - s * pos[1], s * pos[0] + c * pos[1], pos[2])
        vp.camera.look_at((0.0, 0.0, 0.0))
        t = time.perf_counter()
        vp.draw_now()
        if step >= 5:  # warm-up frames excluded
            times.append(1e3 * (time.perf_counter() - t))
    result["orbit_ms_mean"] = round(float(np.mean(times)), 2)
    result["orbit_ms_p95"] = round(float(np.percentile(times, 95)), 2)

    # 3. pick at the centre
    width, height = (int(v) for v in vp.renderer.logical_size)
    t = time.perf_counter()
    picked = vp.pick(width / 2, height / 2)
    result["pick_ms"] = round(1e3 * (time.perf_counter() - t), 2)
    result["picked"] = None if picked is None else [picked.face, list(picked.cell)]

    # 4. UI lateness during a march on the worker
    stl = args.out / "sphere.stl"
    trimesh.creation.icosphere(subdivisions=3).export(str(stl))
    m = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
    m.run_form("load_geometry", {"path": str(stl)})
    m.run_form(
        "build_topology",
        {"topology": "cubed_sphere", "nk": "30", "ds": "0.01", "growth": "1.02",
         "params": "{'ni': 65, 'nj': 128}"},
    )  # fmt: skip
    m.wait_idle()
    late: list[float] = []
    last = [time.perf_counter()]

    def beat() -> None:
        now = time.perf_counter()
        late.append(now - last[0] - 0.010)
        last[0] = now

    timer = QTimer()
    timer.setInterval(10)
    timer.timeout.connect(beat)
    timer.start()
    t = time.perf_counter()
    m.run_form("march", {"march_scheme": "hyperbolic"})
    m.wait_idle()
    timer.stop()
    result["march_s"] = round(time.perf_counter() - t, 2)
    result["ui_lateness_ms_max"] = round(1e3 * max(late), 1)
    result["march_error"] = None if m.last_error is None else str(m.last_error)
    print(json.dumps(result, indent=2))
    app.quit()


if __name__ == "__main__":
    main()
