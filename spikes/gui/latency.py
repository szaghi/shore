"""UI input latency during a long command (ADR 0006, item 4).

A 5 ms ``QTimer`` heartbeat stands in for UI event handling; its lateness is
the input latency a user would feel. A workload runs for ``--seconds`` in a
``QThreadPool`` worker (or a subprocess) while the heartbeat is sampled.

Workloads:

- ``numpy``: Laplacian-style sweeps over a 1e7-node array. Large numpy kernels
  release the GIL, so the UI should stay responsive.
- ``gil``: one long C call that holds the GIL (``list.sort`` on a large
  list). This stands in for a long pybind11 call such as an OCP
  ``BRepMesh_IncrementalMesh`` on a big STEP: pybind11 does not release the GIL
  by default, so the UI freezes for the duration.
- ``gil-process``: the same GIL-holding call in a subprocess (the mitigation).

Usage::

    QT_QPA_PLATFORM=offscreen python latency.py --workload numpy --seconds 30
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import random
import sys
import time

import numpy as np
from PySide6.QtCore import QCoreApplication, QObject, QRunnable, QThreadPool, QTimer, Signal

HEARTBEAT_MS = 5


def numpy_workload(seconds: float) -> int:
    """Smoothing sweeps over a (1e7, 3) float64 array until ``seconds`` elapse."""
    a = np.random.default_rng(0).random((215, 215, 215, 3))
    sweeps = 0
    t_end = time.perf_counter() + seconds
    while time.perf_counter() < t_end:
        a[1:-1] = 0.5 * a[1:-1] + 0.25 * (a[2:] + a[:-2])
        sweeps += 1
    return sweeps


def gil_workload(seconds: float) -> int:
    """Repeated single long GIL-holding C calls (list.sort) for ``seconds``."""
    rng = random.Random(0)
    n = 4_000_000
    calls = 0
    t_end = time.perf_counter() + seconds
    while time.perf_counter() < t_end:
        data = [rng.random() for _ in range(n)]
        data.sort()  # holds the GIL for the whole sort
        calls += 1
    return calls


class _Signals(QObject):
    done = Signal(int)


class Worker(QRunnable):
    """Run a workload in the pool, or in a subprocess for ``gil-process``."""

    def __init__(self, workload: str, seconds: float) -> None:
        super().__init__()
        self.workload = workload
        self.seconds = seconds
        self.signals = _Signals()

    def run(self) -> None:
        if self.workload == "numpy":
            n = numpy_workload(self.seconds)
        elif self.workload == "gil":
            n = gil_workload(self.seconds)
        else:
            with mp.get_context("spawn").Pool(1) as pool:
                n = pool.apply(gil_workload, (self.seconds,))
        self.signals.done.emit(n)


def measure(workload: str, seconds: float) -> dict:
    """Heartbeat lateness statistics while ``workload`` runs."""
    app = QCoreApplication.instance() or QCoreApplication(sys.argv)
    late: list[float] = []
    last = [time.perf_counter()]

    def beat() -> None:
        now = time.perf_counter()
        late.append(max(0.0, 1e3 * (now - last[0]) - HEARTBEAT_MS))
        last[0] = now

    timer = QTimer()
    timer.setInterval(HEARTBEAT_MS)
    timer.timeout.connect(beat)
    result: dict = {}

    def finished(n: int) -> None:
        timer.stop()
        a = np.asarray(late)
        result.update(
            workload=workload,
            seconds=seconds,
            iterations=n,
            beats=int(a.size),
            late_ms_p50=round(float(np.percentile(a, 50)), 2),
            late_ms_p99=round(float(np.percentile(a, 99)), 2),
            late_ms_max=round(float(a.max()), 2),
            budget_ms=100,
            ok=bool(a.max() <= 100.0),
        )
        app.quit()

    worker = Worker(workload, seconds)
    worker.signals.done.connect(finished)
    timer.start()
    QThreadPool.globalInstance().start(worker)
    app.exec()
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workload", choices=("numpy", "gil", "gil-process"), default="numpy")
    ap.add_argument("--seconds", type=float, default=30.0)
    args = ap.parse_args()
    print(json.dumps(measure(args.workload, args.seconds), indent=2))


if __name__ == "__main__":
    main()
