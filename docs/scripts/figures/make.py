"""Draw the docs figures, each in a light and a dark variant.

Run from the project root::

    .venv/bin/python docs/scripts/figures/make.py                 # every figure
    .venv/bin/python docs/scripts/figures/make.py sphere-blocks   # some
    .venv/bin/python docs/scripts/figures/make.py --list

Each module below exposes ``FIGURES``: a dict of figure name -> function of the
mode (``"light"`` or ``"dark"``) that writes ``docs/public/figures/<mode>/<name>.<ext>``
through ``theme.save`` (matplotlib) or ``theme.shot`` (PyVista).
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

MODULES = ("schematics", "charts", "renders_classic", "renders_cad", "gui_shots")


def registry() -> dict:
    out = {}
    for name in MODULES:
        try:
            mod = importlib.import_module(name)
        except ModuleNotFoundError as err:
            if err.name != name:
                raise
            continue
        out.update(mod.FIGURES)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("names", nargs="*", help="figures to draw (default: all)")
    ap.add_argument("--mode", choices=("light", "dark", "both"), default="both")
    ap.add_argument("--list", action="store_true", help="list the figures and exit")
    args = ap.parse_args(argv)
    figures = registry()
    if args.list:
        print("\n".join(sorted(figures)))
        return 0
    unknown = set(args.names) - set(figures)
    if unknown:
        print("unknown figures:", *sorted(unknown))
        return 1
    modes = ("light", "dark") if args.mode == "both" else (args.mode,)
    for name in args.names or sorted(figures):
        for mode in modes:
            t = time.perf_counter()
            figures[name](mode)
            print(f"{name:24s} {mode:5s} {time.perf_counter() - t:6.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
