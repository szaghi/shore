"""SHORE's desktop GUI (ADR 0006, 0037): ``shore gui`` / ``shore-gui``.

Needs the ``gui`` extra (PySide6-Essentials, pygfx, wgpu, rendercanvas).
``shore.gui.session`` and ``shore.gui.forms`` are Qt-free; the window lives
in ``shore.gui.window`` and the renderer behind ``shore.gui.viewport``.
"""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    """``shore-gui [PROJECT_DIR] [--no-replay] [--smoke]``."""
    args = list(sys.argv[1:] if argv is None else argv)
    replay = "--no-replay" not in args
    smoke = "--smoke" in args
    rest = [a for a in args if not a.startswith("--")]
    try:
        from shore.gui.window import main as run
    except ImportError as exc:
        print(f"shore gui needs the gui extra: pip install 'shore-mesh[gui]' ({exc})")
        return 1
    return run(rest[0] if rest else None, replay=replay, smoke=smoke)
