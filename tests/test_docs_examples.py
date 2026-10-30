"""The scripts the docs include run to the end (examples/tutorials, examples/cookbook).

Each script runs in a temporary directory. A script that needs an optional
extra says so on a line of its own: ``# requires: cad`` (build123d),
``# requires: cgns`` (h5py) or ``# requires: vis`` (pyvista); it is skipped
when the extra is missing.
"""

from __future__ import annotations

import re
import runpy
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
SCRIPTS = sorted([*EXAMPLES.glob("tutorials/*.py"), *EXAMPLES.glob("cookbook/*.py")])
EXTRAS = {"cad": "build123d", "cgns": "h5py", "vis": "pyvista", "gui": "PySide6"}


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: f"{p.parent.name}/{p.stem}")
def test_example_runs(script: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for extra in re.findall(r"^# requires: (\w+)\s*$", script.read_text(), re.MULTILINE):
        pytest.importorskip(EXTRAS[extra])
    monkeypatch.chdir(tmp_path)
    runpy.run_path(str(script), run_name="__main__")
