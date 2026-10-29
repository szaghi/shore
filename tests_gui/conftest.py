"""Fixtures of the GUI tests (ADR 0037): an offscreen Qt application.

The GUI tests live outside ``testpaths`` (``tests/``) because they need the
``gui`` extra; the GUI CI job runs ``pytest tests_gui`` with
``QT_QPA_PLATFORM=offscreen`` and a software Vulkan adapter (lavapipe).
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="session")
def qapp():
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="session")
def sphere_stl(tmp_path_factory) -> Path:
    import trimesh

    path = tmp_path_factory.mktemp("stl") / "sphere.stl"
    trimesh.creation.icosphere(subdivisions=3, radius=1.0).export(str(path))
    return path
