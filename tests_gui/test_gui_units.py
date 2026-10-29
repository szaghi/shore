"""Unit tests of the GUI's Qt-free parts and small window behaviours (ADR 0037)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.gui.forms import FORMS, Field, form
from shore.gui.session import Session
from shore.gui.viewport import (
    FACES,
    boundary_values,
    colormap,
    face_grids,
    quad_cell,
    quad_centre,
    quad_indices,
    quad_offsets,
)


class TestForms:
    def test_field_parsing(self):
        assert Field("n", "int").parse(" 7 ") == 7
        assert Field("x", "float").parse("1e-3") == 1e-3
        assert Field("t", "optfloat").parse("") is None
        assert Field("b", "bool").parse("yes") is True
        assert Field("d", "literal").parse("{'ni': 9}") == {"ni": 9}
        assert Field("d", "literal", default=None).parse("") is None

    def test_arguments(self):
        args, kwargs = form("build_topology").arguments(
            {"topology": "ogrid", "nk": "5", "ds": "0.001", "params": "{'ni': 9, 'nj': 16}"}
        )
        assert args == ["ogrid"]
        assert kwargs == {"nk": 5, "ni": 9, "nj": 16}  # defaults left out, params spread

    def test_spread_needs_a_dict(self):
        with pytest.raises(ValueError, match="dict"):
            form("build_primitive").arguments({"kind": "box", "params": "[1, 2]"})

    def test_every_form_has_a_title_and_unique_method(self):
        methods = [f.method for f in FORMS]
        assert len(methods) == len(set(methods))
        assert all(f.title for f in FORMS)
        with pytest.raises(KeyError):
            form("nope")


class TestStructuredArithmetic:
    @pytest.mark.parametrize("shape", [(4, 5, 6), (3, 7, 2), (9, 2, 5)])
    def test_quad_cell_round_trip(self, shape):
        nodes = np.random.default_rng(1).normal(size=(*shape, 3))
        offsets = quad_offsets(shape)
        grids = face_grids(nodes)
        quads = np.concatenate(
            [quad_indices(g.shape[0], g.shape[1]) + 0 for g in grids]
        )  # per-face local numbering, concatenated in FACES order
        assert len(quads) == offsets[-1]
        for q in range(offsets[-1]):
            face, cell = quad_cell(shape, q)
            f = FACES.index(face)
            local = q - offsets[f]
            corners = grids[f].reshape(-1, 3)[quad_indices(*grids[f].shape[:2])[local]]
            assert np.allclose(corners.mean(axis=0), quad_centre(nodes, face, cell))
            assert all(0 <= c <= n - 2 for c, n in zip(cell, shape, strict=True))

    def test_boundary_values_follow_quads(self):
        shape = (4, 5, 6)
        cells = np.arange(np.prod([n - 1 for n in shape]), dtype=float).reshape(3, 4, 5)
        values = boundary_values(cells)
        assert len(values) == quad_offsets(shape)[-1]
        for q in (0, 7, 30, len(values) - 1):
            _, cell = quad_cell(shape, q)
            assert values[q] == cells[cell]

    def test_colormap(self):
        rgba = colormap(np.array([-1.0, 0.0, 1.0, np.nan]), (-1.0, 1.0))
        assert rgba.shape == (4, 4) and rgba.dtype == np.float32
        assert rgba[0, 0] > rgba[0, 2]  # low end red
        assert rgba[2, 2] > rgba[2, 0]  # high end blue
        assert np.allclose(rgba[1, :3], 1.0)  # midpoint white
        assert np.allclose(rgba[3, :3], 0.6)  # NaN grey


class TestSession:
    def test_unknown_command(self):
        with pytest.raises(AttributeError):
            Session().run("_private")
        with pytest.raises(AttributeError):
            Session().run("no_such_command")

    def test_save_needs_a_directory(self):
        with pytest.raises(ValueError, match="directory"):
            Session().save()

    def test_listeners_called_even_on_failure(self):
        s = Session()
        calls = []
        s.subscribe(lambda: calls.append(1))
        with pytest.raises(Exception):  # noqa: B017 - any command error
            s.run("march")
        assert calls == [1]

    def test_tree_and_blocks_of_an_empty_project(self):
        s = Session()
        assert s.blocks() == [] and s.tree() == [("main", True, [])]
        assert s.journal_lines() == []


class TestWindow:
    def test_cancelled_dialogs_do_nothing(self, qapp):
        from shore.gui.window import MainWindow

        w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
        w.open_project()
        w.save_project()
        assert w.last_error is None and w.session.directory is None

    def test_forms_panel_lists_every_form(self, qapp):
        from shore.gui.window import MainWindow

        w = MainWindow(offscreen_viewport=True, ask_path=lambda *a: None)
        assert [w.command.itemData(i) for i in range(w.command.count())] == [
            f.method for f in FORMS
        ]
        w.command.setCurrentIndex(w.command.findData("march"))
        assert set(w._form_values()) == {f.name for f in form("march").fields}
