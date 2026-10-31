"""Command forms of the GUI (ADR 0037): one form per Project command.

A form is a list of fields; submitting it calls ``Session.run(method,
**values)``, so the command and its journal line are the Project's own. The
specs are Qt-free; the window renders them as widgets.

Values are typed by their field: text, int, float, bool, a choice, a path
(file or directory), an optional float, or a Python literal (``ast.
literal_eval``) for structured arguments such as topology options or
primitive parameters.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass, field
from typing import Any, Literal

Kind = Literal["text", "int", "float", "optfloat", "bool", "choice", "file", "dir", "literal"]


@dataclass(frozen=True)
class Field:
    name: str
    kind: Kind
    default: Any = None
    choices: tuple[str, ...] = ()
    label: str = ""
    keyword: bool = True  # passed as name=value; False: positional
    spread: bool = False  # a literal dict spread into keyword arguments
    omit_default: bool = True  # leave the argument out when it equals the default

    def parse(self, raw: Any) -> Any:
        """The typed value of *raw* (a widget's text, a bool, or a value)."""
        if not isinstance(raw, str):
            return raw
        text = raw.strip()
        if self.kind in ("text", "file", "dir", "choice"):
            return text
        if self.kind == "int":
            return int(text)
        if self.kind == "float":
            return float(text)
        if self.kind == "optfloat":
            return None if text == "" else float(text)
        if self.kind == "bool":
            return text.lower() in ("1", "true", "yes", "on")
        if self.kind == "literal":
            return ast.literal_eval(text) if text else self.default
        raise ValueError(f"unknown field kind {self.kind!r}")


@dataclass(frozen=True)
class Form:
    method: str
    title: str
    fields: tuple[Field, ...] = field(default_factory=tuple)
    long: bool = False  # run on a worker thread, with progress and cancel

    def arguments(self, values: dict[str, Any]) -> tuple[list[Any], dict[str, Any]]:
        """Positional and keyword arguments of ``Project.<method>`` from *values*."""
        args: list[Any] = []
        kwargs: dict[str, Any] = {}
        for f in self.fields:
            value = f.parse(values[f.name]) if f.name in values else f.default
            if not f.keyword:
                args.append(value)
            elif f.spread:
                if value:
                    if not isinstance(value, dict):
                        raise ValueError(f"{f.name}: a dict literal is required")
                    kwargs.update(value)
            elif not (f.omit_default and value == f.default):
                kwargs[f.name] = value
        return args, kwargs


TOPOLOGIES = ("cubed_sphere", "ogrid", "cgrid", "ch", "oh")
PRIMITIVES = ("box", "annulus", "flat_caps", "h_grid_rect", "h_grid_tube")

FORMS: tuple[Form, ...] = (
    Form(
        "load_geometry",
        "Load geometry",
        (
            Field("path", "file", "", label="STL / STEP / IGES / BREP", keyword=False),
            Field("name", "text", "body"),
            Field("skip_hygiene", "bool", False),
        ),
    ),
    Form("use_component", "Use component", (Field("name", "text", "main", keyword=False),)),
    Form(
        "build_topology",
        "Build topology",
        (
            Field("topology", "choice", "cubed_sphere", TOPOLOGIES, keyword=False),
            Field("geometry", "text", "body"),
            Field("nk", "int", 40, omit_default=False),
            Field("ds", "float", 1e-3),
            Field("growth", "float", 1.15),
            Field("params", "literal", None, label="ni, nj, ... as a dict", spread=True),
            Field("options", "literal", None, label="topology options (ch / oh)"),
        ),
        long=True,
    ),
    Form(
        "build_primitive",
        "Build primitive",
        (
            Field("kind", "choice", "box", PRIMITIVES, keyword=False),
            Field("params", "literal", None, label="parameters as a dict", spread=True),
        ),
        long=True,
    ),
    Form(
        "march",
        "March",
        (
            Field("smooth_strength", "float", 0.2),
            Field("smooth_iters", "int", 2),
            Field("smooth_mode", "choice", "limited", ("limited", "legacy")),
            Field("march_scheme", "choice", "algebraic", ("algebraic", "hyperbolic")),
        ),
        long=True,
    ),
    Form("snap", "Snap to geometry", (Field("exact", "bool", False),)),
    Form(
        "create_blocking",
        "Create blocking",
        (
            Field("ni", "int", 11, keyword=False),
            Field("nj", "int", 11, keyword=False),
            Field("nk", "int", 11, keyword=False),
            Field("geometry", "text", "body"),
            Field("margin", "float", 0.0),
            Field("vmin", "literal", None, label="(x, y, z), with vmax instead of a geometry"),
            Field("vmax", "literal", None),
        ),
    ),
    Form(
        "split_blocking",
        "Split blocking",
        (
            Field("block", "text", "b0", keyword=False),
            Field("axis", "choice", "i", ("i", "j", "k"), keyword=False),
            Field("fraction", "float", 0.5, keyword=False),
        ),
    ),
    Form(
        "delete_blocks",
        "Delete blocking blocks",
        (Field("blocks", "literal", (), label="labels, e.g. ['b13']", keyword=False),),
    ),
    Form(
        "ogrid",
        "O-grid",
        (
            Field("blocks", "literal", (), label="labels, e.g. ['b0'] (empty: all)", keyword=False),
            Field("faces", "literal", (), label="layered faces, e.g. [('b0', 'i_lo')]"),
            Field("exclude", "literal", (), label="open faces, e.g. [('b0', 'k_lo')]"),
            Field("factor", "float", 0.5),
            Field("nodes", "int", 5),
        ),
    ),
    Form(
        "move_vertex",
        "Move blocking vertex",
        (
            Field("block", "text", "b0", keyword=False),
            Field("vertex", "text", "v0", label="v0..v7", keyword=False),
            Field("xyz", "literal", (0.0, 0.0, 0.0), label="(x, y, z)", keyword=False),
        ),
    ),
    Form(
        "associate_blocking",
        "Associate blocking entity",
        (
            Field("block", "text", "b0", keyword=False),
            Field("entity", "text", "k_lo", label="face, edge or v0..v7", keyword=False),
            Field("geometry", "text", "body", keyword=False),
            Field(
                "kind",
                "choice",
                "body",
                ("body", "surface", "curve", "point", "none"),
                keyword=False,
            ),
            Field("indices", "literal", (), label="entity indices, e.g. [3, 4]", keyword=False),
        ),
    ),
    Form(
        "symmetry_plane",
        "Symmetry plane",
        (
            Field(
                "point", "literal", (0.0, 0.0, 0.0), label="(x, y, z) on the plane", keyword=False
            ),
            Field("normal", "literal", (0.0, 0.0, 1.0), label="(nx, ny, nz)", keyword=False),
        ),
    ),
    Form("premesh", "Premesh blocking", (Field("exact", "bool", False),), long=True),
    Form("orthogonalize", "Orthogonal wall lines", (Field("strength", "float", 1.0),), long=True),
    Form(
        "wall_layer",
        "Wall layer",
        (Field("first_step", "float", 1e-3, keyword=False),),
        long=True,
    ),
    Form(
        "write_xall", "Write Xall run", (Field("directory", "dir", "", keyword=False),), long=True
    ),
    Form("write_cgns", "Write CGNS", (Field("path", "file", "", keyword=False),)),
    Form("write_plot3d", "Write Plot3D", (Field("path", "file", "", keyword=False),)),
)


def form(method: str) -> Form:
    for f in FORMS:
        if f.method == method:
            return f
    raise KeyError(method)
