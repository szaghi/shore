"""Command forms of the GUI (ADR 0037): one form per Project command.

A form is a list of fields; submitting it calls ``Session.run(method,
**values)``, so the command and its journal line are the Project's own. The
specs are Qt-free; the window renders them as widgets.

Values are typed by their field: text, int, float, bool, a choice, a path
(file or directory), an optional int or float, a list (a Python list literal,
or space-separated words), or a Python literal (``ast.literal_eval``) for
structured arguments such as topology options or primitive parameters.

``FORMS`` holds the forms over Project methods: the hand-written ones and,
since ADR 0062, one generated from the signature of every other Project
command. ``ALL_FORMS`` adds the session's own forms (``run_journal``,
``view_file``) and a form per CLI command (``shore.gui.cli_forms``), run
in-process on the session's project.
"""

from __future__ import annotations

import ast
import inspect
import shlex
from dataclasses import dataclass, field
from typing import Any, Literal

from shore.gui.cli_forms import EXCLUDED, MAPPED, cli_forms

Kind = Literal[
    "text", "int", "optint", "float", "optfloat", "bool", "choice", "file", "dir", "list", "literal"
]
Target = Literal["project", "session", "cli"]

__all__ = ["ALL_FORMS", "EXCLUDED", "FORMS", "MAPPED", "Field", "Form", "form"]


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
    help: str = ""  # the tooltip (the CLI option's help)
    convert: str = ""  # a named converter of the value (``CONVERTERS``)

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
        if self.kind == "optint":
            return None if text == "" else int(text)
        if self.kind == "list":
            if text.startswith(("[", "(")):
                return list(ast.literal_eval(text))
            return shlex.split(text)
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
    target: Target = "project"  # Project.<method>, Session.<method>, or a CLI command
    cli: tuple[str, ...] = ()  # the CLI command path of a "cli" form
    group: str = "Project"  # where the form sits in the Commands menu
    help: str = ""

    def arguments(self, values: dict[str, Any]) -> tuple[list[Any], dict[str, Any]]:
        """Positional and keyword arguments of ``Project.<method>`` from *values*."""
        args: list[Any] = []
        kwargs: dict[str, Any] = {}
        for f in self.fields:
            value = f.parse(values[f.name]) if f.name in values else f.default
            if f.convert:
                value = CONVERTERS[f.convert](value)
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
        "write_xall",
        "Write Xall run",
        (
            Field("directory", "dir", "", keyword=False),
            Field(
                "options",
                "file",
                "",
                label="options TOML ([run], [components.<name>])",
                spread=True,
                convert="xall_options",
                help="As shore xall --options.",
            ),
        ),
        long=True,
        group="Export",
    ),
    Form(
        "write_cgns",
        "Write CGNS",
        (
            Field("path", "file", "", keyword=False),
            Field("components", "list", [], label="components (empty: all)"),
        ),
        group="Export",
    ),
    Form(
        "write_plot3d",
        "Write Plot3D",
        (
            Field("path", "file", "", keyword=False),
            Field("components", "list", [], label="components (empty: all)"),
        ),
        group="Export",
    ),
)


def _xall_options(path: Any) -> dict[str, Any]:
    from pathlib import Path

    from shore.cli import _xall_options as read

    return read(Path(path)) if path else {}


#: Value converters named by ``Field.convert``.
CONVERTERS: dict[str, Any] = {"xall_options": _xall_options}

#: Field names whose values are paths, for the generated Project forms.
_PATHS = frozenset(
    {"path", "stem", "body_geo", "other", "adjacency", "grd", "directory", "geo_paths"}
)


def _kind(name: str, annotation: Any) -> Kind:
    text = str(annotation)
    if "bool" in text:
        return "bool"
    if text.startswith("int | None"):
        return "optint"
    if text == "int":
        return "int"
    if text.startswith("float | None"):
        return "optfloat"
    if text == "float":
        return "float"
    if name in _PATHS and ("list" in text or "Sequence" in text or "tuple" in text):
        return "list"
    if name in _PATHS:
        return "file"
    if text in ("str", "str | None"):
        return "text"
    return "literal"


def _project_form(method: str) -> Form:
    """A form from ``Project.<method>``'s signature (ADR 0062)."""
    from shore.project import Project

    call = getattr(Project, method)
    fields: list[Field] = []
    for p in list(inspect.signature(call).parameters.values())[1:]:
        if p.kind is inspect.Parameter.VAR_KEYWORD:
            fields.append(Field(p.name, "literal", None, label=f"{p.name} as a dict", spread=True))
            continue
        kind = _kind(p.name, p.annotation)
        required = p.default is inspect.Parameter.empty
        default = None if required else p.default
        if kind == "list":
            default = list(default) if default else []
        elif kind == "file" and default is None:
            default = ""
        keyword = not required or p.kind is inspect.Parameter.KEYWORD_ONLY
        fields.append(Field(p.name, kind, default, keyword=keyword, omit_default=not required))
    doc = (call.__doc__ or "").strip().splitlines()
    title = method.replace("_", " ").capitalize()
    return Form(method, title, tuple(fields), long=True, help=doc[0] if doc else "")


def _command_methods() -> list[str]:
    import dataclasses

    from shore.project import commands

    return sorted(
        {
            c.method
            for c in vars(commands).values()
            if isinstance(c, type)
            and dataclasses.is_dataclass(c)
            and isinstance(getattr(c, "method", None), str)
        }
    )


_HAND = {f.method for f in FORMS}
SESSION_FORMS: tuple[Form, ...] = (
    Form(
        "run_journal",
        "Run journal",
        (Field("path", "file", "", label="journal (.py)", keyword=False),),
        long=True,
        target="session",
        group="Project",
        help="As shore run: the journal's project becomes the session's.",
    ),
    Form(
        "view_file",
        "View file",
        (Field("path", "file", "", label=".geo blocks or .npz surface", keyword=False),),
        target="session",
        group="Inspect",
        help="As shore view: a .geo file loads as blocks, a .npz surface is shown.",
    ),
)
#: Forms over Project methods: the hand-written ones and one per other
#: Project command (ADR 0037's G12: every one calls ``Project.<method>``).
FORMS = FORMS + tuple(_project_form(m) for m in _command_methods() if m not in _HAND)
#: Every form: Project, session and CLI (ADR 0062).
ALL_FORMS: tuple[Form, ...] = FORMS + SESSION_FORMS + cli_forms()


def form(method: str) -> Form:
    for f in ALL_FORMS:
        if f.method == method:
            return f
    raise KeyError(method)
