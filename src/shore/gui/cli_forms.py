"""GUI forms generated from the CLI (ADR 0062).

Every leaf command of the ``shore`` typer app becomes a form keyed
``"cli:<command path>"``: its options and arguments as typed fields (named as
the CLI's Python parameters, with its defaults, choices and help), less the
ones the GUI owns. Submitting one runs the CLI command in-process on the
session's project (``shore.cli.using_project``), so the GUI has every CLI
capability with the CLI's own defaults and journal lines.

Qt-free: the window renders the fields; ``argv`` turns submitted values back
into the command line.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import click

if TYPE_CHECKING:
    from shore.gui.forms import Field, Form

#: Options the GUI owns: the session is the project, and the GUI journals.
OWNED = frozenset({"project", "component", "no_replay", "journal", "help"})

#: CLI commands that open a project directory or a viewer: the GUI's own forms.
MAPPED: dict[str, str] = {
    "xall": "write_xall",
    "cgns": "write_cgns",
    "plot3d": "write_plot3d",
    "run": "run_journal",
    "view": "view_file",
}

#: CLI commands with no GUI counterpart.
EXCLUDED: frozenset[str] = frozenset({"gui"})

_GROUPS = {
    "mesh": "Mesh",
    "collar": "Mesh",
    "cap-tip": "Mesh",
    "split": "Mesh",
    "primitive": "Primitives",
    "body": "Bodies and profiles",
    "profile": "Bodies and profiles",
    "grd": "Xall files",
    "grd-merge": "Xall files",
    "cc-par": "Xall files",
    "cc-par-merge": "Xall files",
    "proc-input": "Xall files",
    "proc-input-merge": "Xall files",
    "balance": "Xall files",
    "boxes-from-stl": "Xall files",
    "boxes-export": "Xall files",
}


def leaf_commands() -> list[tuple[str, click.Command]]:
    """``(command path, click command)`` of every leaf command of the CLI."""
    from typer.main import get_command

    from shore.cli import app

    out: list[tuple[str, click.Command]] = []

    def walk(group: Any, prefix: str) -> None:
        for name, sub in group.commands.items():
            if isinstance(sub, click.Group):
                walk(sub, f"{prefix}{name} ")
            else:
                out.append((f"{prefix}{name}", sub))

    walk(get_command(app), "")
    return out


def _field(param: click.Parameter) -> Field:
    from shore.gui.forms import Field

    default = param.default if not callable(param.default) else None
    help_text = getattr(param, "help", None) or ""
    positional = isinstance(param, click.Argument)
    kind: str
    choices: tuple[str, ...] = ()
    ptype = param.type
    if param.multiple or param.nargs == -1:
        kind = "list"
        default = list(default) if default else []
    elif param.nargs > 1:
        kind = "literal"
    elif isinstance(param, click.Option) and param.is_flag:
        kind = "bool"
        default = bool(default)
    elif isinstance(ptype, click.Choice):
        kind, choices = "choice", tuple(str(c) for c in ptype.choices)
    elif isinstance(ptype, click.types.IntParamType):
        kind = "int" if default is not None else "optint"
    elif isinstance(ptype, click.types.FloatParamType):
        kind = "float" if default is not None else "optfloat"
    elif isinstance(ptype, click.types.BoolParamType):
        kind = "bool"
        default = bool(default)
    elif isinstance(ptype, click.Path):
        kind = "dir" if (ptype.dir_okay and not ptype.file_okay) else "file"
    else:
        kind = "text"
    if kind in ("file", "dir", "text") and default is not None:
        default = str(default)
    return Field(
        str(param.name),
        kind,  # type: ignore[arg-type]
        default,
        choices,
        keyword=not positional,
        help=help_text,
    )


def cli_forms() -> tuple[Form, ...]:
    """A form per CLI leaf command, ``MAPPED`` and ``EXCLUDED`` aside."""
    from shore.gui.forms import Form

    forms = []
    for path, command in leaf_commands():
        if path in MAPPED or path in EXCLUDED:
            continue
        params = tuple(
            p for p in command.params if p.name not in OWNED and not getattr(p, "hidden", False)
        )
        title = f"shore {path}"
        helptext = (command.help or "").strip().splitlines()
        forms.append(
            Form(
                f"cli:{path}",
                title,
                tuple(_field(p) for p in params),
                long=True,
                target="cli",
                cli=tuple(path.split()),
                group=_GROUPS.get(path.split()[0], "Inspect"),
                help=helptext[0] if helptext else "",
            )
        )
    return tuple(forms)


def argv(form: Form, values: dict[str, Any]) -> list[str]:
    """The command line of *form* (a ``cli:`` form) for submitted *values*."""
    command = dict(leaf_commands())[" ".join(form.cli)]
    params = {p.name: p for p in command.params}
    out: list[str] = [*form.cli]
    options: list[str] = []
    for f in form.fields:
        value = f.parse(values[f.name]) if f.name in values else f.default
        param = params[f.name]
        if isinstance(param, click.Argument):
            if value is None or value == "" or value == []:
                continue
            out.extend(str(v) for v in value) if isinstance(value, list | tuple) else out.append(
                str(value)
            )
            continue
        if value == f.default or value is None or value == "":
            continue
        opt = param.opts[0]
        if f.kind == "bool":
            if value:
                options.append(opt)
            elif param.secondary_opts:
                options.append(param.secondary_opts[0])
        elif f.kind == "list":
            for v in value:
                options.extend([opt, str(v)])
        elif isinstance(value, list | tuple):
            options.append(opt)
            options.extend(str(v) for v in value)
        else:
            options.extend([opt, str(value)])
    return out + options
