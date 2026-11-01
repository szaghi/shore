---
title: shore.gui.cli_forms
outline: [2, 3]
---

# `shore.gui.cli_forms`

GUI forms generated from the CLI (ADR 0062).

Every leaf command of the `shore` typer app becomes a form keyed
`"cli:<command path>"`: its options and arguments as typed fields (named as
the CLI's Python parameters, with its defaults, choices and help), less the
ones the GUI owns. Submitting one runs the CLI command in-process on the
session's project (`shore.cli.using_project`), so the GUI has every CLI
capability with the CLI's own defaults and journal lines.

Qt-free: the window renders the fields; `argv` turns submitted values back
into the command line.

## `leaf_commands` <Badge type="info" text="function" /> {#leaf_commands}

```python
leaf_commands() -> list[tuple[str, click.Command]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/cli_forms.py#L58">source</a>

`(command path, click command)` of every leaf command of the CLI.

## `cli_forms` <Badge type="info" text="function" /> {#cli_forms}

```python
cli_forms() -> tuple[Form, ...]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/cli_forms.py#L119">source</a>

A form per CLI leaf command, `MAPPED` and `EXCLUDED` aside.

## `argv` <Badge type="info" text="function" /> {#argv}

```python
argv(form: Form, values: dict[str, Any]) -> list[str]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/cli_forms.py#L147">source</a>

The command line of *form* (a `cli:` form) for submitted *values*.
