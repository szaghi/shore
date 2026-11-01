---
title: shore.gui.forms
outline: [2, 3]
---

# `shore.gui.forms`

Command forms of the GUI (ADR 0037): one form per Project command.

A form is a list of fields; submitting it calls ``Session.run(method,
**values)``, so the command and its journal line are the Project's own. The
specs are Qt-free; the window renders them as widgets.

Values are typed by their field: text, int, float, bool, a choice, a path
(file or directory), an optional int or float, a list (a Python list literal,
or space-separated words), or a Python literal (`ast.literal_eval`) for
structured arguments such as topology options or primitive parameters.

`FORMS` holds the forms over Project methods: the hand-written ones and,
since ADR 0062, one generated from the signature of every other Project
command. `ALL_FORMS` adds the session's own forms (`run_journal`,
`view_file`) and a form per CLI command (`shore.gui.cli_forms`), run
in-process on the session's project.

**Re-exported**

- `ALL_FORMS`
- [`EXCLUDED`](./shore.gui.cli_forms#EXCLUDED)
- `FORMS`
- [`MAPPED`](./shore.gui.cli_forms#MAPPED)

## `Field` <Badge type="tip" text="dataclass" /> {#Field}

```python
Field(
    name: str,
    kind: Kind,
    default: Any = None,
    choices: tuple[str, ...] = (),
    label: str = '',
    keyword: bool = True,
    spread: bool = False,
    omit_default: bool = True,
    help: str = '',
    convert: str = '',
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L38">source</a>

### `Field.parse` <Badge type="info" text="method" /> {#Field-parse}

```python
parse(raw: Any) -> Any
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L50">source</a>

The typed value of *raw* (a widget's text, a bool, or a value).

## `Form` <Badge type="tip" text="dataclass" /> {#Form}

```python
Form(
    method: str,
    title: str,
    fields: tuple[Field, ...] = field(default_factory=tuple),
    long: bool = False,
    target: Target = 'project',
    cli: tuple[str, ...] = (),
    group: str = 'Project',
    help: str = '',
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L77">source</a>

### `Form.arguments` <Badge type="info" text="method" /> {#Form-arguments}

```python
arguments(values: dict[str, Any]) -> tuple[list[Any], dict[str, Any]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L87">source</a>

Positional and keyword arguments of `Project.<method>` from *values*.

## `form` <Badge type="info" text="function" /> {#form}

```python
form(method: str) -> Form
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L382">source</a>
