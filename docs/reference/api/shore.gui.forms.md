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
(file or directory), an optional float, or a Python literal (``ast.
literal_eval``) for structured arguments such as topology options or
primitive parameters.

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
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L23">source</a>

### `Field.parse` <Badge type="info" text="method" /> {#Field-parse}

```python
parse(raw: Any) -> Any
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L33">source</a>

The typed value of *raw* (a widget's text, a bool, or a value).

## `Form` <Badge type="tip" text="dataclass" /> {#Form}

```python
Form(
    method: str,
    title: str,
    fields: tuple[Field, ...] = field(default_factory=tuple),
    long: bool = False,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L54">source</a>

### `Form.arguments` <Badge type="info" text="method" /> {#Form-arguments}

```python
arguments(values: dict[str, Any]) -> tuple[list[Any], dict[str, Any]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L60">source</a>

Positional and keyword arguments of `Project.<method>` from *values*.

## `form` <Badge type="info" text="function" /> {#form}

```python
form(method: str) -> Form
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/gui/forms.py#L207">source</a>
