---
title: shore.exceptions
outline: [2, 3]
---

# `shore.exceptions`

SHORE-specific exception hierarchy.

## `ShoreError` <Badge type="tip" text="class" /> {#ShoreError}

```python
class ShoreError(Exception)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/exceptions.py#L6">source</a>

Base exception for all SHORE-specific errors.

## `MeshInputError` <Badge type="tip" text="class" /> {#MeshInputError}

```python
class MeshInputError(ShoreError, ValueError)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/exceptions.py#L10">source</a>

Raised when the input geometry is invalid or cannot be read.

## `ParameterError` <Badge type="tip" text="class" /> {#ParameterError}

```python
class ParameterError(ShoreError, ValueError)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/exceptions.py#L14">source</a>

Raised when a function argument is out of range or mutually exclusive.

## `MarchingError` <Badge type="tip" text="class" /> {#MarchingError}

```python
class MarchingError(ShoreError, RuntimeError)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/exceptions.py#L18">source</a>

Raised when hyperbolic marching detects an invalid grid state.
