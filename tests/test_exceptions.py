"""Verify the ShoreError exception hierarchy and dual-inheritance design."""

from __future__ import annotations

import pytest

from shore.exceptions import MarchingError, MeshInputError, ParameterError, ShoreError


def test_exception_hierarchy():
    assert issubclass(MeshInputError, ShoreError)
    assert issubclass(MeshInputError, ValueError)
    assert issubclass(ParameterError, ShoreError)
    assert issubclass(ParameterError, ValueError)
    assert issubclass(MarchingError, ShoreError)
    assert issubclass(MarchingError, RuntimeError)


def test_catch_as_shore_error():
    with pytest.raises(ShoreError):
        raise MeshInputError("bad input")

    with pytest.raises(ShoreError):
        raise ParameterError("bad param")

    with pytest.raises(ShoreError):
        raise MarchingError("bad march")


def test_catch_as_stdlib():
    with pytest.raises(ValueError):
        raise MeshInputError("bad input")

    with pytest.raises(ValueError):
        raise ParameterError("bad param")

    with pytest.raises(RuntimeError):
        raise MarchingError("bad march")
