"""Tests for the hyperbolic marching kernel."""

from __future__ import annotations

import numpy as np
import pytest

from shore.volume.hyperbolic import _compute_normals, march
from tests.conftest import sphere_layer as _sphere_layer


def _flat_layer(ni: int, nj: int) -> np.ndarray:
    """Flat square layer in the xy-plane — normals should point in +z."""
    i = np.linspace(0.0, 1.0, ni)
    j = np.linspace(0.0, 1.0, nj)
    ig, jg = np.meshgrid(i, j, indexing="ij")
    layer = np.stack([ig, jg, np.zeros_like(ig)], axis=-1)
    return layer


def test_normals_flat_layer():
    layer = _flat_layer(10, 12)
    normals = _compute_normals(layer)
    # All normals should point in +z or -z (consistent sign)
    assert np.allclose(np.abs(normals[:, :, 2]), 1.0, atol=1e-6)


def test_normals_sphere_radially_outward():
    layer = _sphere_layer(20, 24, radius=2.0)
    normals = _compute_normals(layer)
    # Dot product of normal with position vector should be ~1 (radially outward)
    pos = layer / np.linalg.norm(layer, axis=-1, keepdims=True)
    dot = np.einsum("ijk,ijk->ij", normals, pos)
    # Interior nodes only (boundary normals are one-sided — slightly less accurate)
    assert np.all(dot[2:-2, :] > 0.95)


def test_march_shape():
    surface = _sphere_layer(16, 20, radius=1.0)
    grid = march(surface, nk=10, ds=0.05, growth=1.1, smooth_strength=0.0)
    assert grid.shape == (16, 20, 10, 3)


def test_march_outer_radius_larger():
    surface = _sphere_layer(16, 20, radius=1.0)
    grid = march(surface, nk=10, ds=0.05, growth=1.1, smooth_strength=0.0)
    r_inner = np.linalg.norm(grid[:, :, 0,  :], axis=-1).mean()
    r_outer = np.linalg.norm(grid[:, :, -1, :], axis=-1).mean()
    assert r_outer > r_inner


def test_min_jacobian_detects_inversion():
    """min_jacobian_two_layers must return negative when layer_k1 is below layer_k."""
    from shore.volume.jacobian import min_jacobian_two_layers

    layer_k = _flat_layer(6, 8)
    # Move k+1 in the -z direction (inward): forces crossing
    layer_k1 = layer_k.copy()
    layer_k1[:, :, 2] = -0.1
    assert min_jacobian_two_layers(layer_k, layer_k1) < 0.0


def test_march_negative_jacobian_raises():
    """A surface with strongly converging normals must trigger the Jacobian check.

    Construct a saddle-like cap layer: rows at the i-boundaries are tilted
    inward so their normals converge toward the centre. With large ds the
    cells at those rows invert before the interior ones do.
    """
    surface = _sphere_layer(8, 10, radius=1.0)
    # Tilt the first and last latitude rows sharply inward by offsetting z
    # so the local normals point toward the axis — they will fold immediately.
    surface[0,  :, 2] *= -5.0
    surface[-1, :, 2] *= -5.0
    with pytest.raises(RuntimeError, match="Negative Jacobian"):
        march(surface, nk=3, ds=0.5, growth=1.0, smooth_strength=0.0)
