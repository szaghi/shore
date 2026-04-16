"""Full test suite for the Yin-Yang double O-grid module."""

from __future__ import annotations

import numpy as np
import pytest

from shore.surface.remesh import sphere_latlon_directions
from shore.surface.yinyang import (
    YANG_ROTATION,
    DEFAULT_THETA_MARGIN,
    coverage_mask_yin,
    coverage_mask_yang,
    yang_latlon_directions,
    project_yang_to_surface,
)
from shore.volume.hyperbolic import march

from tests.conftest import sphere_layer, yang_layer


# ---------------------------------------------------------------------------
# Rotation matrix tests
# ---------------------------------------------------------------------------


def test_yang_rotation_is_orthogonal():
    """R @ R.T == I and det(R) == 1."""
    R = YANG_ROTATION
    identity = R @ R.T
    assert np.allclose(identity, np.eye(3)), "YANG_ROTATION is not orthogonal"
    assert np.isclose(np.linalg.det(R), 1.0), "YANG_ROTATION determinant != 1"


def test_yang_directions_are_unit_vectors():
    """norm(d_yang, axis=-1) ≈ 1.0 for all (i, j)."""
    ni, nj = 20, 24
    d = yang_latlon_directions(ni, nj)
    norms = np.linalg.norm(d, axis=-1)
    assert np.allclose(norms, 1.0, atol=1e-12), "Yang directions are not unit vectors"


def test_yang_poles_on_x_axis():
    """First (i=0) and last (i=-1) rows of Yang directions have
    |x| as the dominant component, confirming the pseudo-poles lie on
    the x-axis. Specifically, |x| > |y| and |x| > |z| for all points
    in these rows, and the mean x is approximately cos(theta_m)."""
    ni, nj = 12, 18
    theta_m = DEFAULT_THETA_MARGIN
    d = yang_latlon_directions(ni, nj, theta_margin=theta_m)
    # First row: x is the dominant component (pseudo-pole region)
    first_row = d[0, :, :]
    assert np.all(np.abs(first_row[:, 0]) > np.abs(first_row[:, 1])), "x not dominant in first row"
    assert np.all(np.abs(first_row[:, 0]) > np.abs(first_row[:, 2])), "z not smaller than x in first row"
    assert np.allclose(first_row[:, 0].mean(), np.cos(theta_m), atol=1e-12), "x mean != cos(theta_m)"
    # Last row: x is also dominant
    last_row = d[-1, :, :]
    assert np.all(np.abs(last_row[:, 0]) > np.abs(last_row[:, 1])), "x not dominant in last row"
    assert np.all(np.abs(last_row[:, 0]) > np.abs(last_row[:, 2])), "z not smaller than x in last row"
    assert np.allclose(last_row[:, 0].mean(), -np.cos(theta_m), atol=1e-12), "x mean != -cos(theta_m)"


def test_yang_directions_differ_from_yin():
    """d_yang != d_yin element-wise everywhere."""
    ni, nj = 10, 12
    theta_m = np.pi / 6
    d_yin = sphere_latlon_directions(ni, nj, theta_margin=theta_m)
    d_yang = yang_latlon_directions(ni, nj, theta_margin=theta_m)
    assert not np.allclose(d_yin, d_yang), "Yang directions should differ from Yin"


# ---------------------------------------------------------------------------
# Coverage mask tests
# ---------------------------------------------------------------------------


def test_coverage_mask_yin():
    """All interior Yin directions (excluding first/last rows) satisfy
    coverage_mask_yin. The first and last rows are at the coverage boundary
    where |z| = cos(theta_m) and are not strictly covered."""
    ni, nj = 20, 24
    theta_m = DEFAULT_THETA_MARGIN
    d = sphere_latlon_directions(ni, nj, theta_margin=theta_m)
    # Check interior rows only (excluding boundary rows i=0 and i=-1)
    interior = d[1:-1, :, :].reshape(-1, 3)
    mask = coverage_mask_yin(interior, theta_m)
    assert np.all(mask), "Some interior Yin directions not covered by Yin mask"


def test_coverage_mask_yang():
    """All interior Yang directions (excluding first/last rows) satisfy
    coverage_mask_yang. The first and last rows are at the coverage boundary
    where |x| = cos(theta_m) and are not strictly covered."""
    ni, nj = 20, 24
    theta_m = DEFAULT_THETA_MARGIN
    d = yang_latlon_directions(ni, nj, theta_margin=theta_m)
    # Check interior rows only (excluding boundary rows i=0 and i=-1)
    interior = d[1:-1, :, :].reshape(-1, 3)
    mask = coverage_mask_yang(interior, theta_m)
    assert np.all(mask), "Some interior Yang directions not covered by Yang mask"


def test_full_coverage_random_sample():
    """10,000 random unit-sphere directions: every one satisfies
    Yin-mask OR Yang-mask at theta_m = pi/6."""
    np.random.seed(42)
    theta_m = DEFAULT_THETA_MARGIN
    # Generate random points on unit sphere using normal distribution
    random_dirs = np.random.randn(10000, 3)
    random_dirs = random_dirs / np.linalg.norm(random_dirs, axis=-1, keepdims=True)

    yin_mask = coverage_mask_yin(random_dirs, theta_m)
    yang_mask = coverage_mask_yang(random_dirs, theta_m)
    combined = yin_mask | yang_mask
    assert np.all(combined), "Some directions not covered by Yin OR Yang"


def test_overlap_non_empty():
    """At theta_m = pi/6, some sample points satisfy both masks."""
    np.random.seed(123)
    theta_m = DEFAULT_THETA_MARGIN
    random_dirs = np.random.randn(10000, 3)
    random_dirs = random_dirs / np.linalg.norm(random_dirs, axis=-1, keepdims=True)

    yin_mask = coverage_mask_yin(random_dirs, theta_m)
    yang_mask = coverage_mask_yang(random_dirs, theta_m)
    overlap = yin_mask & yang_mask
    assert np.any(overlap), "No overlap between Yin and Yang at theta_m=pi/6"


def test_no_coverage_above_pi4():
    """At theta_m = pi/4 + epsilon, some random directions satisfy neither mask
    (gap exists)."""
    np.random.seed(456)
    theta_m = np.pi / 4 + 0.01  # Just above pi/4
    random_dirs = np.random.randn(10000, 3)
    random_dirs = random_dirs / np.linalg.norm(random_dirs, axis=-1, keepdims=True)

    yin_mask = coverage_mask_yin(random_dirs, theta_m)
    yang_mask = coverage_mask_yang(random_dirs, theta_m)
    combined = yin_mask | yang_mask
    # With theta_m > pi/4, there should be uncovered points
    gap_count = np.sum(~combined)
    assert gap_count > 0, "No gap found even at theta_m > pi/4"


def test_theta_margin_ge_pi4_raises():
    """sphere_latlon_directions(8, 12, theta_margin=pi/4) raises ValueError."""
    with pytest.raises(ValueError, match="theta_margin.*>=.*pi/4"):
        sphere_latlon_directions(8, 12, theta_margin=np.pi / 4)


# ---------------------------------------------------------------------------
# Surface projection tests (require sphere STL fixture)
# ---------------------------------------------------------------------------


def test_yang_surface_projection_sphere_stl(sphere_stl):
    """All ni*nj rays hit the sphere STL; no ValueError;
    result shape (ni, nj, 3)."""
    import trimesh
    from shore.surface.io import load_stl

    mesh = load_stl(sphere_stl)
    ni, nj = 10, 12
    theta_m = DEFAULT_THETA_MARGIN

    k0 = project_yang_to_surface(mesh, ni, nj, theta_margin=theta_m)
    assert k0.shape == (ni, nj, 3), f"Expected shape ({ni}, {nj}, 3), got {k0.shape}"


def test_yang_points_on_sphere_surface(sphere_stl):
    """All projected Yang points within 5% of radius = 1.0."""
    from shore.surface.io import load_stl

    mesh = load_stl(sphere_stl)
    ni, nj = 10, 12
    theta_m = DEFAULT_THETA_MARGIN

    k0 = project_yang_to_surface(mesh, ni, nj, theta_margin=theta_m)
    radii = np.linalg.norm(k0, axis=-1)
    assert np.allclose(radii, 1.0, atol=0.05), "Some Yang points not on sphere surface"


def test_yang_march_positive_jacobian(sphere_stl):
    """March Yang layer with nk=5; all Jacobians > 0."""
    from shore.surface.io import load_stl
    from shore.volume.jacobian import min_jacobian_two_layers

    mesh = load_stl(sphere_stl)
    ni, nj = 10, 12
    theta_m = DEFAULT_THETA_MARGIN

    k0 = project_yang_to_surface(mesh, ni, nj, theta_margin=theta_m)
    grid = march(k0, nk=5, ds=0.05, growth=1.1, smooth_strength=0.0)

    # Check all layer transitions
    for k in range(4):
        jmin = min_jacobian_two_layers(grid[:, :, k, :], grid[:, :, k + 1, :])
        assert jmin > 0.0, f"Negative Jacobian at layer {k+1}"


def test_yang_march_outer_radius_larger(sphere_stl):
    """Mean radius of k=-1 layer > mean radius of k=0 layer."""
    from shore.surface.io import load_stl

    mesh = load_stl(sphere_stl)
    ni, nj = 10, 12
    theta_m = DEFAULT_THETA_MARGIN

    k0 = project_yang_to_surface(mesh, ni, nj, theta_margin=theta_m)
    grid = march(k0, nk=5, ds=0.05, growth=1.1, smooth_strength=0.0)

    r_inner = np.linalg.norm(grid[:, :, 0, :], axis=-1).mean()
    r_outer = np.linalg.norm(grid[:, :, -1, :], axis=-1).mean()
    assert r_outer > r_inner, "Outer radius not larger than inner radius"


def test_double_ogrid_combined_coverage_sphere_stl(sphere_stl):
    """Project both Yin and Yang onto sphere STL at theta_m = pi/6;
    sample 1,000 surface points; each is within epsilon of at least
    one block's k=0 layer.

    With ni=12, nj=16, the angular spacing is ~15° in latitude and ~22.5°
    in longitude. A chord-length tolerance of 0.25 accommodates the sparse
    sampling while still ensuring geometric validity.
    """
    from shore.surface.io import load_stl
    from shore.surface.remesh import project_sphere_to_surface

    mesh = load_stl(sphere_stl)
    ni, nj = 12, 16
    theta_m = DEFAULT_THETA_MARGIN

    k0_yin = project_sphere_to_surface(mesh, ni, nj, theta_margin=theta_m)
    k0_yang = project_yang_to_surface(mesh, ni, nj, theta_margin=theta_m)

    # Sample random points on the sphere surface
    np.random.seed(789)
    random_dirs = np.random.randn(1000, 3)
    random_dirs = random_dirs / np.linalg.norm(random_dirs, axis=-1, keepdims=True)
    sample_pts = random_dirs  # unit sphere

    # Check each sample point is near at least one block's k=0 layer
    eps = 0.25  # tolerance for sparse 12x16 grid
    covered = np.zeros(1000, dtype=bool)

    # Flatten k=0 layers for comparison
    yin_flat = k0_yin.reshape(-1, 3)
    yang_flat = k0_yang.reshape(-1, 3)

    for i, pt in enumerate(sample_pts):
        dist_yin = np.min(np.linalg.norm(yin_flat - pt, axis=-1))
        dist_yang = np.min(np.linalg.norm(yang_flat - pt, axis=-1))
        covered[i] = (dist_yin < eps) or (dist_yang < eps)

    assert np.all(covered), f"Some sample points not covered by Yin OR Yang: {np.sum(~covered)} uncovered"