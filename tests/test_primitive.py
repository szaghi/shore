"""Tests for the parametric primitive meshes (box, annulus, flat-caps)."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import ParameterError
from shore.primitive import (
    AnnularCylinderMesh,
    BoxMesh,
    FlatCapsCylinder,
    axis_frame,
    coords_1d,
    resolve_axis,
)
from shore.volume.jacobian import _jacobian_field_open, jacobian_field

# ---------------------------------------------------------------------------
# resolve_axis / axis_frame / coords_1d
# ---------------------------------------------------------------------------


class TestResolveAxis:
    @pytest.mark.parametrize(
        "spec, expected",
        [
            ("x", [1.0, 0.0, 0.0]),
            ("y", [0.0, 1.0, 0.0]),
            ("z", [0.0, 0.0, 1.0]),
            ("X", [1.0, 0.0, 0.0]),
            (" Z ", [0.0, 0.0, 1.0]),
        ],
    )
    def test_label(self, spec, expected):
        assert np.allclose(resolve_axis(spec), expected)

    def test_vector_string(self):
        assert np.allclose(resolve_axis("1 0 0"), [1, 0, 0])

    def test_vector_normalised(self):
        # "2 0 0" → "1 0 0" after normalisation.
        assert np.allclose(resolve_axis((2.0, 0.0, 0.0)), [1, 0, 0])

    def test_vector_oblique(self):
        v = resolve_axis((1.0, 1.0, 1.0))
        assert np.isclose(np.linalg.norm(v), 1.0)
        assert np.allclose(v * np.sqrt(3), [1, 1, 1])

    def test_zero_vector_rejected(self):
        with pytest.raises(ValueError, match="zero-magnitude"):
            resolve_axis((0.0, 0.0, 0.0))

    def test_unparseable_string(self):
        with pytest.raises(ValueError):
            resolve_axis("nope")

    def test_wrong_shape(self):
        with pytest.raises(ValueError):
            resolve_axis((1.0, 0.0))


class TestAxisFrame:
    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_orthonormal_right_handed(self, axis):
        e1, e2, e3 = axis_frame(resolve_axis(axis))
        for v in (e1, e2, e3):
            assert np.isclose(np.linalg.norm(v), 1.0)
        # Pairwise orthogonality
        assert abs(np.dot(e1, e2)) < 1e-12
        assert abs(np.dot(e2, e3)) < 1e-12
        assert abs(np.dot(e1, e3)) < 1e-12
        # Right-handed: e1 x e2 == e3
        assert np.allclose(np.cross(e1, e2), e3)


class TestCoords1D:
    def test_uniform(self):
        c = coords_1d(5, 0.0, 1.0)
        assert np.allclose(c, [0.0, 0.25, 0.5, 0.75, 1.0])

    def test_endpoints_are_pinned_for_every_law(self):
        for law in ("uniform", "tanh", "tanh2"):
            c = coords_1d(20, 2.0, 5.0, law=law, beta=3.0)
            assert np.isclose(c[0], 2.0)
            assert np.isclose(c[-1], 5.0)

    def test_unknown_law_rejected(self):
        with pytest.raises(ValueError, match="unknown law"):
            coords_1d(10, 0.0, 1.0, law="bogus")

    def test_too_few_points(self):
        with pytest.raises(ValueError, match=">= 2"):
            coords_1d(1, 0.0, 1.0)


# ---------------------------------------------------------------------------
# BoxMesh
# ---------------------------------------------------------------------------


class TestBoxMesh:
    def test_shape(self):
        box = BoxMesh.build(vmin=(0, 0, 0), vmax=(10, 5, 3), ni=12, nj=8, nk=6)
        assert box.blocks[0].nodes.shape == (12, 8, 6, 3)

    def test_corners(self):
        box = BoxMesh.build(vmin=(1, 2, 3), vmax=(4, 5, 6), ni=5, nj=4, nk=3)
        n = box.blocks[0].nodes
        assert np.allclose(n[0, 0, 0], [1, 2, 3])
        assert np.allclose(n[-1, -1, -1], [4, 5, 6])

    def test_uniform_spacing_is_uniform(self):
        box = BoxMesh.build(vmin=(0, 0, 0), vmax=(1, 1, 1), ni=11, nj=11, nk=11)
        n = box.blocks[0].nodes
        # x-coord along i should be linspace(0, 1, 11)
        assert np.allclose(n[:, 0, 0, 0], np.linspace(0, 1, 11))
        assert np.allclose(n[0, :, 0, 1], np.linspace(0, 1, 11))
        assert np.allclose(n[0, 0, :, 2], np.linspace(0, 1, 11))

    def test_tanh_clusters_near_start(self):
        box = BoxMesh.build(
            vmin=(0, 0, 0),
            vmax=(1, 1, 1),
            ni=21,
            nj=2,
            nk=2,
            i_law="tanh",
            i_beta=4.0,
        )
        x = box.blocks[0].nodes[:, 0, 0, 0]
        # First cell smaller than last cell when clustering at start.
        d = np.diff(x)
        assert d[0] < d[-1]

    def test_jacobian_positive(self):
        box = BoxMesh.build(vmin=(0, 0, 0), vmax=(2, 3, 4), ni=8, nj=6, nk=5)
        j = _jacobian_field_open(box.blocks[0].nodes)
        assert j.min() > 0

    def test_no_periodicity(self):
        box = BoxMesh.build(vmin=(0, 0, 0), vmax=(1, 1, 1), ni=4, nj=4, nk=4)
        assert box.blocks[0].periodic_i is False
        assert box.blocks[0].periodic_j is False

    def test_min_max_strict_inequality(self):
        with pytest.raises(ParameterError, match="vmax must exceed vmin"):
            BoxMesh.build(vmin=(0, 0, 0), vmax=(0, 5, 3), ni=4, nj=4, nk=4)

    def test_too_few_points(self):
        with pytest.raises(ParameterError, match="ni must be >= 2"):
            BoxMesh.build(vmin=(0, 0, 0), vmax=(1, 1, 1), ni=1, nj=4, nk=4)


# ---------------------------------------------------------------------------
# AnnularCylinderMesh
# ---------------------------------------------------------------------------


class TestAnnulusFull:
    def test_shape_and_periodicity(self):
        ann = AnnularCylinderMesh.build(r_in=0.5, r_out=2.0, z0=0, z1=4, ni=10, nj=32, nk=8)
        assert ann.blocks[0].nodes.shape == (10, 32, 8, 3)
        assert ann.blocks[0].periodic_j is True

    def test_inner_outer_radii_correct(self):
        ann = AnnularCylinderMesh.build(r_in=0.5, r_out=2.0, z0=0, z1=4, ni=10, nj=32, nk=8)
        n = ann.blocks[0].nodes[:, :, 0, :]  # one z slice
        r = np.linalg.norm(n[..., :2], axis=-1)
        assert np.allclose(r[0, :], 0.5)
        assert np.allclose(r[-1, :], 2.0)

    def test_axial_extent_along_z(self):
        ann = AnnularCylinderMesh.build(r_in=0.5, r_out=2.0, z0=1.0, z1=5.0, ni=4, nj=16, nk=10)
        z = ann.blocks[0].nodes[0, 0, :, 2]
        assert np.isclose(z[0], 1.0)
        assert np.isclose(z[-1], 5.0)

    def test_full_annulus_no_duplicate_node(self):
        # Periodic full annulus drops the wrap-around node.
        ann = AnnularCylinderMesh.build(r_in=0.5, r_out=2.0, z0=0, z1=4, ni=4, nj=16, nk=4)
        n = ann.blocks[0].nodes
        # First and last azimuthal nodes should be one cell apart, not coincident.
        gap = np.linalg.norm(n[0, 0, 0] - n[0, -1, 0])
        single_cell = np.linalg.norm(n[0, 0, 0] - n[0, 1, 0])
        assert np.isclose(gap, single_cell)

    def test_jacobian_positive(self):
        ann = AnnularCylinderMesh.build(r_in=0.5, r_out=2.0, z0=0, z1=4, ni=10, nj=32, nk=8)
        j = jacobian_field(ann.blocks[0].nodes)
        assert j.min() > 0

    def test_full_annulus_rejects_nonuniform_j(self):
        with pytest.raises(ParameterError, match="j_law must be 'uniform'"):
            AnnularCylinderMesh.build(
                r_in=0.5,
                r_out=2.0,
                z0=0,
                z1=4,
                ni=10,
                nj=32,
                nk=8,
                j_law="tanh2",
            )

    def test_axis_x_orientation(self):
        ann = AnnularCylinderMesh.build(
            r_in=0.5, r_out=2.0, z0=0, z1=4, ni=4, nj=16, nk=10, axis="x"
        )
        # Axial coordinate (k direction) should now run along x, not z.
        axial = ann.blocks[0].nodes[0, 0, :, :]
        assert np.isclose(axial[0, 0], 0.0)
        assert np.isclose(axial[-1, 0], 4.0)
        # The y / z components along k should be (close to) constant.
        assert np.std(axial[:, 1]) < 1e-12

    def test_axis_arbitrary_vector(self):
        ann = AnnularCylinderMesh.build(
            r_in=0.5,
            r_out=2.0,
            z0=0,
            z1=2,
            ni=4,
            nj=16,
            nk=5,
            axis=(1.0, 1.0, 0.0),
        )
        # Axial direction is along (1,1,0)/sqrt(2); verify the axial
        # spacing between successive k layers is constant magnitude.
        n = ann.blocks[0].nodes
        delta = n[0, 0, 1] - n[0, 0, 0]
        assert np.isclose(np.linalg.norm(delta), 0.5)
        # delta should be parallel to (1,1,0)/sqrt(2)
        unit = delta / np.linalg.norm(delta)
        assert np.allclose(unit, [1 / np.sqrt(2), 1 / np.sqrt(2), 0.0])


class TestAnnulusSector:
    def test_partial_sector_not_periodic(self):
        ann = AnnularCylinderMesh.build(
            r_in=0.5,
            r_out=2.0,
            z0=0,
            z1=4,
            ni=4,
            nj=20,
            nk=4,
            theta_start=0.0,
            theta_end=np.pi,
        )
        assert ann.blocks[0].periodic_j is False

    def test_sector_endpoints(self):
        ann = AnnularCylinderMesh.build(
            r_in=1.0,
            r_out=2.0,
            z0=0,
            z1=1,
            ni=2,
            nj=11,
            nk=2,
            theta_start=0.0,
            theta_end=np.pi / 2,
        )
        n = ann.blocks[0].nodes
        # j=0 lies on the +x axis; j=-1 lies on the +y axis.
        assert np.allclose(n[0, 0, 0, :2], [1.0, 0.0])
        assert np.allclose(n[0, -1, 0, :2], [0.0, 1.0], atol=1e-12)

    def test_sector_rejects_inverted_range(self):
        with pytest.raises(ParameterError, match=r"theta_end .* must exceed"):
            AnnularCylinderMesh.build(
                r_in=0.5,
                r_out=2.0,
                z0=0,
                z1=1,
                ni=2,
                nj=10,
                nk=2,
                theta_start=1.0,
                theta_end=0.5,
            )

    def test_sector_jacobian_positive(self):
        ann = AnnularCylinderMesh.build(
            r_in=0.5,
            r_out=2.0,
            z0=0,
            z1=4,
            ni=10,
            nj=20,
            nk=8,
            theta_end=np.pi,
        )
        j = _jacobian_field_open(ann.blocks[0].nodes)
        assert j.min() > 0


# ---------------------------------------------------------------------------
# FlatCapsCylinder
# ---------------------------------------------------------------------------


class TestFlatCaps:
    def test_block_count_and_labels(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=6, nc=8, nk=3)
        assert len(fc.blocks) == 5
        assert [b.label for b in fc.blocks] == ["center", "sub_e", "sub_n", "sub_w", "sub_s"]

    def test_shapes(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=6, nc=8, nk=3)
        assert fc.blocks[0].nodes.shape == (8, 8, 3, 3)
        for sub in fc.blocks[1:]:
            assert sub.nodes.shape == (6, 8, 3, 3)

    def test_outer_edges_lie_on_circle(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=6, nc=8, nk=3)
        for sub in fc.blocks[1:]:
            outer = sub.nodes[-1, :, 0, :2]  # (nc, 2) at first z slice
            radii = np.linalg.norm(outer, axis=-1)
            assert np.allclose(radii, 2.0)

    def test_center_is_square_and_axis_aligned(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=6, nc=8, nk=3)
        c = fc.blocks[0].nodes[:, :, 0, :2]
        # With inner_frac=0.5 (default), corners sit at radius
        # ``r_inner = inner_frac * r_out = 1.0`` at azimuths
        # ``pi/4 + k*pi/2``, so each side has length
        # ``s = sqrt(2) * r_inner = sqrt(2)``.
        side = np.sqrt(2.0) * 1.0  # r_inner = 0.5 * 2.0
        assert np.isclose(c[-1, 0, 0] - c[0, 0, 0], side)
        assert np.isclose(c[0, -1, 1] - c[0, 0, 1], side)
        # And the four corners of the central square sit at radius r_inner.
        for corner in (c[0, 0], c[-1, 0], c[-1, -1], c[0, -1]):
            assert np.isclose(np.linalg.norm(corner), 1.0)

    def test_inner_edges_match_central_square(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=6, nc=8, nk=3)
        center = fc.blocks[0]
        sub_e, sub_n, sub_w, sub_s = fc.blocks[1:]
        assert np.allclose(sub_e.nodes[0, :, :], center.nodes[-1, :, :])
        assert np.allclose(sub_n.nodes[0, :, :], center.nodes[::-1, -1, :])
        assert np.allclose(sub_w.nodes[0, :, :], center.nodes[0, ::-1, :])
        assert np.allclose(sub_s.nodes[0, :, :], center.nodes[:, 0, :])

    def test_corner_seams_match(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=6, nc=8, nk=3)
        sub_e, sub_n, sub_w, sub_s = fc.blocks[1:]
        # At each square corner the radial line is shared by two sub-blocks.
        assert np.allclose(sub_e.nodes[:, -1, :], sub_n.nodes[:, 0, :])
        assert np.allclose(sub_n.nodes[:, -1, :], sub_w.nodes[:, 0, :])
        assert np.allclose(sub_w.nodes[:, -1, :], sub_s.nodes[:, 0, :])
        assert np.allclose(sub_s.nodes[:, -1, :], sub_e.nodes[:, 0, :])

    def test_jacobian_nonnegative_on_each_block(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=8, nc=10, nk=5)
        for blk in fc.blocks:
            j = _jacobian_field_open(blk.nodes)
            # Sub-block corner cells where i=0 meets the central square
            # corner are near-zero by construction (concentric pole
            # artefact); the bulk of the field must be strictly positive.
            assert j.max() > 0
            # Allow FP-noise negatives at the four corner cells.
            assert j.min() > -1e-12

    def test_axial_extent(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=1.0, z1=5.0, ni=4, nc=4, nk=6)
        for blk in fc.blocks:
            z = blk.nodes[0, 0, :, 2]
            assert np.isclose(z[0], 1.0)
            assert np.isclose(z[-1], 5.0)

    def test_axis_x_orientation(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=4, nc=4, nk=5, axis="x")
        # Axial direction now runs along x.
        n = fc.blocks[0].nodes
        axial = n[0, 0, :, :]
        assert np.isclose(axial[0, 0], 0.0)
        assert np.isclose(axial[-1, 0], 4.0)
        assert np.std(axial[:, 1]) < 1e-12

    def test_too_small_dim_rejected(self):
        with pytest.raises(ParameterError, match="ni must be >= 2"):
            FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=1, nc=4, nk=4)

    def test_inner_frac_default_is_half(self):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=4, nc=6, nk=3)
        c = fc.blocks[0].nodes[:, :, 0, :2]
        # All four corners at radius inner_frac * r_out = 0.5 * 2.0 = 1.0.
        for corner in (c[0, 0], c[-1, 0], c[-1, -1], c[0, -1]):
            assert np.isclose(np.linalg.norm(corner), 1.0)

    @pytest.mark.parametrize("frac", [0.25, 0.5, 0.7])
    def test_inner_frac_controls_square_size(self, frac):
        r_out = 2.0
        fc = FlatCapsCylinder.build(r_out=r_out, z0=0, z1=1, ni=4, nc=6, nk=2, inner_frac=frac)
        c = fc.blocks[0].nodes[:, :, 0, :2]
        for corner in (c[0, 0], c[-1, 0], c[-1, -1], c[0, -1]):
            assert np.isclose(np.linalg.norm(corner), frac * r_out)

    def test_inner_frac_out_of_range_rejected(self):
        with pytest.raises(ParameterError, match="inner_frac must be in"):
            FlatCapsCylinder.build(r_out=2.0, z0=0, z1=1, ni=4, nc=4, nk=2, inner_frac=1.0)
        with pytest.raises(ParameterError, match="inner_frac must be in"):
            FlatCapsCylinder.build(r_out=2.0, z0=0, z1=1, ni=4, nc=4, nk=2, inner_frac=0.0)

    def test_radial_first_step_pinned_per_j_uniform(self):
        """The sub-block first radial cell at column j must equal the
        central square's local tangential cell at column j (aspect-ratio-1
        cells at the seam, even under non-uniform c_law)."""
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=1, ni=12, nc=14, nk=2, inner_frac=0.5)
        center = fc.blocks[0]
        sub_e = fc.blocks[1]
        nc = center.nodes.shape[0]

        for j in range(nc):
            sub_radial = np.linalg.norm(sub_e.nodes[1, j, 0, :] - sub_e.nodes[0, j, 0, :])
            # Local central tangential cell size at column j (centred diff).
            if j == 0:
                tang = np.linalg.norm(center.nodes[-1, 1, 0] - center.nodes[-1, 0, 0])
            elif j == nc - 1:
                tang = np.linalg.norm(center.nodes[-1, -1, 0] - center.nodes[-1, -2, 0])
            else:
                tang = np.linalg.norm(center.nodes[-1, j + 1, 0] - center.nodes[-1, j - 1, 0]) / 2.0
            assert np.isclose(sub_radial, tang, rtol=1e-12), (
                f"j={j}: sub_radial={sub_radial}, tang={tang}"
            )

    def test_radial_first_step_pinned_per_j_tanh2(self):
        """Per-j C1 must hold even when ``c_law`` clusters non-uniformly."""
        fc = FlatCapsCylinder.build(
            r_out=2.0,
            z0=0,
            z1=1,
            ni=12,
            nc=14,
            nk=2,
            inner_frac=0.5,
            c_law="tanh2",
            c_beta=3.0,
        )
        center = fc.blocks[0]
        sub_e = fc.blocks[1]
        nc = center.nodes.shape[0]
        for j in range(1, nc - 1):  # interior columns; centred diff well-defined
            sub_radial = np.linalg.norm(sub_e.nodes[1, j, 0] - sub_e.nodes[0, j, 0])
            tang = np.linalg.norm(center.nodes[-1, j + 1, 0] - center.nodes[-1, j - 1, 0]) / 2.0
            assert np.isclose(sub_radial, tang, rtol=1e-12)

    def test_sub_block_columns_are_straight_radial_lines(self):
        """Each j-column of a sub-block runs from inner[j] outward
        along the line ``inner[j] -> 0``: the unit-vector should be
        constant along i."""
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=1, ni=8, nc=10, nk=2, inner_frac=0.5)
        sub_e = fc.blocks[1]
        for j in range(sub_e.nodes.shape[1]):
            column = sub_e.nodes[:, j, 0, :2]  # (ni, 2)
            radii = np.linalg.norm(column, axis=-1, keepdims=True)
            units = column / radii
            # Every i along this column should have the same unit
            # vector (straight radial line through the origin).
            assert np.allclose(units, units[0:1], atol=1e-12), f"column j={j} not radial"


# ---------------------------------------------------------------------------
# I/O — every primitive uses the inherited write_geo
# ---------------------------------------------------------------------------


class TestPrimitiveIO:
    def test_box_writes_single_geo(self, tmp_path):
        box = BoxMesh.build(vmin=(0, 0, 0), vmax=(1, 1, 1), ni=4, nj=4, nk=4)
        out = tmp_path / "bg.geo"
        box.write_geo(out)
        assert out.exists()

    def test_box_writes_geo_when_extension_omitted(self, tmp_path):
        box = BoxMesh.build(vmin=(0, 0, 0), vmax=(1, 1, 1), ni=4, nj=4, nk=4)
        out = tmp_path / "bg"
        box.write_geo(out)
        assert (tmp_path / "bg.geo").exists()

    def test_annulus_writes_single_geo(self, tmp_path):
        ann = AnnularCylinderMesh.build(r_in=0.5, r_out=2.0, z0=0, z1=4, ni=4, nj=16, nk=4)
        out = tmp_path / "ann.geo"
        ann.write_geo(out)
        assert out.exists()

    def test_flat_caps_writes_five_files(self, tmp_path):
        fc = FlatCapsCylinder.build(r_out=2.0, z0=0, z1=4, ni=4, nc=4, nk=3)
        fc.write_geo(tmp_path / "bg")
        for label in ("center", "sub_e", "sub_n", "sub_w", "sub_s"):
            assert (tmp_path / f"bg_{label}.geo").exists()
