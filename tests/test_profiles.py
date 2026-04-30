"""Tests for shore.profiles — NACA 4-digit / 5-digit + Selig I/O."""

from __future__ import annotations

import numpy as np
import pytest

from shore.exceptions import MeshInputError, ParameterError
from shore.profiles import naca, naca4, naca5, read_section, write_section

# ---------------------------------------------------------------------------
# NACA dispatcher
# ---------------------------------------------------------------------------


class TestNacaDispatcher:
    def test_4digit_dispatches(self):
        c = naca("0012", n_pts=20)
        assert c.shape == (39, 2)

    def test_5digit_dispatches(self):
        c = naca("23012", n_pts=20)
        assert c.shape == (39, 2)

    def test_invalid_designation_rejected(self):
        with pytest.raises(ParameterError, match="digit string"):
            naca("foo", n_pts=20)

    def test_unsupported_digit_count_rejected(self):
        with pytest.raises(ParameterError, match="4-digit and 5-digit"):
            naca("12", n_pts=20)
        with pytest.raises(ParameterError, match="4-digit and 5-digit"):
            naca("123456", n_pts=20)


# ---------------------------------------------------------------------------
# Output shape and ordering
# ---------------------------------------------------------------------------


class TestSeligOrdering:
    @pytest.mark.parametrize("n_pts", [10, 20, 80, 200])
    def test_total_length_is_2n_minus_1(self, n_pts):
        c = naca("0012", n_pts=n_pts)
        assert c.shape == (2 * n_pts - 1, 2)

    def test_le_at_midpoint(self):
        c = naca("0012", n_pts=80)
        # LE sits at index n_pts - 1 (= 79 for n_pts=80).
        le = c[79]
        assert np.allclose(le, [0.0, 0.0], atol=1e-12)

    def test_upper_te_first_lower_te_last(self):
        c = naca("2412", n_pts=80)
        # Upper TE (i=0) sits at x=1, y=cambered upper TE.
        # Lower TE (i=-1) sits at x=1, y=cambered lower TE.
        assert np.isclose(c[0, 0], 1.0)
        assert np.isclose(c[-1, 0], 1.0)
        # For closed_te=True, both TE points coincide.
        assert np.allclose(c[0], c[-1], atol=1e-12)

    def test_upper_surface_has_positive_y_for_symmetric(self):
        c = naca("0012", n_pts=20)
        upper = c[: 20 - 1]  # TE -> LE-1
        # Symmetric airfoil: all upper-surface y are non-negative.
        assert np.all(upper[:, 1] >= -1e-12)

    def test_walks_from_te_to_le_to_te(self):
        c = naca("0012", n_pts=80)
        # Upper half: x decreases monotonically from 1 to 0.
        upper_x = c[:80, 0]
        assert np.all(np.diff(upper_x) <= 1e-12)
        # Lower half: x increases monotonically from 0 to 1.
        lower_x = c[79:, 0]
        assert np.all(np.diff(lower_x) >= -1e-12)


# ---------------------------------------------------------------------------
# 4-digit family
# ---------------------------------------------------------------------------


class TestNaca4:
    def test_naca0012_max_thickness(self):
        # Standard NACA 4-digit max thickness occurs at x = 0.3.
        # Verified against published 12% airfoil ordinate: y ~= 0.060.
        c = naca4("0012", n_pts=200, closed_te=False)
        upper = c[:200]
        i_30 = int(np.argmin(np.abs(upper[:, 0] - 0.30)))
        assert abs(upper[i_30, 1] - 0.060) < 1e-3

    def test_naca0012_le_radius_zero_at_origin(self):
        c = naca4("0012", n_pts=400)
        # LE sits at (0, 0).
        assert np.allclose(c[399], [0.0, 0.0], atol=1e-12)

    def test_naca2412_camber(self):
        # NACA 2412: 2% camber, position p = 0.4.  At x = 0.4 the
        # camber line peaks at y_c = 0.02.
        c = naca4("2412", n_pts=400)
        upper = c[:400]
        lower = c[399:]
        # Find y_c at x = 0.4 by averaging upper + lower y at matching x.
        i_upper_40 = int(np.argmin(np.abs(upper[:, 0] - 0.40)))
        i_lower_40 = int(np.argmin(np.abs(lower[:, 0] - 0.40)))
        y_c = 0.5 * (upper[i_upper_40, 1] + lower[i_lower_40, 1])
        assert abs(y_c - 0.020) < 1e-3

    def test_closed_te_true_closes_exactly(self):
        c = naca4("0012", n_pts=80, closed_te=True)
        assert np.allclose(c[0], c[-1], atol=1e-10)

    def test_closed_te_false_leaves_gap(self):
        c = naca4("0012", n_pts=80, closed_te=False)
        # Open TE: ~0.25% chord gap.
        gap = abs(c[0, 1] - c[-1, 1])
        assert gap > 1e-3

    def test_invalid_designation_rejected(self):
        with pytest.raises(ParameterError, match="4 digits"):
            naca4("00123", n_pts=20)
        with pytest.raises(ParameterError, match="digit string"):
            naca4("abcd", n_pts=20)

    def test_n_pts_too_small_rejected(self):
        with pytest.raises(ParameterError, match="n_pts"):
            naca4("0012", n_pts=4)


# ---------------------------------------------------------------------------
# 5-digit family
# ---------------------------------------------------------------------------


class TestNaca5:
    def test_naca23012_thickness_matches_4digit(self):
        # 5-digit family shares the 4-digit thickness polynomial.
        # NACA 23012 has 12% thickness like NACA 0012.
        c = naca5("23012", n_pts=400, closed_te=False)
        ymax = c[:, 1].max()
        ymin = c[:, 1].min()
        thickness = ymax - ymin
        # Camber adds a small offset to the gross thickness; allow ~1% margin.
        assert 0.115 < thickness < 0.135

    def test_naca23012_camber_positive_forward(self):
        # Standard NACA 230 mean line has positive camber forward of
        # x ~= 0.2 and negative aft of it.  Check the camber peak is
        # forward of the chord midpoint.
        c = naca5("23012", n_pts=400)
        upper = c[:400]
        lower = c[399:]
        # Camber line ~ (y_upper + y_lower) / 2 at matched x.
        # Use upper x as reference; find lower-y at same x via interp.
        x_ref = upper[:, 0]
        y_lower_interp = np.interp(x_ref[::-1], lower[:, 0], lower[:, 1])
        y_camber = 0.5 * (upper[::-1, 1] + y_lower_interp)
        # Peak camber x.
        x_peak = x_ref[::-1][int(np.argmax(y_camber))]
        assert 0.10 < x_peak < 0.30

    def test_unsupported_L_rejected(self):
        # Only L=2 (CL_design = 0.3) is implemented.
        with pytest.raises(ParameterError, match="L=2"):
            naca5("33012", n_pts=20)

    def test_unsupported_S_rejected(self):
        with pytest.raises(ParameterError, match="camber type"):
            naca5("22212", n_pts=20)  # L=2, Q=2, S=2 -> bad S

    def test_unsupported_Q_rejected(self):
        # Q=0 is not in the standard table.
        with pytest.raises(ParameterError, match="mean-line code"):
            naca5("20012", n_pts=20)

    def test_reflexed_family(self):
        # NACA 23112 (reflexed).
        c = naca5("23112", n_pts=80)
        assert c.shape == (159, 2)


# ---------------------------------------------------------------------------
# Selig .dat I/O
# ---------------------------------------------------------------------------


class TestSeligIO:
    def test_roundtrip(self, tmp_path):
        c = naca("0012", n_pts=40)
        path = tmp_path / "naca0012.dat"
        write_section(path, c, name="NACA 0012")
        loaded = read_section(path)
        np.testing.assert_allclose(loaded, c, atol=1e-7)

    def test_no_name_line(self, tmp_path):
        path = tmp_path / "noname.dat"
        path.write_text("0.0 0.0\n0.5 0.1\n1.0 0.0\n0.5 -0.1\n")
        loaded = read_section(path)
        assert loaded.shape == (4, 2)
        assert np.allclose(loaded[1], [0.5, 0.1])

    def test_skips_name_line(self, tmp_path):
        path = tmp_path / "named.dat"
        path.write_text("My Airfoil\n0.0 0.0\n0.5 0.1\n1.0 0.0\n0.5 -0.1\n")
        loaded = read_section(path)
        assert loaded.shape == (4, 2)

    def test_skips_blank_lines(self, tmp_path):
        path = tmp_path / "blanks.dat"
        path.write_text("\n0.0 0.0\n\n0.5 0.1\n1.0 0.0\n\n")
        loaded = read_section(path)
        assert loaded.shape == (3, 2)

    def test_missing_file_rejected(self, tmp_path):
        with pytest.raises(MeshInputError, match="file not found"):
            read_section(tmp_path / "ghost.dat")

    def test_empty_file_rejected(self, tmp_path):
        path = tmp_path / "empty.dat"
        path.write_text("")
        with pytest.raises(MeshInputError, match="no coordinate rows"):
            read_section(path)

    def test_malformed_midfile_rejected(self, tmp_path):
        path = tmp_path / "bad.dat"
        path.write_text("0.0 0.0\nnot a row\n1.0 0.0\n")
        with pytest.raises(MeshInputError, match=r"malformed|non-numeric"):
            read_section(path)

    def test_write_rejects_bad_shape(self, tmp_path):
        with pytest.raises(ParameterError, match="shape"):
            write_section(tmp_path / "bad.dat", np.zeros((10,)))

    def test_write_rejects_short_curve(self, tmp_path):
        with pytest.raises(ParameterError, match="at least 4 rows"):
            write_section(tmp_path / "short.dat", np.zeros((3, 2)))
