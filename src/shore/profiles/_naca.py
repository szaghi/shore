"""NACA 4-digit and 5-digit airfoil section coordinates.

Closed-form polynomial parameterisations.  References:

- 4-digit:   Abbott & von Doenhoff, *Theory of Wing Sections*, §6.5.
- 5-digit:   Jacobs & Pinkerton, NACA TR-537 (1935); reflexed family
             extended in NACA TR-610 (1937).

The 4-digit family is described by three parameters: max camber
``m`` (% chord, first digit / 100), camber position ``p`` (tenths of
chord, second digit / 10), and max thickness ``t`` (% chord, last
two digits / 100).  The 5-digit family uses the same thickness
distribution but a different camber line, parameterised by a design
lift coefficient and a camber line shape.
"""

from __future__ import annotations

import numpy as np

from shore.exceptions import ParameterError

# Standard NACA 4-digit thickness polynomial coefficients.  The last
# coefficient (``-0.1015``) gives an open trailing edge with ~0.1%
# chord thickness.  Replacing it with ``-0.1036`` gives an exactly
# closed TE; we apply that adjustment when ``closed_te=True``.
_T4_COEFS_OPEN = (0.2969, -0.1260, -0.3516, 0.2843, -0.1015)
_T4_COEFS_CLOSED = (0.2969, -0.1260, -0.3516, 0.2843, -0.1036)


def naca(designation: str, n_pts: int = 80, *, closed_te: bool = True) -> np.ndarray:
    """Build a NACA airfoil section by digit-count dispatch.

    Parameters
    ----------
    designation :
        Digit string: 4 chars (4-digit family, e.g. ``"0012"``,
        ``"2412"``) or 5 chars (5-digit family, e.g. ``"23012"``,
        ``"23112"`` reflexed).
    n_pts :
        Number of chord samples per side (LE..TE).  Must be ``>= 5``.
        The returned curve has ``2 * n_pts - 1`` points around the
        full section: ``n_pts`` on the upper surface plus ``n_pts``
        on the lower surface, with the LE node shared.  Cosine-spaced
        along x (denser at LE and TE).
    closed_te :
        If ``True`` (default), apply the small thickness-polynomial
        adjustment that closes the TE exactly.  If ``False``, use
        the standard open-TE polynomial whose residual is ``≈ 0.1%
        chord`` at the TE — useful when you want to detect a blunt
        TE downstream.

    Returns
    -------
    numpy.ndarray, shape ``(2 * n_pts - 1, 2)``
        Selig-ordered coordinates: upper-TE -> LE -> lower-TE.

    Raises
    ------
    ParameterError
        On invalid designation or ``n_pts``.
    """
    if not isinstance(designation, str) or not designation.isdigit():
        raise ParameterError(
            f"naca: designation must be a digit string (e.g. '0012', '23012'), got {designation!r}"
        )
    if len(designation) == 4:
        return naca4(designation, n_pts, closed_te=closed_te)
    if len(designation) == 5:
        return naca5(designation, n_pts, closed_te=closed_te)
    raise ParameterError(
        f"naca: only 4-digit and 5-digit designations are supported "
        f"(got {len(designation)}-digit {designation!r})"
    )


# ---------------------------------------------------------------------------
# 4-digit family
# ---------------------------------------------------------------------------


def naca4(designation: str, n_pts: int = 80, *, closed_te: bool = True) -> np.ndarray:
    """Build a NACA 4-digit airfoil section.

    Parameters
    ----------
    designation :
        4 digits ``MPTT``: ``M`` = max camber as % chord (first digit),
        ``P`` = camber position in tenths of chord (second digit),
        ``TT`` = max thickness as % chord (last two digits).
        Example: ``"2412"`` = 2% camber at 40% chord, 12% thickness.
    n_pts :
        Even total node count, ``>= 8``.
    closed_te :
        Whether to use the closed-TE thickness polynomial.

    Returns
    -------
    numpy.ndarray, shape ``(n_pts, 2)``
        Selig-ordered coordinates.

    Raises
    ------
    ParameterError
        On invalid designation or ``n_pts``.
    """
    digits = _validate_designation(designation, expected_len=4)
    if n_pts < 5:
        raise ParameterError(f"naca4: n_pts must be >= 5, got {n_pts}")

    m = int(digits[0]) / 100.0
    p = int(digits[1]) / 10.0
    t = int(digits[2:4]) / 100.0

    x = _cosine_x(n_pts)  # (n_half,) from 0 to 1
    yt = _thickness4(x, t, closed_te=closed_te)
    yc, dyc_dx = _camber4(x, m, p)

    return _assemble_selig(x, yt, yc, dyc_dx)


def _thickness4(x: np.ndarray, t: float, *, closed_te: bool) -> np.ndarray:
    """NACA 4-digit thickness polynomial ``y_t(x)``."""
    a = _T4_COEFS_CLOSED if closed_te else _T4_COEFS_OPEN
    return 5.0 * t * (a[0] * np.sqrt(x) + a[1] * x + a[2] * x**2 + a[3] * x**3 + a[4] * x**4)


def _camber4(x: np.ndarray, m: float, p: float) -> tuple[np.ndarray, np.ndarray]:
    """NACA 4-digit camber line ``y_c(x)`` and slope ``dy_c/dx``."""
    if m == 0.0 or p == 0.0:
        return np.zeros_like(x), np.zeros_like(x)
    yc = np.where(
        x < p,
        (m / p**2) * (2.0 * p * x - x**2),
        (m / (1.0 - p) ** 2) * ((1.0 - 2.0 * p) + 2.0 * p * x - x**2),
    )
    dyc = np.where(
        x < p,
        (2.0 * m / p**2) * (p - x),
        (2.0 * m / (1.0 - p) ** 2) * (p - x),
    )
    return yc, dyc


# ---------------------------------------------------------------------------
# 5-digit family
# ---------------------------------------------------------------------------


# Tabulated coefficients for the 5-digit camber line from NACA TR-537
# (standard, ``S = 0``) and TR-610 (reflexed, ``S = 1``).  Indexed by
# the second digit Q (mean-line designation: 1..5 mapped to specific
# CL_design / xcm values; with the third digit S indicating standard
# vs reflexed camber).  Source: Abbott & von Doenhoff, App. III.
_NACA5_STANDARD = {
    # Q: (m, k1)   for the camber line  y_c = (k1/6) * (x^3 - 3 m x^2 + m^2 (3-m) x)   x <  m
    #                                       (k1 m^3 / 6) * (1 - x)                     x >= m
    # m is the polynomial root, *not* the design CL position.
    1: (0.0580, 361.40),
    2: (0.1260, 51.640),
    3: (0.2025, 15.957),
    4: (0.2900, 6.643),
    5: (0.3910, 3.230),
}

_NACA5_REFLEXED = {
    # Q: (m, k1, k2_over_k1) for the reflexed camber line
    #   y_c = (k1/6) * ((x-m)^3 - (k2/k1)(1-m)^3 x - m^3 x + m^3)     x <  m
    #         (k1/6) * ((k2/k1)(x-m)^3 - (k2/k1)(1-m)^3 x - m^3 x + m^3)  x >= m
    # Coefficients from NACA TR-610.
    2: (0.1300, 51.990, 0.000764),
    3: (0.2170, 15.793, 0.006770),
    4: (0.3180, 6.520, 0.030300),
    5: (0.4410, 3.191, 0.135500),
}


def naca5(designation: str, n_pts: int = 80, *, closed_te: bool = True) -> np.ndarray:
    """Build a NACA 5-digit airfoil section.

    Parameters
    ----------
    designation :
        5 digits ``LPSTT``: ``L`` = design CL * 20/3 (first digit; for
        the standard family always ``2`` -> CL_design = 0.3), ``P`` =
        camber position in tenths of half-chord-fraction (second
        digit; encodes the mean-line designation Q = 1..5), ``S`` =
        camber type (``0`` = standard, ``1`` = reflexed),
        ``TT`` = max thickness as % chord.  Example: ``"23012"`` =
        standard mean line ``230``, 12% thickness.
    n_pts :
        Even total node count, ``>= 8``.
    closed_te :
        Whether to use the closed-TE thickness polynomial.

    Returns
    -------
    numpy.ndarray, shape ``(n_pts, 2)``
        Selig-ordered coordinates.

    Raises
    ------
    ParameterError
        On invalid designation, unsupported mean-line code, or
        invalid ``n_pts``.
    """
    digits = _validate_designation(designation, expected_len=5)
    if n_pts < 5:
        raise ParameterError(f"naca5: n_pts must be >= 5, got {n_pts}")

    L = int(digits[0])
    Q = int(digits[1])
    S = int(digits[2])
    t = int(digits[3:5]) / 100.0

    if L != 2:
        raise ParameterError(
            f"naca5: only the L=2 (CL_design=0.3) mean-line family is supported, "
            f"got designation {designation!r} with L={L}"
        )
    if S not in (0, 1):
        raise ParameterError(
            f"naca5: third digit (camber type) must be 0 (standard) or 1 (reflexed), got {S}"
        )

    x = _cosine_x(n_pts)
    yt = _thickness4(x, t, closed_te=closed_te)  # 5-digit shares 4-digit thickness
    yc, dyc_dx = _camber5(x, Q, reflexed=bool(S))

    return _assemble_selig(x, yt, yc, dyc_dx)


def _camber5(x: np.ndarray, Q: int, *, reflexed: bool) -> tuple[np.ndarray, np.ndarray]:
    """NACA 5-digit camber line and slope, dispatched by mean-line code Q."""
    if reflexed:
        coefs = _NACA5_REFLEXED.get(Q)
        if coefs is None:
            raise ParameterError(
                f"naca5: unsupported reflexed mean-line code Q={Q} "
                f"(supported: {sorted(_NACA5_REFLEXED)})"
            )
        m, k1, k2_over_k1 = coefs
        # Reflexed camber line (NACA TR-610).
        yc = np.where(
            x < m,
            (k1 / 6.0) * ((x - m) ** 3 - k2_over_k1 * (1.0 - m) ** 3 * x - m**3 * x + m**3),
            (k1 / 6.0)
            * (k2_over_k1 * (x - m) ** 3 - k2_over_k1 * (1.0 - m) ** 3 * x - m**3 * x + m**3),
        )
        dyc = np.where(
            x < m,
            (k1 / 6.0) * (3.0 * (x - m) ** 2 - k2_over_k1 * (1.0 - m) ** 3 - m**3),
            (k1 / 6.0) * (3.0 * k2_over_k1 * (x - m) ** 2 - k2_over_k1 * (1.0 - m) ** 3 - m**3),
        )
        return yc, dyc

    coefs = _NACA5_STANDARD.get(Q)
    if coefs is None:
        raise ParameterError(
            f"naca5: unsupported standard mean-line code Q={Q} "
            f"(supported: {sorted(_NACA5_STANDARD)})"
        )
    m, k1 = coefs
    # Standard camber line (NACA TR-537).
    yc = np.where(
        x < m,
        (k1 / 6.0) * (x**3 - 3.0 * m * x**2 + m**2 * (3.0 - m) * x),
        (k1 * m**3 / 6.0) * (1.0 - x),
    )
    dyc = np.where(
        x < m,
        (k1 / 6.0) * (3.0 * x**2 - 6.0 * m * x + m**2 * (3.0 - m)),
        -(k1 * m**3 / 6.0) * np.ones_like(x),
    )
    return yc, dyc


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _validate_designation(designation: str, *, expected_len: int) -> str:
    if not isinstance(designation, str) or not designation.isdigit():
        raise ParameterError(
            f"NACA designation must be a digit string of length {expected_len}, got {designation!r}"
        )
    if len(designation) != expected_len:
        raise ParameterError(
            f"NACA designation must be {expected_len} digits, "
            f"got {len(designation)}-digit {designation!r}"
        )
    return designation


def _cosine_x(n_pts: int) -> np.ndarray:
    """Cosine-spaced x in [0, 1] with ``n_pts`` samples.

    Cosine spacing clusters points near the LE (x=0) and TE (x=1) —
    the natural NACA parameterisation.  The full Selig curve is
    assembled in :func:`_assemble_selig` from these ``n_pts`` chord
    samples on each side, with the LE node shared.
    """
    beta = np.linspace(0.0, np.pi, n_pts)
    return 0.5 * (1.0 - np.cos(beta))


def _assemble_selig(
    x: np.ndarray,
    yt: np.ndarray,
    yc: np.ndarray,
    dyc_dx: np.ndarray,
) -> np.ndarray:
    """Combine thickness, camber, and slope into Selig-ordered coords.

    Per NACA convention the upper and lower surfaces are offset
    perpendicular to the camber line:

        x_u = x - y_t sin(theta),  y_u = y_c + y_t cos(theta)
        x_l = x + y_t sin(theta),  y_l = y_c - y_t cos(theta)

    where ``theta = atan(dy_c/dx)``.  For symmetric airfoils
    (``y_c = 0``) this collapses to the familiar
    ``(x, ±y_t)``.

    The output is Selig-ordered: upper-TE -> LE -> lower-TE,
    with the LE node shared between the two halves.
    """
    theta = np.arctan(dyc_dx)
    sin_th = np.sin(theta)
    cos_th = np.cos(theta)

    xu = x - yt * sin_th
    yu = yc + yt * cos_th
    xl = x + yt * sin_th
    yl = yc - yt * cos_th

    # Upper surface: TE -> LE
    upper = np.column_stack([xu[::-1], yu[::-1]])
    # Lower surface: LE -> TE (skip LE to avoid duplication)
    lower = np.column_stack([xl[1:], yl[1:]])
    return np.vstack([upper, lower])
