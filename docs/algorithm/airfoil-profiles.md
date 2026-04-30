# Airfoil profile generation

The `shore.profiles` module produces 2D airfoil section coordinates
analytically, with no external dependencies beyond NumPy. This page
documents the polynomial families, the Selig ordering convention,
and the cosine-sampling rule. For the consumer-facing API, see
[`shore profile`](/reference/profile).

## Selig ordering

All sections are returned as ``(N, 2)`` arrays ordered around the
closed curve starting at the upper trailing edge, walking forward
over the suction (upper) surface to the leading edge, then back
over the pressure (lower) surface to the lower trailing edge:

```
out[0]    = upper TE   (x ≈ 1, y ≥ 0)
out[N//2] = leading edge   (x = 0, y = 0)
out[-1]   = lower TE   (x ≈ 1, y ≤ 0)
```

This matches the format used by the
[UIUC airfoil database](https://m-selig.ae.illinois.edu/ads/coord_database.html)
and most CFD pre-processors, so any UIUC `.dat` drops into
[`shore body extrude`](/reference/body) directly.

For a sharp TE (``closed_te=True``, default), ``out[0] == out[-1]``
exactly; the polynomial residual is forced to zero so the curve
closes. For a blunt TE (``closed_te=False``), the two endpoints are
separated by the natural ~0.1% chord polynomial residual.

## Cosine sampling

Points are placed at chord positions ``x_i`` derived from a uniform
angular sampling of a half-cosine:

$$
x_i = \frac{1}{2} \left(1 - \cos\beta_i\right), \quad
\beta_i = \frac{i\pi}{n_{pts} - 1}, \quad i = 0, 1, \dots, n_{pts} - 1
$$

This clusters samples at the leading edge (``x = 0``) and trailing
edge (``x = 1``) — the natural NACA parameterisation, matched to
the rapid curvature changes there. The argument ``n_pts`` is the
**chord-sample count per side**; the closed Selig curve has
``2 * n_pts - 1`` total points (the LE node is shared between
upper and lower).

## NACA 4-digit family

Designation ``MPTT``:

- ``M`` = max camber as % chord (first digit)
- ``P`` = camber position in tenths of chord (second digit)
- ``TT`` = max thickness as % chord (last two digits)

Example: ``2412`` = 2% camber at 40% chord, 12% thickness.

### Thickness distribution

The half-thickness is the standard NACA 4-digit polynomial:

$$
y_t(x) = 5\,t \left(
  a_0 \sqrt{x} + a_1 x + a_2 x^2 + a_3 x^3 + a_4 x^4
\right)
$$

with coefficients

$$
a_0 = 0.2969, \quad a_1 = -0.1260, \quad a_2 = -0.3516, \quad a_3 = 0.2843
$$

and the last coefficient picked to control the TE closure:

| Form | $a_4$ | TE behaviour |
|---|---|---|
| Open TE (original) | $-0.1015$ | $y_t(1) \approx 0.00126$ — leaves a ~0.25% chord gap |
| Closed TE | $-0.1036$ | $y_t(1) = 0$ exactly |

`shore.profiles.naca` uses the closed-TE coefficient by default;
pass ``closed_te=False`` to keep the original open-TE behaviour
(useful when you want to detect the resulting blunt TE downstream
in the [C-grid topology](/topologies/c-grid)).

### Camber line

For ``M = P = 0`` the airfoil is symmetric and the camber line
vanishes. Otherwise:

$$
y_c(x) = \begin{cases}
  \dfrac{m}{p^2} \left(2 p x - x^2\right) & x < p \\
  \dfrac{m}{(1-p)^2} \left((1 - 2p) + 2 p x - x^2\right) & x \geq p
\end{cases}
$$

with $m = M/100$, $p = P/10$.

### Surface coordinates

Upper and lower surfaces are offset perpendicular to the camber
line:

$$
\begin{aligned}
x_u &= x - y_t \sin\theta, & y_u &= y_c + y_t \cos\theta \\
x_l &= x + y_t \sin\theta, & y_l &= y_c - y_t \cos\theta
\end{aligned}
$$

where $\theta = \arctan(\mathrm{d}y_c/\mathrm{d}x)$. For a
symmetric airfoil ($y_c = 0$, $\theta = 0$) this collapses to the
familiar $(x, \pm y_t)$.

## NACA 5-digit family

Designation ``LPSTT``:

- ``L`` = design $C_L \times \tfrac{20}{3}$. Always ``2`` in the
  standard family (giving $C_L^{\rm design} = 0.3$).
- ``P`` = mean-line shape designation $Q \in \{1, 2, 3, 4, 5\}$.
- ``S`` = camber type. ``0`` = standard, ``1`` = reflexed.
- ``TT`` = max thickness as % chord.

Example: ``23012`` = standard mean line ``230`` with 12% thickness;
``23112`` is the reflexed counterpart.

### Thickness distribution

The 5-digit family uses the **same** thickness polynomial as the
4-digit family — only the camber line changes.

### Standard mean line (S = 0)

From NACA TR-537 (Jacobs & Pinkerton, 1935):

$$
y_c(x) = \begin{cases}
  \dfrac{k_1}{6} \left(x^3 - 3 m x^2 + m^2 (3 - m) x\right) & x < m \\[4pt]
  \dfrac{k_1 m^3}{6} (1 - x) & x \geq m
\end{cases}
$$

The coefficients $(m, k_1)$ are tabulated by mean-line code $Q$:

| Q | designation | $m$ | $k_1$ |
|---|---|---|---|
| 1 | ``210`` | 0.0580 | 361.40 |
| 2 | ``220`` | 0.1260 | 51.640 |
| 3 | ``230`` | 0.2025 | 15.957 |
| 4 | ``240`` | 0.2900 | 6.643 |
| 5 | ``250`` | 0.3910 | 3.230 |

### Reflexed mean line (S = 1)

From NACA TR-610. The camber line gains a third segment so the
aft loading is *opposite* to the fore loading — useful for tailless
configurations. Coefficients $(m, k_1, k_2/k_1)$:

| Q | designation | $m$ | $k_1$ | $k_2/k_1$ |
|---|---|---|---|---|
| 2 | ``221`` | 0.1300 | 51.990 | 0.000764 |
| 3 | ``231`` | 0.2170 | 15.793 | 0.006770 |
| 4 | ``241`` | 0.3180 | 6.520 | 0.030300 |
| 5 | ``251`` | 0.4410 | 3.191 | 0.135500 |

The reflexed camber polynomial is

$$
y_c(x) = \begin{cases}
  \dfrac{k_1}{6} \left[(x - m)^3 - \tfrac{k_2}{k_1}(1 - m)^3 x - m^3 x + m^3\right] & x < m \\[6pt]
  \dfrac{k_1}{6} \left[\tfrac{k_2}{k_1}(x - m)^3 - \tfrac{k_2}{k_1}(1 - m)^3 x - m^3 x + m^3\right] & x \geq m
\end{cases}
$$

The surface-offset rule (upper / lower from $y_c$ and $y_t$) is
the same as for the 4-digit family.

## Out of scope

- **Modified 4-digit** (e.g. ``0012-64``) — sharper LE and a
  shifted thickness peak; needs additional polynomial coefficients
  per family.
- **NACA 6-series** (e.g. ``63-012``) — laminar-flow profiles
  whose thickness distribution comes from solving Theodorsen's
  conformal-mapping integral. Not closed-form; needs numerical
  integration.

These can be added later; for the typical airfoil-CFD workflow
(symmetric and cambered 4-digit, the canonical 23xxx family) the
current coverage is enough.

## See also

- [`shore profile`](/reference/profile) — CLI command to generate `.dat` files.
- [`shore body extrude`](/reference/body) — extrude a Selig section into a 3D wing STL.
- [C-grid topology](/topologies/c-grid) — the meshing target for extruded wing STLs.
