---
title: shore.profiles
outline: [2, 3]
---

# `shore.profiles`

2D airfoil profile generators (NACA 4-digit, 5-digit).

Pure NumPy.  No trimesh dependency — :mod:`shore.profiles` produces
`(N, 2)` Selig-ordered coordinates only.  For 3D STL extrusion
into wings, see :mod:`shore.shapes`.

Selig ordering convention
-------------------------
Points are ordered around the closed section starting at the upper
trailing edge, walking forward over the suction (upper) surface to
the leading edge, then back over the pressure (lower) surface to
the lower trailing edge:

```
out[0]    = upper TE        (x ≈ 1, y ≥ 0)
out[N//2] = leading edge    (x ≈ 0, y ≈ 0)
out[-1]   = lower TE        (x ≈ 1, y ≤ 0)

```
For a **sharp** TE (`closed_te=True`, default) `out[0] == out[-1]`
in the limit; the polynomial residual is forced to zero so the curve
is exactly closed.  For a **blunt** TE (`closed_te=False`), the two
endpoints are separated by the analytic NACA polynomial residual
(`≈ 0.1% chord` for the standard 4-digit family).

**Modules**

- [`shore.profiles.io`](./shore.profiles.io) — Selig `.dat` reader and writer for 2D airfoil sections

**Re-exported**

- [`read_section`](./shore.profiles.io#read_section)
- [`write_section`](./shore.profiles.io#write_section)

## `naca` <Badge type="info" text="function" /> {#naca}

```python
naca(designation: str, n_pts: int = 80, *, closed_te: bool = True) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/profiles/_naca.py#L31">source</a>

Build a NACA airfoil section by digit-count dispatch.

**Parameters**

- **`designation`** — Digit string: 4 chars (4-digit family, e.g. `"0012"`, `"2412"`) or 5 chars (5-digit family, e.g. `"23012"`, `"23112"` reflexed).
- **`n_pts`** — Number of chord samples per side (LE..TE).  Must be `>= 5`. The returned curve has `2 * n_pts - 1` points around the full section: `n_pts` on the upper surface plus `n_pts` on the lower surface, with the LE node shared.  Cosine-spaced along x (denser at LE and TE).
- **`closed_te`** — If `True` (default), apply the small thickness-polynomial adjustment that closes the TE exactly.  If `False`, use the standard open-TE polynomial whose residual is `≈ 0.1% chord` at the TE — useful when you want to detect a blunt TE downstream.

**Returns**

- `numpy.ndarray, shape ``(2 * n_pts - 1, 2)``` — Selig-ordered coordinates: upper-TE -&gt; LE -&gt; lower-TE.

**Raises**

- `ParameterError` — On invalid designation or `n_pts`.

## `naca4` <Badge type="info" text="function" /> {#naca4}

```python
naca4(designation: str, n_pts: int = 80, *, closed_te: bool = True) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/profiles/_naca.py#L82">source</a>

Build a NACA 4-digit airfoil section.

**Parameters**

- **`designation`** — 4 digits `MPTT`: `M` = max camber as % chord (first digit), `P` = camber position in tenths of chord (second digit), `TT` = max thickness as % chord (last two digits). Example: `"2412"` = 2% camber at 40% chord, 12% thickness.
- **`n_pts`** — Even total node count, `>= 8`.
- **`closed_te`** — Whether to use the closed-TE thickness polynomial.

**Returns**

- `numpy.ndarray, shape ``(n_pts, 2)``` — Selig-ordered coordinates.

**Raises**

- `ParameterError` — On invalid designation or `n_pts`.

## `naca5` <Badge type="info" text="function" /> {#naca5}

```python
naca5(designation: str, n_pts: int = 80, *, closed_te: bool = True) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/profiles/_naca.py#L178">source</a>

Build a NACA 5-digit airfoil section.

**Parameters**

- **`designation`** — 5 digits `LPSTT`: `L` = design CL \* 20/3 (first digit; for the standard family always `2` -&gt; CL_design = 0.3), `P` = camber position in tenths of half-chord-fraction (second digit; encodes the mean-line designation Q = 1..5), `S` = camber type (`0` = standard, `1` = reflexed), `TT` = max thickness as % chord.  Example: `"23012"` = standard mean line `230`, 12% thickness.
- **`n_pts`** — Even total node count, `>= 8`.
- **`closed_te`** — Whether to use the closed-TE thickness polynomial.

**Returns**

- `numpy.ndarray, shape ``(n_pts, 2)``` — Selig-ordered coordinates.

**Raises**

- `ParameterError` — On invalid designation, unsupported mean-line code, or invalid `n_pts`.
