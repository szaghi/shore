# shore profile

Generate 2D airfoil section coordinates and write them as Selig-format
``.dat`` files (the de-facto standard used by the UIUC airfoil
database and most CFD pre-processors).

`shore profile` is a sub-app with one command for now:

```
shore profile naca DESIGNATION [OPTIONS]
```

Output is one `.dat` file with the section's `(2 * n_pts - 1, 2)`
coordinates ordered upper-TE → LE → lower-TE. Compatible with
[`shore body extrude`](body) and any UIUC-style consumer.

For the algorithmic detail (NACA 4 vs 5-digit polynomial families,
closed vs open TE, cosine sampling), see
[Algorithm — airfoil profiles](/concepts/airfoil-profiles).

## `shore profile naca`

Generate a NACA 4-digit or 5-digit airfoil section.

### Synopsis

```
shore profile naca DESIGNATION [OPTIONS]
```

### Arguments

| Argument | Description |
|---|---|
| `DESIGNATION` | NACA digit string. 4 chars = 4-digit family (`0012`, `2412`); 5 chars = 5-digit family (`23012`, `23112` reflexed). |

### Options

| Option | Default | Description |
|---|---|---|
| `-o` / `--output PATH` | — (required) | Output `.dat` file |
| `--n-pts INTEGER` | 80 | Chord samples per side (LE..TE). The total point count is `2 * n_pts - 1`. |
| `--closed-te` / `--open-te` | `--closed-te` | Apply the small thickness-polynomial adjustment that closes the TE exactly. `--open-te` leaves the natural ~0.1% chord gap. |
| `--name TEXT` | `NACA <designation>` | First-line name comment in the `.dat` file |

### Examples

```bash
# Symmetric 12% thickness:
shore profile naca 0012 -o naca0012.dat

# Cambered 4-digit at higher resolution:
shore profile naca 2412 --n-pts 120 -o naca2412.dat

# 5-digit with explicit name:
shore profile naca 23012 -o naca23012.dat --name "NACA 23012 root section"

# Open TE for downstream blunt-TE detection:
shore profile naca 0012 --open-te -o naca0012-blunt.dat
```

### Selig output format

Each line is two whitespace-separated floats. The first line is an
optional name comment (skipped on read); the rest are coordinate
rows:

```
NACA 0012
   1.0000000   -0.0000000
   0.9996085    0.0000837
   0.9984318    0.0003342
   ...
   0.0003915   -0.0006611
   1.0000000    0.0000000
```

### Supported designations

The 4-digit family is fully supported. The 5-digit family supports
the standard mean lines `21012`, `22012`, `23012`, `24012`, `25012`
(``L = 2``, ``S = 0``) and their reflexed counterparts `22112`,
`23112`, `24112`, `25112`. Other ``L`` values (non-standard design
CL) are rejected with a clear error. The modified 4-digit
(`0012-64`-style) and 6-series families need additional integrals
and aren't implemented in this release.

## See also

- [`shore body extrude`](body) — extrude a Selig section into a 3D wing STL.
- [Algorithm — airfoil profiles](/concepts/airfoil-profiles) — polynomial derivations.
