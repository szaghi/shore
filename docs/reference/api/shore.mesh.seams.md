---
title: shore.mesh.seams
outline: [2, 3]
---

# `shore.mesh.seams`

Point-matched seams between block faces (ADR 0035, 0036).

The face graph names one partner per face, but a face may be shared with
several smaller faces: CH's `body.k_hi` with its three channel blocks,
whose `k_lo` faces all name `body.k_hi`. Every connected pair of faces is
taken once here, and the smaller face is located inside the larger one by
its nodes, which gives each seam its exact index ranges on both sides.

Indices follow the CGNS SIDS: 1-based, ranges ``[[imin, jmin, kmin], [imax,
jmax, kmax]]`, and a `Transform` `T` with `Index2 = T.(Index1 - Begin1)
+ Begin2` (entry `a` is `+-(b + 1)`: index `a`` of this block runs along
index `b` of the donor, reversed when negative).

## `Seam` <Badge type="tip" text="class" /> {#Seam}

```python
class Seam(NamedTuple)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L32">source</a>

One side of a point-matched seam: `face` of this block over `rng`,
onto `donor`'s `donor_face` over `donor_rng` (begin = the image of
`rng`'s begin), with `transform`.

## `seams` <Badge type="info" text="function" /> {#seams}

```python
seams(
    blocks: list[HexBlock],
    tol: float,
    *,
    only: Collection[tuple[str, str]] | None = None,
    bcs: Collection[FaceBC] = CONNECTED,
) -> dict[str, list[Seam]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L45">source</a>

Per block label, its seams; every connected pair of faces once on each side.

*tol* is the largest node distance accepted (absolute). With *only*,
just the pairs with a side in it (`(label, face)`) are located; *bcs*
are the face BCs taken as connections.

**Raises**

- `ParameterError` — A face connected to a block outside *blocks*, or a pair whose nodes do not coincide.

## `shared_faces` <Badge type="info" text="function" /> {#shared_faces}

```python
shared_faces(blocks: list[HexBlock]) -> set[tuple[str, str]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L103">source</a>

`(label, face)` of the faces with more than one partner.

A face is shared with several faces when more than one face names it as
partner, or when a face naming it is not its own partner.

## `transform_vector` <Badge type="info" text="function" /> {#transform_vector}

```python
transform_vector(
    face: str,
    partner: str,
    along: tuple[tuple[str, bool], tuple[str, bool]],
) -> tuple[int, int, int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L124">source</a>

`Transform` of a face connection from its AxisMap.

Along the faces this is the AxisMap; across them a step into this block
is a step out of the donor.

## `inverse_transform` <Badge type="info" text="function" /> {#inverse_transform}

```python
inverse_transform(t: tuple[int, int, int]) -> tuple[int, int, int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L142">source</a>

## `affine` <Badge type="info" text="function" /> {#affine}

```python
affine(
    t: tuple[int, int, int],
    origin: np.ndarray,
    image: np.ndarray,
    index: np.ndarray,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L149">source</a>

`T (index - origin) + image`.

## `point_range` <Badge type="info" text="function" /> {#point_range}

```python
point_range(block: HexBlock, face: str) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L159">source</a>

1-based `[[imin, jmin, kmin], [imax, jmax, kmax]]` of *face*.

## `face_size` <Badge type="info" text="function" /> {#face_size}

```python
face_size(block: HexBlock, face: str) -> int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L169">source</a>

## `Aliases` <Badge type="tip" text="class" /> {#Aliases}

```python
Aliases(blocks: list[HexBlock], tol: float)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L233">source</a>

Nodes that are one point across SHARED / DIRICHLET seams (ADR 0040).

Built topologies share most such nodes as numpy views, so writing one
writes all; DIRICHLET seams and blocks read back from saved grids (ADR
0033) hold copies. After an in-place operation, `sync` gives every
group the value of its last member the operation changed (block order),
as views would; groups it did not change, or changed alike, are left
as they are.

### `Aliases.sync` <Badge type="info" text="method" /> {#Aliases-sync}

```python
sync(changed: list[np.ndarray]) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/seams.py#L292">source</a>

Give each group the value of its last member that *changed* marks.

*changed* holds one boolean `(ni, nj, nk)` mask per block.
