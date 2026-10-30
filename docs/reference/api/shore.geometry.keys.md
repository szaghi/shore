---
title: shore.geometry.keys
outline: [2, 3]
---

# `shore.geometry.keys`

Persistent entity identity: fingerprints and matching (ADR 0003).

OCCT shapes carry no identity that survives re-reading an edited CAD file,
and `EntityRef` indices are only valid for one load. An `EntityKey` is
a geometric fingerprint of an entity, normalised by the model's bounding
box so that a uniform rescale leaves it unchanged:

- `kind` and `geom_type` (`Plane`, `Cylinder`, `BSplineSurface`,
  `Line`, `Circle`, ...) must match exactly;
- `centroid`, `measure` (area, length, 0 for a point) and `extent`
  (bounding-box size) are compared numerically.

`match_keys` re-finds stored keys in a geometry by a global one-to-one
assignment per (kind, type) group and reports each key as `same`
(unchanged within `tol_same`), `changed` (matched, but moved or
resized, as after a parametric edit) or `lost` (no candidate within
`tol_match`).

## `EntityKey` <Badge type="tip" text="dataclass" /> {#EntityKey}

```python
EntityKey(
    kind: EntityKind,
    geom_type: str,
    centroid: tuple[float, float, float],
    measure: float,
    extent: tuple[float, float, float],
    name: str | None = None,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L38">source</a>

Geometric fingerprint of one entity, normalised by the model bbox.

Coordinates are `(x - bbox_min) / diag`; `measure` is area / diag**2
(surface), length / diag (curve), 0 (point, body uses area).

### `EntityKey.to_json` <Badge type="info" text="method" /> {#EntityKey-to_json}

```python
to_json() -> dict[str, Any]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L52">source</a>

### `EntityKey.from_json` <Badge type="info" text="classmethod" /> {#EntityKey-from_json}

```python
from_json(data: dict[str, Any], *, source: str = '<key>') -> EntityKey
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L65">source</a>

### `EntityKey.distance` <Badge type="info" text="method" /> {#EntityKey-distance}

```python
distance(other: EntityKey) -> float
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L78">source</a>

Dissimilarity in normalised units; `inf` across kinds or types.

## `make_key` <Badge type="info" text="function" /> {#make_key}

```python
make_key(
    kind: EntityKind,
    geom_type: str,
    centroid: np.ndarray,
    measure: float,
    bbox: np.ndarray,
    model_bounds: np.ndarray,
    name: str | None = None,
) -> EntityKey
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L94">source</a>

Normalise raw entity properties into an `EntityKey`.

*bbox* and *model_bounds* are `(2, 3)` `[min, max]` boxes; *measure*
is raw area / length (0 for points).

## `KeyMatch` <Badge type="tip" text="dataclass" /> {#KeyMatch}

```python
KeyMatch(key: EntityKey, ref: EntityRef | None, status: MatchStatus, cost: float)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L124">source</a>

Where a stored key landed in a geometry.

## `MatchReport` <Badge type="tip" text="dataclass" /> {#MatchReport}

```python
MatchReport(matches: list[KeyMatch] = field(default_factory=list))
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L134">source</a>

Outcome of `match_keys`: one `KeyMatch` per stored key, in order.

### `MatchReport.count` <Badge type="info" text="method" /> {#MatchReport-count}

```python
count(status: MatchStatus) -> int
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L139">source</a>

### `MatchReport.all_matched` <Badge type="info" text="property" /> {#MatchReport-all_matched}

```python
all_matched: bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L143">source</a>

### `MatchReport.summary` <Badge type="info" text="method" /> {#MatchReport-summary}

```python
summary() -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L146">source</a>

## `match_keys` <Badge type="info" text="function" /> {#match_keys}

```python
match_keys(
    keys: list[EntityKey],
    geometry: Geometry,
    *,
    tol_same: float = 1e-09,
    tol_match: float = 0.25,
) -> MatchReport
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/keys.py#L153">source</a>

Re-find *keys* in *geometry* by one-to-one assignment.

Keys are grouped by `(kind, geom_type)`; within a group the
assignment minimises the total `EntityKey.distance`. A pairing with
cost above *tol_match* is rejected (`lost`); one at or below
*tol_same* is `same`; anything in between is `changed`. Names, when
both sides have one, must agree.
