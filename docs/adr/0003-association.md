# 0003 — Block-to-geometry association replaces `body_projected`

- **Status:** Proposed
- **Date:** 2026-10-25

## Context

ICEM CFD Hexa, Salome HexaBlock and Pointwise share one core abstraction.
Each topological entity of the block structure is **associated** with a
geometric entity, and its grid nodes are projected onto that entity:

| Block entity | Geometry entity |
|---|---|
| vertex | CAD vertex, or a point on a curve or surface |
| edge | CAD curve (or a chain of curves) |
| face | CAD surface (or a set of surfaces) |
| block | none (interior) |

SHORE's block model already has the shape for this. `Vertex`, `Edge` and
`Face` exist, and each carries `body_projected: bool`. The flag can only say
"projected onto the one body". It cannot say *which* curve or surface, so
nothing can be meshed with edges pinned to CAD feature lines (wing leading
edge, trailing edge, wing–fuselage intersection).

The hard part is not the data structure but **persistent naming**. OCCT gives
shapes no identity that survives re-import of an edited STEP file. FreeCAD
calls this the "topological naming problem". EGADS solves it with persistent
attributes.

## Decision

1. **`Association` value type** in `shore.mesh.association`:

   ```python
   @dataclass(frozen=True)
   class Association:
       geometry: str                # key of the Geometry in the project
       kind: EntityKind             # BODY | SURFACE | CURVE | POINT
       keys: tuple[EntityKey, ...]  # one or more entities (curve chain, surface set)
   ```

   `Vertex.association`, `Edge.association` and `Face.association` replace
   `body_projected`. `body_projected=True` becomes
   `Association(geometry, BODY, ())`, so existing behaviour is preserved and
   the property `body_projected` stays as a read-only shim.

2. **Projection follows association.** After any operation that moves nodes
   (TFI, smoothing, marching at k=0, interactive drag), the nodes of an
   associated entity are projected onto its geometry, with priority
   **vertex → edge → face**. Vertices pin edges, and edges pin face
   boundaries. This is ICEM's projection order.

3. **`EntityKey` is a fingerprint, not an OCCT index:**
   - entity type and, if present, the XCAF label name from the STEP file;
   - bounding box, length or area, and centroid, rounded to a tolerance
     relative to the model size.

   On re-import, keys are matched by name first, then by fingerprint nearest
   neighbour within tolerance. Unmatched associations are reported, never
   silently dropped. A topology that loses an association fails validation,
   in the spirit of the `AxisMap` wiring invariant.

4. **Serialisation.** Associations are written next to the face graph in
   `.adjacency.json` (schema version bump). They are part of the project
   file (ADR 0004).

5. **Templates instantiate against geometry.** Each existing topology gains
   an association step:
   - the C-grid associates its TE column edges to the detected TE curve;
   - the cubed-sphere associates all k=0 faces to the body.

   New interactive blockings (ADR 0006) are built from the same primitives.

## Alternatives considered

- **Keep `body_projected` and add a side table.** Rejected: two sources of
  truth for one fact, and the side table would not follow blocks through
  `split_blocks`.
- **Store OCCT `TopoDS_Shape` handles directly.** Rejected: they are not
  serialisable and do not survive re-import.
- **Adopt EGADS attributes.** Rejected as a dependency (ADR 0002). Its
  persistence semantics are the reference for the fingerprint matching.

## Consequences

- `split_blocks` must propagate associations to the child entities. An edge
  split yields two edges with the same curve association, and the new
  internal faces are unassociated.
- `.adjacency.json` gains a schema version. Readers accept the old schema and
  map `body_projected` as above.
- Fingerprint matching fails on large topological edits (a face split in CAD).
  The tool must then ask the user to re-associate, which the GUI provides.

## Verification

- **Round-trip:** associations serialise, reload and match after STEP
  re-export of an unchanged model (100 % match).
- **Robustness:** after a small CAD edit (wing span +5 %), names or
  fingerprints match every surviving entity, and the report lists the
  changed ones.
- **Regression:** with `StlGeometry`, all existing topology tests pass
  unchanged through the shim.
