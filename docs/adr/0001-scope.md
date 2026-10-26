# 0001 — Scope: a CAD-associative structured/overset mesher

- **Status:** Proposed
- **Date:** 2026-10-25

## Context

The goal "SHORE becomes a full CFD mesher" does not constrain any design
decision. ANSYS ICEM CFD and Cadence Fidelity Pointwise are *general* meshers:
structured, unstructured, hybrid, polyhedral, dozens of solver exports. They
represent decades of commercial engineering. A general unstructured/hybrid
mesher is also a commodity: Gmsh (OCCT kernel, Python API) is free and good.

The open-source landscape (survey, 2026-10) has a precise gap:

- **Salome HexaBlock** has an ICEM-style block model with CAD association, but
  no hyperbolic marching and no overset packaging; it lives inside Salome.
- **Gmsh** transfinite meshing needs every sub-volume to have 5–6 faces of 3–4
  edges, has no elliptic or hyperbolic grid smoothing, and drops transfinite
  constraints on OCC transforms.
- **MDOlab pyHyp** marches hyperbolic volume grids, but its own documentation
  states the input surface meshes "generally come from ICEM".
- **NASA Chimera Grid Tools** (OVERGRID, SURGRD, HYPGEN) covers the overset
  workflow, but distribution is restricted to US persons.
- **ICEM CFD** blocking is documented in 2026 R1 but frozen: dropped from the
  Student edition, features migrating to SpaceClaim/Fluent Meshing, now under
  Synopsys.

No open tool combines CAD-associative multiblock blocking, hyperbolic
near-body marching and Chimera solver packaging. SHORE already has the last
two.

## Decision

SHORE is a **CAD-associative structured multiblock and overset mesher for
Chimera CFD solvers**.

In scope:

1. Geometry input from STL (existing) and STEP/IGES/BREP (new, ADR 0002).
2. Block topology authored by templates (existing topologies), by script, or
   interactively, with every block entity associated to geometry (ADR 0003).
3. Structured grid algorithms: 1-D distributions, algebraic TFI, elliptic
   smoothing, algebraic and hyperbolic marching for volume grids, hyperbolic
   marching for surface (collar) grids (ADR 0005).
4. Overset packaging for Xall/Xnavis (existing), plus CGNS and Plot3D export
   so other structured/overset solvers can use SHORE.
5. A desktop GUI over the same command layer as the CLI (ADRs 0004, 0006).

Non-goals:

- unstructured, hybrid, prismatic-layer or polyhedral meshes (use Gmsh);
- CAD authoring beyond the parametric bodies SHORE already generates
  (build123d is the authoring tool);
- domain connectivity, hole cutting and donor search (the solver chain's
  overset preprocessor owns these; SHORE provides marker boxes only).

## Alternatives considered

- **General mesher (ICEM/Pointwise parity).** Rejected: unbounded scope,
  competes with free Gmsh on unstructured meshes, and dilutes the structured
  quality that is SHORE's reason to exist.
- **Stay an STL-only pipeline.** Rejected: STL facets cap surface accuracy at
  the tessellation error, and centroid ray-casting limits SHORE to star-shaped
  bodies. Neither limit is acceptable for a production mesher.

## Consequences

- Every later decision is judged by whether it serves structured/overset
  quality. Feature requests for unstructured output are declined by
  reference to this record.
- CGNS export becomes a requirement, not a nice-to-have, because it is the
  path to solvers other than Xnavis.
- The existing topologies (cubed-sphere, O-, C-, OH-, CH-grid) stay. They
  become *templates* instantiated against a geometry (ADR 0003).

## Verification

Revisit this record if a funded use case needs hybrid prism layers on the
SHORE surface grid. That would be an export path to an external hybrid mesher,
not a new SHORE algorithm.
