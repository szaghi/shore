---
title: Python API
---

# Python API

Generated from the docstrings by `docs/scripts/api_docs.py`; do not edit
these pages by hand. The command-line interface is documented under
[CLI commands](/reference/), the project journal under
[Project scripts](/reference/run).

## Package

| Module | Summary |
|---|---|
| [`shore`](./shore) | SHORE — Structured Hyperbolic Overlapping Remeshing Engine |
| [`shore.config`](./shore.config) | TOML-driven configuration for the `shore mesh` pipeline |
| [`shore.exceptions`](./shore.exceptions) | SHORE-specific exception hierarchy |
| [`shore.pipeline`](./shore.pipeline) | Meshing pipeline steps shared by the CLI and the command layer  |
| [`shore.quality`](./shore.quality) | Structured-grid quality metrics (roadmap phase 4) |

## Geometry

| Module | Summary |
|---|---|
| [`shore.geometry`](./shore.geometry) | Geometry backends behind one protocol  |
| [`shore.geometry.base`](./shore.geometry.base) | Geometry protocol |
| [`shore.geometry.keys`](./shore.geometry.keys) | Persistent entity identity |
| [`shore.geometry.occt`](./shore.geometry.occt) | OCCT backend of the `Geometry` protocol  |
| [`shore.geometry.stl`](./shore.geometry.stl) | STL backend of the `Geometry` protocol  |
| [`shore.profiles`](./shore.profiles) | 2D airfoil profile generators (NACA 4-digit, 5-digit) |
| [`shore.profiles.io`](./shore.profiles.io) | Selig `.dat` reader and writer for 2D airfoil sections |
| [`shore.shapes`](./shore.shapes) | 3D body STL builders backed by trimesh |

## Block model

| Module | Summary |
|---|---|
| [`shore.mesh`](./shore.mesh) | OOP mesh layer for SHORE |
| [`shore.mesh.association`](./shore.mesh.association) | Block-to-geometry association  |
| [`shore.mesh.block`](./shore.mesh.block) | Block ABC and HexBlock |
| [`shore.mesh.blocking`](./shore.mesh.blocking) | Blocking |
| [`shore.mesh.body`](./shore.mesh.body) | Body |
| [`shore.mesh.edge`](./shore.mesh.edge) | Edge |
| [`shore.mesh.face`](./shore.mesh.face) | Face and FaceBC |
| [`shore.mesh.inspect`](./shore.mesh.inspect) | File-based inspection utilities |
| [`shore.mesh.mesh`](./shore.mesh.mesh) | Mesh |
| [`shore.mesh.orthogonal`](./shore.mesh.orthogonal) | Grid lines leaving the walls along their normal  and wall layers  |
| [`shore.mesh.redistribute`](./shore.mesh.redistribute) | Redistributing the nodes of a block edge  |
| [`shore.mesh.seams`](./shore.mesh.seams) | Point-matched seams between block faces  |
| [`shore.mesh.spacing`](./shore.mesh.spacing) | Spacing1D |
| [`shore.mesh.topology`](./shore.mesh.topology) | Topology ABC and concrete topologies (cubed-sphere, O-grid, C-grid, OH-/CH-hybrid) |
| [`shore.mesh.vertex`](./shore.mesh.vertex) | Vertex |

## Surface

| Module | Summary |
|---|---|
| [`shore.surface`](./shore.surface) | Surface I/O and unstructured-to-structured remeshing |
| [`shore.surface.anchor`](./shore.surface.anchor) | Anchor-point resolution strategies for ray-cast projections |
| [`shore.surface.cap`](./shore.surface.cap) | Polar cap k=0 seed face via gnomonic projection of a flat tangent square |
| [`shore.surface.cgrid`](./shore.surface.cgrid) | C-grid surface construction |
| [`shore.surface.equator`](./shore.surface.equator) | Equator band construction from cap seam boundaries |
| [`shore.surface.intersect`](./shore.surface.intersect) | Intersection loop of two bodies, the start curve of a collar grid  |
| [`shore.surface.io`](./shore.surface.io) | STL reader |
| [`shore.surface.march`](./shore.surface.march) | Hyperbolic surface marching |
| [`shore.surface.projector`](./shore.surface.projector) | Ray-mesh projector with cached BVH and per-ray status tracking |
| [`shore.surface.remesh`](./shore.surface.remesh) | Unstructured STL surface → structured (ni, nj) quad grid via sphere projection |
| [`shore.surface.tfi`](./shore.surface.tfi) | Surface grids on one parametric face |

## Volume and topologies

| Module | Summary |
|---|---|
| [`shore.primitive`](./shore.primitive) | Parametric structured-grid primitives for Chimera/overset assembly |
| [`shore.primitive.annulus`](./shore.primitive.annulus) | Annular cylinder primitive — a single block with optional theta-periodicity |
| [`shore.primitive.base`](./shore.primitive.base) | Primitive mesh base class and shared utilities |
| [`shore.primitive.box`](./shore.primitive.box) | Cartesian box primitive — a single axis-aligned hexahedral block |
| [`shore.primitive.flat_caps`](./shore.primitive.flat_caps) | Flat-caps cylinder primitive — 5-block conforming alternative to a hollow Chimera background |
| [`shore.primitive.h_grid_rect`](./shore.primitive.h_grid_rect) | Rectangular-duct H-grid primitive |
| [`shore.primitive.h_grid_tube`](./shore.primitive.h_grid_tube) | Circular-bore H-grid primitive |
| [`shore.topology`](./shore.topology) | Topology decomposition utilities for multi-block structured grids |
| [`shore.topology.cgrid_cap`](./shore.topology.cgrid_cap) | Cap-fill builder for the C-grid topology's wing-tip closures |
| [`shore.topology.ch_channel`](./shore.topology.ch_channel) | Outer-quad helpers for the CH (C-grid + H-channel) topology |
| [`shore.topology.oh_channel`](./shore.topology.oh_channel) | Channel-block frustum builder for the OH (cubed-sphere + H-channel) topology |
| [`shore.topology.sphere_cap`](./shore.topology.sphere_cap) | Polar cap block generation for sphere O-grids |
| [`shore.volume`](./shore.volume) | Hyperbolic normal-extrusion |
| [`shore.volume.hyperbolic`](./shore.volume.hyperbolic) | Hyperbolic marching primitives |
| [`shore.volume.hyperbolic_scheme`](./shore.volume.hyperbolic_scheme) | Hyperbolic volume marching |
| [`shore.volume.jacobian`](./shore.volume.jacobian) | Jacobian field computation for structured volume grids |
| [`shore.volume.spacing`](./shore.volume.spacing) | 1D grid spacing laws for surface angle redistribution and wall-normal extrusion |

## Collars

| Module | Summary |
|---|---|
| [`shore.collar`](./shore.collar) | The collar block of two intersecting bodies  |

## Projects

| Module | Summary |
|---|---|
| [`shore.project`](./shore.project) | Command layer |
| [`shore.project.commands`](./shore.project.commands) | Commands |
| [`shore.project.project`](./shore.project.project) | Project |

## Chimera and Xall

| Module | Summary |
|---|---|
| [`shore.balance`](./shore.balance) | Multi-way equal-partition split planner for Xall MPI load balancing |
| [`shore.boxes`](./shore.boxes) | Generate :class:`shore.io.boxes.Box` instances from an STL surface |
| [`shore.io`](./shore.io) | I/O |
| [`shore.io.adjacency`](./shore.io.adjacency) | Adjacency JSON sidecar for split block lists |
| [`shore.io.boxes`](./shore.io.boxes) | Hexahedral marker boxes for Xall `cc.par` |
| [`shore.io.cc_par`](./shore.io.cc_par) | Xall `cc.par` writer |
| [`shore.io.cgns`](./shore.io.cgns) | CGNS export of structured block lists  |
| [`shore.io.geo`](./shore.io.geo) | Xall geo format writer/reader |
| [`shore.io.grd`](./shore.io.grd) | Xall binary `.grd` writer |
| [`shore.io.plot3d`](./shore.io.plot3d) | Plot3D export of structured block lists  |
| [`shore.io.proc_input`](./shore.io.proc_input) | Xall `proc.input` writer + merger |
| [`shore.io.vtk`](./shore.io.vtk) | VTK export for ParaView |
| [`shore.split`](./shore.split) | Split structured blocks along ijk vertex indices, preserving adjacency |

## GUI and visualisation

| Module | Summary |
|---|---|
| [`shore.gui`](./shore.gui) | SHORE's desktop GUI  |
| [`shore.gui.forms`](./shore.gui.forms) | Command forms of the GUI  |
| [`shore.gui.session`](./shore.gui.session) | The GUI's owner of a project  |
| [`shore.gui.viewport`](./shore.gui.viewport) | The GUI's renderer interface  |
| [`shore.gui.viewport.pygfx_backend`](./shore.gui.viewport.pygfx_backend) | pygfx/wgpu implementation of the GUI viewport  |
| [`shore.gui.window`](./shore.gui.window) | SHORE's main window  |
| [`shore.vis`](./shore.vis) | Visualisation helpers for SHORE grids |
