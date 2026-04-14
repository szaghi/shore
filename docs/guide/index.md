# About SHORE

SHORE (**S**tructured **H**yperbolic **O**verlapping **R**emeshing **E**ngine) is a Python tool that converts raw triangulated surface meshes (STL) into structured body-fitted volume grids suitable for Chimera (overset) CFD solvers.

## Motivation

Chimera mesh generation requires two classes of grids:

- **Near-body grids** — structured, body-fitted, high near-wall resolution
- **Background grids** — Cartesian or simple structured blocks that fill the far field

SHORE generates the near-body grids. Starting from the surface of a body (e.g. a fuselage, wing section, or bluff body), it builds a structured O-grid that wraps tightly around it by marching layers of cells outward along the surface normal direction. The result is directly usable as input to the Xall_GPU CFD solver.

## Design philosophy

**Stages are explicit.** The pipeline has two clearly separated stages — surface projection and hyperbolic marching — each with its own CLI command. You can run them independently, inspect intermediate results, and rerun only the stage that changed.

**Failures are loud.** Every layer is checked for cell inversion before the grid advances. A negative Jacobian aborts marching immediately with a precise error message. There are no silent bad meshes.

**No topology guessing.** SHORE does not attempt to decompose arbitrary surfaces into structured patches automatically. The current O-grid approach assumes the body is star-shaped with respect to its centroid. Complex multi-block topologies are a future extension.

## Scope

| In scope | Out of scope (roadmap) |
|---|---|
| Raw triangulated STL input | STEP / parametric input |
| Single O-grid block | Multi-block structured decomposition |
| Xall_GPU geo output | Other CFD formats (Plot3D, CGNS) |
| Sphere-projected surface | Advanced surface remeshing |
| Hyperbolic normal extrusion | Elliptic volume smoothing |

## Related tools

- **Xall_GPU** — the target Chimera CFD solver; SHORE generates its near-body geo blocks
- **trimesh** — STL I/O and BVH-accelerated ray casting
- **NumPy / SciPy** — all numerical kernels
