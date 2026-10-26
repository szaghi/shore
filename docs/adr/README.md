# Architecture decision records

These records fix the decisions that take SHORE from a semi-automatic
STL → Chimera pipeline to a CAD-associative structured/overset mesher with an
interactive GUI. Each record states the decision, the forces behind it, the
alternatives rejected and the consequences accepted. A record is changed only
by a later record that supersedes it.

| # | Title | Status |
|---|---|---|
| [0001](./0001-scope.md) | Scope: a CAD-associative structured/overset mesher, not a general mesher | Proposed |
| [0002](./0002-geometry-backend.md) | Geometry protocol with STL and OCCT (build123d/OCP) backends | Proposed |
| [0003](./0003-association.md) | Block-to-geometry association replaces `body_projected` | Proposed |
| [0004](./0004-command-layer.md) | Command layer and replayable journal under CLI, API and GUI | Proposed |
| [0005](./0005-marching-schemes.md) | Name the current marcher "algebraic"; add a true hyperbolic scheme | Proposed |
| [0006](./0006-gui-stack.md) | GUI: PySide6 + pygfx/wgpu, desktop, Python + numpy only | Proposed |
| [0007](./0007-surface-grids.md) | Surface grids: TFI + elliptic smoothing in (u, v) per face | Proposed |

The phased plan that implements them is in [roadmap.md](./roadmap.md).

## Template

```markdown
# NNNN — Title

- **Status:** Proposed | Accepted | Superseded by NNNN
- **Date:** YYYY-MM-DD

## Context
## Decision
## Alternatives considered
## Consequences
## Verification
```
