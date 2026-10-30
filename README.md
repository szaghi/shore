<div align="center">

# SHORE

**Structured grids for overset CFD**

[![Latest Version](https://img.shields.io/pypi/v/shore-mesh.svg)](https://pypi.org/project/shore-mesh/)
[![Build Status](https://github.com/szaghi/shore/actions/workflows/python-package.yml/badge.svg)](https://github.com/szaghi/shore/actions)
[![Coverage](https://img.shields.io/endpoint?url=https://szaghi.github.io/shore/coverage-badge.json)](https://szaghi.github.io/shore/)
[![Supported Python versions](https://img.shields.io/badge/Py-%203.11,%203.12-blue.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/szaghi/shore/main/docs/public/figures/dark/hero-wing-body.webp">
  <img alt="A half wing-body blocked and premeshed with SHORE" src="https://raw.githubusercontent.com/szaghi/shore/main/docs/public/figures/light/hero-wing-body.webp" width="640">
</picture>

**[Documentation](https://szaghi.github.io/shore/)** ·
[Quick start](https://szaghi.github.io/shore/guide/quickstart) ·
[Tutorials](https://szaghi.github.io/shore/tutorials/) ·
[Cookbook](https://szaghi.github.io/shore/cookbook/) ·
[Python API](https://szaghi.github.io/shore/reference/api/)

</div>

SHORE (Structured Hyperbolic Overlapping Remeshing Engine) makes the
structured multi-block grids of an overset (Chimera) CFD case and writes the
files the solver reads. It is a Python package with a CLI, a scriptable
`Project` API and a desktop GUI; all three record the same replayable
journal.

- **Marched near-body grids** — O-grid, cubed sphere, C-grid with wing-tip
  caps, OH and CH hybrids, extruded from a body with seam-aware normals and
  a Jacobian guard (algebraic or Steger–Chaussee hyperbolic scheme).
- **Interactive blocking** — split, hole, O-grid, associate with geometry,
  premesh with seams matched bit for bit, wall-normal lines and a wall layer
  with an exact first cell.
- **Collars** — intersection curves, hyperbolic surface marching and the
  collar volume at body junctions, half models and sharp corners included.
- **STL and CAD** — STEP / IGES / BREP through OpenCascade with exact
  projection; fingerprinted entities so journals replay on edited models.
- **Solver output** — the Xall/Xnavis bundle (`.grd`, `cc.par`,
  `proc.input`, marker boxes) for multi-component runs, CGNS, Plot3D, VTK.

## Install

```bash
pipx install shore-mesh                        # CLI
pipx install 'shore-mesh[cad,gui,cgns,vis]'    # with CAD, GUI, CGNS and PyVista
```

Python 3.11+. For library use, `pip install shore-mesh` in a virtual
environment. [Installation](https://szaghi.github.io/shore/guide/installation)
has the details and the one-line installer.

## A first grid

```bash
shore body sphere -o sphere.stl --subdivisions 4
shore mesh sphere.stl --ni 33 --nj 64 --nk 25 --ds 0.01 --growth 1.15 \
    --theta-cap-deg 35 -o sphere --journal sphere.py
shore check sphere_sub0.geo
shore run sphere.py              # replays the journal, bit for bit
```

The same in Python:

```python
from shore.project import Project

p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=25, ni=33, nj=64, ds=0.01, growth=1.15, theta_cap=0.61)
p.march()
p.write_geo("sphere")
p.save_journal("sphere.py")
```

Then follow the [tutorials](https://szaghi.github.io/shore/tutorials/): a
wing with a C-grid, a Chimera assembly and its Xall bundle, a collar, projects
and replay, and a CAD wing–body blocked in Python and in the GUI.

## Development

```bash
git clone https://github.com/szaghi/shore.git && cd shore
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
make test    # ~2000 tests
make lint
```

New algorithms are accepted against specs frozen before they are written
(`tests/spec/`, [ADR 0005](https://szaghi.github.io/shore/adr/0005-marching-schemes));
[Contributing](https://szaghi.github.io/shore/guide/contributing) explains the
workflow, and the [design records](https://szaghi.github.io/shore/adr/) what
was decided and why.

## Using SHORE with Claude Code

This repo ships a [Claude Code](https://claude.com/claude-code) skill at
[`.claude/skills/shore/`](.claude/skills/shore/) that teaches the agent how to
*operate* the `shore` CLI end-to-end — topology choice, the spacing knobs, and
the full STL → mesh → Xall-bundle (`grd` / `cc-par` / `proc-input`) workflow,
including the argument-order gotchas. Claude Code discovers it automatically
when you open this repository.

To make it available in **every** project (so you can drive `shore` from any
directory), copy it into your user skills directory:

```bash
cp -r .claude/skills/shore ~/.claude/skills/shore
```


## Authors

**[Stefano Zaghi](https://github.com/szaghi)** · stefano.zaghi@cnr.it
> *Chief Yak Shaver, Accidental Research Scientist, and HPC Farmer* — CFD researcher who decided that one more day debugging Fortran build systems was one day too many, opened a Python REPL "just to prototype," and now finds himself maintaining a meshing library, a chimera assembler, an MPI load balancer, and the seven blog tabs he keeps meaning to read.

**Andrea Di Mascio** · andrea.dimascio@univaq.it
> *Sommo Vate & Resident Gandalf* — responsible for the deep CFD intuitions SHORE pretends to know.

**[Claude](https://claude.ai)** (Anthropic)
> *Omniscient Code Oracle & Tireless Rubber Duck* — AI pair programmer, responsible for writing the boring parts so humans don't have to.

## License

SHORE is released under the [MIT License](LICENSE).
