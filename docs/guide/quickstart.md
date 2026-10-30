# Quick start

Ten minutes from an empty directory to a checked near-body grid and the
journal that rebuilds it. You need SHORE installed
([installation](/guide/installation)); `shore export` also needs the `vis`
extra.

## 1. Make a body

Any closed STL will do. SHORE can make test bodies itself:

```bash
shore body sphere -o sphere.stl --subdivisions 4
```

## 2. Mesh it

```bash
shore mesh sphere.stl --ni 33 --nj 64 --nk 25 --ds 0.01 --growth 1.15 \
    --theta-cap-deg 35 -o sphere --journal sphere.py
```

This builds the default topology, a six-block **cubed sphere** (four
equatorial blocks and two polar caps), and marches 24 layers out from the
wall: first cell `--ds 0.01`, each next one `--growth 1.15` times thicker.
`--theta-cap-deg` sets the size of the polar caps. It writes
`sphere_sub0.geo` … `sphere_cap_south.geo` and the journal `sphere.py`.

::: tip Nodes, not cells
`--ni`, `--nj` and `--nk` count **nodes**: a block is a `(ni, nj, nk, 3)`
array. Under the cubed sphere `--nj` is the node count of the whole equator
ring, shared by the four equatorial blocks. Only the binary `.grd` stores
cell counts, and SHORE converts them for you.
:::

## 3. Check it

```bash
shore check sphere_sub0.geo
shore info sphere_sub0.geo
```

`check` reports the minimum cell Jacobian layer by layer; a grid SHORE
writes never has an inverted cell, because marching stops at the first one.
[Grid quality](/concepts/quality) explains the metrics.

## 4. Look at it

```bash
shore export sphere_*.geo -o sphere.vtm   # all six blocks, for ParaView
shore view sphere_sub0.geo                # or one block in an interactive window
```

## 5. Replay it

The journal is a Python script of `Project` calls:

```python
# SHORE journal: replay with `shore run <this file>`.
from shore.project import Project

p = Project()
p.load_geometry('sphere.stl')
p.build_topology('cubed_sphere', nk=25, ds=0.01, ni=33, nj=64, theta_cap=0.6108652381980153, ...)
p.march()
p.write_geo('sphere')
```

`shore run sphere.py` rewrites the same blocks, bit for bit. Change a
number in it, point `load_geometry` at another body, and replay: that is
how SHORE is scripted. The same calls work in your own Python:

```python
from shore.project import Project

p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=25, ni=33, nj=64, ds=0.01, growth=1.15, theta_cap=0.61)
p.march()
p.write_geo("sphere")
```

## Next

- [Tutorial 1](/tutorials/sphere) takes this grid further: spacing laws,
  quality, VTK and the marching schemes.
- [Choosing a workflow](/guide/workflows) maps your body to the right tool.
- The [cookbook](/cookbook/) has one recipe per task.
