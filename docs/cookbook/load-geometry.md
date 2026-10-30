# Load STL and CAD geometry

Every SHORE command reads its body through one `Geometry` interface, whether
the file is a triangulated STL or a STEP, IGES or BREP model. These recipes
load both, list what a CAD model contains, and show when the exact CAD
projection matters.

The code is
[`examples/cookbook/load_geometry.py`](https://github.com/szaghi/shore/blob/main/examples/cookbook/load_geometry.py)
(about 4 s). CAD input needs the `cad` extra: `pip install "shore-mesh[cad]"`
(build123d and OCP); STL needs nothing extra.

## Load a file and list its entities

The CAD test model is a half wing–body built with build123d and written in
the three CAD formats:

<<< @/../examples/cookbook/load_geometry.py#make_cad{python}

`load_geometry` picks the backend from the suffix: `.stl` gives an
`StlGeometry`, `.step`/`.stp`, `.iges`/`.igs` and `.brep`/`.brp` give an
`OcctGeometry`.

<<< @/../examples/cookbook/load_geometry.py#entities{python}

```text
wing_body.step  {'body': 1, 'surface': 8, 'curve': 19, 'point': 12}
wing_body.brep  {'body': 1, 'surface': 8, 'curve': 19, 'point': 12}
wing_body.igs   {'body': 1, 'surface': 8, 'curve': 19, 'point': 12}
bounds [[-3.0, -0.5, -0.0], [3.0, 0.5, 3.0]]
tessellation: 5078 triangles over 8 faces
```

<Figure name="cad-entities" caption="The CAD wing–body: faces coloured per face, edges as dark curves." />

- Entities have four kinds: `BODY` (always one), `SURFACE` (CAD faces),
  `CURVE` (CAD edges) and `POINT` (CAD vertices). An `EntityRef(kind,
  index)` names one; the index is what `Project.associate(...,
  indices=...)` and the GUI's Association dock use.
- An STL has one `BODY` and no faces. Its `CURVE` and `POINT` entities are
  the feature curves found at load time: see
  [Feature curves of an STL](/cookbook/stl-curves).
- IGES stores faces as independent trimmed surfaces; the loader sews them
  (tolerance `1e-6`) so adjacent faces share their edges, which association
  needs. The counts above match the STEP file after sewing.
- A CAD file without the `cad` extra raises `MeshInputError` with the
  `pip install` hint; an unknown suffix raises too.

## Choose between tessellated and exact projection

A CAD geometry answers `project` from a tessellation by default (tier 1) and
from the exact B-rep with `exact=True` (tier 2). The recipe projects 200
points onto a unit sphere three ways:

<<< @/../examples/cookbook/load_geometry.py#exact{python}

```text
STL        max |r - 1| = 4.3e-03
CAD tier 1 max |r - 1| = 2.6e-03
CAD exact  max |r - 1| = 1.1e-16
(u, v) of the first point on face 0: [5.4731 1.2934]
```

- The tier-1 tessellation has a linear deflection of `1e-3` times the
  bounding-box diagonal (`OcctGeometry(rel_deflection=...)` or
  `deflection=` changes it). Its error is the same order as a fine STL's.
- Use the exact tier when the wall spacing is near the faceting error: a
  first cell of `1e-5` on a body of unit size cannot sit on a surface that
  is wrong by `1e-3`. `Project.snap(exact=True)` and
  `Project.premesh(exact=True)` put the nodes on the exact surface.
- Projecting onto one `SURFACE` returns its `(u, v)` parameters in
  `Projection.params`; projecting onto a `CURVE` returns the curve
  parameter. An STL has no parametrisation (`params` is `None`).
- On a trimmed face, when the foot point on the underlying surface falls
  outside the trim, the exact projection returns the closest point of the
  face itself, which lies on its boundary.

## Name geometries in a project

A project holds several geometries under names; the commands refer to them by
name (`geometry="wing_body"`). `Body.from_file` is the single-body entry
point the marching topologies use.

<<< @/../examples/cookbook/load_geometry.py#named{python}

```text
wing_body: OcctGeometry, CAD=True
ball: StlGeometry, CAD=False
anchor [0. 0. 0.] CAD True
```

- `Project.load_geometry(path, name=...)` defaults to the name `"body"`.
  Loading a second file under the same name replaces the first and
  discards the meshes of the components built on it (undo restores them).
- `Body.from_file` also resolves the anchor, the interior point from which
  the cubed-sphere and O-grid topologies cast rays. On CAD the anchor is
  found on the merged tessellation. A body without a valid anchor (an open
  half model) still loads; the error is raised only when something needs
  the ray caster.
- The CLI takes CAD wherever it takes STL; for example a cubed sphere on
  the STEP sphere:

  ```bash
  shore mesh sphere.step --ni 33 --nj 64 --nk 21 --theta-cap-deg 35 -o sphere_cad
  ```

## See also

- [Geometry](/concepts/geometry): the backends, the two projection tiers,
  entity fingerprints.
- [Associate with geometry](/cookbook/associate) and
  [Rebuild a face on a CAD surface](/cookbook/cad-face): what the entity
  indices are for.
- [ADR 0002](/adr/0002-geometry-backend): one geometry protocol, two
  backends.
- API: [`load_geometry`](/reference/api/shore.geometry#load_geometry),
  [`OcctGeometry`](/reference/api/shore.geometry.occt#OcctGeometry),
  [`StlGeometry`](/reference/api/shore.geometry.stl#StlGeometry),
  [`Body`](/reference/api/shore.mesh.body#Body).
