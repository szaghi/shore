"""Projects and journals: two components, undo/redo, save, load, replay after an input edit.

Writes ``sphere.stl``, ``journal.py``, the project directory ``proj/``, the Xall
run ``run/``, ``proj.cgns`` and ``edited/sphere.stl`` in the current directory.
Takes about 10 s.
"""

# requires: cgns

from __future__ import annotations

from pathlib import Path

import numpy as np

from shore import shapes
from shore.project import Project
from shore.project.project import run_journal

# #region components
shapes.sphere(radius=0.5, subdivisions=3).export("sphere.stl")

p = Project()
p.load_geometry("sphere.stl", name="sphere")  # a named geometry, shared by all components
# component "main" (the default): an O-grid on the sphere, marched
p.build_topology("ogrid", nk=17, ds=2e-3, growth=1.2, geometry="sphere", ni=17, nj=32)
p.march()
# component "background": a Cartesian box round it
p.use_component("background")
p.build_primitive("box", vmin=(-3.0, -3.0, -3.0), vmax=(3.0, 3.0, 3.0), ni=13, nj=13, nk=13)
print("components:", {name: len(m.blocks) for name, m in p.components.items() if m})
# #endregion components

# #region undo
p.undo()  # drops the box
print("after undo:", p.mesh)
p.redo()  # builds it again
print("after redo:", [b.label for b in p.mesh.blocks])
# #endregion undo

# #region journal
print(p.journal())
p.save_journal("journal.py")  # `shore run journal.py` replays it
# #endregion journal

# #region save
p.save("proj")  # project.toml, journal.py, inputs/, grids/<component>/
q = Project.load("proj")  # replays the journal and checks every grid bit for bit
r = Project.load("proj", replay=False)  # reads the saved grids: no march
same = np.array_equal(q.components["main"].blocks[0].nodes, r.components["main"].blocks[0].nodes)
print("replay == saved grids:", same)
# #endregion save

# #region outputs
p.write_xall("run", np=4)  # per component under run/<component>/, then merged files
p.write_cgns("proj.cgns")  # a CGNS base per component
print(sorted(f.name for f in Path("run").iterdir()))
# #endregion outputs

# #region replay
# edit the input: a larger sphere, same file name, in another directory
Path("edited").mkdir(exist_ok=True)
shapes.sphere(radius=0.6, subdivisions=3).export("edited/sphere.stl")
e = run_journal("journal.py", input_dir="edited")  # 'sphere.stl' resolves in edited/
wall = e.components["main"].blocks[0].nodes[:, :, 0]
print(f"edited wall radius: {np.linalg.norm(wall, axis=-1).mean():.3f}")
# #endregion replay
