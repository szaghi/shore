"""Save, load and replay a project; load without replay; undo, redo; edit a journal.

Writes ``sphere.stl``, the project directory ``sphere.shore/``, the journals
``sphere_journal.py`` and ``sphere_fine.py``, and ``out/sphere_*.geo`` in the
current directory.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from shore import shapes
from shore.exceptions import ParameterError
from shore.project import Project, run_journal

shapes.sphere(radius=1.0, subdivisions=3).export("sphere.stl")

# #region undo_redo
p = Project()
p.load_geometry("sphere.stl")
p.build_topology(
    "cubed_sphere", nk=9, ds=0.02, growth=1.15, ni=33, nj=64, theta_cap=np.deg2rad(35.0)
)
p.march()
marched = [b.nodes.copy() for b in p.mesh.blocks]
p.redistribute("sub0", "j_lo_k_lo", "tanh2", beta=2.0)

p.undo()  # back to the marched grid, bit for bit
print(
    "undo restores:",
    all(np.array_equal(b.nodes, n) for b, n in zip(p.mesh.blocks, marched, strict=True)),
)
p.redo()  # and forward again
print([type(c).__name__ for c in p.history])
# #endregion undo_redo

# #region save_load
p.write_geo("out/sphere")
p.save("sphere.shore")  # project.toml, journal.py, inputs/, adjacency.json, grids/
print(sorted(f.name for f in Path("sphere.shore").iterdir()))

t0 = time.perf_counter()
q = Project.load("sphere.shore")  # replays the journal, checks every grid bit for bit
t1 = time.perf_counter()
r = Project.load("sphere.shore", replay=False)  # reads the saved grids instead
t2 = time.perf_counter()
print(f"replayed load {t1 - t0:.2f} s, load without replay {t2 - t1:.2f} s")
same = all(
    np.array_equal(a.nodes, b.nodes) for a, b in zip(q.mesh.blocks, r.mesh.blocks, strict=True)
)
print("same blocks:", same, "same journal:", q.journal() == r.journal())
# #endregion save_load

# #region no_replay_limits
try:
    r.undo()  # the loaded history was recorded, not run: it cannot be undone
except ParameterError as err:
    print("refused:", err)
try:
    r.march()  # nor can a loaded mesh be marched again
except ParameterError as err:
    print("refused:", err)
r.snap()  # new commands work, and can be undone
r.undo()
# #endregion no_replay_limits

# #region edit_journal
p.save_journal("sphere_journal.py")
print(Path("sphere_journal.py").read_text())

# a finer wall-normal schedule: edit the script, replay it
text = Path("sphere_journal.py").read_text()
Path("sphere_fine.py").write_text(text.replace("nk=9, ds=0.02", "nk=17, ds=0.005"))
fine = run_journal("sphere_fine.py")  # or: shore run sphere_fine.py
print("fine:", fine.mesh.blocks[0].nodes.shape)
# #endregion edit_journal
