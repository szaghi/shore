# Cookbook

Short recipes for one task each. Every recipe's code is a script in
[`examples/cookbook/`](https://github.com/szaghi/shore/tree/main/examples/cookbook)
that runs as it stands and writes its output to the current directory;
the pages show the output you should get. Recipes marked CAD need the
`cad` extra (`pip install "shore-mesh[cad]"`).

For the ideas behind them see [Concepts](/concepts/); for every flag and
argument see the [CLI reference](/reference/) and the
[Python API](/reference/api/).

## Geometry

- [Write an analytic test body as STL](/cookbook/test-bodies#write-an-analytic-body-as-stl)
- [Generate a NACA section, sharp or blunt TE](/cookbook/test-bodies#generate-a-naca-section)
- [Extrude a wing with taper, sweep and twist](/cookbook/test-bodies#extrude-a-wing-with-taper-sweep-and-twist)
- [Load STEP, IGES or BREP and list its faces, edges and vertices](/cookbook/load-geometry#load-a-file-and-list-its-entities) (CAD)
- [Choose between tessellated and exact CAD projection](/cookbook/load-geometry#choose-between-tessellated-and-exact-projection) (CAD)
- [Load several geometries under names](/cookbook/load-geometry#name-geometries-in-a-project) (CAD)
- [List the feature curves and points of an STL](/cookbook/stl-curves#list-the-curves-and-points)
- [Choose the crease angle for STL features](/cookbook/stl-curves#choose-the-crease-angle)
- [Snap a grid edge onto an STL feature curve](/cookbook/stl-curves#snap-a-grid-edge-onto-a-feature-curve)

## Spacing

- [Reach a boundary-layer thickness with geometric layers](/cookbook/wall-spacing#reach-a-thickness-with-geometric-layers)
- [Hit a target first cell with a tanh or tanh2 law](/cookbook/wall-spacing#fix-the-thickness-with-a-tanh-law)
- [Size the first cell for a target y+](/cookbook/wall-spacing#choose-the-first-cell-from-a-target-y)
- [Know which k-schedule each block marches with](/cookbook/wall-spacing#know-which-schedule-a-block-uses)
- [Cluster one edge of a box primitive](/cookbook/per-edge-spacing#cluster-one-edge-of-a-box)
- [Set the k-schedule of a cubed sphere or O-grid edge by edge](/cookbook/per-edge-spacing#set-the-k-schedule-of-a-cubed-sphere)
- [Cluster the nodes along one block edge of an existing grid](/cookbook/redistribute#cluster-the-nodes-along-one-block-edge)
- [Redistribute an edge through the seams of a marched grid](/cookbook/redistribute#redistribute-through-seams-on-a-marched-grid)
- [Know what a redistribution refuses](/cookbook/redistribute#know-what-is-refused)

## Geometry association and surfaces

- [Find the CAD entities to associate with](/cookbook/associate#find-the-geometry-entities) (CAD)
- [Associate a face, an edge and a vertex, then snap](/cookbook/associate#associate-a-face-an-edge-and-a-vertex-then-snap) (CAD)
- [Clear or undo an association](/cookbook/associate#clear-or-undo-an-association) (CAD)
- [Find the CAD surface under a block face](/cookbook/cad-face#find-the-surface-under-a-block-face) (CAD)
- [Rebuild a block face on a CAD surface](/cookbook/cad-face#fill-the-face-interior-on-the-surface) (CAD)

## Half models

- [Cut a body at its symmetry plane](/cookbook/half-models#cut-a-body-into-a-half-body)
- [Build a half-collar at a junction on the symmetry plane](/cookbook/half-models#build-a-half-collar-at-the-junction)
- [Write SYMMETRY faces to `cc.par`](/cookbook/half-models#write-the-symmetry-faces-to-cc-par)

## Blocking

- [Split a blocking block at a node](/cookbook/blocking#split-a-block-at-a-node)
- [Delete a block to make a hole round the body](/cookbook/blocking#delete-a-block-to-make-a-hole-round-the-body)
- [Insert an O-grid around a body in a box](/cookbook/blocking#insert-an-o-grid-around-the-body)
- [Make grid lines leave the wall along its normal](/cookbook/blocking#make-grid-lines-leave-the-wall-along-its-normal)
- [Cluster a wall layer with an absolute first cell](/cookbook/blocking#cluster-a-wall-layer-with-an-absolute-first-cell)
- [Move a blocking vertex, freely or onto the geometry](/cookbook/blocking#move-a-vertex-freely-or-onto-the-geometry)

## Quality

- [Check a `.geo` file for folded cells](/cookbook/quality#check-a-geo-file-for-folded-cells)
- [Summarise the quality of a whole grid](/cookbook/quality#summarise-a-whole-grid)
- [Compare blocks and find the worst cell](/cookbook/quality#compare-blocks-and-find-the-worst-cell)
- [Compare march settings on wall orthogonality](/cookbook/quality#compare-march-settings)

## Solver files

- [Export a project as CGNS](/cookbook/exports#write-a-cgns-file)
- [Export a project as Plot3D](/cookbook/exports#write-a-plot3d-grid)
- [Export to VTK for ParaView](/cookbook/exports#write-vtk-for-paraview)
- [Write one Xall run from several components](/cookbook/multi-component#write-the-run-with-an-options-file)
- [Merge separately built `.grd`, `cc.par` and `proc.input` files](/cookbook/multi-component#merge-separate-component-files)
- [Split blocks by hand, keeping the seams](/cookbook/split-balance#split-blocks-by-hand)
- [Plan a balanced split for N ranks](/cookbook/split-balance#plan-a-balanced-split)
- [Keep enough cells for the multigrid levels](/cookbook/split-balance#keep-the-multigrid-floor)
- [Split after the overset pre-processor with `bsplit`](/cookbook/after-overset-split)

## Projects

- [Undo and redo a command](/cookbook/projects#undo-and-redo-a-command)
- [Save a project and load it back](/cookbook/projects#save-a-project-and-load-it-back)
- [Load a large project fast, without replay](/cookbook/projects#load-a-project-fast-without-replay)
- [Edit a journal and replay it](/cookbook/projects#edit-a-journal-and-replay-it)
- [Build a project from the CLI](/cookbook/projects#build-a-project-from-the-cli)
