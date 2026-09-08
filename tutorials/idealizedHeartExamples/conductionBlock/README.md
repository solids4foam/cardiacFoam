# Bundle branch block on the idealized biventricular heart

Monodomain tissue simulation modelling left or right bundle branch block via
a modified Purkinje graph with one bundle severed. No `ionicConstantOverrides`
are used — the conduction delay is entirely structural, same approach as
`tutorials/PATHOS/LBBB`/`RBBB`. Two variants of the same case, selected at
run time; no default — see Execution below.

## Stack

- electro model: `monodomainSolver`
- ionic model: `BuenoOrovio`
- conduction system: `purkinjeGraphModel` with `monodomain1DSolver`, graph
  `purkinjeGraph` (a symlink to `purkinjeGraph.lbbb` or `purkinjeGraph.rbbb`,
  set by `Allrun`)
- ECG: `ecgDomains.ECG.ecgSolver pseudoECG`

`constant/electroProperties` itself never changes between variants — it
always reads `graphFile purkinjeGraph;`, same as `tutorials/PATHOS/RBBB`'s
approach. Only the graph file differs.

## Mesh and shared anatomy fields

Shared with `electrophysiologyHeart`/`electromechanicalHeart`/
`ionicPathology` via `../mesh/` — see `../README.md`. `Allrun` copies in
`fiber`/`sheet`/`sheetNormal`/`tm`/`tv`/`apicobasal`/`Conductivity`/
`polyMesh`, and the baseline `constant/purkinjeGraph` (as
`purkinjeGraph.healthy`, the source both variants are derived from).

## Tissue heterogeneity

`monodomainSolverCoeffs.ionicHeterogeneity` classifies cells into
`endocardialCells`/`mCells`/`epicardialCells` bands (`endoMInterface 0.3`,
`mEpiInterface 0.7`) from `field t;`, `system/setExprFieldsDict`'s
`t = 1 - tm` (`tm`: 0 at epicardium, 1 at endocardium — the opposite
orientation `ionicHeterogeneity` requires).

## Bundle-branch-block graph derivation

`Allrun` derives both `constant/purkinjeGraph.lbbb` and
`constant/purkinjeGraph.rbbb` from `purkinjeGraph.healthy` every run (a
single `awk` pass, no Python dependency), then symlinks
`constant/purkinjeGraph` to whichever was requested.

`purkinjeGraph`'s `rootNode` has exactly two direct branches — the left-
and right-bundle roots. Zeroing either branch's root-adjacent edge
conductance fully disconnects that ventricle's tree from the root (it's a
tree — no alternate path exists). Same structural technique
`tutorials/PATHOS/*/make_bundle_block_graphs.py` uses; the edge indices
below were re-derived for this graph's own node numbering (by
BFS-partitioning the graph from each of `rootNode`'s two neighbours and
checking which partition contains the node coincident with the LV/RV
growth seed), not copied from that script's Strocchi-mesh-specific values.

- LBB bridge = `conductionEdges[0]` → disconnects the LV subtree (LBBB)
- RBB bridge = `conductionEdges[22]` → disconnects the RV subtree (RBBB)

`severBundleBranch()` in `Allrun` verifies the target edge connects the
expected node pair before zeroing it, and aborts instead of zeroing the
wrong edge if the baseline graph is ever regenerated with different node
numbering.

## ECG electrodes

`ecgDomains.ECG.electrodePositions` are transferred to this anatomy from a
reference heart's validated V1-V6 placement, not measured on this mesh — a
normalized approximation, not patient-specific placement. 24-46mm from the
epicardium, comparable to the reference case's own spread.

## Execution

```bash
./Allrun lbbb
./Allrun rbbb
./Allrun lbbb parallel
```
