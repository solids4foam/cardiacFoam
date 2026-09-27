# Bundle branch block on the idealized biventricular heart

Monodomain tissue simulation modelling left or right bundle branch block via
a modified Purkinje graph with one bundle severed. No `ionicConstantOverrides`
are used — the conduction delay is entirely structural. Two variants of the
same case, selected at run time; no default — see Execution below.

## Stack

- electro model: `monodomainSolver`
- ionic model: `BuenoOrovio`
- conduction system: `purkinjeGraphModel` with `monodomain1DSolver`, graph
  `purkinjeGraph` (a symlink to `purkinjeGraph.lbbb` or `purkinjeGraph.rbbb`,
  set by `Allrun`)
- ECG: `ecgDomains.ECG.ecgSolver pseudoECG`

`constant/electroProperties` itself never changes between variants — it
always reads `graphFile purkinjeGraph;`. Only the graph file differs.

Purkinje conduction and junction coupling: `purkinjeConductivity 0.4`
(~3.3 m/s along the tree) and `pvjRadius 1.65e-3`, the smallest junction
radius valid on this mesh.

## Mesh and anatomy fields

`Allrun` copies from `../../mesh/`:
`fiber`/`sheet`/`sheetNormal`/`tm`/`tv`/`apicobasal`/`Conductivity`/
`polyMesh`, and the baseline `constant/purkinjeGraph` (as
`purkinjeGraph.healthy`, the source both variants are derived from).

## Tissue heterogeneity

`monodomainSolverCoeffs.ionicHeterogeneity` (`mode namedRegions;`) classifies
cells into `endocardialCells`/`mCells`/`epicardialCells` regions
(`range (0 0.3)`/`range (0.3 0.7)`/`range (0.7 1)`) from `field t;`,
`system/setExprFieldsDict`'s
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
tree — no alternate path exists). The edge indices below were derived by
BFS-partitioning the graph from each of `rootNode`'s two neighbours and
checking which partition contains the node coincident with the LV/RV
growth seed.

- LBB bridge = `conductionEdges[0]` → disconnects the LV subtree (LBBB)
- RBB bridge = `conductionEdges[22]` → disconnects the RV subtree (RBBB)

`severBundleBranch()` in `Allrun` verifies the target edge connects the
expected node pair before zeroing it, and aborts instead of zeroing the
wrong edge if the baseline graph is ever regenerated with different node
numbering.

## ECG electrodes

`ecgDomains.ECG.electrodePositions` are placed by angle around the LV long
axis in the LV frame (`L` apex-to-base, `S` LV centre to RV centre, anterior
`A = L x S`, here `-z`): V1..V6 at 35, 65, 100, 135, 170, 205 deg from `S`
toward `A`, each at its original apex-base height and 25 mm from the nearest
tissue.

## Execution

```bash
./Allrun lbbb
./Allrun rbbb
./Allrun lbbb parallel
```
