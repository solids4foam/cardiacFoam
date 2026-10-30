# monodomain1D3D

Coupled 1D-3D monodomain manufactured solution: a reaction-diffusion Purkinje
graph and the 3D myocardium on the unit cube, joined at two Purkinje-ventricular
junctions (PVJs) by `reactionDiffusionPvjCoupler`. It checks the node-edge graph
operator, the 3D monodomain operator and the PVJ source (sampling the tissue
voltage, depositing the junction current) under joint refinement, for
unidirectional transfer with explicit and implicit tissue-side assembly and for
bidirectional transfer, which is always explicit on the tissue side. The exact
junction source is built with the same sampling and deposition map as the
solver, so the case checks the coupled workflow against that map; it does not check the map against an independent continuum coupling,
a calibrated junction resistance, graph branch points with more than two edges,
or a wall treatment (the 3D field has a zero normal gradient on every face).

Origin: the 3D field is the monodomain problem of the FDA regulatory science
tool (Pathmanathan & Gray 2014); the graph field and the PVJ coupling problem
are this work.

## Problem

The tissue field, the graph field and the junction are defined in the
construction document
[monodomain_mms.pdf](../generators/docs/monodomain/monodomain_mms.pdf), sections
"Coupled 1D--3D junction", "Purkinje graph" and "FDA cube", with the generators
[`coupled_1d3d_mms.py`](../generators/coupled_1d3d_mms.py),
[`graph_monodomain_mms.py`](../generators/graph_monodomain_mms.py) and
[`fda_monodomain_mms.py`](../generators/fda_monodomain_mms.py) and the NumPy
references [`coupled_1d3d.py`](../generators/reference/coupled_1d3d.py),
[`graph_monodomain.py`](../generators/reference/graph_monodomain.py) and
[`fda_monodomain.py`](../generators/reference/fda_monodomain.py). Neither the
graph nor the tissue receives a bulk manufactured source; the applied terms are
the exact junction source and, for bidirectional coupling, the graph terminal
term.

The junction current at a terminal is deposited as `S_pvj` with a linear kernel
of radius `pvjRadius 0.11` and `R_pvj = 1` (`-rPvj`). `unidirectional` keeps
only the tissue source; `bidirectional` also applies the negative junction
current at the graph terminal, which the network loses from the volume of its
terminal node, so the case sets `purkinjeFibreRadius` to `1/sqrt(pi)` m (the
unit cross-section the manufactured junction current assumes) instead of the
20 um default. `explicit` deposits the terminal current; `implicit` assembles a
cell-local sink in the tissue voltage (unidirectional only; bidirectional
transfer is always explicit). The verifiers add the exact junction source,
computed from the exact fields with the same map and kernel, and write the
difference between the numerical and the exact junction source.

### Stack

- myocardium solver: `monodomainSolver`, `solutionAlgorithm implicit`,
  `timeCouplingScheme godunov`, `ddt(Vm) backward`, `Gauss linear` gradient,
  `staggeredElectrophysicsAdvanceScheme`
- graph: `purkinjeGraphModel` with `monodomain1DSolver` (reaction step, then
  implicit node-edge diffusion solved by Hines elimination)
- coupler: `reactionDiffusionPvjCoupler`, `couplingMode` (`-coupling`),
  `pvjCouplingScheme` (`-pvj`), `rPvj` (`-rPvj`), `pvjRadius 0.11`,
  `pvjKernel linear`
- ionic model: `monodomainFDAManufactured` in both domains
- verifiers: `manufacturedFDAMonodomainVerifier` (3D),
  `manufacturedGraphVerifier` (graph), `coupled1D3DMonodomainVerifier` (PVJ
  source)

## Mesh

- 3D (`system/blockMeshDict.3D`): `N` cells per unit length in each direction.
- Graph (`constant/purkinjeGraph.nodes<k>`, `k` = 3, 11, 21, 41, 81, 161):
  `(k - 1)/2` equal edges per branch, edge length `1/(k - 1)`. The coupled
  studies refine both together, `k = N + 1`. `setup/generate_purkinje_graphs.py`
  writes these files to `constant/` and VTK views to `constant/graphFiles/`.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

Every option writes a file in `local/`, which the dictionaries include;
`Allclean` removes `local/`, the mesh and the verifier outputs. `Allrun` stops
with `run ./Allclean first` if a mesh log exists.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` (6 subdomains, scotch), parallel `cardiacFoam`, `reconstructPar -latestTime` | serial |
| `-N n` | cells per unit length of the cube | 20 |
| `-graph nodes003\|...\|nodes161` | Purkinje graph | `nodes<N+1>` (edge length `1/N`) |
| `-domain coupled\|graph` | `cardiacFoam` on the coupled problem, or `runPurkinjeGraph` on the graph alone (no junction, serial) | coupled |
| `-coupling unidirectional\|bidirectional` | PVJ `couplingMode` | unidirectional |
| `-pvj explicit\|implicit` | PVJ `pvjCouplingScheme` (unidirectional only) | explicit |
| `-rPvj x` | junction resistance | 1 |
| `-dt x` | time step | `endTime/ceil(endTime/(0.008971136 (10/N)^2))` |
| `-endTime x` | end time | 0.1 |

Outputs:

- `postProcessing/3D_<n>_cells.dat`: 3D errors of `Vm`, `u1`, `u2` and the
  bulk, boundary and gradient rows of
  [Bulk, boundary and gradient rows](../README.md#bulk-boundary-and-gradient-rows)
- `postProcessing/graph_1D_<k>_nodes.dat`: graph errors of `Vm1D`, `u1`, `u2`
- `verification/coupled1D3DMonodomain_diagnostics.csv`: at the last step, the
  numerical and exact junction sources (`sourceL1`, `exactSourceL1`), their
  difference (`sourceErrorL1`, `sourceErrorL2`) and the deposited currents
- `postProcessing/purkinjeNetwork.dat`: graph node voltages and terminal
  sources at the write times

Regression (the defaults: `N = 20`, the 21-node graph, unidirectional explicit
transfer, serial): the 3D errors with their boundary and gradient rows, the
junction source, its exact value and its error, and the graph errors.

```bash
./regression/regressionTest.sh
```

## Studies

`setup/studies/` holds one study per question, each run in copies of the case
with the commands in its README; `setup/studies/summarize.py` prints the
tables.

| study | question |
|---|---|
| [graphConvergence](setup/studies/graphConvergence/README.md) | spatial order of the graph operator alone, at a fixed time step |
| [coupledConvergence](setup/studies/coupledConvergence/README.md) | joint 1D-3D order of the fields and the junction source, for unidirectional transfer with both assemblies and for bidirectional transfer |
| [weakCouplingControl](setup/studies/weakCouplingControl/README.md) | whether the tissue and graph errors depend on the junction |
