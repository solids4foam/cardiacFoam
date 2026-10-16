# monodomain1D3D

Coupled 1D-3D monodomain manufactured solution: a reaction-diffusion Purkinje
graph and the 3D myocardium on the unit cube, joined at two Purkinje-ventricular
junctions (PVJs) by `reactionDiffusionPvjCoupler`. It checks the node-edge graph
operator, the 3D monodomain operator and the PVJ source (sampling the tissue
voltage, depositing the junction current) under joint refinement, for
unidirectional and bidirectional transfer and explicit and implicit tissue-side
assembly. The exact junction source is built with the same sampling and
deposition map as the solver, so the case checks the coupled workflow against
that map; it does not check the map against an independent continuum coupling,
a calibrated junction resistance, graph branch points with more than two edges,
or a wall treatment (the 3D field has a zero normal gradient on every face).

Origin: the 3D field is the monodomain problem of the FDA regulatory science
tool (Pathmanathan & Gray 2014); the graph field and the PVJ coupling problem
are this work.

## Problem

```text
chi Cm dVm/dt = div(G grad Vm) - chi Cm Iion + I_s + S_pvj          in [0,1]^3
n . grad Vm = 0                                                      on every face
chi Cm dV1D/dt = d/ds(g dV1D/ds) - chi Cm Iion1D + I_s1D (+ I_app1D) on each graph edge
```

`Vm = sqrt(1 + t) cos(pi x) cos(2 pi y) cos(3 pi z)` and `V1D = sqrt(1 + t)
cos(pi x)` at graph position `x`, with `G = diag(0.111453302, 0.121585420,
0.030396355)`, graph conductivity `g = 0.111453302`, `chi = 3`, `Cm = 2` and
the FDA ionic closure (`u1`, `u2`, `u3`) in both domains; `I_s` and `I_s1D`
are the manufactured sources added by the verifiers.

The graph runs from the root `(1/2, 1/6, 1/3)` along two straight branches to
the terminals `(0, 1/6, 1/3)` and `(1, 1/6, 1/3)`. There `d Vm/dx = 0`, so the
3D no-flux condition holds, and `Vm = -V1D/2`, so the junction current does not
vanish: with `R_pvj = 1` the exact terminal currents are `+-(3/2) sqrt(1 + t)`.

At terminal `k` the junction current `I_pvj,k = (V1D_k - M_k[Vm])/R_pvj`
samples the tissue voltage with the volume-weighted map `M_k` and is deposited
as `S_pvj` with a linear kernel of radius 0.11 normalised to unit integral.
`unidirectional` keeps only the tissue source; `bidirectional` also applies
`I_app1D = -I_pvj,k` at the graph terminal, which the network loses from the
volume `π ρ² L` of its terminal node, so the case sets `purkinjeFibreRadius`
to `1/√π` m (the unit cross-section the manufactured junction current assumes)
instead of the 20 µm default. `explicit` deposits the terminal current;
`implicit` assembles a cell-local sink in the tissue voltage. The
verifiers add the exact junction source, computed from the exact fields with
the same map and kernel, to the manufactured sources and write the difference
between the numerical and the exact junction source.

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
| `-pvj explicit\|implicit` | PVJ `pvjCouplingScheme` | explicit |
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
| [coupledConvergence](setup/studies/coupledConvergence/README.md) | joint 1D-3D order of the fields and the junction source, for the four transfer and assembly combinations |
| [weakCouplingControl](setup/studies/weakCouplingControl/README.md) | whether the tissue and graph errors depend on the junction |
