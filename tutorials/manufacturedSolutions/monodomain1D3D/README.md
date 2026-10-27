# monodomain1D3D

Coupled 1D-3D monodomain manufactured solution: a reaction-diffusion Purkinje
graph and the 3D myocardium on the unit cube, joined at two Purkinje-ventricular
junctions (PVJs) by `reactionDiffusionPvjCoupler`. It checks the node-edge graph
operator, the 3D monodomain operator and the PVJ source (sampling the tissue
voltage, depositing the junction current) under joint refinement, for
unidirectional transfer with explicit and implicit tissue-side assembly and for
bidirectional transfer, which is always explicit on the tissue side.

The junction forcing has two references (`-reference`). `discrete` (the
default) builds the exact junction source with the same sampling and
deposition map as the solver, so map errors cancel and the case checks the
coupled workflow against that map. `continuum` builds it from the continuum
kernel normalised over the tissue part of the support and averaged over each
cell (`pvjContinuumReference`), so the discrete map's own error reaches the
solution and the junction is verified against its continuum definition. The
case does not check a calibrated junction resistance, graph branch points with
more than two edges, or a wall treatment (the 3D field has a zero normal
gradient on every face); the sealed curved wall is
[monodomain1D3DShell](../monodomain1D3DShell/README.md).

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
as `S_pvj` with a kernel of radius 0.11 (linear by default, `-kernel`)
normalised to unit integral.
`unidirectional` keeps only the tissue source; `bidirectional` also applies
`I_app1D = -I_pvj,k` at the graph terminal, which the network loses from the
volume `π ρ² L` of its terminal node, so the case sets `purkinjeFibreRadius`
to `1/√π` m (the unit cross-section the manufactured junction current assumes)
instead of the 20 µm default. `explicit` deposits the terminal current;
`implicit` assembles a cell-local sink in the tissue voltage (unidirectional
only; bidirectional transfer is always explicit). With `-reference discrete` the
verifiers add the exact junction source, computed from the exact fields with
the same map and kernel, to the manufactured sources and write the difference
between the numerical and the exact junction source. With `-reference
continuum` they add the continuum junction source instead (below).

### Continuum reference

For a terminal `t` with radius `R` and kernel `w`, over the tissue part of the
ball `Omega cap B`:

```text
W   = int_{Omega cap B} w dx,     <F> = int_{Omega cap B} w F dx / W
I_c = sqrt(1 + t) (F1D(t) - <F>)/R_pvj
S_c = I_c w/W                       (explicit: the gathered operator)
S_c = w (V1D(t) - Vm(x))/(R_pvj W)  (implicit: the local operator)
```

`W` and `<F>` come from Gauss quadrature in spherical coordinates about the
terminal (exact for a terminal on a flat face and for one on the inner sphere
of the shell); the tissue forcing is the average of `S_c` over each cell,
integrated over `K cap B` by tetrahedral subdivision with the sphere clipped
at a cut size of `R/96` (uniform, gaussian), `R/48` (linear) or `R/24`
(cosine). Both forcings are evaluated at the new time level `t0 + dt` of each
step, so the staggered lag of the graph behind the tissue is part of the
error. The kernels are `uniform` (`w = 1`), `linear` (`1 - r/R`), `gaussian`
(`exp(-4.5 r^2/R^2)` cut at `R`) and `cosine` (`(1 + cos(pi r/R))/2`), written
in `pvjContinuumReference` from these definitions and not shared with
`pvjMapper`.

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
| `-kernel uniform\|linear\|gaussian\|cosine` | PVJ `pvjKernel` | linear |
| `-radius x` | PVJ `pvjRadius` | 0.11 |
| `-reference discrete\|continuum` | junction forcing of `coupled1D3DMonodomainVerifier` | discrete |
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
- under both references (the diagnostics file above is `-reference discrete`
  only, the steps file `-reference continuum` only):
  - `verification/coupled1D3DMonodomain_map.csv`: for every kernel, the
    discrete map against the continuum one at `t = 0`: normalisation
    `Wh/Wc - 1`, gathered exact field `<F>h - <F>c`, deposited source per
    unit current against the cell averages of the continuum source (L1, L2,
    relative L2), and `sum omega_K/Wc - 1` (the mesh-integrated against the
    exact normalisation: geometry and quadrature of the reference)
  - `verification/coupled1D3DMonodomain_steps.csv` (continuum): every step, the junction
    charge leaving the graph (from the graph's own charge budget over the
    step) and entering the tissue (from the deposited source), their
    difference relative to the charge crossing all junctions, and each
    terminal's deposited and continuum currents
  - `verification/coupled1D3DMonodomain_local.csv`: final-time tissue errors
    over the cells within `R` and `2R` of each terminal, and the errors of the
    graph terminal voltage, the gathered tissue voltage and the junction
    current against the continuum ones (the step maxima are continuum only)

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
| [continuumCoupling](setup/studies/continuumCoupling/README.md) | the junction against its continuum definition: map convergence of every kernel, support-local orders, a junction-suppression control, charge balance and the time-step order |
