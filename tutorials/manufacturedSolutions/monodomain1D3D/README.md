# monodomain1D3D

Manufactured 1D-3D monodomain coupling test case.

The case couples a 3D monodomain myocardium domain to a small 1D Purkinje
graph through `reactionDiffusionPvjCoupler` and verifies the coupled fields with
`coupled1D3DMonodomainVerifier`.

The active 1D graph input is `constant/purkinjeGraph`.  Refined graph inputs and
matching VTK geometry files can be generated with:

```bash
setup/generate_purkinje_graphs.py
```

By default the graph lies on `y=1/6, z=1/3`, with PVJ terminals at
`(0, 1/6, 1/3)` and `(1, 1/6, 1/3)`.  The terminal faces still satisfy the
homogeneous Neumann boundary condition for the 3D manufactured solution because
the terminals are on `x=0` and `x=1`, where `d cos(pi x)/dx = 0`.  Unlike the
older `y=0.5, z=1/3` placement, the terminal values do not cancel:
`F_3D = -0.5 F_1D` at the PVJs.

This creates:

- `constant/purkinjeGraph.nodes003`
- `constant/purkinjeGraph.nodes011`
- `constant/purkinjeGraph.nodes021`
- `constant/purkinjeGraph.nodes041`
- `constant/purkinjeGraph.nodes081`
- `constant/purkinjeGraph.nodes161`
- `constant/graphFiles/purkinjeGraph.nodes*.vtk`

Select one graph as the active input with:

```bash
setup/select_purkinje_graph.sh nodes041
```

Run the case from this directory:

```bash
./Allrun
```

`Allrun` regenerates `constant/polyMesh` from `system/blockMeshDict.3D` before
launching `cardiacFoam`.

For a graph-only diagnostic, run:

```bash
blockMesh -dict system/blockMeshDict.3D
runPurkinjeGraph -case .
```

For graph-only manufactured convergence rates, run:

```bash
setup/run_coupling1D3D_hex.sh
```

The sweep selects each `constant/purkinjeGraph.nodes*` input, runs
`runPurkinjeGraph`, copies the graph verifier summaries, and writes:

- `outputs/1dGraphConvergence/graph_convergence_summary.csv`
- `outputs/1dGraphConvergence/graph_convergence_rates.csv`

For coupled 1D-3D manufactured convergence rates, run:

```bash
setup/run_coupling1D3D_hex.sh
```

The coupled sweep can exercise the explicit and implicit PVJ source split without
manual dictionary edits:

```bash
PVJ_COUPLING_SCHEME=explicit OUTPUT_SUFFIX=_pvjExplicit setup/run_coupling1D3D_hex.sh
PVJ_COUPLING_SCHEME=implicit OUTPUT_SUFFIX=_pvjImplicit setup/run_coupling1D3D_hex.sh
```

The tissue diffusion algorithm can also be switched with
`SOLUTION_ALGORITHM=explicit` or `SOLUTION_ALGORITHM=implicit`.

To run the canonical coupling-scheme matrix for the new 1D-3D tests, use:

```bash
setup/run_coupling1D3D_hex.sh
```

This runs four full convergence sweeps:

- `uni_pvjExplicit`: `couplingMode unidirectional`, `pvjCouplingScheme explicit`
- `uni_pvjImplicit`: `couplingMode unidirectional`, `pvjCouplingScheme implicit`
- `bi_pvjExplicit`: `couplingMode bidirectional`, `pvjCouplingScheme explicit`
- `bi_pvjImplicit`: `couplingMode bidirectional`, `pvjCouplingScheme implicit`

All four canonical runs keep `solutionAlgorithm implicit`, because
`pvjCouplingScheme implicit` only becomes a true matrix split when the myocardium
diffusion solve is implicit.  For a quick smoke check of the same matrix, run:

```bash
COUPLED_1D3D_PAIRS="10:nodes011" ENDTIME=0.02 setup/run_coupling1D3D_hex.sh
```

The suite writes `outputs/coupled1D3DSchemeSuite/manifest.csv`.  Individual
convergence tables and plots are written under
`outputs/coupled1D3DConvergence_<case_id>/`.

The coupled sweep runs `cardiacFoam` under joint 1D/3D refinement, copies the
myocardium, graph, and PVJ coupling verifier summaries, and writes:

- `outputs/coupled1D3DConvergence/coupled_convergence_summary.csv`
- `outputs/coupled1D3DConvergence/coupled_convergence_rates.csv`
- `outputs/coupled1D3DConvergence/coupled_1D3D_convergence.png`
- `outputs/coupled1D3DConvergence/coupled_1D3D_convergence.pdf`

The coupled post-processing script also writes the same PNG/PDF convergence plot
for any alternate output directory passed via `--output-dir`, including the
bidirectional sweep directory.

## Convergence verification

**Why this study exists.** Before the coupled 1D-3D monodomain can be trusted in
the paper, the implementation has to be shown to converge at the expected rate
against a known solution. The Method of Manufactured Solutions (MMS) supplies that
known solution: an analytical `V_exact` is imposed on both domains, the matching
forcing terms are injected, and the solver error is measured under joint
`h_1D ≈ h_3D` refinement with `dt ~ h²` so temporal error never masks the spatial
order. The full setup, mesh/time-step pairing, per-sweep tables and the
root-cause analysis are recorded in
[`MMS_CONVERGENCE_NOTES.md`](MMS_CONVERGENCE_NOTES.md).

**Current interpretation.**

- The 1D Purkinje and 3D myocardium solvers each converge at **O(h²)** on their own
  (1D standalone, 3D standalone, and the negligible-coupling sweep with
  `rPvj=1e6`).
- Active-coupling MMS runs use the production PVJ operator directly: the coupled
  verifier does not overwrite terminal currents, 3D source fields, or implicit
  source coefficients.
- The coupled verifier adds the exact PVJ source to the manufactured ionic
  residual, so the analytical 1D/3D fields solve the coupled MMS equations while
  the numerical PVJ source remains raw.
- The coupled verifier also reports the diagnostic residual
  `S_pvj(V_num) - S_pvj(V_exact)`.
- The final N=10→20→40→80 baseline scheme matrix confirms **O(h²)** convergence
  for unidirectional/bidirectional and explicit/implicit PVJ assembly.
- The production PVJ coupler does not contain FDA-specific manufactured-solution
  formulas. The FDA reference is used inside `coupled1D3DMonodomainVerifier` to
  form the MMS source and diagnostics.
- The default graph uses the non-cancelling `y=1/6, z=1/3` terminal placement.
  Boundary fluxes remain zero when terminals lie on x-boundary faces, but the
  PVJ terms are no longer hidden by a cancelling placement.

A first-step bug (V_1D uninitialised at `t=0`, fixed via `preInitialize()`) found
during this work is documented in the notes file.
