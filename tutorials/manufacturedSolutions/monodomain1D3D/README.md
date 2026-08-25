# manufacturedSolutions/monodomain1D3D tutorial

Manufactured 1D-3D monodomain coupling test case. Couples a 3D monodomain myocardium domain to a small 1D Purkinje graph through `reactionDiffusionPvjCoupler` and verifies coupled fields with `coupled1D3DMonodomainVerifier`.

## Overview

Default graph lies on `y=1/6, z=1/3`, with PVJ terminals at `(0, 1/6, 1/3)` and `(1, 1/6, 1/3)`. Terminal faces satisfy homogeneous Neumann boundary condition for 3D manufactured solution because terminals are on `x=0` and `x=1` where `d cos(π x)/dx = 0`. Terminal values do not cancel: `F_3D = -0.5 F_1D` at PVJs (unlike older `y=0.5, z=1/3` placement).

## Graph Configuration

### Generating Refined Graphs

```bash
setup/generate_purkinje_graphs.py
```

Creates:

- `constant/purkinjeGraph.nodes003` through `.nodes161`
- `constant/graphFiles/purkinjeGraph.nodes*.vtk`

### Selecting Active Graph

```bash
setup/select_purkinje_graph.sh nodes041
```

The active 1D graph input is `constant/purkinjeGraph`.

## Usage

### Manual Execution

```bash
./Allrun
```

`Allrun` regenerates `constant/polyMesh` from `system/blockMeshDict.3D` before launching `cardiacFoam`.

### Graph-Only Diagnostic

```bash
blockMesh -dict system/blockMeshDict.3D
runPurkinjeGraph -case .
```

### Graph-Only Convergence Rates

Example of running a graph convergence study via driverFOAM (selects each `constant/purkinjeGraph.nodes*` input, runs `runPurkinjeGraph`):

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomain1D3D/setup/studies/coupledConvergence/sweep_active.json
```

Writes:
- `outputs/1dGraphConvergence/graph_convergence_summary.csv`
- `outputs/1dGraphConvergence/graph_convergence_rates.csv`

### Coupled 1D-3D Convergence Sweeps (Suggested)

Active coupling (bidirectional PVJ):

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomain1D3D/setup/studies/coupledConvergence/sweep_active.json
```

Bidirectional coupling:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomain1D3D/setup/studies/coupledConvergence/sweep_bidirectional.json
```

Decoupled (negligible PVJ coupling with `rPvj=1e6`):

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomain1D3D/setup/studies/coupledConvergence/sweep_decoupled.json
```

Coupled sweeps run `cardiacFoam` under joint 1D/3D refinement, copy verifier summaries, and write:

- `outputs/coupled1D3DConvergence/coupled_convergence_summary.csv`
- `outputs/coupled1D3DConvergence/coupled_convergence_rates.csv`
- `outputs/coupled1D3DConvergence/coupled_1D3D_convergence.png`
- `outputs/coupled1D3DConvergence/coupled_1D3D_convergence.pdf`

## Convergence verification

**Why this study exists.** The coupled 1D-3D monodomain implementation has to be shown to converge at the expected rate
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
