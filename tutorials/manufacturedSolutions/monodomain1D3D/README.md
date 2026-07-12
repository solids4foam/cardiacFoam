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
setup/run_purkinje_graph_sweep.sh
```

The sweep selects each `constant/purkinjeGraph.nodes*` input, runs
`runPurkinjeGraph`, copies the graph verifier summaries, and writes:

- `outputs/1dGraphConvergence/graph_convergence_summary.csv`
- `outputs/1dGraphConvergence/graph_convergence_rates.csv`

For coupled 1D-3D manufactured convergence rates, run:

```bash
setup/run_coupled_1D3D_sweep.sh
```

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

**Conclusions.**

- The 1D Purkinje and 3D myocardium solvers each converge at **O(h²)** on their own
  (1D standalone, 3D standalone, and the negligible-coupling sweep with
  `rPvj=1e6`).
- The original active-coupling MMS floor came from a missing verification forcing
  term: the standalone manufactured monodomain forcing did not cancel the PVJ
  coupling operator evaluated at the manufactured reference.
- `coupled1D3DMonodomainVerifier` now subtracts that manufactured coupling residual
  from the PVJ current/source buffers, using the actual staggered time levels
  (`t0/t0` for secondary preparation and `t0/t0+dt` for primary preparation).
- With active coupling at `rPvj=1`, `pvjRadius=0.11`, and boundary-face terminals,
  both unidirectional and bidirectional coupled sweeps recover **O(h²)** convergence.
  The verified N=10→20→40→80 3D Vm rates are approximately
  `1.84, 1.99, 1.99`; bidirectional 1D Vm rates are approximately
  `1.94, 2.00, 1.99`.
- The corrected MMS verifies that the solvers integrate the source produced by
  `pvjMapper` at the expected order. It does not independently verify the spatial
  accuracy of `pvjMapper`'s fixed-radius sphere-average geometry, because the
  verifier intentionally duplicates that same mapper to cancel the manufactured
  residual.

A first-step bug (V_1D uninitialised at `t=0`, fixed via `preInitialize()`) found
during this work is documented in the notes file.
