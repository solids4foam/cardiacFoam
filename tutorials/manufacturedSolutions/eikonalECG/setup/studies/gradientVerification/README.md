# gradientVerification - eikonalECG

## Purpose

This study exercises the gradient reconstruction operator in isolation
(`gradientReconstructionOrder`, `applications/test/gradientReconstructionOrder/`)
against an exact analytic field on the tet mesh, comparing `gaussLinear` vs
`leastSquares` across `N = 10, 20, 40, 80` -- the full comparison used in the
paper's qualitative discussion. For the registered `eikonal_gradient_tet`
Paper I table (`leastSquares` only), see `../gradient_reconstruction/`.

`gradientReconstructionOrder` needs only the mesh and `system/fvSchemes` --
no cardiacFoam solve is strictly required for it -- but the driverFOAM
tutorial (`manufactured_eikonal_ecg.py`, `gradient_reconstruction=True`)
still runs the full mesh -> solve pipeline and appends it as a workflow_dag
step after the solve, so every case in this study is a complete, auditable
driverFOAM run rather than a bespoke bash loop. Its stdout is captured under
`postProcessing/workflow_logs/` like any workflow step and archived by the
sweep's own generic snapshot/diff collector.

`analyse_error_localisation.py` (the single-case, retained-time-directory
spatial-correlation deep dive over the SOLVED field's cellwise error -- not
the gradient-operator error above) is a read-only analysis script, run
manually after a case with `error_localisation_analysis=True` and
`writeErrorField` enabled has completed; see `../errorLocalisation/README.md`
for the equivalent driverFOAM-driven N-ladder, and this script's own
docstring for direct invocation.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradientVerification/sweep_gradient_tet.json
python3 tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradientVerification/aggregate_gradient_verification.py
```

This writes `results/eikonal_gradient_tet.csv` (full gaussLinear/leastSquares
comparison; distinct from the canonical `setup/results/eikonal_gradient_tet.csv`
the registered experiment writes).

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, which is explicitly ignored by git. Do not commit
generated OpenFOAM data.
