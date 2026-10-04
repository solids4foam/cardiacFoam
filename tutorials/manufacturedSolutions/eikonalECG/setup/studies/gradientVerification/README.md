# gradientVerification - eikonalECG

## Status (2026-09-17)

Configured (spec, mesh template and aggregator in place) but not rerun in
this verification pass. It meshes with `../tetConvergence/box.geo.template`.
Its one existing case, `leastSquares_axis_40_errorField`, supplies the
`@tbl-eikonal-localisation` data.

## Purpose

This study exercises the gradient reconstruction operator in isolation
(`gradientReconstructionOrder`, `applications/test/gradientReconstructionOrder/`)
against an exact analytic field on the tet mesh, comparing `gaussLinear` vs
`leastSquares` across `N = 10, 20, 40, 80`.

Each case runs the full mesh -> solve pipeline, with
`gradientReconstructionOrder` appended as a workflow_dag step after the solve
(`gradient_reconstruction=True`). Its stdout is captured under
`postProcessing/workflow_logs/` and archived by the sweep's snapshot/diff
collector.

`analyse_error_localisation.py` (the single-case, retained-time-directory
spatial-correlation deep dive over the SOLVED field's cellwise error -- not
the gradient-operator error above) is a read-only analysis script, run
manually after a case with `error_localisation_analysis=True` and
`writeErrorField` enabled has completed; see the script's docstring for
direct invocation.

## Execution

Spec: `tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradientVerification/sweep_gradient_tet.json`.

```bash
[omnidriver command to run]
python3 tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradientVerification/aggregate_gradient_verification.py <sweep output dir>
```

This writes `eikonal_gradient_tet.csv` into the sweep output directory (full
gaussLinear/leastSquares comparison; the registered experiment's leastSquares-only
table is `../gradient_reconstruction/`'s).

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, which is explicitly ignored by git. Do not commit
generated OpenFOAM data.
