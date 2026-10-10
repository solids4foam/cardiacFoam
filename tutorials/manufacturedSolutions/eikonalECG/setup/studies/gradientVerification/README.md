# gradientVerification - eikonalECG

## Status (2026-09-17)

Configured (spec, mesh template and aggregator in place) but not rerun in
this verification pass. It meshes with `../tetConvergence/box.geo.template`.
Its one existing case, `leastSquares_axis_40_errorField`, supplies the
`@tbl-eikonal-localisation` data.

## Purpose

This study exercises the gradient reconstruction operator in isolation
against an exact analytic field on the tet mesh, comparing `gaussLinear` vs
`leastSquares` across `N = 10, 20, 40, 80`. Every run's eikonal verifier
reports it as the `gradReconActivationTime` rows (exact cell and boundary-face
values through the configured scheme, bulk and boundary parts).

The spec also appends a `gradient_reconstruction=True` workflow_dag step after
the solve, which runs a `gradientReconstructionOrder` executable that this
repository does not build; its stdout would be captured under
`postProcessing/workflow_logs/`.

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
python3 tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradientVerification/aggregate_gradient_verification.py
```

This writes `results/eikonal_gradient_tet.csv` (full gaussLinear/leastSquares
comparison; distinct from the canonical `setup/results/eikonal_gradient_tet.csv`
the registered experiment writes).

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, which is explicitly ignored by git. Do not commit
generated OpenFOAM data.
