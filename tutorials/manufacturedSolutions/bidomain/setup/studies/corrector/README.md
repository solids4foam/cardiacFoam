# corrector - bidomain

## Purpose

This study validates the predictor-corrector inner loop convergence and stability for the bidomain solver.

## Execution

```bash
driverFoam sweep-run --spec tutorials/manufacturedSolutions/bidomain/setup/studies/corrector/sweep_corrector_study.json
```

Sweeps `N = 10, 20, 40` across the four reported variants (`baseline`,
`outer2`, `nonorth1`, `combined`), each a `(nOuterCorrectors,
nNonOrthogonalCorrectors)` pair applied via driverFOAM's own
`n_outer_correctors`/`n_nonorthogonal_correctors` overrides
(`manufactured_monodomain_pseudo_ecg._apply_case`); the short, fixed
step-count screening window (2/9/36 steps) is set via `control_dict_overrides`
on `writeControl`/`writeInterval`/`writeFormat`. Each of the 12 cases gets
its own output directory (named from its `(N, nOuterCorrectors,
nNonOrthogonalCorrectors)` values), and its raw `postProcessing/` output is
archived into that same directory's own `sweepCases/` subfolder via
`archive_dir_name` -- identical mechanism to every other sweep in this
tutorial, nothing corrector-specific.

Comparing the 12 results across variants/resolutions is left to the
postprocessing module, not a bespoke script.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
