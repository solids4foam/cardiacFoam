# gradient_reconstruction - eikonalECG

## Purpose

Canonical entry point for the registered `eikonal_gradient_tet` verification
experiment (`applications/scripts/driverFoam/verification_experiments.json`):
isolated `leastSquares` gradient reconstruction on the tet mesh across
`N = 10, 20, 40, 80`.

`sweep_gradient_tet.json` drives `gradientReconstructionOrder` as a
`gradient_reconstruction=True` workflow_dag step (see
`../gradientVerification/README.md` for how that step is wired) -- the
`leastSquares`-only subset of `../gradientVerification/`'s own
`sweep_gradient_tet.json`, which covers both gradient schemes for the
paper's qualitative gaussLinear-vs-leastSquares discussion.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradient_reconstruction/sweep_gradient_tet.json
python3 tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradient_reconstruction/aggregate_gradient_reconstruction.py
```

This writes the canonical `setup/results/eikonal_gradient_tet.csv`.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, which is explicitly ignored by git. Do not commit
generated OpenFOAM data.
