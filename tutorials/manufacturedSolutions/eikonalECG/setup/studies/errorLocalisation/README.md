# errorLocalisation - eikonalECG

## Purpose

This is the canonical driverFOAM path for eikonalECG's solved-field
bulk/boundary error decomposition across the tet `N` ladder
(`N = 10, 20, 40, 80`) on both gradient schemes (`GaussLinear`,
`leastSquares`). It sweeps the `manufacturedEikonalECG` entry with
`eikonalSolverCoeffs.verificationModel.writeErrorField` enabled to get the
`L2_bulk`/`L2_boundary`/`L2_total` trend the paper reports, using the same
bulk/boundary split convention (`manufacturedEikonalVerifier.C`'s
`computeBoundaryBulkNorms`) as the standalone `gradientReconstructionOrder`
utility's own decomposition (see `../gradientVerification/`).

Complementary to `../gradientVerification/run_error_localisation.sh`, which
is a single-case (default `N=40`, `leastSquares`) spatial-correlation deep
dive kept for direct/manual use -- see that script's header.

## Execution

```bash
cd tutorials/manufacturedSolutions/eikonalECG
rm -rf setup/studies/errorLocalisation/results/sweepCases setup/studies/errorLocalisation/results/sweepRun
mkdir -p setup/studies/errorLocalisation/results
../../../applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec setup/studies/errorLocalisation/sweep_tet_error_localisation.json \
    --output-dir setup/studies/errorLocalisation/results/sweepRun
python3 setup/studies/errorLocalisation/aggregate_bulk_boundary.py
```

Current status: `sweep-plan` reaches case materialisation but fails before
OpenFOAM because the spec asks the driver to write
`eikonalSolverCoeffs.verificationModel.writeErrorField` and the current case
dictionary does not contain that key. This is a stale spec/dictionary
contract, not a solver runtime result; do not interpret the command above as
verified until the key is added or the override is moved to the current
dictionary scope.

This writes `setup/results/eikonal_bulk_boundary_tet.csv`.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, which is explicitly ignored by git. Do not commit
generated OpenFOAM data.
