# errorLocalisation - eikonalECG

## Purpose

This is the canonical omnidriver path for eikonalECG's solved-field
bulk/boundary error decomposition across the tet `N` ladder
(`N = 10, 20, 40, 80`) on both gradient schemes (`GaussLinear`,
`leastSquares`). It sweeps the `manufacturedEikonalECG` entry with
`eikonalSolverCoeffs.verificationModel.writeErrorField` enabled to get the
`L2_bulk`/`L2_boundary`/`L2_total` trend the paper reports, using the same
bulk/boundary split convention (`manufacturedEikonalVerifier.C`'s
`computeBoundaryBulkNorms`) as the standalone `gradientReconstructionOrder`
utility's own decomposition (see `../gradientVerification/`).

The workflow also writes the cellwise `activationTimeError` field, mesh-quality
fields (`checkMesh -writeAllFields`), and cell centres after the solve. These
are required inputs to the retained spatial-localisation analysis.

## Execution

First materialize and inspect the eight omnidriver cases. Spec: `tutorials/manufacturedSolutions/eikonalECG/setup/studies/errorLocalisation/sweep_tet_error_localisation.json`.

    [omnidriver command to run]

Then run the same manifest:

    [omnidriver command to run]

After a successful sweep, run the in-repository bulk/boundary aggregation over
the archived outputs:

    python3 tutorials/manufacturedSolutions/eikonalECG/setup/studies/errorLocalisation/aggregate_bulk_boundary.py <sweep output dir>

For a selected completed case, the retained coordinate-based localisation
analysis can then read the driver-generated `activationTimeError`, `Cx`, `Cy`,
and `Cz` fields:

    python3 tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradientVerification/analyse_error_localisation.py

This writes `eikonal_bulk_boundary_tet.csv` into the sweep output directory.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the
local `results/` folder, which is explicitly ignored by git. Do not commit
generated OpenFOAM data.
