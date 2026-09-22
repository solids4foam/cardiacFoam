# cartesianConvergence - eikonalECG

## Purpose

This study validates spatial convergence on structured hexahedral (cartesian) meshes using the eikonalECG exact solution.

## Execution

First materialize and inspect the 12 Cartesian cases:

    driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/cartesianConvergence/sweep_hex_convergence.json

Then run the same manifest with driverFOAM:

    driverFoam sweep-run \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/cartesianConvergence/sweep_hex_convergence.json

Each case runs both the manufactured activation-time verifier and the
manufactured eikonal-ECG verifier. The ECG reference uses Gauss-Legendre order 48 with no check orders.
Against order 96, the worst order-48 difference over the closest eight
electrodes (standoff 0.051-0.082) is 9.1e-7.
The former aggregate.py command has been removed; use the archived driverFOAM
outputs as the run record until an in-repository canonical collector is
added.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
