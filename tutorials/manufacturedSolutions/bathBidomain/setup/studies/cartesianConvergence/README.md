# cartesianConvergence - bathBidomain

## Purpose

This study validates spatial convergence on structured hexahedral (cartesian) meshes using the bathBidomain exact solution.

## Execution

Execute the sweep using `foamctl run --sweep sweep_hex_convergence.json` or by invoking the local runner script.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
