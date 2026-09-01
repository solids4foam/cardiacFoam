# tetConvergence - monodomainPseudoECG

## Purpose

This study validates spatial convergence on unstructured tetrahedral meshes using the monodomainPseudoECG exact solution.

It is the retained generic-Delaunay family: four resolutions for each of the
axis-aligned/rotated conductivity tensors and Gauss--linear/least-squares
gradient reconstructions.  The pseudo-ECG verifier is enabled at every level
with the same reference quadrature as the Cartesian sweep.  No Frontal/Netgen
duplicate is scheduled for this rerun.

## Execution

From the repository root, use the current wrapper:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetConvergence/sweep_tet_generic.json \
    --output-dir .tmp/driverfoam/monodomainPseudoECG-tet
applications/scripts/driverFoam/bin/driverFoam sweep-run \
    --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetConvergence/sweep_tet_generic.json \
    --output-dir .tmp/driverfoam/monodomainPseudoECG-tet
```

Resolve the runtime preflight before expecting OpenFOAM execution.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
