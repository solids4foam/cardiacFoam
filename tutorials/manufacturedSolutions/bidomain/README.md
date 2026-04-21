# manufacturedSolutions/bidomain tutorial

This tutorial is the one-mesh manufactured-solution bidomain check.

This is the manufactured-solution verification workflow for the one-mesh
bidomain stack.

## Stack

- electro model: `bidomainSolver`
- ionic model: `bidomainFDAManufactured`
- field verification:
  `manufacturedFDABidomainVerifier` from `libverificationModels`
- shared analytical oracle:
  `verificationModels`

## Purpose

This case verifies:

- field convergence against an analytical oracle
- manufactured ionic export variables (`u1`, `u2`, `u3`)
- bidomain potentials (`phiE`, `phiI`) on a single mesh

## Key configuration idea

For this workflow, the ionic model exposes manufactured verification metadata,
and the exact manufactured reference no longer lives inside the ionic-model
folder.

That split is now:

- ionic-model-side behavior: `bidomainFDAManufactured`
- analytical oracle: `verificationModels`
- generic field verification hook: `modelPrePostProcessors`
- concrete manufactured field verifier:
  `verificationModels/bidomainVerification`

## Outputs

Typical outputs include:

- manufactured field summaries in `postProcessing/`
- `Vm`
- `phiE`
- `phiI`

## Execution

Manual:

```bash
blockMesh -dict system/blockMeshDict.1D
./Allrun
./runRegressionTest.sh
```

Driver-managed sweeps:

```bash
applications/scripts/driverFoam/bin/driverFoam all --tutorial manufacturedFDABidomain --config tutorials/manufacturedSolutions/bidomain/setupManufacturedFDA/driver_config.json
```

The driver is the intended entrypoint for large dimension / cell-count / time-step
sweeps.
