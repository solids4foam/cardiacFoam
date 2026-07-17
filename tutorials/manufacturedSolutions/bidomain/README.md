# manufacturedSolutions/bidomain tutorial

This tutorial is the one-mesh manufactured-solution bidomain check.

This is the manufactured-solution verification workflow for the one-mesh
bidomain stack.

## Stack

- myocardium solver: `bidomainSolver`
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

- ionic-model-side behavior: `bidomainFDAManufactured`
- analytical oracle: `verificationModels`
- field verification hook: `modelPrePostProcessors`
- field verifier: `verificationModels/bidomainVerification`

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
./regressionTest.sh
```

Driver-managed sweeps:

```bash
applications/scripts/driverFoam/bin/driverFoam all --entry manufacturedFDABidomain --config tutorials/manufacturedSolutions/bidomain/setup/driver_config.json
```

After the sweep completes, persist the canonical Paper I convergence table:

```bash
python3 applications/scripts/paperI_results/aggregate.py bidomain
```
