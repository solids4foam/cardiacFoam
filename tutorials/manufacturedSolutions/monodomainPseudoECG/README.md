# manufacturedSolutions/monodomainPseudoECG tutorial

This is the manufactured-solution verification workflow for the monodomain
stack with pseudo-ECG verification.

## Stack

- electro model: `MonoDomainSolver`
- ionic model: `monodomainFDAManufactured`
- field verification:
  `manufacturedFDAMonodomainVerifier` from `libverificationModels`
- shared analytical oracle:
  `verificationModels`
- optional manufactured pseudo-ECG verification:
  `pseudoECGManufacturedVerifier` from `libverificationModels`

## Purpose

This case verifies:

- field convergence against an analytical oracle
- manufactured ionic export variables (`u1`, `u2`, `u3`)
- manufactured pseudo-ECG reference output

## Key configuration idea

For this workflow, the ionic model exposes manufactured verification metadata,
and the exact manufactured reference no longer lives inside the ionic-model
folder.

That split is now:

- ionic-model-side behavior: `monodomainFDAManufactured`
- analytical oracle: `verificationModels`
- generic field verification hook: `modelPrePostProcessors`
- concrete manufactured field verifier:
  `verificationModels/monodomainVerification`
- manufactured ECG verification:
  `verificationModels/pseudoECGVerification`

## Outputs

Typical outputs include:

- manufactured field summaries in `postProcessing/`
- `postProcessing/pseudoECG.dat`
- `postProcessing/manufacturedPseudoECG.dat`
- `postProcessing/manufacturedPseudoECGSummary.dat`

## Execution

Manual:

```bash
./Allrun
./runRegressionTest.sh
```

Driver-managed sweeps:

```bash
applications/scripts/driverFoam/bin/driverFoam all --tutorial manufacturedFDA --config tutorials/manufacturedSolutions/monodomainPseudoECG/setupManufacturedFDA/driver_config.json
```

The driver is the intended entrypoint for large dimension / cell-count / time-step
sweeps.
