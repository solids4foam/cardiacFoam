# manufacturedSolutions/monodomainPseudoECG tutorial

This is the manufactured-solution verification workflow for the monodomain
stack with pseudo-ECG verification.

## Stack

- myocardium solver: `monodomainSolver`
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

- ionic-model-side behavior: `monodomainFDAManufactured`
- analytical oracle: `verificationModels`
- field verification hook: `modelPrePostProcessors`
- field verifier: `verificationModels/monodomainVerification`
- ECG verifier: `verificationModels/ecgVerification`

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
./regressionTest.sh
```

Driver-managed sweeps:

```bash
applications/scripts/driverFoam/bin/driverFoam all --entry manufacturedFDA --config tutorials/manufacturedSolutions/monodomainPseudoECG/setup/driver_config.json
```
