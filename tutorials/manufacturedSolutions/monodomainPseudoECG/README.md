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

After the sweep completes, persist the canonical Paper I convergence tables (reads
`driverPostProcessingArchive_postProcessing/`, writes `setup/results/*.csv`, never
touches the sweep's own output):

```bash
python3 applications/scripts/paperI_results/aggregate.py mono_spatial
python3 applications/scripts/paperI_results/aggregate.py pseudo_ecg_spatial
```

Temporal-discretization sweep:

```bash
applications/scripts/driverFoam/bin/driverFoam all --entry manufacturedFDA --config tutorials/manufacturedSolutions/monodomainPseudoECG/setup/driver_config_temporal_3d.json
```

Recommended temporal MMS settings:

- spatial baseline: `N = 80`, so `dx = 1/80 = 0.0125`
- temporal sweep: `dt = [2.24215e-3, 1.121075e-3, 5.60538e-4, 2.80269e-4]`

These values keep the existing spatial sweep unchanged while providing a fixed-`dx`
3D run for temporal-order measurements. The temporal config writes to
`postProcessingTemporal/` and disables ECG post-processing so the field
convergence results are isolated from the pseudo-ECG workflow.
