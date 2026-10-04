# cartesianConvergence - bathBidomain

## Purpose

This study validates spatial convergence on structured hexahedral (cartesian) meshes using the bathBidomain exact solution.

## Execution

From the repository root:

```bash
[omnidriver command to run]
```

Resolve the runtime preflight before expecting OpenFOAM execution.

`sweep_hex_convergence.json` and `sweep_hex_groundElectrode.json` run the
`groundElectrode` variant; `sweep_hex_electrodePair.json` runs the case's own
`electrodePair`. The first two are the same study.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
