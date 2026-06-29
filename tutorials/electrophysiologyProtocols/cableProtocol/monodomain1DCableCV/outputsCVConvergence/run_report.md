# driverFoam Run Report

## Summary

- Entry: `monodomainAndEikonal1DCableCVConvergence`
- Entry kind: `registered_tutorial`
- Entry path: `coreProtocols/cableProtocol/monodomain1DCableCV`
- Requested action: `all`
- Run ID: `1782674663-5afa9cba`
- Status: `failed`
- Post-process status: `not_started`
- Dry run: `False`
- Continue on error: `False`
- Case root: `/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/coreProtocols/cableProtocol/monodomain1DCableCV`
- Setup root: `/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/coreProtocols/cableProtocol/monodomain1DCableCV/setupMonodomain1DCableCV`
- Output dir: `/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/coreProtocols/cableProtocol/monodomain1DCableCV/outputsCVConvergence`
- Started at (UTC): `2026-06-28T19:24:23.539198+00:00`
- Updated at (UTC): `2026-06-28T19:24:23.602165+00:00`
- Finished at (UTC): `None`
- Current case: `None`
- Total cases: `5`
- Planned cases: `0`
- Completed cases: `0`
- Failed cases: `1`

## Run Error

```text
"Scope 'eikonalSolverCoeffs' not found"
```

## Case Results

### eikonalSolver_eikonal_none_DT1000_DX0.5_COND01
- Status: `failed`
- Index: `1` / `5`
- Duration (s): `0.006162`
- Started at (UTC): `2026-06-28T19:24:23.584246+00:00`
- Finished at (UTC): `2026-06-28T19:24:23.590415+00:00`
- Parameters: `{"conductivity": "[-1 -3 3 0 0 2 0] (0.1334 0 0 0.1334 0 0.1334)", "conductivity_id": 1, "dt_ms": 1000.0, "dx_mm": 0.5, "ionicModel": "eikonal", "solver": "eikonalSolver", "tissue": "none"}`
- Error:

```text
"Scope 'eikonalSolverCoeffs' not found"
```

## Artifacts

- Plots manifest: `/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/coreProtocols/cableProtocol/monodomain1DCableCV/outputsCVConvergence/plots.json`
