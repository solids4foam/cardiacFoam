# eikonalHeart

Standalone 1D-3D eikonal activation tutorial using the heart mesh, Purkinje
graph, PVJ locations, transmural field, and ECG electrodes prepared from the
PATHOS AcuteIschemia case.

## Stack

- myocardium solver: `eikonalSolver`
- conduction system solver: `eikonalSolver1D`
- coupling: `eikonalPvjCoupler`
- ECG: `eikonalECG`

The Purkinje graph activates from `rootStimulus.node 0`; terminal activation
times are deposited into the 3D heart mesh as constraints for `psi`. The 3D
solve uses the same PIMPLE style as
`tutorials/manufacturedSolutions/eikonalECG`: `nOuterCorrectors 5000`,
strict `psi` residual control, and a configured `asymmetric_psi` solver for the
advection-diffusion eikonal form. The heart mesh uses `PCG` with no
preconditioner for `psi`; the manufactured case's diagonal/DILU preconditioners
trap on this assembled heart-mesh matrix.

## Execution

```bash
./Allrun
./Allrun parallel
```

Expected primary outputs:

- `1/psi`
- `postProcessing/purkinjeNetwork.dat`
- `postProcessing/eikonalECG.dat`
