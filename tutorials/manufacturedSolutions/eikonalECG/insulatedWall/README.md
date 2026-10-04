# eikonalECG/insulatedWall

Eikonal manufactured solution with insulated walls, `n.M.grad(psi) = 0`, on
the faces `x = 0` and `x = 1`, while `n.grad(psi) != 0` and the conductivity is
oblique to the wall. The parent case prescribes the exact activation time on
every face, so it never exercises the insulated wall.

Origin: this work.

## Stack

- myocardium solver: `eikonalSolver` (advection-diffusion form)
- field verifier: `manufacturedEikonalVerifier`
- ECG verifier: `manufacturedEikonalECGVerifier`

## Problem

- Domain: unit cube. `walls` = the faces `x = 0` and `x = 1`. `sides` = the
  remaining faces.
- Exact solution: `psi = exp(k . x)` with
  `k = (0.5, 0.78097651867925399, 2.0745795101815179)`, chosen so that
  `(M k) . e_x = 0`. `M` is the eikonal conductivity tensor.
- `verificationModel/productionPatches (walls)` leaves `walls` to the solver's
  wall treatment. The verifier writes the exact solution on every other face.
  The verifier stops if `n.(M k) != 0` on a listed patch.
- The ECG verifier uses the same `k`.

## Wall treatment

Set in `constant/electroProperties` under `eikonalSolverCoeffs`. `walls` is
`zeroGradient` in `0/activationTime`.

| label | keys | wall |
|---|---|---|
| exact | no `productionPatches`, `fixedValue` on `walls` | exact values on every face (parent case) |
| 0 | `sealedHeartBoundary false;` | `zeroGradient` |
| A | `sealedHeartBoundary true;` | zero wall face conductivity |
| A+B | `sealedHeartBoundary true; sealedWallTrace conormal;` | plus `conormalZeroFlux` wall value |

The case ships A+B.

## Usage

```bash
./Allrun parallel
./regression/regressionTest.sh
```

Tetrahedral meshes: [setup/studies/tetConvergence](setup/studies/tetConvergence/README.md).
