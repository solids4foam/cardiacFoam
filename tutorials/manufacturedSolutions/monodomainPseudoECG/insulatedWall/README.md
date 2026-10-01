# monodomainPseudoECG/insulatedWall

Monodomain manufactured solution whose walls are insulated in the conormal
sense, `n.G.grad(Vm) = 0`, while `n.grad(Vm) != 0` and the tangential gradient
is non-zero. The parent case uses `sin^2` profiles whose gradient vanishes on
every wall, so it cannot distinguish `n.grad(Vm) = 0` from `n.G.grad(Vm) = 0`.
This case can.

Origin: this work.

## Stack

- myocardium solver: `monodomainSolver`
- ionic model: `monodomainFDAManufactured`
- field verifier: `manufacturedConormalMonodomainVerifier`
- analytical profile: `verificationModels/conormalManufacturedProfile.H`

## Problem

- Domain: unit cube. `walls` = the faces `x = 0` and `x = 1`. `y0/y1` and
  `z0/z1` are cyclic pairs.
- Conductivity: the rotated tensor of the parent case.
- Exact solution: `Vm = sqrt(1 + t) Re[v(x) exp(i 2 pi (y + z))]`. The quadratic `v`
  satisfies the conormal condition at `x = 0` and `x = 1`.

## Wall treatment

Set in `constant/electroProperties` under `monodomainSolverCoeffs`.

| label | keys | wall |
|---|---|---|
| 0 | `sealedHeartBoundary false;` | `zeroGradient`, tensor face flux keeps `mag(S) (G n)_t . grad_t(Vm)` |
| A | `sealedHeartBoundary true;` | zero wall face conductivity |
| A+B | `sealedHeartBoundary true; sealedWallTrace conormal;` | plus `conormalZeroFlux` wall value, `dVm/dn = -(G n)_t . grad_t(Vm) / (n . G n)` |

The case ships A+B.

## Usage

```bash
blockMesh -dict system/blockMeshDict.3D
./Allrun parallel
./regression/regressionTest.sh
```

Tetrahedral meshes: [setup/studies/tetConvergence](setup/studies/tetConvergence/README.md).
