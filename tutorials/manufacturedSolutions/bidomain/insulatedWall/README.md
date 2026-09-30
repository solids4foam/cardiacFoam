# bidomain/insulatedWall

Heart-only bidomain manufactured solution with insulated walls,
`n.G_i.grad(phi_i) = 0` and `n.G_e.grad(phi_e) = 0`, and unequal anisotropy
(`G_e` is not proportional to `G_i`). On the walls `n.grad(Vm)`, `n.grad(phi_e)`
and the tangential gradients are non-zero, and the intracellular wall terms
`n.G_i.grad(Vm)` and `n.G_i.grad(phi_e)` cancel. The parent case uses
diagonal conductivities on axis-aligned walls, where `n.grad(u) = 0` and
`n.G.grad(u) = 0` coincide.

Origin: this work.

## Stack

- myocardium solver: `bidomainSolver`
- ionic model: `bidomainFDAManufactured`
- field verifier: `manufacturedUnequalConormalBidomainVerifier`
- analytical reference: `verificationModels/bidomainVerification/manufacturedUnequalConormalBidomainReference.H`

## Problem

- Domain: unit cube. `walls` = the faces `x = 0` and `x = 1`. `y0/y1` and
  `z0/z1` are cyclic pairs.
- `G_i = (0.1 0.02 0; 0.035 0; 0.025)`, `G_e = (0.05 -0.01 0; 0.035 0; 0.045)`
  (symmetric, upper triangle).
- The exact solution is closed form and needs no `phi_e` source.
- `phi_e` and `phi_i` are compared after removing the gauge at `phiERefPoint`.

## Wall treatment

Set in `constant/electroProperties` under `bidomainSolverCoeffs`.

| label | keys | wall |
|---|---|---|
| 0 | `sealedHeartBoundary false;` | `zeroGradient` for `Vm` and `phiE` |
| A | `sealedHeartBoundary true;` | zero wall face conductivity |
| A+B | `sealedHeartBoundary true; sealedWallTrace conormal;` | plus `conormalZeroFlux` wall values. `phiE` uses `G_e`. `Vm` uses `G_i` with offset `phiE`: `dVm/dn = dphi_i/dn - dphi_e/dn` |

The case ships A+B.

## Usage

```bash
blockMesh -dict system/blockMeshDict.3D
./Allrun parallel
./regression/regressionTest.sh
```

Tetrahedral meshes: [setup/studies/tetConvergence](setup/studies/tetConvergence/README.md).
