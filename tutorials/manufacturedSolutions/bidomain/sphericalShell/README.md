# bidomain/sphericalShell

Heart-only bidomain manufactured solution on a spherical heart shell, meshed
as a cubed sphere (`-mesh hex`, default) or with gmsh tetrahedra
(`-mesh tet`): sealed curved boundaries approximated by planar faces, two
unequal rotated anisotropic conductivity tensors with all six components that
turn through the wall, and potentials whose conormal fluxes vanish on both
spheres.

Origin: this work.

## Problem

    chi Cm dVm/dt = div(Gi grad(Vm + phiE)) - chi Cm Iion + I_s     in r1 < |x| < r2
    div((Gi + Ge) grad phiE) = -div(Gi grad Vm) + s_e                in r1 < |x| < r2
    n . Gi grad(phiI) = 0,  n . Ge grad(phiE) = 0                   on |x| = r1, |x| = r2

with `phiI = Vm + phiE`. The shell has `r1 = 0.5` and `r2 = 1`; the inner
sphere is the patch `endocardium` and the outer sphere the patch
`epicardium`. The exact potentials are `phiI = sqrt(1 + t) Fi` and
`phiE = sqrt(1 + t) Fe`, with `n . Gi grad(Fi) = 0` and
`n . Ge grad(Fe) = 0` on each sphere; `Fi` is the field of
[monodomainPseudoECG/sphericalShell](../../monodomainPseudoECG/sphericalShell/README.md)
and `Gi`, `Ge` are the heart tensors of
[bathBidomain/sphericalShell](../../bathBidomain/sphericalShell/README.md),
which this case shares its geometry, mesh generators and `N` convention with.
The ionic model closes the `Vm` equation with `I_s = sqrt(1 + t) (beta F -
div(Gi grad Fi))`, `F = Fi - Fe`; the `phiE` equation takes
`s_e = sqrt(1 + t) (div(Gi grad Fi) + div(Ge grad Fe))`, whose integral over the
shell is zero. The generator is
[`shell_bidomain_mms.py`](../../../../src/verificationModels/bathBidomainVerification/generator/README.md).

### Stack

- myocardium solver: `bidomainSolver`, `solutionAlgorithm implicit`,
  `timeCouplingScheme sbdf2`
- ionic model: `monodomainFDAManufactured` (`dimension "3D"`)
- verifier: `manufacturedShellBidomainVerifier`: adds `s_e` to the `phiE`
  right-hand side and requests the compatibility projection of the
  pure-Neumann `phiE` system
- conductivity: `conductivitySource field`; `0/ConductivityIntracellular` and
  `0/ConductivityExtracellular` written by `setManufacturedFields`
- walls: `sealedHeartBoundary true`, `sealedWallTrace zeroGradient`
- reference: `phiERefPoint` inside the shell pins `phiE`; the error rows
  compare `phiE` with the exact field after removing the volume-weighted mean
  error (`phiE_gauge`)
- gradient: `sealedWallLeastSquares` for `grad(Vm)`, `grad(Vm.component(0))`,
  `grad(phiE)` and `grad(phiE.component(0))` (`system/fvSchemes`): wall cells
  take no value from the sealed walls; `leastSquares` elsewhere
- `phiE` operator: the tensor Laplacian with `conormalStabilisation`
  (`numericsModels`): the face-conormal magnitude `|G n|` is its implicit
  coefficient, so the correctors converge where the tensors are oblique to all
  three face directions of the cubed sphere; the term vanishes at convergence

`log.cardiacFoam` prints `phiE compatibility projection: imbalance <value>`
at every `phiE` solve: the sum of the equation's source and boundary
coefficients that the projection removes.

## Mesh

- `-mesh hex` (`system/blockMeshDict`): six projected blocks of a cubed
  sphere. `N`, a multiple of 4, sets `k = N/4`: `3k` cells along a cube-face
  edge and `2k` radially, so the shell has `N` radial cells per unit length.
- `-mesh tet` (`shell.geo.template`, copied to `local/shell.geo`): Delaunay
  tetrahedra of size `1/N`, a fixed random seed and one thread, so a gmsh
  version gives the same mesh for the same `N` and seed.

The mesh has exactly two patches, `endocardium` and `epicardium`, and one
cellZone, `myocardium`.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

Every option writes a file in `local/`, which the dictionaries include;
`Allclean` removes `local/` and the mesh. `Allrun` stops with
`run ./Allclean first` if a mesh log exists, and, for `-mesh tet`, when `gmsh`
or `gmshToFoam` is missing or gmsh reports an error.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` (4 subdomains, scotch), parallel `cardiacFoam`, `reconstructPar -latestTime` | serial |
| `-mesh hex\|tet` | cubed sphere or gmsh tetrahedra | hex |
| `-N n` | hex: a multiple of 4, radial cells per unit length; tet: element size `1/n` | 12 |
| `-seed n` | gmsh random seed (`-mesh tet`) | 1 |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.2 |

`postProcessing/3D_<hex|tet>_<N>_cells.dat` holds the error summary
(`field L1 L2 Linf`), with the `Vm`, `u1`, `u2`, `phiE_gauge`, bulk, boundary
and gradient rows of
[Bulk, boundary and gradient rows](../../README.md#bulk-boundary-and-gradient-rows).

Regression (hex cubed sphere, `-N 12`, short end time, parallel on four
ranks; no gmsh needed):

```bash
./regression/regressionTest.sh
```

## Studies

`setup/studies/` holds one study per mesh family, each run in a copy of the
case: `hexConvergence` (cubed sphere) and `tetConvergence` (gmsh
tetrahedra).
