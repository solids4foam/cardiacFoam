# monodomainPseudoECG/sphericalShell

Monodomain manufactured solution on a spherical heart shell, meshed as a
cubed sphere (`-mesh hex`, default) or with gmsh tetrahedra (`-mesh tet`):
sealed curved boundaries approximated by planar faces, a rotated anisotropic
conductivity with all six tensor components that turns through the wall, and
a solution whose conormal flux vanishes on both spheres.

Origin: this work.

## Problem

    chi Cm dVm/dt = div(G grad Vm) - chi Cm Iion + I_s    in r1 < |x| < r2
    n . G grad Vm = 0                                       on |x| = r1, |x| = r2

The shell has `r1 = 0.5` and `r2 = 1`; the inner sphere is the patch
`endocardium` and the outer sphere the patch `epicardium`. The case is the
heart-only counterpart of
[bathBidomain/sphericalShell](../../bathBidomain/sphericalShell/README.md) and
shares its geometry, its mesh generators and its `N` convention.

### Stack

- myocardium solver: `monodomainSolver`, `solutionAlgorithm implicit`,
  `timeCouplingScheme godunov`
- ionic model: `monodomainFDAManufactured` (`dimension "3D"`)
- verifier: `manufacturedShellMonodomainVerifier`
- conductivity: `conductivitySource field`; `0/Conductivity` written by
  `setManufacturedFields`
- walls: `sealedHeartBoundary true`, `sealedWallTrace zeroGradient` or
  `conormal` (`-trace`)
- gradient: `sealedWallLeastSquares` for `grad(Vm)` and `grad(Vm.component(0))`
  (`system/fvSchemes`): wall cells take no value from the sealed walls;
  `leastSquares` elsewhere

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
| `-trace zeroGradient\|conormal` | `sealedWallTrace` | zeroGradient |

`postProcessing/3D_<hex|tet>_<N>_cells.dat` holds the error summary
(`field L1 L2 Linf`), with the `Vm`, `u1`, `u2`, bulk, boundary and gradient
rows of
[Bulk, boundary and gradient rows](../../README.md#bulk-boundary-and-gradient-rows).

Regression (hex cubed sphere, `-N 12`, short end time, parallel on four
ranks; no gmsh needed):

```bash
./regression/regressionTest.sh
```

## Studies

`setup/studies/` holds one study per mesh family, each run in a copy of the
case: `hexConvergence` (cubed sphere) and `tetConvergence` (gmsh
tetrahedra). Each README gives the commands, the ladder and the measured
behaviour of the shipped gradient and of the stock `leastSquares` control.
