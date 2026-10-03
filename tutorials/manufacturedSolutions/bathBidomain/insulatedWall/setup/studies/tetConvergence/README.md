# tetConvergence - bathBidomain/insulatedWall

## Purpose

Spatial convergence of the conormal bath-bidomain exact solution on
unstructured tetrahedra. The domain is the tutorial slab
`[-1,2] x [0,1] x [0,0.25]`: heart `[0,1]` (cellZone `myocardium`), bath
`[-1,0]` and `[1,2]` (cellZone `bath`), periodic in `y` and `z`, insulated
outer walls `xMin` and `xMax`. The exact solution does not depend on `z`, so
the verifier keeps `dimension "2D"`.

| file | role |
|---|---|
| `slab.geo.template` | gmsh Delaunay tetrahedra, `lc = 1/N`, conforming heart/bath interfaces, periodic `y` and `z` surfaces |
| `createPatchDict` | turns `y0/y1` and `z0/z1` into the cyclic pairs `cy0/cy1`, `cz0/cz1` |

Cells with gmsh 4.15.2: 4053 (N = 10), 30223 (N = 20), about 2.4e5 (N = 40).

## Ladder

| N | lc | deltaT |
|---|---|---|
| 10 | 0.1 | 1e-3 |
| 20 | 0.05 | 2.5e-4 |
| 40 | 0.025 | 6.25e-5 |

`endTime 0.2` for every level.

Configurations, set in `constant/electroProperties` under `bidomainSolverCoeffs`:

| name | keys |
|---|---|
| 0 | `sealedHeartBoundary false; bathHeartPhiETrace zeroGradient;`, no `sealedWallTrace`, `interfaceConductivityInterpolation distanceWeightedHarmonic;` |
| A | `sealedHeartBoundary true; bathHeartPhiETrace zeroGradient;`, no `sealedWallTrace`, `interfaceConductivityInterpolation distanceWeightedHarmonic;` |
| A+T+B+C | tutorial defaults: `sealedHeartBoundary true; bathHeartPhiETrace global; sealedWallTrace conormal;`, `interfaceConductivityInterpolation conormalHarmonic;` |

## Execution

In a copy of the tutorial, for each `N` and configuration:

```bash
N=20
DT=2.5e-4
cp -r tutorials/manufacturedSolutions/bathBidomain/insulatedWall run_N${N}
cd run_N${N}
sed "s/__LC__/$(python3 -c "print(1/${N})")/" setup/studies/tetConvergence/slab.geo.template > slab.geo
gmsh -3 slab.geo -o slab.msh
gmshToFoam slab.msh
cp setup/studies/tetConvergence/createPatchDict system/
createPatch -overwrite
checkMesh
foamDictionary -precision 17 -entry deltaT -set ${DT} system/controlDict
setTorsoOrganConductivityField
decomposePar
mpirun -np 8 cardiacFoam -parallel
reconstructPar -latestTime
```

- `gmshToFoam` builds the `myocardium` and `bath` cellZones from the gmsh physical volumes, so `blockMesh` and `topoSet` are not run.
- The log reports `bathHeartPhiETrace global: <n> exposed heart faces` for the global trace (138 at N = 10).
- The error summary is `postProcessing/2D_<n>_cells.dat` (`Vm`, `phiE`, `phiI` rows: L1 L2 Linf). `phiE` has its volume mean removed.
- Observed order: `p = log2(L2(N)/L2(2N))`.

## Outputs

Run copies and meshes stay outside the repository.
