# tetConvergence - bidomain/insulatedWall

Spatial convergence of the insulated-wall exact solution on Delaunay
tetrahedra that are periodic in `y` and `z`.

| file | role |
|---|---|
| `box.geo.template` | gmsh unit cube, `lc = 1/N`, periodic `y` and `z` surfaces, `walls` on `x = 0, 1` |
| `createPatchDict` | turns `y0/y1` and `z0/z1` into the cyclic pairs `cy0/cy1`, `cz0/cz1` |

## Ladder

| N | lc | deltaT |
|---|---|---|
| 10 | 0.1 | 0.00892857 |
| 20 | 0.05 | 0.00224215 |
| 40 | 0.025 | 0.000560538 |

`endTime 0.2` for every level.

## Execution

In a copy of the case, for each `N`:

```bash
N=20
DT=0.00224215
cp -r tutorials/manufacturedSolutions/bidomain/insulatedWall run_N${N}
cd run_N${N}
sed "s/__LC__/$(python3 -c "print(1/${N})")/" setup/studies/tetConvergence/box.geo.template > box.geo
gmsh -3 box.geo -o box.msh
gmshToFoam box.msh
cp setup/studies/tetConvergence/createPatchDict system/
createPatch -overwrite
checkMesh
foamDictionary -precision 17 -entry PIMPLE/nNonOrthogonalCorrectors -add 1 system/fvSolution
foamDictionary -precision 17 -entry deltaT -set ${DT} system/controlDict
./Allrun parallel
```

- The error summary is `postProcessing/*_cells.dat` (`Vm`, `phiE_gauge`, `phiI_gauge` rows: L1 L2 Linf).
- Observed order: `p = log2(L2(N)/L2(2N))`.

## Outputs

Run copies and meshes stay outside the repository.
