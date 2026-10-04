# tetConvergence - eikonalECG/insulatedWall

Spatial convergence of the insulated-wall eikonal solution on Delaunay
tetrahedra.

| file | role |
|---|---|
| `box.geo.template` | gmsh unit cube, `lc = 1/N`, `walls` on `x = 0, 1`, `sides` elsewhere |

## Execution

In a copy of the case, for each `N` in 10, 20, 40:

```bash
N=20
cp -r tutorials/manufacturedSolutions/eikonalECG/insulatedWall run_N${N}
cd run_N${N}
sed "s/__LC__/$(python3 -c "print(1/${N})")/" setup/studies/tetConvergence/box.geo.template > box.geo
gmsh -3 box.geo -o box.msh
gmshToFoam box.msh
checkMesh
decomposePar
mpirun -np 6 cardiacFoam -parallel
reconstructPar -latestTime
```

- `gmshToFoam` replaces `blockMesh`, so do not run `Allrun`, which runs `blockMesh`.
- The activation-time summary is `postProcessing/*3D_*_cells_*.dat`.
- Observed order: `p = log2(L2(N)/L2(2N))`.

## Outputs

Run copies and meshes stay outside the repository.
