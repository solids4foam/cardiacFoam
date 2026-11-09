# Monodomain manufactured solutions

Verification of the monodomain solver against manufactured solutions, one case
per geometry. The solutions are defined in
[`../generators/`](../generators/README.md); each case README describes the
problem, the options and the checks.

| Case | Geometry | What it verifies |
| --- | --- | --- |
| [`cube`](cube/README.md) | unit cube on hexahedra (1D, 2D, 3D) or gmsh tetrahedra (3D) | the interior operators (diffusion, gradient reconstruction, non-orthogonal correction, ionic coupling, time integration) with the FDA cosine field and a rotated anisotropic field, and the 3D pseudo-ECG against its reference integral |
| [`sphericalShell`](sphericalShell/README.md) | spherical heart shell on a cubed sphere or gmsh tetrahedra | sealed curved walls crossed by a rotated full conductivity tensor; with `-ecg`, the pseudo-ECG of that field |
| [`coupled1D3DCube`](coupled1D3DCube/README.md) | Purkinje graph joined to the unit cube | the graph operator, the tissue operator and the Purkinje-myocardial junction source under joint refinement |
| [`coupled1D3DShell`](coupled1D3DShell/README.md) | Purkinje graph joined to the spherical shell on tetrahedra | the junction against its continuum definition, with supports cut by the sealed curved endocardium |

## Run a case

```bash
cd tutorials/manufacturedSolutions/monodomain/cube
./Allrun -h        # options of this case
./Allrun           # build the mesh, run, post-process
```

The regression of one case, from the repository root:

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight \
    ./tutorials/Alltest-regression --only 'monodomain/cube'
```

`coupled1D3DShell` needs gmsh and `gmshToFoam`; `sphericalShell` and `cube`
need them only with `-mesh tet`.
