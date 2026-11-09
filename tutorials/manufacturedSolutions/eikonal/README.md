# Eikonal manufactured solutions

Verification of the eikonal solver against manufactured and exact solutions,
one case per geometry. The manufactured solutions are defined in
[`../generators/`](../generators/README.md); each case README describes the
problem, the options and the checks.

| Case | Geometry | What it verifies |
| --- | --- | --- |
| [`cube`](cube/README.md) | unit cube on hexahedra (1D, 2D, 3D) or gmsh tetrahedra (3D) | the interior discretisation of the eikonal-diffusion equation and the eikonal-driven ECG against its reference integral |
| [`sphericalShell`](sphericalShell/README.md) | spherical heart shell on a cubed sphere or gmsh tetrahedra | sealed curved walls crossed by a rotated full conductivity tensor, with wavefronts that collide on the walls; with `-ecg`, the eikonal ECG of that field |
| [`coupled1D3DSphere`](coupled1D3DSphere/README.md) | one-dimensional graph into a sealed tissue sphere, gmsh tetrahedra | the path from a graph into tissue and the curved front, against an exact solution with no manufactured source |

## Run a case

```bash
cd tutorials/manufacturedSolutions/eikonal/cube
./Allrun -h        # options of this case
./Allrun           # build the mesh, run, post-process
```

The regression of one case, from the repository root:

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight \
    ./tutorials/Alltest-regression --only 'eikonal/cube'
```

`coupled1D3DSphere` needs gmsh and `gmshToFoam`; `cube` and `sphericalShell`
need them only with `-mesh tet`.
