# Bidomain with bath manufactured solutions

Verification of the bidomain solver coupled to a conductive bath: the global
heart-and-bath `phiE` equation and the heart-bath interface. The case README
describes the problem, the options and the checks.

| Case | Geometry | What it verifies |
| --- | --- | --- |
| [`slab`](slab/README.md) | myocardium slab between two bath slabs with planar interfaces, on hexahedra (1D, 2D, 3D) or conformal gmsh tetrahedra (3D) | the heart-and-bath `phiE` equation, current continuity across the interfaces, the bath boundary conditions, the ionic coupling and the predictor-corrector bath coupling |

The heart-only bidomain is in [`../bidomain`](../bidomain/README.md).

## Run a case

```bash
cd tutorials/manufacturedSolutions/bathBidomain/slab
./Allrun -h        # options of this case
./Allrun           # build the mesh, run, post-process
```

The regression, from the repository root:

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight \
    ./tutorials/Alltest-regression --only 'bathBidomain/slab'
```

`-mesh tet` needs gmsh and `gmshToFoam`.
