# Bidomain manufactured solutions

Verification of the heart-only bidomain solver (the elliptic `phiE` solve
coupled to the parabolic `Vm` solve) against manufactured solutions, one case
per geometry. Each case README describes the problem, the options and the
checks.

| Case | Geometry | What it verifies |
| --- | --- | --- |
| [`cube`](cube/README.md) | unit cube on hexahedra (1D, 2D, 3D) or gmsh tetrahedra (3D) | the interior operators of the one-mesh bidomain stack with diagonal conductivities and a zero normal gradient on every face |
| [`sphericalShell`](sphericalShell/README.md) | spherical heart shell on a cubed sphere or gmsh tetrahedra | sealed curved walls crossed by two unequal rotated conductivity tensors, and the pure-Neumann `phiE` system with its compatibility projection |

The bidomain with a conductive bath is in
[`../bathBidomain`](../bathBidomain/README.md).

## Run a case

```bash
cd tutorials/manufacturedSolutions/bidomain/cube
./Allrun -h        # options of this case
./Allrun           # build the mesh, run, post-process
```

The regression of one case, from the repository root:

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight \
    ./tutorials/Alltest-regression --only 'bidomain/cube'
```

`-mesh tet` needs gmsh and `gmshToFoam`.
