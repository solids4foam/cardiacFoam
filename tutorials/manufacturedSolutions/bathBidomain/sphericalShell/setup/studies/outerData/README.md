# outerData - bathBidomain/sphericalShell

Exact current (Neumann, default) against exact potential (Dirichlet) on the
outer box, at fixed configuration.

## Execution

```bash
./Allrun parallel -mesh tet -N <n> -decomposition scotch -np <ranks> -elliptic -nOuter 10 -nNonOrth 10
./Allrun parallel -mesh tet -N <n> -decomposition scotch -np <ranks> -elliptic -nOuter 10 -nNonOrth 10 -outer dirichlet
```

## Expected behaviour

On tetrahedra `phiE_boundaryBand` converges at about 2.05 with the exact
current and at about 1.1 with the exact potential: a fixedValue face value is
placed at the normal foot of a skewed boundary cell, whereas a fixedGradient
face value is defined there. The interior potentials converge at about 1.9 in
both.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
