# configurations - bathBidomain/sphericalShell

The interface configurations 0, A, AI, ATI and ATBI at fixed resolution, on
the cubed sphere and on tetrahedra.

## Execution

```bash
./Allrun parallel -N 24 -decomposition scotch -np <ranks> -config <config>
./Allrun parallel -mesh tet -N 10 -decomposition scotch -np <ranks> -config <config>
```

On the cubed sphere ATBI converges and ATI holds `Vm` at the wall (see
`hexConvergence`); on tetrahedra ATBI diverges (see `tetConvergence`).

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
