# patchTest - bathBidomain/obliqueSlab

The G0 affine solution: piecewise affine potentials in the interface normal
coordinate with constant conductivities, reproduced exactly by a consistent
discretisation. See [G0: three roles](../../../README.md#g0-three-roles).

## Execution

```bash
./Allrun -solution affine -elliptic -nOuter 10 -nNonOrth 10 -shear 0 -N 10
./Allrun -solution affine -elliptic -nOuter 10 -nNonOrth 10 -N 10
./Allrun -solution affine -elliptic -nOuter 10 -nNonOrth 10 -N 10 -grad GaussLinear
./Allrun -solution affine -elliptic -nOuter 10 -nNonOrth 10 -mesh tet -N 10
./Allrun -solution affine -elliptic -nOuter 10 -nNonOrth 10 -mesh tet -N 10 -grad GaussLinear
```

## Expected behaviour

| mesh | gradient | `phiE_heart` Linf |
|---|---|---|
| orthogonal hex | leastSquares | round-off |
| sheared hex | leastSquares | round-off |
| sheared hex | GaussLinear | round-off |
| tet (N = 10, 20) | leastSquares | 1.4e-5, 1.2e-5 |
| tet (N = 10) | GaussLinear | 5.6e-3 |

On tetrahedra the jump-aware interface gradient does not reproduce the
piecewise-linear field (`gradJumpAwareTangentialError` Linf 1.8e-3 at every
`N`), which leaves the leastSquares floor of about 1e-5.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
