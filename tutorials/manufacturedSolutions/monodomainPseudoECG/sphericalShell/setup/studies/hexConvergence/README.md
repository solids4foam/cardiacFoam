# hexConvergence - monodomainPseudoECG/sphericalShell

Spatial convergence on the cubed sphere (`Allrun -mesh hex`, default): sealed
curved walls approximated by planar faces, an oblique full conductivity tensor
and nearly orthogonal wall cells.

## Ladder

| N | k | cells |
|---|---|---|
| 8 | 2 | 864 |
| 12 | 3 | 2916 |
| 16 | 4 | 6912 |
| 24 | 6 | 23328 |
| 32 | 8 | 55296 |

`N` is a multiple of 4 and `k = N/4`. The verifier writes
`postProcessing/3D_hex_<N>_cells.dat`: rows `Vm`, `u1`, `u2`, the bulk,
boundary and gradient rows, and the header lines `# h`, `# time` and
`# cellsPerDirection`.

## Execution

In a copy of the case, for each `N`, once per wall trace:

```bash
./Allrun parallel -mesh hex -N <n> -trace zeroGradient
./Allrun parallel -mesh hex -N <n> -trace conormal
```

The default end time is `T = 0.2` (the short window) and the time step is
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`, which scales as `h^2`.

### Long window

The long window runs to `T = 2` with a control at `T = 4`, and repeats the
run with the default time step halved:

```bash
./Allrun parallel -mesh hex -N <n> -trace <trace> -endTime 2
./Allrun parallel -mesh hex -N <n> -trace <trace> -endTime 4
./Allrun parallel -mesh hex -N <n> -trace <trace> -endTime 2 -dt <half the default>
```

## Measured behaviour (T = 0.2)

Errors are L2 norms against the exact fields, except the `Vm` Linf column;
`h` is the verifier's `# h`. Orders are between successive levels.

### `-trace zeroGradient`

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 8 | 0.1603 | 4.468e-03 | 2.699e-02 | 3.940e-03 | 2.402e-01 | 1.382e-01 |
| 12 | 0.1074 | 2.475e-03 | 1.750e-02 | 2.052e-03 | 1.846e-01 | 7.605e-02 |
| 16 | 0.0807 | 1.654e-03 | 1.132e-02 | 1.387e-03 | 1.624e-01 | 4.491e-02 |
| 24 | 0.0539 | 1.008e-03 | 6.762e-03 | 8.793e-04 | 1.385e-01 | 1.928e-02 |
| 32 | 0.0404 | 7.452e-04 | 6.242e-03 | 6.526e-04 | 1.235e-01 | 1.009e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 8 -> 12 | 1.48 | 1.08 | 1.63 | 0.66 | 1.49 |
| 12 -> 16 | 1.41 | 1.52 | 1.37 | 0.45 | 1.84 |
| 16 -> 24 | 1.22 | 1.28 | 1.13 | 0.39 | 2.09 |
| 24 -> 32 | 1.05 | 0.28 | 1.04 | 0.40 | 2.25 |

The `Vm` L2 order falls from 1.48 to 1.05 over the ladder.

### `-trace conormal`

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 8 | 0.1603 | 4.139e-03 | 2.601e-02 | 3.574e-03 | 2.774e-01 | 1.382e-01 |
| 12 | 0.1074 | 2.007e-03 | 1.587e-02 | 1.454e-03 | 1.386e-01 | 7.605e-02 |
| 16 | 0.0807 | 1.156e-03 | 9.221e-03 | 7.241e-04 | 8.106e-02 | 4.491e-02 |
| 24 | 0.0539 | 5.224e-04 | 4.200e-03 | 2.671e-04 | 3.818e-02 | 1.928e-02 |
| 32 | 0.0404 | 2.972e-04 | 2.639e-03 | 1.356e-04 | 2.314e-02 | 1.009e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 8 -> 12 | 1.81 | 1.24 | 2.25 | 1.73 | 1.49 |
| 12 -> 16 | 1.93 | 1.90 | 2.44 | 1.88 | 1.84 |
| 16 -> 24 | 1.97 | 1.95 | 2.47 | 1.86 | 2.09 |
| 24 -> 32 | 1.96 | 1.62 | 2.36 | 1.74 | 2.25 |

The `Vm` L2 order is 1.96 between `N = 24` and `32`.

## Serial and decomposed runs

A serial run and a four-rank scotch run (`parallel`) with `-N 12 -endTime 0.02`
agree to a maximum relative difference of 9.1e-12 over the `Vm`, `u1` and `u2`
rows.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
