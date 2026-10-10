# tetConvergence - monodomainPseudoECG/sphericalShell

Spatial convergence on gmsh tetrahedra (`Allrun -mesh tet`): sealed curved
walls approximated by planar faces, an oblique full conductivity tensor and
skewed wall cells.

## Ladder

| N | lc | cells (three seeds) |
|---|---|---|
| 6 | 0.167 | 4430 - 4500 |
| 8 | 0.125 | 9815 - 9963 |
| 10 | 0.1 | 18059 - 18177 |
| 14 | 0.0714 | 46846 - 47129 |

Meshes at different `N` are generated independently. Every level runs with
`-seed 1`, `-seed 2` and `-seed 3`; the tables give the mean over the seeds.
The verifier writes `postProcessing/3D_tet_<N>_cells.dat` (rows and header
lines as in [hexConvergence](../hexConvergence/README.md)).

## Execution

In a copy of the case, for each `N` and `s`, once per wall trace:

```bash
./Allrun parallel -mesh tet -N <n> -seed <s> -trace zeroGradient
./Allrun parallel -mesh tet -N <n> -seed <s> -trace conormal
```

The default end time is `T = 0.2` (the short window) and the time step is
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`.

### Long window

The long window runs to `T = 2` with a control at `T = 4`, and repeats the
run with the default time step halved:

```bash
./Allrun parallel -mesh tet -N <n> -seed <s> -trace <trace> -endTime 2
./Allrun parallel -mesh tet -N <n> -seed <s> -trace <trace> -endTime 4
./Allrun parallel -mesh tet -N <n> -seed <s> -trace <trace> -endTime 2 -dt <half the default>
```

## Measured behaviour (T = 0.2)

Errors are L2 norms against the exact fields, except the `Vm` Linf column;
`h` is the verifier's `# h`. Values are means over the three seeds; orders
are between successive levels of the means.

### `-trace zeroGradient`

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 6 | 0.0934 | 6.111e-03 | 4.217e-02 | 4.530e-03 | 2.057e-01 | 1.543e-01 |
| 8 | 0.0718 | 4.364e-03 | 3.679e-02 | 3.153e-03 | 1.743e-01 | 1.143e-01 |
| 10 | 0.0587 | 3.528e-03 | 3.387e-02 | 2.483e-03 | 1.606e-01 | 9.547e-02 |
| 14 | 0.0427 | 2.306e-03 | 2.475e-02 | 1.678e-03 | 1.410e-01 | 7.222e-02 |

| N | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|
| 6 -> 8 | 1.28 | 0.52 | 1.38 | 0.63 | 1.14 |
| 8 -> 10 | 1.05 | 0.41 | 1.18 | 0.40 | 0.89 |
| 10 -> 14 | 1.34 | 0.99 | 1.23 | 0.41 | 0.88 |

The `Vm` L2 order is between 1.05 and 1.34. Over the seeds at one `N`, the
`Vm` L2 error varies by at most 3 % (4.308e-03 to 4.431e-03 at `N = 8`;
2.275e-03 to 2.334e-03 at `N = 14`).

### `-trace conormal`

| N | h | Vm L2 | Vm Linf | Vm_boundary L2 | gradVm_boundary L2 | gradReconVm_boundary L2 |
|---|---|---|---|---|---|---|
| 6 | 0.0934 | 1.533e-02 | 1.435e+00 | 1.362e-02 | 8.401e+00 | 1.543e-01 |
| 8 | 0.0718 | 5.793e-03 | 9.625e-02 | 4.470e-03 | 2.252e-01 | 1.143e-01 |
| 10 | 0.0587 | 4.496e-03 | 7.338e-02 | 3.365e-03 | 1.816e-01 | 9.547e-02 |
| 14 | 0.0427 | 5.417e-03 | 8.927e-01 | 2.492e-03 | 8.663e+00 | 7.222e-02 |

These means include three meshes with a local wall-gradient blow-up, so the
orders between levels are not meaningful (they range from -12.17 to +13.77
for `gradVm_boundary`). Per seed:

| run | Vm L2 | Vm Linf | gradVm_boundary L2 | gradVm_boundary Linf |
|---|---|---|---|---|
| N 6, seed 1 | 8.36e-03 | 9.76e-02 | 3.42e-01 | 1.29e+01 |
| N 6, seed 2 | 8.05e-03 | 7.19e-02 | 3.03e-01 | 6.85e+00 |
| N 6, seed 3 | 2.96e-02 | 4.14e+00 | 2.46e+01 | 3.75e+03 |
| N 8, seed 1 | 6.04e-03 | 7.48e-02 | 2.25e-01 | 4.10e+00 |
| N 8, seed 2 | 5.65e-03 | 9.15e-02 | 2.17e-01 | 3.45e+00 |
| N 8, seed 3 | 5.68e-03 | 1.22e-01 | 2.34e-01 | 1.12e+01 |
| N 10, seed 1 | 4.48e-03 | 7.30e-02 | 1.71e-01 | 3.88e+00 |
| N 10, seed 2 | 4.49e-03 | 9.04e-02 | 2.02e-01 | 1.64e+01 |
| N 10, seed 3 | 4.52e-03 | 5.67e-02 | 1.71e-01 | 3.44e+00 |
| N 14, seed 1 | 2.54e-03 | 4.51e-02 | 1.09e-01 | 3.22e+00 |
| N 14, seed 2 | 4.66e-03 | 7.62e-01 | 7.89e+00 | 2.26e+03 |
| N 14, seed 3 | 9.05e-03 | 1.87e+00 | 1.80e+01 | 5.38e+03 |

Three of the twelve meshes (`N = 6` seed 3, `N = 14` seeds 2 and 3) show a
local wall-gradient blow-up, with `gradVm_boundary` Linf of order 1e3; on the
other nine `gradVm_boundary` Linf stays between 3.2 and 16.4.
`gradReconVm_boundary`, computed from the exact field, is unaffected.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
