# cartesianConvergence - bidomain

Spatial convergence on Cartesian hexahedra (`Allrun`, `-mesh hex`) in 1D, 2D
and 3D, with the case defaults (`godunov` coupling, `backward`, two outer
sweeps, `leastSquares`).

## Ladder

| N | dt | hex cells (3D) |
|---|---|---|
| 10 | 0.00892857 | 1000 |
| 20 | 0.00224215 | 8000 |
| 40 | 0.000560538 | 64000 |
| 80 | 0.000140174 | 512000 |

The time step falls as h^2, so the slopes are space-time slopes. These steps
do not divide `endTime 0.2`, so each run stops at the step nearest to 0.2
(0.196429, 0.199551, 0.200112, 0.200028; the verifier writes the final time
in its summary and evaluates the exact field there). Without `-dt` the case uses
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`, which ends exactly at 0.2.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<d>D_N<n>`, the names that `../summarize.py` reads. From this folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.00892857 0.00224215 0.000560538 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    for d in 1 2 3
    do
        rsync -a --exclude setup ../../../ results/${d}D_N$n/
        results/${d}D_N$n/Allrun -dim $d -N $n -dt $dt
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_hex_convergence.json`.

## Measured behaviour

Results of the last run of this study (`results/`), gauge-shifted potentials.

`Vm`:

| dim | N | L2 | Linf | order L2 | order Linf |
|---|---|---|---|---|---|
| 1D | 10 | 2.7785e-04 | 5.5740e-04 |  |  |
| 1D | 20 | 6.9473e-05 | 1.4106e-04 | 2.00 | 1.98 |
| 1D | 40 | 1.7327e-05 | 3.5274e-05 | 2.00 | 2.00 |
| 1D | 80 | 4.3262e-06 | 8.8128e-06 | 2.00 | 2.00 |
| 2D | 10 | 5.3690e-04 | 1.2740e-03 |  |  |
| 2D | 20 | 1.4051e-04 | 3.4990e-04 | 1.93 | 1.86 |
| 2D | 40 | 3.5507e-05 | 8.9439e-05 | 1.98 | 1.97 |
| 2D | 80 | 8.8881e-06 | 2.2456e-05 | 2.00 | 1.99 |
| 3D | 10 | 1.0461e-03 | 2.9627e-03 |  |  |
| 3D | 20 | 2.7395e-04 | 8.1882e-04 | 1.93 | 1.86 |
| 3D | 40 | 6.9223e-05 | 2.0963e-04 | 1.98 | 1.97 |
| 3D | 80 | 1.7330e-05 | 5.2654e-05 | 2.00 | 1.99 |

`phiE_gauge`:

| dim | N | L2 | Linf | order L2 | order Linf |
|---|---|---|---|---|---|
| 1D | 10 | 1.7221e-04 | 2.9918e-04 |  |  |
| 1D | 20 | 4.2442e-05 | 7.4991e-05 | 2.02 | 2.00 |
| 1D | 40 | 1.0545e-05 | 1.8703e-05 | 2.01 | 2.00 |
| 1D | 80 | 2.6298e-06 | 4.6689e-06 | 2.00 | 2.00 |
| 2D | 10 | 3.7001e-04 | 8.3581e-04 |  |  |
| 2D | 20 | 9.8042e-05 | 2.3303e-04 | 1.92 | 1.84 |
| 2D | 40 | 2.4851e-05 | 5.9785e-05 | 1.98 | 1.96 |
| 2D | 80 | 6.2253e-06 | 1.5024e-05 | 2.00 | 1.99 |
| 3D | 10 | 7.3326e-04 | 2.0549e-03 |  |  |
| 3D | 20 | 1.9324e-04 | 5.7204e-04 | 1.92 | 1.84 |
| 3D | 40 | 4.8903e-05 | 1.4671e-04 | 1.98 | 1.96 |
| 3D | 80 | 1.2248e-05 | 3.6868e-05 | 2.00 | 1.99 |

`phiI_gauge`:

| dim | N | L2 | Linf | order L2 | order Linf |
|---|---|---|---|---|---|
| 1D | 10 | 1.5386e-04 | 2.5822e-04 |  |  |
| 1D | 20 | 3.9289e-05 | 6.6066e-05 | 1.97 | 1.97 |
| 1D | 40 | 9.8517e-06 | 1.6570e-05 | 2.00 | 2.00 |
| 1D | 80 | 2.4633e-06 | 4.1440e-06 | 2.00 | 2.00 |
| 2D | 10 | 1.7860e-04 | 4.3818e-04 |  |  |
| 2D | 20 | 4.5437e-05 | 1.1687e-04 | 1.97 | 1.91 |
| 2D | 40 | 1.1398e-05 | 2.9655e-05 | 2.00 | 1.98 |
| 2D | 80 | 2.8484e-06 | 7.4327e-06 | 2.00 | 2.00 |
| 3D | 10 | 3.1412e-04 | 9.0772e-04 |  |  |
| 3D | 20 | 8.1013e-05 | 2.4678e-04 | 1.96 | 1.88 |
| 3D | 40 | 2.0393e-05 | 6.2919e-05 | 1.99 | 1.97 |
| 3D | 80 | 5.0999e-06 | 1.5786e-05 | 2.00 | 1.99 |

The auxiliary states `u1` and `u2` converge at second order in every
dimension (L2 orders 2.00-2.01 at the finest pair). `Vm`, `phiE` and `phiI`
reach second order in L2 and Linf in 1D, 2D and 3D.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
