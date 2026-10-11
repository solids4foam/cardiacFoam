# temporalConvergence - bidomain

Fixed-grid temporal convergence on Cartesian hexahedra for the four
combinations of ionic-PDE time coupling (`-coupling godunov|sbdf2`) and time
scheme of `Vm` (`-ddt Euler|backward`). The 2D `N = 640` mesh is fine enough
that the time step dominates the error over the ladder; the 3D `N = 80` mesh
shows where the spatial error takes over. The error is measured against the
exact field at `endTime 0.2`.

## Ladder

| dimension | N | cells | dt |
|---|---|---|---|
| 2D | 640 | 409600 | 0.025, 0.0125, 0.00625, 0.003125, 0.0015625 |
| 3D | 80 | 512000 | as above |

Every time step divides `endTime 0.2`, so every run ends at 0.2. The runs
use `-gradient GaussLinear`, which on these orthogonal meshes with diagonal
conductivities gives the same operators as `leastSquares`. At 2D `N = 640`
the `phiE` solve stops at the PCG limit of 1000 iterations, with a final
residual near 1e-6, so the last digits of these runs depend on the
decomposition.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<d>D_N<n>_<coupling>_<ddt>_dt<dt>`, the names that `../summarize.py` reads
(each `<d>D_N<n>_<coupling>_<ddt>` is one temporal ladder). From this folder:

```bash
mkdir -p results
for spec in "2 640" "3 80"
do
    set -- $spec; d=$1; n=$2
    for c in godunov sbdf2
    do
        for t in Euler backward
        do
            for dt in 0.025 0.0125 0.00625 0.003125 0.0015625
            do
                r=results/${d}D_N${n}_${c}_${t}_dt$dt
                rsync -a --exclude setup ../../../ $r/
                $r/Allrun parallel -dim $d -N $n -dt $dt -coupling $c -ddt $t -gradient GaussLinear
            done
        done
    done
done
python3 ../summarize.py results
```

## Measured behaviour

Results of the last run of this study (`results/`).

2D `N = 640`, `Vm` L2:

| coupling | ddt | 0.025 | 0.0125 | 0.00625 | 0.003125 | 0.0015625 | orders |
|---|---|---|---|---|---|---|---|
| godunov | Euler | 1.2593e-03 | 6.3443e-04 | 3.1833e-04 | 1.5936e-04 | 7.9653e-05 | 0.99, 0.99, 1.00, 1.00 |
| godunov | backward | 1.0475e-03 | 5.1405e-04 | 2.5439e-04 | 1.2645e-04 | 6.2957e-05 | 1.03, 1.01, 1.01, 1.01 |
| sbdf2 | Euler | 4.4622e-04 | 1.8376e-04 | 8.1555e-05 | 3.8071e-05 | 1.8272e-05 | 1.28, 1.17, 1.10, 1.06 |
| sbdf2 | backward | 2.8633e-04 | 7.3635e-05 | 1.8533e-05 | 4.5129e-06 | 9.7853e-07 | 1.96, 1.99, 2.04, 2.21 |

Godunov splitting is first order whichever time scheme `Vm` uses, SBDF2 with
`Euler` is first order, and SBDF2 with `backward` is second order over this
window (2.21 at the finest pair, where the spatial error starts to enter).

3D `N = 80`, `Vm` L2:

| coupling | ddt | 0.025 | 0.0125 | 0.00625 | 0.003125 | 0.0015625 |
|---|---|---|---|---|---|---|
| godunov | Euler | 9.8437e-04 | 4.8373e-04 | 2.3182e-04 | 1.0557e-04 | 4.2510e-05 |
| godunov | backward | 8.3035e-04 | 3.9687e-04 | 1.8605e-04 | 8.2292e-05 | 3.1140e-05 |
| sbdf2 | Euler | 3.1911e-04 | 1.1485e-04 | 3.7765e-05 | 5.8156e-06 | 8.5106e-06 |
| sbdf2 | backward | 2.0997e-04 | 3.8640e-05 | 7.9868e-06 | 1.8030e-05 | 2.0754e-05 |

On the 3D mesh the SBDF2 errors stop decreasing below `dt = 0.00625`
(`backward`) and `dt = 0.003125` (`Euler`): the error has reached the spatial
error of the mesh, and the 3D ladder gives no temporal order there.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
