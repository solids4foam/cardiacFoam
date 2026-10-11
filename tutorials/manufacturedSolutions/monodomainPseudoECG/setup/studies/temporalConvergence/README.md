# temporalConvergence - monodomainPseudoECG

Fixed-grid temporal convergence of the FDA cosine field on Cartesian
hexahedra, for the four combinations of ionic-PDE time coupling
(`-coupling godunov|sbdf2`) and time scheme of `Vm` (`-ddt Euler|backward`).
The mesh is fine enough that the time step dominates the error over the
ladder; the error is measured against the exact field at `endTime 0.2`.

## Ladder

| dimension | N | dt |
|---|---|---|
| 1D, 2D | 160, 320, 640 | 0.025, 0.0125, 0.00625, 0.003125, 0.0015625 |
| 1D, 2D (`sbdf2`, `backward` only) | as above | also 0.00078125, 0.000390625 |
| 1D (`sbdf2`, `backward` only) | 1280 | 0.025 down to 0.000390625 |
| 3D (`sbdf2`, `backward` only) | 80 | 0.025 down to 0.000390625 |

Every time step divides `endTime 0.2`, so every run ends at 0.2.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<d>D_N<n>_<coupling>_<ddt>_dt<dt>`, the names that `../summarize.py` reads
(each `<d>D_N<n>_<coupling>_<ddt>` is one temporal ladder). From this folder,
the main matrix:

```bash
mkdir -p results
for d in 1 2
do
    for n in 160 320 640
    do
        for c in godunov sbdf2
        do
            for t in Euler backward
            do
                for dt in 0.025 0.0125 0.00625 0.003125 0.0015625
                do
                    r=results/${d}D_N${n}_${c}_${t}_dt$dt
                    rsync -a --exclude setup ../../../ $r/
                    $r/Allrun -dim $d -field fda -N $n -dt $dt -coupling $c -ddt $t
                done
            done
        done
    done
done
python3 ../summarize.py results
```

The extended `sbdf2`/`backward` ladders add `-dt 0.00078125` and
`-dt 0.000390625`, `-dim 1 -N 1280` and `-dim 3 -N 80` with the same command.

## Measured behaviour

Results of the last run of this study (`results/`).

`Vm` L2 order between `dt = 0.003125` and `0.0015625` at `N = 640`:

| coupling | ddt | 1D | 2D |
|---|---|---|---|
| godunov | Euler | 1.00 | 1.00 |
| godunov | backward | 1.00 | 1.00 |
| sbdf2 | Euler | 1.04 | 1.11 |
| sbdf2 | backward | 2.04 | 2.58 |

For `sbdf2`/`backward` in 2D the orders of the coarser pairs are 1.97, 2.01
and 2.11 (`dt = 0.025` to `0.003125`); the finest pair already contains
spatial error.

Godunov splitting is first order whichever time scheme `Vm` uses; SBDF2
coupling with `backward` is second order. With SBDF2 and `backward` the error
reaches the spatial error of the mesh below `dt` of about 0.003 in 2D (and on
the 3D `N = 80` mesh), after which the observed order is no longer temporal.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
