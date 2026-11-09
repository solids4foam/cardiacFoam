# tetConvergence - bidomain

Spatial convergence on gmsh Delaunay tetrahedra of the unit cube
(`Allrun -mesh tet`, template `box.geo.template`) for both cell-gradient
reconstructions: the same mesh, problem and solver path with only
`-gradient` changed. The other settings are the case defaults (`godunov`
coupling, `backward`, two outer sweeps, `phiE` tolerance `1e-15`).

## Ladder

| N | lc | tetrahedra, seed 1 | seed 2 | seed 3 |
|---|---|---|---|---|
| 10 | 0.1 | 4874 | 4874 | 4874 |
| 20 | 0.05 | 37040 | 37038 | 37038 |
| 40 | 0.025 | 289388 | 289540 | 289327 |
| 80 | 0.0125 | 2278889 | | |

Cell counts are those of gmsh 4.15.2. `leastSquares` runs `N = 10-80` with
`-seed 1` and `N = 10-40` with `-seed 2` and `-seed 3`; `Gauss linear` runs
`N = 10-40` with `-seed 1`. Every run uses the `Allrun` default time step,
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`: it scales as `h^2` and every
level ends exactly at `endTime 0.2`, so the slopes are space-time slopes
([tetTemporalControl](../tetTemporalControl/README.md) shows they are
spatial). On tetrahedra the verifier labels its summary `3D_<n>_cells.dat`
with `n = round(cbrt(cells))` (17 at `N = 10`) and reports `h = 1/n`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_tet_N<n>` (seed 1) and `leastSquares_s<seed>_tet_N<n>`, the names
that `../summarize.py` reads. From this folder:

```bash
mkdir -p results
for n in 10 20 40 80
do
    r=results/leastSquares_tet_N$n
    rsync -a --exclude setup ../../../ $r/
    $r/Allrun parallel -mesh tet -gradient leastSquares -N $n
done
for n in 10 20 40
do
    r=results/GaussLinear_tet_N$n
    rsync -a --exclude setup ../../../ $r/
    $r/Allrun parallel -mesh tet -gradient GaussLinear -N $n
    for s in 2 3
    do
        r=results/leastSquares_s${s}_tet_N$n
        rsync -a --exclude setup ../../../ $r/
        $r/Allrun parallel -mesh tet -seed $s -gradient leastSquares -N $n
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_generic.json` (the seed-1
runs).

## Measured behaviour

Results of the last run of this study (`results/`), as printed by
`../summarize.py`; gauge-shifted potentials, orders between consecutive
levels of the verifier's `h`. `gradReconPhiE` is the error of the cell
gradient reconstructed from the exact `phiE`.

`leastSquares`, seed 1:

| N | h | Vm L2 | Vm Linf | phiE L2 | phiE Linf | phiI L2 | u1 L2 | gradReconPhiE L2 | order Vm L2 | order Vm Linf | order phiE L2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.0588235 | 6.9530e-03 | 3.8951e-02 | 4.9149e-03 | 2.7222e-02 | 2.0686e-03 | 1.1690e-03 | 6.9638e-01 |  |  |  |
| 20 | 0.030303 | 2.4118e-03 | 1.6799e-02 | 1.7096e-03 | 1.1683e-02 | 7.0792e-04 | 5.0241e-04 | 4.2948e-01 | 1.60 | 1.27 | 1.59 |
| 40 | 0.0151515 | 6.3151e-04 | 5.7885e-03 | 4.4747e-04 | 4.1461e-03 | 1.8553e-04 | 1.5532e-04 | 2.7603e-01 | 1.93 | 1.54 | 1.93 |
| 80 | 0.00757576 | 1.6045e-04 | 1.7178e-03 | 1.1365e-04 | 1.2290e-03 | 4.7180e-05 | 4.2243e-05 | 1.8355e-01 | 1.98 | 1.75 | 1.98 |

`Gauss linear`, seed 1:

| N | h | Vm L2 | Vm Linf | phiE L2 | phiE Linf | phiI L2 | u1 L2 | gradReconPhiE L2 | order Vm L2 | order Vm Linf | order phiE L2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.0588235 | 1.3921e-02 | 7.9104e-02 | 9.8631e-03 | 5.6086e-02 | 4.0658e-03 | 2.3607e-03 | 1.3723e+00 |  |  |  |
| 20 | 0.030303 | 8.0524e-03 | 6.2525e-02 | 5.7110e-03 | 4.4227e-02 | 2.3425e-03 | 1.6621e-03 | 1.3793e+00 | 0.83 | 0.35 | 0.82 |
| 40 | 0.0151515 | 5.0344e-03 | 4.5213e-02 | 3.5670e-03 | 3.1940e-02 | 1.4680e-03 | 1.0675e-03 | 1.3934e+00 | 0.68 | 0.47 | 0.68 |

With `leastSquares` the L2 orders rise to 1.98 at the finest pair and the
`Vm` Linf order to 1.75. With `Gauss linear` the L2 orders fall from 0.83 to
0.68 and the exact-field gradient error does not decrease (1.37 to 1.39), so
the ladder stops at `N = 40`.

`leastSquares` over three gmsh seeds:

| seed | N | Vm L2 | Vm Linf | phiE L2 | order Vm L2 | order Vm Linf | order phiE L2 |
|---|---|---|---|---|---|---|---|
| 1 | 10 | 6.9530e-03 | 3.8951e-02 | 4.9149e-03 |  |  |  |
| 1 | 20 | 2.4118e-03 | 1.6799e-02 | 1.7096e-03 | 1.60 | 1.27 | 1.59 |
| 1 | 40 | 6.3151e-04 | 5.7885e-03 | 4.4747e-04 | 1.93 | 1.54 | 1.93 |
| 2 | 10 | 6.9549e-03 | 3.8951e-02 | 4.9160e-03 |  |  |  |
| 2 | 20 | 2.4071e-03 | 1.6799e-02 | 1.7063e-03 | 1.60 | 1.27 | 1.60 |
| 2 | 40 | 6.3260e-04 | 5.7288e-03 | 4.4823e-04 | 1.93 | 1.55 | 1.93 |
| 3 | 10 | 6.9521e-03 | 3.8951e-02 | 4.9146e-03 |  |  |  |
| 3 | 20 | 2.4071e-03 | 1.6740e-02 | 1.7063e-03 | 1.60 | 1.27 | 1.59 |
| 3 | 40 | 6.3264e-04 | 5.7885e-03 | 4.4826e-04 | 1.93 | 1.53 | 1.93 |

The three seeds agree to within 0.2 % in the L2 errors and give the same
orders to within 0.02.

Every linear solve reaches its absolute tolerance `1e-15` below the explicit
`maxIter 5000`; the largest iteration count is 963 for `phiE`
(`leastSquares`, `N = 80`) and 13 for `Vm`.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
