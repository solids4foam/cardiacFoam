# tetConvergence - bidomain

Spatial convergence on gmsh Delaunay tetrahedra of the unit cube
(`Allrun -mesh tet`, template `box.geo.template`) for both cell-gradient
reconstructions: the same mesh, problem and solver path with only
`-gradient` changed. The other settings are the case defaults (`godunov`
coupling, `backward`, two outer sweeps, `phiE` tolerance `1e-15`).

## Ladder

| N | lc | dt | tetrahedra (gmsh 4.15.2, seed 1) |
|---|---|---|---|
| 10 | 0.1 | 0.00892857 | 4874 |
| 20 | 0.05 | 0.00224215 | 37040 |
| 40 | 0.025 | 0.000560538 | 289388 |
| 80 | 0.0125 | 0.000140174 | |

The time steps are those of
[cartesianConvergence](../cartesianConvergence/README.md), so the slopes are
space-time slopes ([tetTemporalControl](../tetTemporalControl/README.md)).
On tetrahedra the verifier labels its summary `3D_<n>_cells.dat` with
`n = round(cbrt(cells))` (17 at `N = 10`) and reports `h = 1/n`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_tet_N<n>`, the names that `../summarize.py` reads. From this
folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.00892857 0.00224215 0.000560538 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    for g in leastSquares GaussLinear
    do
        rsync -a --exclude setup ../../../ results/${g}_tet_N$n/
        results/${g}_tet_N$n/Allrun parallel -mesh tet -gradient $g -N $n -dt $dt
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_generic.json`.

## Measured behaviour

Results of the last run of this study (`results/`), on tetrahedra from gmsh
4.14.0 (4979, 37142, 289960 and 2276675 cells). Orders use
`h = cells^(-1/3)`; gauge-shifted potentials.

`leastSquares`:

| N | h | Vm L2 | order | Vm Linf | order | phiE L2 | order | phiE Linf | order | phiI L2 | order |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.0586 | 6.8568e-03 |  | 4.5165e-02 |  | 4.8427e-03 |  | 3.1871e-02 |  | 2.0293e-03 |  |
| 20 | 0.0300 | 2.4031e-03 | 1.565 | 1.5819e-02 | 1.566 | 1.7021e-03 | 1.561 | 1.1279e-02 | 1.551 | 7.0384e-04 | 1.581 |
| 40 | 0.0151 | 6.3413e-04 | 1.945 | 5.4785e-03 | 1.548 | 4.4894e-04 | 1.946 | 3.8757e-03 | 1.559 | 1.8589e-04 | 1.944 |
| 80 | 0.0076 | 1.6012e-04 | 2.004 | 1.8968e-03 | 1.544 | 1.1333e-04 | 2.004 | 1.3494e-03 | 1.536 | 4.6976e-05 | 2.002 |

`GaussLinear`:

| N | h | Vm L2 | order | Vm Linf | order | phiE L2 | order | phiE Linf | order | phiI L2 | order |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.0586 | 1.3238e-02 |  | 7.7293e-02 |  | 9.3680e-03 |  | 5.4511e-02 |  | 3.8734e-03 |  |
| 20 | 0.0300 | 8.2035e-03 | 0.714 | 7.8342e-02 | -0.020 | 5.8095e-03 | 0.713 | 5.5398e-02 | -0.024 | 2.3946e-03 | 0.718 |
| 40 | 0.0151 | 5.1002e-03 | 0.694 | 4.1856e-02 | 0.915 | 3.6100e-03 | 0.695 | 2.9640e-02 | 0.913 | 1.4908e-03 | 0.692 |
| 80 | 0.0076 | 3.7548e-03 | 0.446 | 2.4883e-02 | 0.757 | 2.6561e-03 | 0.447 | 1.7565e-02 | 0.762 | 1.0993e-03 | 0.443 |

With `leastSquares` the potentials reach second order in L2 (2.00 at the
finest pair) and about 1.55 in Linf; `u1` and `u2` approach second order (1.90
at the finest pair). With `Gauss linear` the L2 orders fall from 0.7 to 0.45
and the Linf error does not decrease from `N = 10` to `N = 20`.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
