# tetConvergence - monodomainPseudoECG

Spatial convergence on gmsh Delaunay tetrahedra of the unit cube
(`Allrun -mesh tet`, template `box.geo.template`), for both manufactured fields
and both cell-gradient reconstructions: the same mesh, problem and solver path
with only `-gradient` changed. The pseudo-ECG is evaluated at all 161
electrodes.

## Ladder

| N | lc | dt | tetrahedra (gmsh 4.15.2, seed 1) |
|---|---|---|---|
| 10 | 0.1 | 0.00892857 | 4874 |
| 20 | 0.05 | 0.00224215 | 37040 |
| 40 | 0.025 | 0.000560538 | |
| 80 | 0.0125 | 0.000140174 | |

The time steps are those of
[cartesianConvergence](../cartesianConvergence/README.md). On tetrahedra the
verifier labels its summary `3D_<n>_cells.dat` with `n = round(cbrt(cells))`
(17 at `N = 10`) and reports `h = 1/n`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_<field>_tet_N<n>`, the names that `../summarize.py` reads. From
this folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.00892857 0.00224215 0.000560538 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    for g in leastSquares GaussLinear
    do
        for f in fda anisotropic
        do
            rsync -a --exclude setup ../../../ results/${g}_${f}_tet_N$n/
            results/${g}_${f}_tet_N$n/Allrun parallel -mesh tet -gradient $g -field $f -N $n -dt $dt
        done
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_generic.json`.

## Measured behaviour

Results of the last run of this study (`results/`): volume-weighted `Vm` L2
and the 161-electrode mean of each electrode's pseudo-ECG L2 error against
the reference integral.

| gradient | field | N | h | Vm L2 | order | ECG mean L2 | order |
|---|---|---|---|---|---|---|---|
| leastSquares | FDA (axis-aligned) | 10 | 0.0588 | 1.0844e-02 | | 1.0170e-03 | |
| | | 20 | 0.0303 | 2.8380e-03 | 2.02 | 2.7944e-04 | 1.94 |
| | | 40 | 0.0152 | 6.8162e-04 | 2.05 | 7.5553e-05 | 1.88 |
| | | 80 | 0.0076 | 1.6659e-04 | 2.03 | 2.0499e-05 | 1.88 |
| leastSquares | anisotropic | 10 | 0.0588 | 2.8001e-03 | | 4.3475e-04 | |
| | | 20 | 0.0303 | 7.0690e-04 | 2.07 | 1.0259e-04 | 2.17 |
| | | 40 | 0.0152 | 1.7462e-04 | 2.01 | 2.4409e-05 | 2.07 |
| | | 80 | 0.0076 | 4.3828e-05 | 1.99 | 5.8151e-06 | 2.06 |
| GaussLinear | FDA (axis-aligned) | 10 | 0.0588 | 2.2650e-02 | | 1.3276e-03 | |
| | | 20 | 0.0303 | 1.3201e-02 | 0.81 | 4.0780e-04 | 1.77 |
| | | 40 | 0.0152 | 1.1193e-02 | 0.23 | 2.2066e-04 | 0.88 |
| | | 80 | 0.0076 | 1.0780e-02 | 0.05 | 1.9313e-04 | 0.19 |
| GaussLinear | anisotropic | 10 | 0.0588 | 6.0136e-03 | | 2.9845e-04 | |
| | | 20 | 0.0303 | 3.4220e-03 | 0.84 | 9.5035e-05 | 1.72 |
| | | 40 | 0.0152 | 2.5821e-03 | 0.40 | 5.4393e-05 | 0.80 |
| | | 80 | 0.0076 | 2.2730e-03 | 0.18 | 4.4284e-05 | 0.29 |

With `leastSquares` the field and the pseudo-ECG converge at second order on
this family; with `Gauss linear` the field stalls and the pseudo-ECG
decelerates with it.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
