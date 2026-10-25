# errorLocalisation - eikonalECG

Where the tetrahedral activation-time error lies, for both cell-gradient
reconstructions with the case defaults otherwise (rotated conductivity,
advection-diffusion formulation): the verifier's bulk/boundary split of the
volume-weighted L2 (`activationTime_bulk`, `activationTime_boundary`; the
boundary cells own a face on a physical patch, both parts are normalised by
the total volume, so `total^2 = bulk^2 + boundary^2`). The runs also write
the cellwise error and the cell centres (`-errorField`) for
`analyse_error_localisation.py`.

## Ladder

| N | lc | tetrahedra (gmsh 4.15.2, seed 1) |
|---|---|---|
| 10 | 0.1 | 4874 |
| 20 | 0.05 | 37040 |
| 30 | 0.0333 | 123211 |
| 40 | 0.025 | 289388 |

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_N<n>`, the names that `../summarize.py` reads. From this folder:

```bash
mkdir -p results
for n in 10 20 30 40
do
    for g in GaussLinear leastSquares
    do
        rsync -a --exclude setup ../../../ results/${g}_N$n/
        results/${g}_N$n/Allrun parallel -mesh tet -N $n -gradient $g -errorField
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_error_localisation.json`.

`analyse_error_localisation.py` bins the cellwise error of one run by the
distance from the planes `x, y, z = 0` next to the imposed cells and by the
distance from the nearest face, and reports the correlations:

```bash
(cd results/leastSquares_N40 && python3 ../../analyse_error_localisation.py)
```

## Measured behaviour

Results of the last run of this study (`results/`): activation-time L2 split
into bulk and boundary parts, with `h = 1/round(cbrt(cells))`, and the share
of the boundary part in the error energy (`boundary^2/total^2`).

| gradient | N | bulk | boundary | total | order bulk | order boundary | order total | boundary share |
|---|---|---|---|---|---|---|---|---|
| GaussLinear | 10 | 8.7532e-02 | 2.7964e-02 | 9.1891e-02 | | | | 9.3 % |
| | 20 | 8.0490e-02 | 1.0636e-02 | 8.1190e-02 | 0.13 | 1.46 | 0.19 | 1.7 % |
| | 30 | 8.1759e-02 | 5.5525e-03 | 8.1948e-02 | -0.04 | 1.56 | -0.02 | 0.46 % |
| | 40 | 7.9666e-02 | 3.6477e-03 | 7.9750e-02 | 0.09 | 1.51 | 0.10 | 0.21 % |
| leastSquares | 10 | 1.1315e-02 | 1.8572e-02 | 2.1747e-02 | | | | 73 % |
| | 20 | 4.9045e-03 | 6.9271e-03 | 8.4875e-03 | 1.26 | 1.49 | 1.42 | 67 % |
| | 30 | 2.6262e-03 | 3.6836e-03 | 4.5239e-03 | 1.50 | 1.52 | 1.51 | 66 % |
| | 40 | 1.8187e-03 | 2.4195e-03 | 3.0268e-03 | 1.32 | 1.51 | 1.45 | 64 % |

With `Gauss linear` the bulk error stalls near 8e-2 while the boundary part
converges at 1.5, so the error becomes bulk-dominated. With `leastSquares`
both parts converge at 1.3-1.5 and the boundary-adjacent cells keep about two
thirds of the error energy. The reconstruction of the exact gradient shows the
same split (`gradReconActivationTime` rows, see
[tetConvergence](../tetConvergence/README.md)): its `Gauss linear` bulk part
does not converge.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
