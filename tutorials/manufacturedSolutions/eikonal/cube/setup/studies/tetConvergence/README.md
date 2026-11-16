# tetConvergence - eikonalECG

Spatial convergence on gmsh Delaunay tetrahedra of the unit cube
(`Allrun -mesh tet`, template `box.geo.template`), with one setting changed
at a time on the same meshes: the cell-gradient reconstruction
(`Gauss linear` against `leastSquares`, axis conductivity), the conductivity
(axis against rotated, `leastSquares`) and the formulation of the metric term
(explicit against advection-diffusion, rotated conductivity). Every run also
reports the reconstruction of the exact gradient alone
(`gradReconActivationTime` rows), which depends only on the mesh and the
gradient scheme.

## Ladder

| N | lc | tetrahedra (gmsh 4.15.2, seed 1) | h = 1/round(cbrt(cells)) |
|---|---|---|---|
| 10 | 0.1 | 4874 | 1/17 |
| 20 | 0.05 | 37040 | 1/33 |
| 30 | 0.0333 | 123211 | 1/50 |
| 40 | 0.025 | 289388 | 1/66 |

The levels are not nested; orders use the realised `h`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_<conductivity>_<formulation>[_sine]_N<n>`, the names that
`../summarize.py` reads. From this folder:

```bash
mkdir -p results
for n in 10 20 30 40
do
    for c in "GaussLinear axis explicit" "leastSquares axis explicit" \
             "leastSquares rotated explicit" "leastSquares rotated advectionDiffusion"
    do
        set -- $c
        for t in singleEvent sine
        do
            run=results/$1_$2_$3_N$n
            [ $t = sine ] && run=results/$1_$2_$3_sine_N$n
            rsync -a --exclude setup ../../../ $run/
            $run/Allrun parallel -mesh tet -N $n -gradient $1 -conductivity $2 -formulation $3 -template $t
        done
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_generic.json`.

## Measured behaviour

Results of the last run of this study (`results/`): activation-time L2 and
Linf, the largest L2 error over E1-E5 of the ECG with `-template sine`, and
the outer iterations to the residual tolerance `1e-8`.

| gradient | conductivity | formulation | N | h | psi L2 | order | psi Linf | order | ECG max L2 | order | outer iterations |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GaussLinear | axis | explicit | 10 | 0.0588 | 7.2899e-02 | | 2.8343e-01 | | 4.4699e-02 | | 110 |
| | | | 20 | 0.0303 | 5.8120e-02 | 0.34 | 2.5080e-01 | 0.18 | 7.2922e-03 | 2.73 | 380 |
| | | | 30 | 0.0200 | 5.7760e-02 | 0.01 | 2.2262e-01 | 0.29 | 3.1978e-03 | 1.98 | 811 |
| | | | 40 | 0.0152 | 5.5746e-02 | 0.13 | 2.0430e-01 | 0.31 | 1.1751e-03 | 3.61 | 1393 |
| leastSquares | axis | explicit | 10 | 0.0588 | 2.2472e-02 | | 2.3915e-01 | | 3.7517e-02 | | 114 |
| | | | 20 | 0.0303 | 8.7251e-03 | 1.43 | 1.5229e-01 | 0.68 | 8.9773e-03 | 2.16 | 394 |
| | | | 30 | 0.0200 | 4.4245e-03 | 1.63 | 1.0902e-01 | 0.80 | 6.0407e-03 | 0.95 | 837 |
| | | | 40 | 0.0152 | 2.8968e-03 | 1.53 | 8.7813e-02 | 0.78 | 1.6613e-03 | 4.65 | 1440 |
| leastSquares | rotated | explicit | 10 | 0.0588 | 2.2572e-02 | | 2.4153e-01 | | 2.5775e-02 | | 115 |
| | | | 20 | 0.0303 | 8.8130e-03 | 1.42 | 1.5080e-01 | 0.71 | 6.9730e-03 | 1.97 | 399 |
| | | | 30 | 0.0200 | 4.5838e-03 | 1.57 | 1.0936e-01 | 0.77 | 4.5146e-03 | 1.05 | 851 |
| | | | 40 | 0.0152 | 3.0698e-03 | 1.44 | 8.6561e-02 | 0.84 | 1.2466e-03 | 4.64 | 1463 |
| leastSquares | rotated | advectionDiffusion | 10 | 0.0588 | 2.1747e-02 | | 2.4532e-01 | | 2.3917e-02 | | 129 |
| | | | 20 | 0.0303 | 8.4877e-03 | 1.42 | 1.5163e-01 | 0.73 | 6.9006e-03 | 1.87 | 419 |
| | | | 30 | 0.0200 | 4.5239e-03 | 1.51 | 1.0947e-01 | 0.78 | 4.5710e-03 | 0.99 | 876 |
| | | | 40 | 0.0152 | 3.0268e-03 | 1.45 | 8.8702e-02 | 0.76 | 1.2583e-03 | 4.65 | 1492 |

With `leastSquares` the activation time converges at about 1.5 in L2 and 0.7-0.8
in Linf, for both conductivities and both formulations; the two formulations
differ by at most 4 %. With `Gauss linear` the activation time stalls near
5.6e-2. The ECG at E1-E5 decreases faster than the field under both
gradients, without a stable pairwise order on this non-nested family, and
under `Gauss linear` it keeps decreasing while the field it is computed from
has stalled.

With the single-event template the activation time is the same, and the
largest L2 error over E1-E5 of the ECG is:

| gradient | conductivity | formulation | N | ECG max L2 | order |
|---|---|---|---|---|---|
| GaussLinear | axis | explicit | 10 | 3.4933e-04 | |
| | | | 20 | 2.8349e-04 | 0.31 |
| | | | 30 | 3.2797e-04 | -0.35 |
| | | | 40 | 3.2441e-04 | 0.04 |
| leastSquares | axis | explicit | 10 | 2.7274e-04 | |
| | | | 20 | 1.1155e-04 | 1.35 |
| | | | 30 | 6.6450e-05 | 1.25 |
| | | | 40 | 4.8864e-05 | 1.11 |
| leastSquares | rotated | explicit | 10 | 2.1296e-04 | |
| | | | 20 | 8.3085e-05 | 1.42 |
| | | | 30 | 4.8243e-05 | 1.31 |
| | | | 40 | 3.6291e-05 | 1.03 |
| leastSquares | rotated | advectionDiffusion | 10 | 2.3475e-04 | |
| | | | 20 | 7.2249e-05 | 1.78 |
| | | | 30 | 4.1792e-05 | 1.32 |
| | | | 40 | 3.2227e-05 | 0.94 |

Under `Gauss linear` the single-event ECG stalls with the field; under
`leastSquares` it decreases with the field. The sine ECG converges at orders
2 to 3 while the activation time it is computed from stalls; its integral
cancels.

Reconstruction of the exact gradient (`gradReconActivationTime`, L2 split as
in [errorLocalisation](../errorLocalisation/README.md)):

| gradient | N | bulk L2 | boundary L2 | total L2 | Linf | order bulk | order boundary | order total |
|---|---|---|---|---|---|---|---|---|
| leastSquares | 10 | 2.1819e-01 | 1.2487e+00 | 1.2676e+00 | 1.9887e+01 | | | |
| | 20 | 1.1573e-01 | 8.7954e-01 | 8.8712e-01 | 2.7381e+01 | 0.96 | 0.53 | 0.54 |
| | 30 | 7.8510e-02 | 7.0664e-01 | 7.1099e-01 | 2.6716e+01 | 0.93 | 0.53 | 0.53 |
| | 40 | 5.9348e-02 | 6.0811e-01 | 6.1099e-01 | 2.6400e+01 | 1.01 | 0.54 | 0.55 |
| GaussLinear | 10 | 3.7998e+00 | 1.6485e+00 | 4.1420e+00 | 3.3776e+01 | | | |
| | 20 | 4.1890e+00 | 1.2622e+00 | 4.3750e+00 | 4.6729e+01 | -0.15 | 0.40 | -0.08 |
| | 30 | 4.2284e+00 | 1.0548e+00 | 4.3580e+00 | 4.5087e+01 | -0.02 | 0.43 | 0.01 |
| | 40 | 4.2599e+00 | 9.1351e-01 | 4.3567e+00 | 5.0801e+01 | -0.03 | 0.52 | 0.00 |

The `leastSquares` reconstruction is first order in the bulk and about half
order next to the boundary, which sets its total; the `Gauss linear`
reconstruction does not converge in the bulk, the part that stalls the field.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
