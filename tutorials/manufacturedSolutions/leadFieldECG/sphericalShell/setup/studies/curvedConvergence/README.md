# curvedConvergence - leadFieldECG/sphericalShell

At what order do the lead field and the direct forward solve converge to the
exact torso potential when the current generated in a curved anisotropic
heart crosses both spheres into a torso of variable conductivity?

The rows are the static rows of the verifier: the lead field and the direct
solve applied to the exact `Vm` at the heart cell centres, so the error is
that of the passive solve and of the recording alone. The body-surface
electrodes are the ECG and carry the convergence result; the interior
electrodes are probes that show where an error comes from. Three groups of
rows:

- body surface: the eight electrodes `surface1` to `surface8` at
  `|x| = 1.47`, 0.03 below the outer sphere, as the maximum and the RMS over
  the eight of the lead-field error against the exact potential at the
  electrode cell centroid, and `phiStatic_surface`, the RMS and maximum error
  of the direct static potential over the cells of the outer boundary;
- interior probes `outerNear`, `outerSide`, `outerFar` (outer bath) and
  `pool` (bath ball): the lead-field error at the electrode cell centroid;
- `phiStatic`: the direct static potential against the exact one at every
  cell centre, over all cells, the heart, the torso and the band of cells
  next to the spheres.

The gap between the lead field and the direct solve is the largest over all
twelve electrodes. A single electrode can sit near a zero of the error and
give erratic orders between levels; the group maximum and RMS and the
surface field row do not depend on one point.

`h` is the cube root of the mean cell volume, the volume of the ball
`|x| < 1.5` divided by the number of cells. The orders in brackets are
between consecutive levels.

## Ladder

| mesh | N | cells |
|---|---|---|
| hex (cubed sphere) | 8, 12, 16, 24, 32, 40 | 2808, 9477, 22464, 75816, 179712, 351000 |
| tet (gmsh 4.15.2, seed 1) | 8, 12, 16, 20, 24, 32 | 36754, 118535, 278906, 526677, 908200, 2133385 |
| tet (gmsh 4.15.2, seeds 2 and 3) | 8, 12, 16 | see [Tetrahedra: seeds](#tetrahedra-seeds) |

Every run is one time step of the case defaults; the static rows do not
depend on it. Each run builds 12 lead fields, one per electrode. They need
29 to 50 passes of their adjoint fixed point on the hexahedra and 61 to 200
on the tetrahedra, each to the `leadFieldFixedPoint` tolerance `1e-10` of
the case. The largest `phiE` iteration count is 1032 (`tet_N32`), below the
`maxIter 5000` of `system/fvSolution`.

The seed 2 tetrahedra at `N = 12` stop in `buildLeadFields`: the fixed point
of the first electrode, `surface1`, diverges. Its `phiE` initial residual
falls to `2e-5` by pass 100, then grows by a factor of about 2.7 per 100
passes and levels near `5e-3`, and the run stops at `maxPasses 1000` with a
relative change of `1e-2`. More passes do not converge it. The same mesh
needed 496 to 620 passes with the bidomain wall conductivity `Ge`; with the
lead-field wall `Gi + Ge` the fixed point on this mesh is no longer a
contraction. The levels and seeds that ran converge as before.

All levels and seeds run in parallel on 6 ranks (`parallel`, the case's
`decomposeParDict`).

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<mesh>_N<n>`, the names that `../summarize.py` reads. From this folder:

```bash
mkdir -p results
run()
{
    name=$1; shift
    rsync -a --exclude setup --exclude regression ../../../ results/$name/
    bash results/$name/Allrun parallel "$@"
}
for n in 8 12 16 24 32 40
do
    run hex_N$n -mesh hex -N $n
done
for n in 8 12 16 20 24 32
do
    run tet_N$n -mesh tet -N $n
done
python3 ../summarize.py results
```

The tetrahedral seeds 2 and 3 go to `results/seeds/tet_N<n>_seed<s>`:

```bash
mkdir -p results/seeds
for s in 2 3
do
    for n in 8 12 16
    do
        name=seeds/tet_N${n}_seed$s
        rsync -a --exclude setup --exclude regression ../../../ results/$name/
        bash results/$name/Allrun parallel -mesh tet -N $n -seed $s
    done
done
```

`../summarize.py` reads one seed at a time from a folder of
`tet_N<n>` entries, for example links to `results/seeds/tet_N<n>_seed2`.

The omnidriver specs for this study are `sweep_hex_convergence.json`,
`sweep_tet_convergence.json` and `sweep_tet_seeds.json`, for the record
`manufacturedLeadFieldECGShell`, which omnidriver does not have. They name the
`Allrun` options as axes (`mesh`, `seed`, `parallel`; `-N` is `numberCells` on hexahedra and `tetNumberCells` on tetrahedra).

## Measured behaviour

Errors against the exact potential; the orders in brackets are between
consecutive levels. `max` and `RMS` are over the eight body-surface
electrodes; `lead field - direct` is the largest gap over all twelve
electrodes.

### Hexahedra: body surface

| N | h | max | RMS | `phiStatic_surface` RMS | `phiStatic_surface` max | lead field - direct |
|---|---|---|---|---|---|---|
| 8 | 0.1714 | 2.392e-01 | 1.213e-01 | 1.219e-01 | 2.396e-01 | 2.6e-10 |
| 12 | 0.1143 | 1.013e-01 (2.12) | 5.410e-02 (1.99) | 5.495e-02 (1.97) | 1.057e-01 (2.02) | 3.3e-10 |
| 16 | 0.0857 | 5.390e-02 (2.19) | 2.922e-02 (2.14) | 2.984e-02 (2.12) | 5.924e-02 (2.01) | 4.2e-10 |
| 24 | 0.0571 | 2.301e-02 (2.10) | 1.296e-02 (2.00) | 1.309e-02 (2.03) | 2.639e-02 (1.99) | 4.2e-10 |
| 32 | 0.0428 | 1.278e-02 (2.04) | 7.241e-03 (2.02) | 7.302e-03 (2.03) | 1.475e-02 (2.02) | 5.0e-10 |
| 40 | 0.0343 | 8.041e-03 (2.08) | 4.599e-03 (2.03) | 4.687e-03 (1.99) | 9.389e-03 (2.02) | 5.7e-10 |

Single body-surface electrodes:

| N | surface1 | surface2 | surface3 | surface4 | surface5 | surface6 | surface7 | surface8 |
|---|---|---|---|---|---|---|---|---|
| 8 | 8.191e-04 | 1.940e-01 | 3.972e-03 | 1.110e-01 | 5.856e-02 | 4.502e-02 | 2.392e-01 | 7.145e-02 |
| 12 | 8.712e-03 (-5.83) | 9.392e-02 (1.79) | 6.800e-05 (10.03) | 5.263e-02 (1.84) | 2.759e-02 (1.86) | 1.592e-02 (2.56) | 1.013e-01 (2.12) | 2.190e-02 (2.92) |
| 16 | 3.466e-03 (3.20) | 5.098e-02 (2.12) | 1.653e-03 (-11.09) | 2.903e-02 (2.07) | 1.733e-02 (1.62) | 2.706e-03 (6.16) | 5.390e-02 (2.19) | 1.274e-02 (1.88) |
| 24 | 2.253e-03 (1.06) | 2.301e-02 (1.96) | 7.299e-04 (2.02) | 1.417e-02 (1.77) | 7.437e-03 (2.09) | 1.379e-03 (1.66) | 2.295e-02 (2.11) | 4.985e-03 (2.31) |
| 32 | 1.413e-03 (1.62) | 1.278e-02 (2.04) | 5.011e-04 (1.31) | 7.991e-03 (1.99) | 4.466e-03 (1.77) | 3.232e-04 (5.04) | 1.264e-02 (2.07) | 3.205e-03 (1.54) |
| 40 | 9.615e-04 (1.73) | 8.041e-03 (2.08) | 3.278e-04 (1.90) | 5.262e-03 (1.87) | 2.906e-03 (1.93) | 7.571e-05 (6.50) | 7.976e-03 (2.06) | 1.944e-03 (2.24) |

### Hexahedra: field

| N | all RMS | all max | heart RMS | heart max | torso RMS | torso max | band RMS | band max |
|---|---|---|---|---|---|---|---|---|
| 8 | 1.214e-01 | 2.902e-01 | 1.202e-01 | 2.902e-01 | 1.218e-01 | 2.681e-01 | 1.250e-01 | 2.902e-01 |
| 12 | 5.478e-02 (1.96) | 1.244e-01 (2.09) | 5.439e-02 (1.96) | 1.244e-01 (2.09) | 5.492e-02 (1.96) | 1.144e-01 (2.10) | 5.625e-02 (1.97) | 1.244e-01 (2.09) |
| 16 | 2.976e-02 (2.12) | 6.942e-02 (2.03) | 2.956e-02 (2.12) | 6.942e-02 (2.03) | 2.983e-02 (2.12) | 6.500e-02 (1.96) | 3.058e-02 (2.12) | 6.942e-02 (2.03) |
| 24 | 1.306e-02 (2.03) | 3.093e-02 (1.99) | 1.298e-02 (2.03) | 3.093e-02 (1.99) | 1.309e-02 (2.03) | 2.941e-02 (1.96) | 1.342e-02 (2.03) | 3.093e-02 (1.99) |
| 32 | 7.280e-03 (2.03) | 1.738e-02 (2.00) | 7.234e-03 (2.03) | 1.738e-02 (2.00) | 7.297e-03 (2.03) | 1.662e-02 (1.98) | 7.496e-03 (2.02) | 1.732e-02 (2.02) |
| 40 | 4.671e-03 (1.99) | 1.112e-02 (2.00) | 4.638e-03 (1.99) | 1.112e-02 (2.00) | 4.683e-03 (1.99) | 1.069e-02 (1.98) | 4.819e-03 (1.98) | 1.107e-02 (2.01) |

### Hexahedra: interior probes

| N | h | outerNear | outerSide | outerFar | pool |
|---|---|---|---|---|---|
| 8 | 0.1714 | 2.048e-01 | 1.623e-01 | 8.302e-03 | 7.854e-02 |
| 12 | 0.1143 | 9.352e-02 (1.93) | 7.147e-02 (2.02) | 2.120e-03 (3.37) | 3.262e-02 (2.17) |
| 16 | 0.0857 | 5.379e-02 (1.92) | 3.916e-02 (2.09) | 1.643e-03 (0.89) | 1.689e-02 (2.29) |
| 24 | 0.0571 | 2.402e-02 (1.99) | 1.719e-02 (2.03) | 5.846e-04 (2.55) | 6.791e-03 (2.25) |
| 32 | 0.0428 | 1.366e-02 (1.96) | 9.557e-03 (2.04) | 1.640e-04 (4.42) | 3.585e-03 (2.22) |
| 40 | 0.0343 | 8.631e-03 (2.06) | 5.796e-03 (2.24) | 2.315e-04 (-1.55) | 2.156e-03 (2.28) |

### Tetrahedra: body surface

| N | h | max | RMS | `phiStatic_surface` RMS | `phiStatic_surface` max | lead field - direct |
|---|---|---|---|---|---|---|
| 8 | 0.0727 | 1.066e-01 | 6.122e-02 | 5.813e-02 | 1.188e-01 | 7.0e-11 |
| 12 | 0.0492 | 4.341e-02 (2.30) | 2.561e-02 (2.23) | 2.376e-02 (2.29) | 5.655e-02 (1.90) | 2.4e-11 |
| 16 | 0.0370 | 2.297e-02 (2.23) | 1.268e-02 (2.47) | 1.200e-02 (2.39) | 3.208e-02 (1.99) | 2.2e-11 |
| 20 | 0.0299 | 1.386e-02 (2.38) | 7.969e-03 (2.19) | 7.311e-03 (2.34) | 2.084e-02 (2.04) | 4.6e-11 |
| 24 | 0.0250 | 9.428e-03 (2.12) | 5.534e-03 (2.01) | 4.843e-03 (2.27) | 1.419e-02 (2.12) | 5.7e-12 |
| 32 | 0.0188 | 4.945e-03 (2.27) | 2.781e-03 (2.42) | 2.683e-03 (2.07) | 1.022e-02 (1.15) | 1.5e-11 |

Single body-surface electrodes:

| N | surface1 | surface2 | surface3 | surface4 | surface5 | surface6 | surface7 | surface8 |
|---|---|---|---|---|---|---|---|---|
| 8 | 1.770e-02 | 1.066e-01 | 2.665e-03 | 7.133e-02 | 4.668e-02 | 3.656e-03 | 1.017e-01 | 2.643e-02 |
| 12 | 9.244e-03 (1.66) | 4.341e-02 (2.30) | 1.558e-03 (1.37) | 3.141e-02 (2.10) | 1.770e-02 (2.48) | 5.891e-03 (-1.22) | 4.292e-02 (2.21) | 9.858e-03 (2.53) |
| 16 | 4.301e-03 (2.68) | 2.297e-02 (2.23) | 2.416e-03 (-1.54) | 1.464e-02 (2.67) | 1.061e-02 (1.79) | 3.415e-03 (1.91) | 1.961e-02 (2.75) | 3.262e-03 (3.88) |
| 20 | 9.187e-04 (7.28) | 1.386e-02 (2.38) | 3.355e-04 (9.32) | 9.787e-03 (1.90) | 6.647e-03 (2.21) | 2.626e-03 (1.24) | 1.256e-02 (2.10) | 3.225e-03 (0.05) |
| 24 | 3.904e-04 (4.71) | 9.428e-03 (2.12) | 7.418e-04 (-4.37) | 6.876e-03 (1.94) | 6.581e-03 (0.06) | 1.072e-03 (4.93) | 7.776e-03 (2.64) | 1.782e-03 (3.26) |
| 32 | 2.776e-04 (1.20) | 4.945e-03 (2.27) | 1.818e-04 (4.94) | 3.820e-03 (2.07) | 2.069e-03 (4.06) | 6.317e-04 (1.86) | 4.120e-03 (2.23) | 1.028e-03 (1.93) |

### Tetrahedra: field

| N | all RMS | all max | heart RMS | heart max | torso RMS | torso max | band RMS | band max |
|---|---|---|---|---|---|---|---|---|
| 8 | 5.991e-02 | 1.458e-01 | 6.095e-02 | 1.458e-01 | 5.954e-02 | 1.324e-01 | 6.266e-02 | 1.450e-01 |
| 12 | 2.444e-02 (2.30) | 6.738e-02 (1.98) | 2.458e-02 (2.33) | 6.738e-02 (1.98) | 2.439e-02 (2.29) | 6.000e-02 (2.03) | 2.551e-02 (2.30) | 6.570e-02 (2.03) |
| 16 | 1.244e-02 (2.37) | 3.471e-02 (2.33) | 1.262e-02 (2.34) | 3.471e-02 (2.33) | 1.237e-02 (2.38) | 3.208e-02 (2.19) | 1.297e-02 (2.37) | 3.222e-02 (2.50) |
| 20 | 7.577e-03 (2.34) | 2.221e-02 (2.11) | 7.692e-03 (2.34) | 2.221e-02 (2.11) | 7.536e-03 (2.34) | 2.084e-02 (2.04) | 7.949e-03 (2.31) | 2.020e-02 (2.20) |
| 24 | 5.004e-03 (2.28) | 1.452e-02 (2.34) | 5.106e-03 (2.26) | 1.452e-02 (2.34) | 4.968e-03 (2.29) | 1.419e-02 (2.12) | 5.298e-03 (2.23) | 1.373e-02 (2.13) |
| 32 | 2.746e-03 (2.11) | 1.182e-02 (0.72) | 2.808e-03 (2.10) | 8.982e-03 (1.69) | 2.724e-03 (2.11) | 1.182e-02 (0.64) | 2.951e-03 (2.06) | 1.182e-02 (0.53) |

### Tetrahedra: interior probes

| N | h | outerNear | outerSide | outerFar | pool |
|---|---|---|---|---|---|
| 8 | 0.0727 | 1.143e-01 | 7.587e-02 | 6.608e-04 | 3.934e-02 |
| 12 | 0.0492 | 4.152e-02 (2.59) | 3.018e-02 (2.36) | 6.336e-04 (0.11) | 1.404e-02 (2.64) |
| 16 | 0.0370 | 2.167e-02 (2.28) | 1.486e-02 (2.48) | 3.569e-05 (10.09) | 7.654e-03 (2.13) |
| 20 | 0.0299 | 1.374e-02 (2.15) | 9.037e-03 (2.35) | 1.557e-04 (-6.95) | 4.218e-03 (2.81) |
| 24 | 0.0250 | 8.411e-03 (2.70) | 5.795e-03 (2.45) | 1.501e-04 (0.20) | 2.507e-03 (2.86) |
| 32 | 0.0188 | 4.697e-03 (2.05) | 3.017e-03 (2.29) | 3.788e-04 (-3.25) | 1.274e-03 (2.38) |

### Tetrahedra: seeds

The same levels meshed with three gmsh seeds. Body surface and field: the
maximum and RMS over the eight surface electrodes, the RMS and maximum error
of the direct static potential over the surface cells and over all cells,
and the RMS over the interface band. The `max/min` rows are the largest over
the smallest value of the seeds that ran. The lead fields of seed 2 at
`N = 12` do not converge (see [Ladder](#ladder)).

| N | seed | cells | surface max | surface RMS | phiStatic_surface RMS | phiStatic_surface max | all RMS | all max | band RMS |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 1 | 36754 | 1.066e-01 | 6.122e-02 | 5.813e-02 | 1.188e-01 | 5.991e-02 | 1.458e-01 | 6.266e-02 |
| 8 | 2 | 36660 | 1.145e-01 | 6.546e-02 | 6.294e-02 | 1.346e-01 | 6.445e-02 | 1.621e-01 | 6.741e-02 |
| 8 | 3 | 36853 | 1.229e-01 | 6.834e-02 | 6.362e-02 | 1.332e-01 | 6.513e-02 | 1.587e-01 | 6.809e-02 |
| 8 | max/min |  | 1.15 | 1.12 | 1.09 | 1.13 | 1.09 | 1.11 | 1.09 |
| 12 | 1 | 118535 | 4.341e-02 | 2.561e-02 | 2.376e-02 | 5.655e-02 | 2.444e-02 | 6.738e-02 | 2.551e-02 |
| 12 | 2 | | not converged | | | | | | |
| 12 | 3 | 119209 | 4.555e-02 | 2.511e-02 | 2.328e-02 | 5.191e-02 | 2.392e-02 | 6.585e-02 | 2.493e-02 |
| 12 | max/min |  | 1.05 | 1.02 | 1.02 | 1.09 | 1.02 | 1.02 | 1.02 |
| 16 | 1 | 278906 | 2.297e-02 | 1.268e-02 | 1.200e-02 | 3.208e-02 | 1.244e-02 | 3.471e-02 | 1.297e-02 |
| 16 | 2 | 278490 | 2.158e-02 | 1.246e-02 | 1.152e-02 | 2.991e-02 | 1.191e-02 | 3.374e-02 | 1.250e-02 |
| 16 | 3 | 278405 | 2.281e-02 | 1.323e-02 | 1.204e-02 | 3.022e-02 | 1.243e-02 | 3.465e-02 | 1.302e-02 |
| 16 | max/min |  | 1.06 | 1.06 | 1.05 | 1.07 | 1.04 | 1.03 | 1.04 |

Interior probes, static lead-field error:

| N | seed | outerNear | outerSide | outerFar | pool |
|---|---|---|---|---|---|
| 8 | 1 | 1.143e-01 | 7.587e-02 | 6.608e-04 | 3.934e-02 |
| 8 | 2 | 1.006e-01 | 8.582e-02 | 2.268e-03 | 4.315e-02 |
| 8 | 3 | 1.069e-01 | 8.335e-02 | 9.249e-03 | 4.252e-02 |
| 8 | max/min | 1.14 | 1.13 | 14.00 | 1.10 |
| 12 | 1 | 4.152e-02 | 3.018e-02 | 6.336e-04 | 1.404e-02 |
| 12 | 2 | not converged | | | |
| 12 | 3 | 3.786e-02 | 3.045e-02 | 2.779e-03 | 1.457e-02 |
| 12 | max/min | 1.10 | 1.01 | 4.39 | 1.04 |
| 16 | 1 | 2.167e-02 | 1.486e-02 | 3.569e-05 | 7.654e-03 |
| 16 | 2 | 1.990e-02 | 1.475e-02 | 1.365e-03 | 6.891e-03 |
| 16 | 3 | 2.216e-02 | 1.474e-02 | 5.048e-04 | 7.026e-03 |
| 16 | max/min | 1.11 | 1.01 | 38.26 | 1.11 |

### Reading

- The body-surface ECG converges at second order on both meshes, in the
  maximum and in the RMS over the eight electrodes, and the error over the
  surface cells converges at second order in the RMS. Its maximum is second
  order on the hexahedra and on the tetrahedra up to `N = 24`; see the next
  point for `N = 32`. On the tetrahedra the surface group converges slightly
  above second order.
- Single surface electrodes give erratic orders between levels. On the
  hexahedra the signed error of `surface3` changes sign between levels, that
  of `surface1` grows from `N = 8` to `N = 12` before it falls, and that of
  `surface6` falls to more than an order of magnitude below the group
  maximum; on the tetrahedra the signed error of
  `surface1` changes sign between levels, `surface3` and `surface6` stay well
  below the group maximum, and `surface8` scatters between levels. The group
  maximum and RMS do not follow these points.
- The field error is second order on both meshes, in the RMS and in the
  maximum, in the heart, in the torso and in the band of cells next to the
  spheres. The band carries the same error as the rest of the field, so the
  heart-torso face conductivity adds no error of lower order.
- The exception is the maximum on the tetrahedra at `N = 32`: over all
  cells, the torso and the band it falls only from `1.45e-2` to `1.18e-2`
  (order 0.5 to 0.7), and over the surface cells at order 1.2, while every
  RMS stays second order. This gmsh mesh has a worse cell than the coarser
  levels: `checkMesh` gives a minimum cell volume of `4.4e-7`, a seventh of
  the `3.2e-6` at `N = 24` where the size scaling gives about `1.3e-6`, a
  maximum skewness of 1.00 against 0.88, and a maximum non-orthogonality of
  75.5 degrees against 69.8. The maximum is a single-cell value; the cell
  holding it is not located.
- The probes `outerNear`, `outerSide` and `pool` converge at second order on
  both meshes. At `outerFar` the signed error changes sign between the
  levels on the tetrahedra, and on the hexahedra it stays one to two orders
  of magnitude below the other probes with erratic orders: the electrode
  lies near a zero of the discretisation error, so its order is not defined.
  Its error stays below the field maximum of the same level.
- The convergence trend is seed-independent. On the tetrahedra the surface
  group, the surface field row, the field RMS, the band RMS and the errors
  at `outerNear`, `outerSide` and `pool` differ by less than a fifth between
  seeds at each level, and fall at about second order or above for every
  seed (for seed 2 from `N = 8` to `N = 16`). Only `outerFar` depends on the
  seed, by up to more than an order of magnitude, as expected for an
  electrode at a zero of the discretisation error.
- The lead field and the direct solve agree to `6e-10` at every level and
  every electrode, at the level of the fixed-point and linear-solver
  tolerances.
