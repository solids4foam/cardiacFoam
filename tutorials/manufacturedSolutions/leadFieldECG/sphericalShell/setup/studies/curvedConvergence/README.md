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

- body surface: the eight electrodes `surface1` to `surface8` on the outer
  sphere at `|x| = 1.495`, as the maximum and the RMS over the eight of the
  lead-field error against the exact potential at the electrode cell
  centroid, and `phiStatic_surface`, the RMS and maximum error of the direct
  static potential over the cells of the outer boundary;
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
| hex (cubed sphere) | 8, 12, 16, 24, 32 | 2808, 9477, 22464, 75816, 179712 |
| tet (gmsh 4.15.2, seed 1) | 8, 12, 16, 20 | 36754, 118535, 278906, 526677 |
| tet (gmsh 4.15.2, seeds 2 and 3) | 8, 12, 16 | see [Tetrahedra: seeds](#tetrahedra-seeds) |

Every run is one time step of the case defaults; the static rows do not
depend on it. Each run builds 12 lead fields, one per electrode. They need
29 to 49 passes of their adjoint fixed point on the hexahedra and 61 to 208
on the tetrahedra, each to the `leadFieldFixedPoint` tolerance `1e-10` of
the case. The exception is the seed 2 tetrahedra at `N = 12`, which need
482 to 620 passes, within the `maxPasses 1000` of `system/fvSolution`.

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
for n in 8 12 16 24 32
do
    run hex_N$n -mesh hex -N $n
done
for n in 8 12 16 20
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

## Measured behaviour

Errors against the exact potential; the orders in brackets are between
consecutive levels. `max` and `RMS` are over the eight body-surface
electrodes; `lead field - direct` is the largest gap over all twelve
electrodes.

### Hexahedra: body surface

| N | h | max | RMS | `phiStatic_surface` RMS | `phiStatic_surface` max | lead field - direct |
|---|---|---|---|---|---|---|
| 8 | 0.1714 | 2.473e-01 | 1.253e-01 | 1.259e-01 | 2.477e-01 | 2.9e-10 |
| 12 | 0.1143 | 1.061e-01 (2.09) | 5.660e-02 (1.96) | 5.742e-02 (1.94) | 1.097e-01 (2.01) | 3.7e-10 |
| 16 | 0.0857 | 5.675e-02 (2.18) | 3.059e-02 (2.14) | 3.115e-02 (2.13) | 6.121e-02 (2.03) | 4.6e-10 |
| 24 | 0.0571 | 2.441e-02 (2.08) | 1.366e-02 (1.99) | 1.377e-02 (2.01) | 2.737e-02 (1.98) | 4.8e-10 |
| 32 | 0.0428 | 1.343e-02 (2.08) | 7.611e-03 (2.03) | 7.661e-03 (2.04) | 1.527e-02 (2.03) | 5.5e-10 |

Single body-surface electrodes:

| N | surface1 | surface2 | surface3 | surface4 | surface5 | surface6 | surface7 | surface8 |
|---|---|---|---|---|---|---|---|---|
| 8 | 1.530e-03 | 1.994e-01 | 7.401e-03 | 1.158e-01 | 5.857e-02 | 4.695e-02 | 2.473e-01 | 7.442e-02 |
| 12 | 7.286e-03 (-3.85) | 9.764e-02 (1.76) | 3.558e-04 (7.48) | 5.477e-02 (1.85) | 2.919e-02 (1.72) | 1.790e-02 (2.38) | 1.061e-01 (2.09) | 2.451e-02 (2.74) |
| 16 | 2.524e-03 (3.69) | 5.298e-02 (2.12) | 1.747e-03 (-5.53) | 3.015e-02 (2.08) | 1.815e-02 (1.65) | 3.686e-03 (5.49) | 5.675e-02 (2.18) | 1.406e-02 (1.93) |
| 24 | 1.780e-03 (0.86) | 2.402e-02 (1.95) | 5.859e-04 (2.69) | 1.483e-02 (1.75) | 7.791e-03 (2.09) | 1.909e-03 (1.62) | 2.441e-02 (2.08) | 5.617e-03 (2.26) |
| 32 | 1.127e-03 (1.59) | 1.332e-02 (2.05) | 4.414e-04 (0.98) | 8.355e-03 (2.00) | 4.631e-03 (1.81) | 6.100e-04 (3.97) | 1.343e-02 (2.08) | 3.517e-03 (1.63) |

### Hexahedra: field

| N | all RMS | all max | heart RMS | heart max | torso RMS | torso max | band RMS | band max |
|---|---|---|---|---|---|---|---|---|
| 8 | 1.257e-01 | 2.958e-01 | 1.248e-01 | 2.958e-01 | 1.260e-01 | 2.781e-01 | 1.295e-01 | 2.958e-01 |
| 12 | 5.728e-02 (1.94) | 1.287e-01 (2.05) | 5.679e-02 (1.94) | 1.287e-01 (2.05) | 5.745e-02 (1.94) | 1.200e-01 (2.07) | 5.875e-02 (1.95) | 1.287e-01 (2.05) |
| 16 | 3.107e-02 (2.13) | 7.161e-02 (2.04) | 3.078e-02 (2.13) | 7.161e-02 (2.04) | 3.117e-02 (2.12) | 6.710e-02 (2.02) | 3.189e-02 (2.12) | 7.161e-02 (2.04) |
| 24 | 1.372e-02 (2.02) | 3.206e-02 (1.98) | 1.357e-02 (2.02) | 3.206e-02 (1.98) | 1.377e-02 (2.02) | 3.047e-02 (1.95) | 1.406e-02 (2.02) | 3.206e-02 (1.98) |
| 32 | 7.627e-03 (2.04) | 1.796e-02 (2.01) | 7.532e-03 (2.05) | 1.796e-02 (2.01) | 7.660e-03 (2.04) | 1.717e-02 (1.99) | 7.828e-03 (2.04) | 1.790e-02 (2.03) |

### Hexahedra: interior probes

| N | h | outerNear | outerSide | outerFar | pool |
|---|---|---|---|---|---|
| 8 | 0.1714 | 2.137e-01 | 1.679e-01 | 6.959e-03 | 8.319e-02 |
| 12 | 0.1143 | 9.645e-02 (1.96) | 7.550e-02 (1.97) | 9.105e-04 (5.02) | 3.546e-02 (2.10) |
| 16 | 0.0857 | 5.531e-02 (1.93) | 4.148e-02 (2.08) | 1.116e-03 (-0.71) | 1.853e-02 (2.26) |
| 24 | 0.0571 | 2.480e-02 (1.98) | 1.833e-02 (2.01) | 2.293e-04 (3.90) | 7.647e-03 (2.18) |
| 32 | 0.0428 | 1.409e-02 (1.97) | 1.015e-02 (2.05) | 4.979e-05 (5.31) | 4.077e-03 (2.19) |

### Tetrahedra: body surface

| N | h | max | RMS | `phiStatic_surface` RMS | `phiStatic_surface` max | lead field - direct |
|---|---|---|---|---|---|---|
| 8 | 0.0727 | 1.122e-01 | 6.467e-02 | 6.126e-02 | 1.237e-01 | 6.8e-11 |
| 12 | 0.0492 | 4.623e-02 (2.27) | 2.760e-02 (2.18) | 2.549e-02 (2.25) | 5.924e-02 (1.89) | 2.0e-11 |
| 16 | 0.0370 | 2.444e-02 (2.23) | 1.403e-02 (2.37) | 1.326e-02 (2.29) | 3.421e-02 (1.93) | 3.0e-11 |
| 20 | 0.0299 | 1.467e-02 (2.41) | 8.855e-03 (2.17) | 8.082e-03 (2.34) | 2.203e-02 (2.08) | 5.6e-11 |

Single body-surface electrodes:

| N | surface1 | surface2 | surface3 | surface4 | surface5 | surface6 | surface7 | surface8 |
|---|---|---|---|---|---|---|---|---|
| 8 | 1.966e-02 | 1.122e-01 | 5.843e-03 | 7.586e-02 | 4.889e-02 | 4.396e-03 | 1.077e-01 | 2.609e-02 |
| 12 | 8.099e-03 (2.27) | 4.623e-02 (2.27) | 3.421e-03 (1.37) | 3.285e-02 (2.14) | 2.390e-02 (1.83) | 4.382e-03 (0.01) | 4.548e-02 (2.21) | 1.193e-02 (2.01) |
| 16 | 5.231e-03 (1.53) | 2.444e-02 (2.23) | 3.179e-03 (0.26) | 1.693e-02 (2.32) | 1.135e-02 (2.61) | 2.066e-03 (2.64) | 2.229e-02 (2.50) | 4.797e-03 (3.19) |
| 20 | 3.497e-03 (1.90) | 1.467e-02 (2.41) | 9.403e-04 (5.75) | 1.056e-02 (2.23) | 8.914e-03 (1.14) | 8.888e-04 (3.98) | 1.379e-02 (2.27) | 4.136e-03 (0.70) |

### Tetrahedra: field

| N | all RMS | all max | heart RMS | heart max | torso RMS | torso max | band RMS | band max |
|---|---|---|---|---|---|---|---|---|
| 8 | 6.334e-02 | 1.525e-01 | 6.476e-02 | 1.525e-01 | 6.283e-02 | 1.389e-01 | 6.622e-02 | 1.525e-01 |
| 12 | 2.628e-02 (2.25) | 6.934e-02 (2.02) | 2.653e-02 (2.29) | 6.934e-02 (2.02) | 2.620e-02 (2.24) | 6.315e-02 (2.02) | 2.735e-02 (2.27) | 6.863e-02 (2.05) |
| 16 | 1.372e-02 (2.28) | 3.643e-02 (2.26) | 1.389e-02 (2.27) | 3.643e-02 (2.26) | 1.366e-02 (2.28) | 3.421e-02 (2.15) | 1.420e-02 (2.30) | 3.404e-02 (2.46) |
| 20 | 8.357e-03 (2.34) | 2.337e-02 (2.09) | 8.440e-03 (2.35) | 2.309e-02 (2.15) | 8.327e-03 (2.34) | 2.337e-02 (1.80) | 8.688e-03 (2.32) | 2.337e-02 (1.77) |

### Tetrahedra: interior probes

| N | h | outerNear | outerSide | outerFar | pool |
|---|---|---|---|---|---|
| 8 | 0.0727 | 1.178e-01 | 8.144e-02 | 1.497e-03 | 4.435e-02 |
| 12 | 0.0492 | 4.399e-02 (2.52) | 3.317e-02 (2.30) | 2.272e-04 (4.83) | 1.674e-02 (2.50) |
| 16 | 0.0370 | 2.326e-02 (2.23) | 1.693e-02 (2.36) | 7.157e-04 (-4.02) | 9.308e-03 (2.06) |
| 20 | 0.0299 | 1.469e-02 (2.17) | 1.040e-02 (2.30) | 5.239e-04 (1.47) | 5.302e-03 (2.66) |

### Tetrahedra: seeds

The same levels meshed with three gmsh seeds. Body surface and field: the
maximum and RMS over the eight surface electrodes, the RMS and maximum error
of the direct static potential over the surface cells and over all cells,
and the RMS over the interface band. The `max/min` rows are the largest over
the smallest value of the three seeds.

| N | seed | cells | surface max | surface RMS | phiStatic_surface RMS | phiStatic_surface max | all RMS | all max | band RMS |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 1 | 36754 | 1.122e-01 | 6.467e-02 | 6.126e-02 | 1.237e-01 | 6.334e-02 | 1.525e-01 | 6.622e-02 |
| 8 | 2 | 36660 | 1.191e-01 | 6.756e-02 | 6.510e-02 | 1.380e-01 | 6.677e-02 | 1.665e-01 | 6.984e-02 |
| 8 | 3 | 36853 | 1.238e-01 | 6.962e-02 | 6.572e-02 | 1.367e-01 | 6.728e-02 | 1.619e-01 | 7.028e-02 |
| 8 | max/min |  | 1.10 | 1.08 | 1.07 | 1.12 | 1.06 | 1.09 | 1.06 |
| 12 | 1 | 118535 | 4.623e-02 | 2.760e-02 | 2.549e-02 | 5.924e-02 | 2.628e-02 | 6.934e-02 | 2.735e-02 |
| 12 | 2 | 118663 | 4.840e-02 | 2.776e-02 | 2.646e-02 | 5.823e-02 | 2.734e-02 | 6.986e-02 | 2.845e-02 |
| 12 | 3 | 119209 | 4.818e-02 | 2.704e-02 | 2.520e-02 | 5.470e-02 | 2.593e-02 | 6.902e-02 | 2.692e-02 |
| 12 | max/min |  | 1.05 | 1.03 | 1.05 | 1.08 | 1.05 | 1.01 | 1.06 |
| 16 | 1 | 278906 | 2.444e-02 | 1.403e-02 | 1.326e-02 | 3.421e-02 | 1.372e-02 | 3.643e-02 | 1.420e-02 |
| 16 | 2 | 278490 | 2.272e-02 | 1.346e-02 | 1.257e-02 | 3.154e-02 | 1.298e-02 | 3.530e-02 | 1.353e-02 |
| 16 | 3 | 278405 | 2.455e-02 | 1.449e-02 | 1.314e-02 | 3.160e-02 | 1.357e-02 | 4.202e-02 | 1.412e-02 |
| 16 | max/min |  | 1.08 | 1.08 | 1.05 | 1.08 | 1.06 | 1.19 | 1.05 |

Interior probes, static lead-field error:

| N | seed | outerNear | outerSide | outerFar | pool |
|---|---|---|---|---|---|
| 8 | 1 | 1.178e-01 | 8.144e-02 | 1.497e-03 | 4.435e-02 |
| 8 | 2 | 1.028e-01 | 8.999e-02 | 1.350e-03 | 4.717e-02 |
| 8 | 3 | 1.090e-01 | 8.744e-02 | 1.002e-02 | 4.542e-02 |
| 8 | max/min | 1.15 | 1.11 | 7.42 | 1.06 |
| 12 | 1 | 4.399e-02 | 3.317e-02 | 2.272e-04 | 1.674e-02 |
| 12 | 2 | 4.507e-02 | 3.443e-02 | 6.011e-04 | 1.797e-02 |
| 12 | 3 | 4.088e-02 | 3.365e-02 | 1.716e-03 | 1.728e-02 |
| 12 | max/min | 1.10 | 1.04 | 7.55 | 1.07 |
| 16 | 1 | 2.326e-02 | 1.693e-02 | 7.157e-04 | 9.308e-03 |
| 16 | 2 | 2.115e-02 | 1.670e-02 | 8.796e-04 | 8.502e-03 |
| 16 | 3 | 2.343e-02 | 1.673e-02 | 2.316e-05 | 8.810e-03 |
| 16 | max/min | 1.11 | 1.01 | 37.99 | 1.09 |

### Reading

- The body-surface ECG converges at second order on both meshes, in the
  maximum and in the RMS over the eight electrodes, and the error over the
  surface cells converges at second order in the RMS and the maximum. On the
  tetrahedra the surface group converges slightly above second order.
- Single surface electrodes give erratic orders between levels. On the
  hexahedra the signed errors of `surface1` and `surface3` change sign
  between levels, and that of `surface6` falls to more than an order of
  magnitude below the group maximum; on the tetrahedra `surface1` changes
  sign, `surface3` and `surface6` stay well below the group maximum, and
  `surface5` and `surface8` scatter between levels. The group maximum and
  RMS do not follow these points.
- The field error is second order on both meshes, in the RMS and in the
  maximum, in the heart, in the torso and in the band of cells next to the
  spheres. The band carries the same error as the rest of the field, so the
  heart-torso face conductivity adds no error of lower order.
- The probes `outerNear`, `outerSide` and `pool` converge at second order on
  both meshes. At `outerFar` the signed error changes sign between the
  levels on the hexahedra: the electrode lies near a zero of the
  discretisation error, so its order is not defined. Its error stays below
  the field maximum of the same level.
- The convergence trend is seed-independent. On the tetrahedra the surface
  group, the surface field row, the field RMS, the band RMS and the errors
  at `outerNear`, `outerSide` and `pool` differ by less than a fifth between
  seeds at each level, and fall at about second order or above for every
  seed. The field maximum over all cells, a single-cell value, spreads more
  at the finest level. Only `outerFar` depends on the seed, by up to more
  than an order of magnitude, as expected for an electrode at a zero of the
  discretisation error.
- The lead field and the direct solve agree to `6e-10` at every level and
  every electrode, at the level of the fixed-point and linear-solver
  tolerances.
