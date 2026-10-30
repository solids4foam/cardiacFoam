# curvedConvergence - leadFieldECG/sphericalShell

At what order do the lead field and the direct forward solve converge to the
exact torso potential when the current generated in a curved anisotropic
heart crosses both spheres into a torso of variable conductivity?

The rows are the static rows of the verifier: the lead field and the direct
solve applied to the exact `Vm` at the heart cell centres, so the error is
that of the passive solve and of the recording alone. Two groups of rows:

- electrodes `outerFar`, `outerNear`, `outerSide` (outer bath) and `pool`
  (bath ball): the lead-field potential against the exact potential at the
  electrode cell centroid, and the gap between the lead field and the direct
  solve;
- `phiStatic`: the direct static potential against the exact one at every
  cell centre, over all cells, the heart, the torso and the band of cells
  next to the spheres.

`h` is the cube root of the mean cell volume, the volume of the ball
`|x| < 1.5` divided by the number of cells. The orders in brackets are
between consecutive levels.

## Ladder

| mesh | N | cells |
|---|---|---|
| hex (cubed sphere) | 8, 12, 16, 24, 32 | 2808, 9477, 22464, 75816, 179712 |
| tet (gmsh 4.15.2, seed 1) | 8, 12, 16, 20 | 36754, 118535, 278906, 526677 |

Every run is one time step of the case defaults; the static rows do not
depend on it. The lead fields need 31 to 49 passes of their adjoint fixed
point on the hexahedra and 71 to 203 on the tetrahedra, each to the
`leadFieldFixedPoint` tolerance `1e-10` of the case.

All levels were run in parallel on 6 ranks (`parallel`, the case's
`decomposeParDict`).

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<mesh>_N<n>`, the names that `../summarize.py` reads. From this folder:

```bash
mkdir -p results
run()
{
    name=$1; shift
    rsync -a --exclude setup ../../../ results/$name/
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

## Measured behaviour

### Hexahedra: electrodes

| N | h | outerFar | outerNear | outerSide | pool | lead field - direct |
|---|---|---|---|---|---|---|
| 8 | 0.1714 | 6.959e-03 | 2.137e-01 | 1.679e-01 | 8.319e-02 | 3.0e-10 |
| 12 | 0.1143 | 9.105e-04 (5.02) | 9.645e-02 (1.96) | 7.550e-02 (1.97) | 3.546e-02 (2.10) | 2.8e-10 |
| 16 | 0.0857 | 1.116e-03 (-0.71) | 5.531e-02 (1.93) | 4.148e-02 (2.08) | 1.853e-02 (2.26) | 3.8e-10 |
| 24 | 0.0571 | 2.293e-04 (3.90) | 2.480e-02 (1.98) | 1.833e-02 (2.01) | 7.647e-03 (2.18) | 4.3e-10 |
| 32 | 0.0428 | 4.979e-05 (5.31) | 1.409e-02 (1.97) | 1.015e-02 (2.05) | 4.077e-03 (2.19) | 4.8e-10 |

### Hexahedra: field

| N | all RMS | all max | heart RMS | heart max | torso RMS | torso max | band RMS | band max |
|---|---|---|---|---|---|---|---|---|
| 8 | 1.257e-01 | 2.958e-01 | 1.248e-01 | 2.958e-01 | 1.260e-01 | 2.781e-01 | 1.295e-01 | 2.958e-01 |
| 12 | 5.728e-02 (1.94) | 1.287e-01 (2.05) | 5.679e-02 (1.94) | 1.287e-01 (2.05) | 5.745e-02 (1.94) | 1.200e-01 (2.07) | 5.875e-02 (1.95) | 1.287e-01 (2.05) |
| 16 | 3.107e-02 (2.13) | 7.161e-02 (2.04) | 3.078e-02 (2.13) | 7.161e-02 (2.04) | 3.117e-02 (2.12) | 6.710e-02 (2.02) | 3.189e-02 (2.12) | 7.161e-02 (2.04) |
| 24 | 1.372e-02 (2.02) | 3.206e-02 (1.98) | 1.357e-02 (2.02) | 3.206e-02 (1.98) | 1.377e-02 (2.02) | 3.047e-02 (1.95) | 1.406e-02 (2.02) | 3.206e-02 (1.98) |
| 32 | 7.627e-03 (2.04) | 1.796e-02 (2.01) | 7.532e-03 (2.05) | 1.796e-02 (2.01) | 7.660e-03 (2.04) | 1.717e-02 (1.99) | 7.828e-03 (2.04) | 1.790e-02 (2.03) |

### Tetrahedra: electrodes

| N | h | outerFar | outerNear | outerSide | pool | lead field - direct |
|---|---|---|---|---|---|---|
| 8 | 0.0727 | 1.497e-03 | 1.178e-01 | 8.144e-02 | 4.435e-02 | 3.0e-11 |
| 12 | 0.0492 | 2.272e-04 (4.83) | 4.399e-02 (2.52) | 3.317e-02 (2.30) | 1.674e-02 (2.50) | 1.4e-11 |
| 16 | 0.0370 | 7.157e-04 (-4.02) | 2.326e-02 (2.23) | 1.693e-02 (2.36) | 9.308e-03 (2.06) | 2.7e-11 |
| 20 | 0.0299 | 5.239e-04 (1.47) | 1.469e-02 (2.17) | 1.040e-02 (2.30) | 5.302e-03 (2.66) | 9.5e-12 |

### Tetrahedra: field

| N | all RMS | all max | heart RMS | heart max | torso RMS | torso max | band RMS | band max |
|---|---|---|---|---|---|---|---|---|
| 8 | 6.334e-02 | 1.525e-01 | 6.476e-02 | 1.525e-01 | 6.283e-02 | 1.389e-01 | 6.622e-02 | 1.525e-01 |
| 12 | 2.628e-02 (2.25) | 6.934e-02 (2.02) | 2.653e-02 (2.29) | 6.934e-02 (2.02) | 2.620e-02 (2.24) | 6.315e-02 (2.02) | 2.735e-02 (2.27) | 6.863e-02 (2.05) |
| 16 | 1.372e-02 (2.28) | 3.643e-02 (2.26) | 1.389e-02 (2.27) | 3.643e-02 (2.26) | 1.366e-02 (2.28) | 3.421e-02 (2.15) | 1.420e-02 (2.30) | 3.404e-02 (2.46) |
| 20 | 8.357e-03 (2.34) | 2.337e-02 (2.09) | 8.440e-03 (2.35) | 2.309e-02 (2.15) | 8.327e-03 (2.34) | 2.337e-02 (1.80) | 8.688e-03 (2.32) | 2.337e-02 (1.77) |

### Reading

- The field error is second order on both meshes, in the RMS and in the
  maximum, in the heart, in the torso and in the band of cells next to the
  spheres. The band carries the same error as the rest of the field, so the
  heart-torso face conductivity adds no error of lower order.
- The electrodes `outerNear`, `outerSide` and `pool` converge at second
  order on both meshes.
- At `outerFar` the signed error changes sign between the levels (hexahedra:
  `-9.1e-4` at N = 12, `-2.3e-4` at N = 24, `+5.0e-5` at N = 32): the
  electrode lies near a zero of the discretisation error, so its order is not
  defined. Its error stays below the field maximum of the same level.
- The lead field and the direct solve agree to `5e-10` at every level, at
  the level of the fixed-point and linear-solver tolerances.
