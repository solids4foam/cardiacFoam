# cartesianConvergence - bidomain

Spatial convergence on Cartesian hexahedra (`Allrun`, `-mesh hex`) in 1D, 2D
and 3D, with the case defaults (`godunov` coupling, `backward`, two outer
sweeps, `leastSquares`).

## Ladder

| N | hex cells (3D) |
|---|---|
| 10 | 1000 |
| 20 | 8000 |
| 40 | 64000 |
| 80 | 512000 |

Every run uses the `Allrun` default time step,
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`: it scales as `h^2` and every
level ends exactly at `endTime 0.2`, so the slopes are space-time slopes.
Two controls on the 3D ladder: `-gradient GaussLinear` at `N = 10, 20, 40`
(parallel), and `N = 40` at half the default time step.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<d>D_N<n>`, `GaussLinear_3D_N<n>` and `3D_N40_dt<dt>`, the names that
`../summarize.py` reads (`3D_N40_dt<dt>` is a two-level temporal ladder; its
default-step level is a link to `3D_N40`). From this folder:

```bash
mkdir -p results
for n in 10 20 40 80
do
    for d in 1 2 3
    do
        rsync -a --exclude setup ../../../ results/${d}D_N$n/
        results/${d}D_N$n/Allrun -dim $d -N $n
    done
done
for n in 10 20 40
do
    rsync -a --exclude setup ../../../ results/GaussLinear_3D_N$n/
    results/GaussLinear_3D_N$n/Allrun parallel -dim 3 -N $n -gradient GaussLinear
done
set -- $(awk -v N=40 'BEGIN { h = 0.00892857*(10/N)^2; k = int(0.2/h); if (k*h < 0.2) k++; printf "%.17g %.17g", 0.2/k, 0.1/k }')
ln -s 3D_N40 results/3D_N40_dt$1
rsync -a --exclude setup ../../../ results/3D_N40_dt$2/
results/3D_N40_dt$2/Allrun -dim 3 -N 40 -dt $2
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_hex_convergence.json` (the
ladder and the `Gauss linear` control).

## Measured behaviour

Results of the last run of this study (`results/`), as printed by
`../summarize.py`; gauge-shifted potentials, orders between consecutive
levels of `h = 1/N`.

| dim | N | Vm L2 | Vm Linf | phiE L2 | phiE Linf | phiI L2 | u1 L2 | u2 L2 | order Vm L2 | order Vm Linf | order phiE L2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1D | 10 | 2.7326e-04 | 5.4861e-04 | 1.6896e-04 | 2.9385e-04 | 1.5204e-04 | 4.8571e-04 | 1.3464e-04 |  |  |  |
| 1D | 20 | 6.8869e-05 | 1.3986e-04 | 4.2047e-05 | 7.4314e-05 | 3.8997e-05 | 1.2401e-04 | 3.4423e-05 | 1.99 | 1.97 | 2.01 |
| 1D | 40 | 1.7191e-05 | 3.5002e-05 | 1.0459e-05 | 1.8554e-05 | 9.7826e-06 | 3.1078e-05 | 8.6293e-06 | 2.00 | 2.00 | 2.01 |
| 1D | 80 | 4.2997e-06 | 8.7597e-06 | 2.6136e-06 | 4.6407e-06 | 2.4495e-06 | 7.7796e-06 | 2.1603e-06 | 2.00 | 2.00 | 2.00 |
| 2D | 10 | 5.5650e-04 | 1.3069e-03 | 3.8443e-04 | 8.6064e-04 | 1.8308e-04 | 2.3732e-04 | 9.1707e-05 |  |  |  |
| 2D | 20 | 1.4155e-04 | 3.5136e-04 | 9.8808e-05 | 2.3420e-04 | 4.5644e-05 | 5.9943e-05 | 2.3186e-05 | 1.98 | 1.90 | 1.96 |
| 2D | 40 | 3.5600e-05 | 8.9493e-05 | 2.4921e-05 | 5.9848e-05 | 1.1409e-05 | 1.4958e-05 | 5.7876e-06 | 1.99 | 1.97 | 1.99 |
| 2D | 80 | 8.9102e-06 | 2.2474e-05 | 6.2412e-06 | 1.5040e-05 | 2.8520e-06 | 3.7417e-06 | 1.4478e-06 | 2.00 | 1.99 | 2.00 |
| 3D | 10 | 1.0750e-03 | 3.0313e-03 | 7.5407e-04 | 2.1051e-03 | 3.2213e-04 | 7.8365e-05 | 3.4127e-05 |  |  |  |
| 3D | 20 | 2.7526e-04 | 8.2166e-04 | 1.9417e-04 | 5.7415e-04 | 8.1378e-05 | 1.8572e-05 | 8.0989e-06 | 1.97 | 1.88 | 1.96 |
| 3D | 40 | 6.9285e-05 | 2.0965e-04 | 4.8947e-05 | 1.4674e-04 | 2.0410e-05 | 4.5270e-06 | 1.9751e-06 | 1.99 | 1.97 | 1.99 |
| 3D | 80 | 1.7348e-05 | 5.2677e-05 | 1.2261e-05 | 3.6885e-05 | 5.1053e-06 | 1.1272e-06 | 4.9183e-07 | 2.00 | 1.99 | 2.00 |

The `Vm` L2, `Vm` Linf and `phiE` L2 orders are 1.99-2.00 at the finest pair
in 1D, 2D and 3D; the `phiE` Linf, `phiI`, `u1` and `u2` errors fall about
fourfold per level.

`Gauss linear` against `leastSquares`, 3D:

| gradient | N | Vm L2 | Vm Linf | phiE L2 | phiE Linf | gradReconPhiE_boundary L2 | phiE iterations |
|---|---|---|---|---|---|---|---|
| `leastSquares` | 10 | 1.0750e-03 | 3.0313e-03 | 7.5407e-04 | 2.1051e-03 | 2.9385e-01 | 45 |
| `leastSquares` | 20 | 2.7526e-04 | 8.2166e-04 | 1.9417e-04 | 5.7415e-04 | 5.8779e-02 | 95 |
| `leastSquares` | 40 | 6.9285e-05 | 2.0965e-04 | 4.8947e-05 | 1.4674e-04 | 1.0822e-02 | 192 |
| `Gauss linear` | 10 | 1.0750e-03 | 3.0313e-03 | 7.5407e-04 | 2.1051e-03 | 3.0546e-01 | 60 |
| `Gauss linear` | 20 | 2.7526e-04 | 8.2166e-04 | 1.9417e-04 | 5.7415e-04 | 8.7224e-02 | 113 |
| `Gauss linear` | 40 | 6.9285e-05 | 2.0965e-04 | 4.8947e-05 | 1.4674e-04 | 2.7731e-02 | 223 |

On these hexahedra `Gauss linear` gives the `leastSquares` field errors to
every printed digit (`Vm`, `phiE`, `phiI`, `u1`, `u2`, L2 and Linf, at every
level); only the gradient reconstruction of the exact fields at the boundary
and the `phiE` iteration count differ.

Half the default time step, 3D `N = 40`:

| dt | Vm L2 | Vm Linf | phiE L2 | phiI L2 | u1 L2 | u2 L2 |
|---|---|---|---|---|---|---|
| 0.000557103 | 6.9285e-05 | 2.0965e-04 | 4.8947e-05 | 2.0410e-05 | 4.5270e-06 | 1.9751e-06 |
| 0.000278552 | 7.7895e-05 | 2.2195e-04 | 5.5074e-05 | 2.2824e-05 | 3.6538e-06 | 1.5571e-06 |
| change | +12.43 % | +5.87 % | +12.52 % | +11.83 % | -19.29 % | -21.16 % |

Halving the time step raises the `Vm`, `phiE` and `phiI` L2 errors by about
12 %: the Lie-Godunov splitting error is opposite in sign to the spatial
error and both scale as `h^2` along the ladder, so at the default step they
partly cancel. The `u1` and `u2` errors fall by about 20 %.

Every linear solve reaches its absolute tolerance `1e-15` below the explicit
`maxIter 5000`; the largest iteration count is 414 for `phiE`
(3D `N = 80`) and 7 for `Vm`.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
