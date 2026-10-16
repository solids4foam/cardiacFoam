# tetConvergence - bathBidomain

Spatial convergence and interface diagnostics on conformal gmsh tetrahedra
(`Allrun -mesh tet`, template `three_domain_box.geo.template`): the
ground-electrode problem with the predictor-corrector (the case defaults on
`-mesh tet`), the assembled interface current of the exact `phiE` field on the
same meshes, and the electrode-pair variant. The interfaces `x = 0` and
`x = 1` are internal faces shared by heart and bath cells, so the heart side
of the interface is a non-orthogonal reconstructed stencil.

## Ladder

| N | lc | dt | cells (gmsh 4.15.2, seed 1) |
|---|---|---|---|
| 10 | 0.1 | 0.00892857 | 14792 |
| 20 | 0.05 | 0.00224215 | 111482 |
| 40 | 0.025 | 0.000560538 | |
| 80 | 0.0125 | 0.000140174 | |

The time steps are those of
[cartesianConvergence](../cartesianConvergence/README.md) and do not divide
`endTime`; each run stops at the last step before it. The verifier labels its
summary `3D_<n>_cells.dat` with `n = round(cbrt(myocardium cells))` (17 at
`N = 10`); `../summarize.py` computes the orders from the whole mesh,
`h = (cells/3)^(-1/3)`. The electrode pair stops at `N = 40`.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<configuration>_tet_N<n>`, the names that `../summarize.py` reads. From this
folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.00892857 0.00224215 0.000560538 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    rsync -a --exclude setup ../../../ results/predictor_tet_N$n/
    results/predictor_tet_N$n/Allrun parallel -mesh tet -N $n -dt $dt
    rsync -a --exclude setup ../../../ results/exactInterface_tet_N$n/
    results/exactInterface_tet_N$n/Allrun parallel -mesh tet -N $n -dt $dt -endTime 0.02
    if [ $n -le 40 ]
    then
        rsync -a --exclude setup ../../../ results/electrodePair_tet_N$n/
        results/electrodePair_tet_N$n/Allrun parallel -mesh tet -N $n -dt $dt -endTime 0.02 -variant electrodePair
    fi
done
python3 ../summarize.py results
```

The omnidriver specs for this study are `sweep_tet_electrodePair.json`, `sweep_tet_exactInterface.json` and `sweep_tet_generic.json`.

## Measured behaviour

Results of the last run of this study (`results/`), on gmsh 4.14.0 meshes
(14838, 111654, 868416 and 6830547 cells) and with a limit of 1000 PCG
iterations for `phiE`, which the `N = 80` solves reached (`system/fvSolution`
sets `maxIter 5000`).

Ground electrode with the predictor-corrector, field errors (L2):

| N | h | Vm | order | phiE | order | phiI | order |
|---|---|---|---|---|---|---|---|
| 10 | 0.0587 | 1.1106e-03 | | 4.9622e-03 | | 3.5857e-03 | |
| 20 | 0.0300 | 3.7398e-04 | 1.62 | 1.4852e-03 | 1.79 | 1.0834e-03 | 1.78 |
| 40 | 0.0151 | 1.0375e-04 | 1.88 | 4.4474e-04 | 1.76 | 3.2855e-04 | 1.74 |
| 80 | 0.0076 | 2.7316e-05 | 1.94 | 1.5344e-04 | 1.55 | 1.1463e-04 | 1.53 |

Interface rows of the same runs (L2; potentials volume weighted, current
densities face-area weighted, in A/m^2 against the manufactured current
`alpha = 0.01`):

| N | heart phiE | bath phiE | x=0 current jump | x=1 current jump | x=0 leakage | x=1 leakage | x=0 assembled | x=1 assembled |
|---|---|---|---|---|---|---|---|---|
| 10 | 3.8238e-03 | 5.4431e-03 | 4.0366e-03 | 3.9818e-03 | 3.9685e-03 | 3.8765e-03 | 5.9904e-04 | 5.9373e-04 |
| 20 | 1.1528e-03 | 1.6260e-03 | 6.0793e-04 | 6.2093e-04 | 5.8971e-04 | 6.0458e-04 | 3.8649e-04 | 3.9304e-04 |
| 40 | 3.4651e-04 | 4.8643e-04 | 1.5267e-04 | 1.5793e-04 | 1.5064e-04 | 1.5509e-04 | 3.5142e-04 | 3.9299e-04 |
| 80 | 1.1887e-04 | 1.6792e-04 | 9.5046e-05 | 9.9408e-05 | 5.7130e-05 | 5.6690e-05 | 4.1201e-04 | 4.2078e-04 |

Orders between consecutive levels:

| row | 10-20 | 20-40 | 40-80 |
|---|---|---|---|
| heart phiE | 1.78 | 1.76 | 1.56 |
| bath phiE | 1.80 | 1.76 | 1.55 |
| x=0 / x=1 current jump | 2.81 / 2.76 | 2.02 / 2.00 | 0.69 / 0.67 |
| x=0 / x=1 leakage | 2.83 / 2.76 | 2.00 / 1.99 | 1.41 / 1.46 |
| x=0 / x=1 assembled | 0.65 / 0.61 | 0.14 / 0.00 | -0.23 / -0.10 |

The regional potentials converge at 1.8 to 1.55 and the current jump and
leakage at about 2 up to `N = 40`, then flatten. The assembled current density
of the solved field stops falling after `N = 20` and rises slightly at
`N = 80`: the shared-face flux is conservative, but its constitutive error
does not converge on this family. These runs evaluated the assembled current
with the `unweightedHarmonic` interface conductivity of the utility while the
solver used `distanceWeightedHarmonic`; `Allrun` gives both the method of
`-interface`.

Assembled interface current density of the exact `phiE` field
(`-endTime 0.02`, gmsh 4.15.2 meshes):

| N | final time | x=0 | x=1 | order (x=0, x=1) |
|---|---|---|---|---|
| 10 | 0.0178571 | 1.6071e-03 | 1.6404e-03 | |
| 20 | 0.0201793 | 8.8496e-04 | 8.6846e-04 | 0.88, 0.94 |
| 40 | 0.0201794 | 4.4223e-04 | 4.5193e-04 | 1.01, 0.95 |
| 80 | 0.0200392 | 2.5222e-04 | 2.5152e-04 | 0.82, 0.85 |

(orders of `h = (cells/3)^(-1/3)`; the `N = 80` run used `dt = 0.0001401345`.)
The interface flux construction is first order on the exact field; the solved
assembled current above lies below it up to `N = 40` and above it at
`N = 80`.

Electrode pair (`-endTime 0.02`), L2 errors after removing the volume means
of `phiE` and `phiI`, and interface rows:

| N | Vm | phiE | phiI | heart phiE | bath phiE | x=0 leakage | x=1 leakage | x=0 assembled | x=1 assembled |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 1.9133e-04 | 2.6915e-03 | 1.4039e-04 | 1.9858e-04 | 3.2933e-03 | 3.3822e-03 | 3.3243e-03 | 4.4416e-04 | 4.6743e-04 |
| 20 | 1.5430e-04 | 8.7756e-04 | 5.3710e-05 | 1.1348e-04 | 1.0718e-03 | 5.3634e-04 | 5.4046e-04 | 2.8141e-04 | 3.0035e-04 |
| 40 | 6.9834e-05 | 2.7428e-04 | 2.1288e-05 | 5.1132e-05 | 3.3399e-04 | 1.3524e-04 | 1.3487e-04 | 2.2833e-04 | 2.3953e-04 |

The interface rows follow the ground-electrode pattern: the leakage falls at
about 2.8 and 2.0, the assembled current flattens, and the two interfaces
agree to two or three digits.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
