# tetTemporalControl - bidomain

The [tetConvergence](../tetConvergence/README.md) ladder refines the time
step with the mesh (`dt ~ h^2`) under first-order `godunov` coupling, so its
slopes are space-time slopes. This control keeps the `N = 40` mesh and halves
the time step; a change much smaller than the change between successive
levels shows that the tetrahedral errors are spatial.

## Ladder

| N | dt | dt/2 |
|---|---|---|
| 40 | `Allrun` default | half the default |

for both gradient reconstructions, with the other settings of
tetConvergence. The `Allrun` default time step is
`endTime/ceil(endTime/(0.00892857 (10/N)^2))`, which ends exactly at
`endTime 0.2`; half of it ends there too.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<gradient>_tet_N40_dt<dt>`, the names that `../summarize.py` reads; the
default-step level is a link to the tetConvergence run. From this folder,
after [tetConvergence](../tetConvergence/README.md):

```bash
mkdir -p results
set -- $(awk -v N=40 'BEGIN { h = 0.00892857*(10/N)^2; k = int(0.2/h); if (k*h < 0.2) k++; printf "%.17g %.17g", 0.2/k, 0.1/k }')
for g in leastSquares GaussLinear
do
    ln -s ../../tetConvergence/results/${g}_tet_N40 results/${g}_tet_N40_dt$1
    r=results/${g}_tet_N40_dt$2
    rsync -a --exclude setup ../../../ $r/
    $r/Allrun parallel -mesh tet -gradient $g -N 40 -dt $2
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_dt_half.json`.

## Measured behaviour

Results of the last run of this study (`results/`), as printed by
`../summarize.py`; gauge-shifted potentials.

| gradient | dt | Vm L2 | Vm Linf | phiE L2 | phiE Linf | phiI L2 |
|---|---|---|---|---|---|---|
| `leastSquares` | 0.000557103 | 6.3151e-04 | 5.7885e-03 | 4.4747e-04 | 4.1461e-03 | 1.8553e-04 |
| `leastSquares` | 0.000278552 | 6.3130e-04 | 5.7677e-03 | 4.4656e-04 | 4.1048e-03 | 1.8511e-04 |
| `leastSquares` | change | -0.03 % | -0.36 % | -0.20 % | -1.00 % | -0.23 % |
| `Gauss linear` | 0.000557103 | 5.0344e-03 | 4.5213e-02 | 3.5670e-03 | 3.1940e-02 | 1.4680e-03 |
| `Gauss linear` | 0.000278552 | 5.0281e-03 | 4.5211e-02 | 3.5589e-03 | 3.1938e-02 | 1.4697e-03 |
| `Gauss linear` | change | -0.13 % | -0.00 % | -0.23 % | -0.01 % | +0.12 % |

Change of the L2 errors from `N = 20` to `N = 40` in
[tetConvergence](../tetConvergence/README.md), for comparison:

| gradient | Vm L2 | phiE L2 |
|---|---|---|
| `leastSquares` | -73.82 % | -73.83 % |
| `Gauss linear` | -37.48 % | -37.54 % |

Halving the time step changes the L2 errors by at most 0.23 %, against 74 %
(`leastSquares`) and 37 % (`Gauss linear`) between the two finest levels: the
tetrahedral slopes are spatial.

Every linear solve reaches its absolute tolerance `1e-15` below the explicit
`maxIter 5000`; the largest iteration count is 494 for `phiE` and
13 for `Vm`.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
