# interfaceCurrentConvergence - bathBidomain

Interface diagnostics on the conformal tetrahedral meshes of
[tetConvergence](../tetConvergence/README.md) for the two heart-bath
interface options of the global `phiE` equation
(`interfaceConductivityInterpolation`, `Allrun -interface`). Both take the
flux-continuous face conductivity on heart-bath faces; they differ in the
cell gradient of the explicit corrections: the standard one
(`distanceWeightedHarmonic`, the case default) or the jump-aware one of
`conormalHarmonic`, with the `nInterfaceGradientCorrectors 2` passes of the
case `fvSolution`. The assembled current density of
`bathBidomainInterfaceMetrics` is the flux of the `distanceWeightedHarmonic`
operator in both runs.

## Ladder

| N | dt | endTime |
|---|---|---|
| 10 | 0.00892857 | 0.02 |
| 20 | 0.00224215 | 0.02 |
| 40 | 0.000560538 | 0.02 |
| 80 | 0.000140174 | 0.02 |

for both methods. The steps do not divide `endTime`; each run stops at the
last step before it.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<method>_tet_N<n>`, the names that `../summarize.py` reads. From this
folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.00892857 0.00224215 0.000560538 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    for m in distanceWeightedHarmonic conormalHarmonic
    do
        rsync -a --exclude setup ../../../ results/${m}_tet_N$n/
        results/${m}_tet_N$n/Allrun parallel -mesh tet -N $n -dt $dt -endTime 0.02 -interface $m
    done
done
python3 ../summarize.py results
python3 ../summarize.py results --relative-to distanceWeightedHarmonic_tet
```

The omnidriver specs for this study are `sweep_tet_distanceWeightedHarmonic.json` and `sweep_tet_conormalHarmonic.json`.

## Measured quantities

For each method and level `../summarize.py` reports the L2 and Linf errors of
`Vm`, `phiE` and `phiI` and the interface rows (L2): heart and bath `phiE`,
current jump, intracellular leakage and assembled current density at `x = 0`
and `x = 1`, with the orders between levels. With
`--relative-to distanceWeightedHarmonic_tet` it reports the relative change
of `conormalHarmonic` against `distanceWeightedHarmonic` at each level.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
