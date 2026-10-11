# interfaceCurrentConvergence - bathBidomain

Interface diagnostics on the conformal tetrahedral meshes of
[tetConvergence](../tetConvergence/README.md) for the two interface
conductivities of the global `phiE` equation
(`interfaceConductivityInterpolation`, `Allrun -interface`): the harmonic mean
of the two cell conductivities (`unweightedHarmonic`) and the harmonic mean
weighted by the face-to-centre distances (`distanceWeightedHarmonic`, the
case default). The solver and `bathBidomainInterfaceMetrics` use the same
method.

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
    for m in distanceWeightedHarmonic unweightedHarmonic
    do
        rsync -a --exclude setup ../../../ results/${m}_tet_N$n/
        results/${m}_tet_N$n/Allrun parallel -mesh tet -N $n -dt $dt -endTime 0.02 -interface $m
    done
done
python3 ../summarize.py results
python3 ../summarize.py results --relative-to distanceWeightedHarmonic_tet
```

## Measured behaviour

Results of the last run of this study (`results/`, `N = 10` only, gmsh
4.15.2, 14792 cells), L2:

| method | Vm | phiE | heart phiE | bath phiE | x=0 current jump | x=0 leakage | x=0 assembled | x=1 assembled |
|---|---|---|---|---|---|---|---|---|
| distanceWeightedHarmonic | 1.9567e-04 | 3.7496e-03 | 2.8942e-03 | 4.1110e-03 | 3.3805e-03 | 3.6315e-03 | 3.5649e-04 | 3.6961e-04 |
| unweightedHarmonic | 1.9558e-04 | 3.7817e-03 | 2.9427e-03 | 4.1379e-03 | 3.4014e-03 | 3.6188e-03 | 3.3727e-04 | 3.5564e-04 |

At `N = 10` the method changes the potentials and the current jump by less
than 2 % and the assembled current density by 4 to 5 %; the finer levels
have not been run.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
