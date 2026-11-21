# benchmarkLadder - NiedererEtAl2011verification/bidomain

The Niederer et al. (2011) resolutions with the bidomain model: Δx = 0.5,
0.2 and 0.1 mm at Δt = 0.005 ms, for the published conductivities and for
the equal-anisotropy set. The equal set solves the discrete monodomain
problem of the [monodomain case](../../../../README.md) (see
[What the equal set reproduces](../../../README.md#what-the-equal-set-reproduces)),
so its rows are compared with the monodomain case at the same Δx and Δt; the
published rows give the bidomain activation times of the benchmark.

## Ladder

| Δx (mm) | `-dx` | cells | `-endTime` (s) | run |
|---|---|---|---|---|
| 0.5 | `5e-4` | 40 x 6 x 14 = 3360 | 0.2 | serial |
| 0.2 | `2e-4` | 100 x 15 x 35 = 52500 | 0.08 | serial |
| 0.1 | `1e-4` | 200 x 30 x 70 = 420000 | 0.055 | `parallel` |

Every run uses `-dt 5e-6`. The end times are those of the monodomain
`cartesianConvergence` study, after the last probe activates at each Δx.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<conductivity>_dx<dx>`. From this folder:

```bash
mkdir -p results
for conductivity in equal published
do
    for level in "5e-4 0.2" "2e-4 0.08" "1e-4 0.055 parallel"
    do
        set -- $level
        rsync -a --exclude setup ../../../ results/${conductivity}_dx$1/
        results/${conductivity}_dx$1/Allrun $3 -dx $1 -dt 5e-6 -endTime $2 -conductivity $conductivity
    done
done
for run in results/*_dx*
do
    echo "$run $(tail -n 1 $run/postProcessing/Niedererpoints/0/activationTime)"
done
```

The 0.2 mm runs can also take `parallel`. `system/decomposeParDict` sets 6
subdomains.

The omnidriver sweep format of this ladder is `sweep_hex_ladder.json`: its
axes are the `Allrun` options (`dx`, `dt`, `endTime`, `conductivity`,
`parallel`) and its entry `niederer2011Bidomain` names this case.

## Comparison with the monodomain case

The monodomain rows with the same Δx, Δt and end time are those of the
monodomain [`cartesianConvergence`](../../../../setup/studies/cartesianConvergence/sweep_hex_convergence.json)
study at `deltaT 5e-6`. Equal against monodomain is a consistency check of
the bidomain implementation: the two differ by the tolerances of the linear
solves and of the outer `phiE`-`Vm` iterations. Published against equal is
the bidomain effect of unequal anisotropy ratios at fixed monodomain tensor.

## Measured behaviour

Equal set against the monodomain case, serial unless stated, on the last
runs on one machine (OpenFOAM v2412, arm64):

| Δx (mm) | Δt (ms) | end (s) | activation-time difference | bidomain / monodomain wall time |
|---|---|---|---|---|
| 0.5 | 0.005 | 0.2 | none in the 6 written digits: all 3360 cells, the nine probes and the 101 line points | 1490 s / 310 s |
| 0.5 | 0.01 | 0.015 | at most 4.5e-15 s over the 129 activated cells, `Vm` at most 4.9e-13 V, with `writePrecision 15` | 70 s / 13 s |
| 0.2 | 0.01 | 0.015 | none against the monodomain regression reference: probe 0 `0.00119338` s, probe 4 `0.0112285` s | 590 s on 6 ranks |

Probe activation times at Δx 0.5 mm, Δt 0.005 ms (s), both models:

| probe | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| activation | 0.0011922 | 0.132193 | 0.0465417 | 0.142004 | 0.0358999 | 0.133508 | 0.0559116 | 0.142914 | 0.0657129 |

On the 0.5 mm mesh the `PIMPLE` loop converges in 2 or 3 outer correctors
per step and a `phiE` solve takes at most 78 PCG iterations, far below the
5000 cap.
