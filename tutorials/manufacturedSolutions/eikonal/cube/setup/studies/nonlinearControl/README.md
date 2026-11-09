# nonlinearControl - eikonalECG

Share of the outer fixed-point stopping tolerance in the tetrahedral errors:
the three `leastSquares` configurations of
[tetConvergence](../tetConvergence/README.md) (axis conductivity with the
explicit formulation, rotated conductivity with each formulation) rerun with
the outer residual tolerance tightened from `1e-8` to `1e-10`
(`-tolerance 1e-10`), on the same meshes.

## Ladder

| N | lc | tetrahedra (gmsh 4.15.2, seed 1) |
|---|---|---|
| 10 | 0.1 | 4874 |
| 20 | 0.05 | 37040 |
| 30 | 0.0333 | 123211 |
| 40 | 0.025 | 289388 |

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<conductivity>_<formulation>_tolerance1e-10_N<n>`, the names that
`../summarize.py` reads. From this folder:

```bash
mkdir -p results
for n in 10 20 30 40
do
    for c in "axis explicit" "rotated explicit" "rotated advectionDiffusion"
    do
        set -- $c
        run=results/$1_$2_tolerance1e-10_N$n
        rsync -a --exclude setup ../../../ $run/
        $run/Allrun parallel -mesh tet -N $n -conductivity $1 -formulation $2 -tolerance 1e-10
    done
done
python3 ../summarize.py results
```

The omnidriver spec for this study is `sweep_tet_outer_tolerance.json`.

## Measured behaviour

Results of the last run of this study (`results/`), against the matching
`tolerance 1e-8` runs of tetConvergence: activation-time L2 and Linf, outer
iterations at each tolerance, and the largest absolute change of the L1, L2
and Linf errors.

| configuration | N | L2 | Linf | outer iterations (1e-8 / 1e-10) | largest change |
|---|---|---|---|---|---|
| axis, explicit | 10 | 2.2472e-02 | 2.3915e-01 | 114 / 147 | 2.0e-07 |
| | 20 | 8.7256e-03 | 1.5229e-01 | 394 / 518 | 1.2e-06 |
| | 30 | 4.4251e-03 | 1.0902e-01 | 837 / 1119 | 2.3e-06 |
| | 40 | 2.8975e-03 | 8.7812e-02 | 1440 / 1945 | 3.5e-06 |
| rotated, explicit | 10 | 2.2572e-02 | 2.4153e-01 | 115 / 148 | 2.0e-07 |
| | 20 | 8.8134e-03 | 1.5080e-01 | 399 / 525 | 1.1e-06 |
| | 30 | 4.5845e-03 | 1.0937e-01 | 851 / 1136 | 2.0e-06 |
| | 40 | 3.0706e-03 | 8.6560e-02 | 1463 / 1975 | 2.6e-06 |
| rotated, advectionDiffusion | 10 | 2.1746e-02 | 2.4532e-01 | 129 / 166 | 1.1e-06 |
| | 20 | 8.4877e-03 | 1.5163e-01 | 419 / 551 | 6.0e-08 |
| | 30 | 4.5240e-03 | 1.0947e-01 | 876 / 1170 | 1.0e-06 |
| | 40 | 3.0268e-03 | 8.8702e-02 | 1492 / 2013 | 2.9e-07 |

Tightening the tolerance by two orders changes the activation-time errors by
at most 3.5e-6 (0.3 % of L1, 0.03 % of L2) and the orders by at most 0.01:
the tolerance `1e-8` does not limit the tetrahedral errors over this ladder.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
