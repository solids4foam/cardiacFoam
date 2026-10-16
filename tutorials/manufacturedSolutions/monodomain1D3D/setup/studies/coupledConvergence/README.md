# coupledConvergence - monodomain1D3D

Joint refinement of the 3D mesh and the graph (`k = N + 1` nodes, edge length
`1/N`) with `dt` proportional to `h^2`, with `R_pvj = 1`, for unidirectional
transfer with explicit and implicit tissue-side assembly (`-pvj
explicit|implicit`) and for bidirectional transfer (`-coupling bidirectional`),
which the coupler always assembles explicitly on the tissue side, so `-pvj`
has no effect on it. It measures the order of the 3D and
graph voltages and of the junction source error, the difference between the
numerical junction source and the exact one built with the same map.

## Ladder

| N | graph | dt | final time |
|---|---|---|---|
| 10 | nodes011 | 0.008971136 | 0.0986825 |
| 20 | nodes021 | 0.002242784 | 0.100925 |
| 40 | nodes041 | 0.000560696 | 0.0998039 |
| 80 | nodes081 | 0.000140174 | 0.0999441 |

These steps do not divide `endTime 0.1`, so each run stops at the step
nearest 0.1 (the verifiers write the final time and evaluate the exact fields
there). Without `-dt` the case uses
`endTime/ceil(endTime/(0.008971136 (10/N)^2))`, which ends exactly at 0.1.

## Execution

Each run is a copy of the case in `results/` (ignored by git), named
`<coupling>_<pvj>_N<n>`, the names that `../summarize.py` reads. From this
folder:

```bash
mkdir -p results
N=(10 20 40 80)
DT=(0.008971136 0.002242784 0.000560696 0.000140174)
for i in 0 1 2 3
do
    n=${N[$i]}; dt=${DT[$i]}
    for arm in "unidirectional explicit" "unidirectional implicit" \
               "bidirectional explicit"
    do
        set -- $arm; c=$1; p=$2
        rsync -a --exclude setup ../../../ results/${c}_${p}_N$n/
        results/${c}_${p}_N$n/Allrun -N $n -dt $dt -coupling $c -pvj $p
    done
done
python3 ../summarize.py results
```

The omnidriver specs for this study are `sweep_active.json` and `sweep_bidirectional.json`.

## Measured behaviour

Results of the last run of this study (`results/`).

Unidirectional, explicit:

| N | h | 3D Vm L2 | 3D Vm Linf | 1D Vm L2 | 1D Vm Linf | source error L1 | order 3D L2 | order 3D Linf | order 1D L2 | order 1D Linf | order source |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.1 | 1.8499e-03 | 4.9623e-03 | 1.8097e-04 | 3.4317e-04 | 3.3123e-03 | | | | | |
| 20 | 0.05 | 4.8860e-04 | 1.3812e-03 | 4.5652e-05 | 8.8475e-05 | 7.8308e-04 | 1.92 | 1.85 | 1.99 | 1.96 | 2.08 |
| 40 | 0.025 | 1.2168e-04 | 3.4859e-04 | 1.1195e-05 | 2.1956e-05 | 1.9263e-04 | 2.01 | 1.99 | 2.03 | 2.01 | 2.02 |
| 80 | 0.0125 | 3.0525e-05 | 8.7734e-05 | 2.7873e-06 | 5.5004e-06 | 4.8107e-05 | 2.00 | 1.99 | 2.01 | 2.00 | 2.00 |

Unidirectional, implicit:

| N | h | 3D Vm L2 | 3D Vm Linf | 1D Vm L2 | 1D Vm Linf | source error L1 | order 3D L2 | order 3D Linf | order 1D L2 | order 1D Linf | order source |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.1 | 1.8483e-03 | 4.9600e-03 | 1.8097e-04 | 3.4317e-04 | 7.9656e-04 | | | | | |
| 20 | 0.05 | 4.8737e-04 | 1.3812e-03 | 4.5652e-05 | 8.8475e-05 | 2.2710e-04 | 1.92 | 1.84 | 1.99 | 1.96 | 1.81 |
| 40 | 0.025 | 1.2138e-04 | 3.4859e-04 | 1.1195e-05 | 2.1956e-05 | 5.9484e-05 | 2.01 | 1.99 | 2.03 | 2.01 | 1.93 |
| 80 | 0.0125 | 3.0448e-05 | 8.7734e-05 | 2.7873e-06 | 5.5004e-06 | 1.4978e-05 | 2.00 | 1.99 | 2.01 | 2.00 | 1.99 |

Bidirectional (explicit on the tissue side; a run with `-pvj implicit` gives
the same numbers to every printed digit):

| N | h | 3D Vm L2 | 3D Vm Linf | 1D Vm L2 | 1D Vm Linf | source error L1 | order 3D L2 | order 3D Linf | order 1D L2 | order 1D Linf | order source |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.1 | 1.8501e-03 | 4.9622e-03 | 1.9770e-04 | 4.1325e-04 | 3.3370e-03 | | | | | |
| 20 | 0.05 | 4.8865e-04 | 1.3812e-03 | 4.9769e-05 | 1.1425e-04 | 7.9166e-04 | 1.92 | 1.85 | 1.99 | 1.85 | 2.08 |
| 40 | 0.025 | 1.2170e-04 | 3.4859e-04 | 1.2074e-05 | 2.9211e-05 | 1.9496e-04 | 2.01 | 1.99 | 2.04 | 1.97 | 2.02 |
| 80 | 0.0125 | 3.0529e-05 | 8.7734e-05 | 2.9816e-06 | 7.3769e-06 | 4.8702e-05 | 2.00 | 1.99 | 2.02 | 1.99 | 2.00 |

All three quantities converge at second order in every arm. The 3D voltage
hardly depends on the transfer direction (bidirectional transfer changes its
L2 error by under 0.02%); the implicit assembly changes its L2 error by under
0.3% and its Linf error by under 0.05%. Under unidirectional transfer the
graph receives no current, so its error is the same for both assemblies.
Bidirectional transfer, which takes each junction current from the volume of
the terminal node, raises the graph L2 error by 7 to 9% and its Linf error by
20 to 34%, with a Linf order of 1.85 on the coarsest pair. Under
unidirectional transfer the implicit assembly lowers the junction source error
by a factor of 3.2 to 4.2, with a lower order on the coarse pairs (1.81, 1.93)
that reaches 1.99 on the finest pair; the bidirectional source error is within
1.2% of the unidirectional explicit one.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
