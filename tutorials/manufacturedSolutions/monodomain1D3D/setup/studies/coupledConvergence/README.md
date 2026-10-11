# coupledConvergence - monodomain1D3D

Joint refinement of the 3D mesh and the graph (`k = N + 1` nodes, edge length
`1/N`) with `dt` proportional to `h^2`, for the four combinations of transfer
(`-coupling unidirectional|bidirectional`) and tissue-side assembly (`-pvj
explicit|implicit`), with `R_pvj = 1`. It measures the order of the 3D and
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
    for c in unidirectional bidirectional
    do
        for p in explicit implicit
        do
            rsync -a --exclude setup ../../../ results/${c}_${p}_N$n/
            results/${c}_${p}_N$n/Allrun -N $n -dt $dt -coupling $c -pvj $p
        done
    done
done
python3 ../summarize.py results
```

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

Bidirectional, explicit:

| N | h | 3D Vm L2 | 3D Vm Linf | 1D Vm L2 | 1D Vm Linf | source error L1 | order 3D L2 | order 3D Linf | order 1D L2 | order 1D Linf | order source |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.1 | 1.8499e-03 | 4.9623e-03 | 1.8183e-04 | 3.4724e-04 | 3.3137e-03 | | | | | |
| 20 | 0.05 | 4.8860e-04 | 1.3812e-03 | 4.5747e-05 | 8.9172e-05 | 7.8330e-04 | 1.92 | 1.85 | 1.99 | 1.96 | 2.08 |
| 40 | 0.025 | 1.2168e-04 | 3.4859e-04 | 1.1205e-05 | 2.2053e-05 | 1.9266e-04 | 2.01 | 1.99 | 2.03 | 2.02 | 2.02 |
| 80 | 0.0125 | 3.0525e-05 | 8.7734e-05 | 2.7884e-06 | 5.5128e-06 | 4.8111e-05 | 2.00 | 1.99 | 2.01 | 2.00 | 2.00 |

Bidirectional, implicit:

| N | h | 3D Vm L2 | 3D Vm Linf | 1D Vm L2 | 1D Vm Linf | source error L1 | order 3D L2 | order 3D Linf | order 1D L2 | order 1D Linf | order source |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.1 | 1.8483e-03 | 4.9600e-03 | 1.8163e-04 | 3.4627e-04 | 7.9590e-04 | | | | | |
| 20 | 0.05 | 4.8737e-04 | 1.3812e-03 | 4.5759e-05 | 8.9267e-05 | 2.2681e-04 | 1.92 | 1.84 | 1.99 | 1.96 | 1.81 |
| 40 | 0.025 | 1.2138e-04 | 3.4859e-04 | 1.1206e-05 | 2.2068e-05 | 5.9448e-05 | 2.01 | 1.99 | 2.03 | 2.02 | 1.93 |
| 80 | 0.0125 | 3.0448e-05 | 8.7734e-05 | 2.7886e-06 | 5.5150e-06 | 1.4973e-05 | 2.00 | 1.99 | 2.01 | 2.00 | 1.99 |

All three quantities converge at second order in every combination. The 3D
voltage does not depend on the transfer direction; the implicit assembly
changes its L2 error by under 0.3% and its Linf error by under 0.05%. Under
unidirectional transfer the graph receives no current, so its error is the
same for both assemblies; bidirectional transfer changes it by under 1.2%. The
implicit assembly lowers the junction source error by a factor of 3.2 to 4.2,
with a lower order on the coarse pairs (1.81, 1.93) that reaches 1.99 on the
finest pair.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
