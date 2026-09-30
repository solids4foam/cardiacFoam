# SBDF2 manufactured temporal check

This check adopts the coarse 2D temporal template from the historical GPU
backup. It holds the 640×640 mesh fixed, runs to `t=0.2`, and uses
`dt = 0.025, 0.0125, 0.00625, 0.003125, 0.0015625`. The test uses the scalar
`monodomainFDAManufactured` model with `RKF45` and a fixed initial ODE step;
that manufactured model has no Rush–Larsen setting. It is not a cardiac
action-potential test. Its analytic solution lets the verifier measure the
tissue time-integration error directly.

The tissue scheme is an implicit monodomain solve, `timeCouplingScheme
sbdf2`, and `ddt(Vm) backward`. In the batched case, the model dictionary also
selects `batchedIntegrator euler`, which advances the manufactured ionic
states `u1` and `u2` and changes the computed ionic current. It does not
replace the backward SBDF2 update of tissue voltage `Vm`.

For the batched ladder, `batchedSubsteps` doubles whenever tissue `deltaT`
halves: 25, 50, 100, 200, 400. This makes the ionic Euler substep decrease
quadratically with tissue `deltaT` (`deltaT_ion = deltaT/batchedSubsteps`),
so its first-order error should not dominate the second-order tissue scheme.
All five substep counts and actual ionic steps are written into `summary.csv`.
The scalar reference retains its adaptive RKF45 ODE solver.

## Acceptance rule

The four coarse levels through `dt=0.003125` must have monotonically
decreasing L2 errors for `Vm`, `u1`, and `u2`, with each observed order
between 1.8 and 2.3. This window brackets the backup's measured `Vm` orders
(1.97, 2.01, 2.12). The `dt=0.0015625` error and order are reported but are
excluded from the gate because the backup indicates that this level is
beginning to approach the fixed-mesh spatial floor. CPU-batched and
CUDA-batched fields and states must also match within `1e-10`. A missing
case or verifier norm fails the run. Results include per-case logs,
`summary.csv`, and `acceptance.txt`.

## Run

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
python3 tutorials/ioniGPUKernel/run_sbdf2_manufactured_temporal.py \
    /tmp/cardiac_sbdf2_mms_temporal_manual
```

The five-case scalar temporal ladder completed on xenosim as job 9359. Its
output is `/tmp/cardiac_sbdf2_mms_temporal_20260928/`.

| `deltaT` | `Vm L2` | Observed order | Runtime (s) |
|---:|---:|---:|---:|
| 0.025 | 4.36024e-4 | — | 17.29 |
| 0.0125 | 1.11309e-4 | 1.96983 | 21.97 |
| 0.00625 | 2.76440e-5 | 2.00954 | 30.07 |
| 0.003125 | 6.41272e-6 | 2.10796 | 44.62 |
| 0.0015625 | 1.07423e-6 | 2.57763* | 69.70 |

The acceptance rule passed. The first three halving orders are the declared
second-order window; the last order is reported but excluded because that
level approaches the fixed-mesh spatial floor. This validates the scalar
SBDF2/backward tissue method at this setup.

The batched SBDF2 ladder completed as job 9362 using the same mesh,
parameters, zero stimulus, end time, and tissue `deltaT` ladder, with batched
Euler substeps 25, 50, 100, 200, and 400 as `deltaT` halves. CPU and CUDA
both passed monotonic L2-error reduction and second-order gates for `Vm`,
`u1`, and `u2` over the four coarsest levels:

| Field | Observed orders (`dt` .025 to .003125) |
|---|---|
| `Vm` | 1.97109, 2.00964, 2.10618 |
| `u1` | 1.91077, 1.96285, 2.00616 |
| `u2` | 1.96748, 1.98566, 2.00410 |

CPU-batched and CUDA-batched fields and saved states matched at all five
levels within `1e-10`; the largest reported absolute difference was below
`2e-14`. Batched `Vm L2` errors were about 1.8–1.9% above the scalar SBDF2
reference on the four coarse levels and 2.8% at the finest reported level.
The finest level is excluded from the order gate due to the spatial floor.
Job 9355 shared the node, so job 9362 timings are not performance results.

The manufactured model is not used to test cardiac Rush–Larsen because it
does not have Rush–Larsen gate equations. Repeat the batched ladder with
`sbatch tutorials/ioniGPUKernel/run_manufactured_batched_sbdf2.slurm` after
moving or renaming the existing output directory.

## Fixed-substep Euler study

GPU job 9450 completed all 30 fixed-substep SBDF2 runs for Euler substep
counts 1, 5, 25, 100, 400, and 800 across the five tissue steps. At fixed
substep counts, the first three observed `Vm` orders (`deltaT` 0.025 to
0.003125) were 1.601, 1.451, 1.307 for 1; 1.861, 1.800, 1.718 for 5;
1.946, 1.960, 2.001 for 25; 1.964, 1.997, 2.079 for 100; 1.968, 2.006,
2.101 for 400; and 1.969, 2.008, 2.104 for 800 substeps. The low-substep
results show the first-order Euler contribution lowering the apparent
SBDF2 order. With 25 or more substeps, the runs recover a near-second-order
trend over these levels. Orders at the finest tissue step were excluded from
that comparison because they rise as the fixed-mesh spatial floor is
approached.

Against the 800-substep run at each `deltaT`, the 400-substep `Vm` RMSE was
at most 1.45% of the 800-substep analytic `Vm L2` error (maximum absolute
RMSE `2.62e-7`). This suggests that 400 substeps are sufficient for this
manufactured setup at the tested tissue steps, while 800 remains a numerical
reference rather than an exact solution. The outputs are correctness data;
job 9355 shared the node, so no timing comparison is made. Results are in
`/tmp/cardiac_manufactured_euler_substeps_20260929/` (`summary.csv` and
`substep_effect_vs_finest.csv`).
