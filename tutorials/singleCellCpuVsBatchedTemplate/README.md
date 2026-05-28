# singleCellCpuVsBatchedBuenoOrovio

Clean benchmark harness for the first migrated batched ionic model.

It compares the adaptive CPU `BuenoOrovio` single-cell path against these
`BuenoOrovioBatched` modes:

- `batched_euler`: cell-major Euler, `batchedSubsteps 100`
- `batched_rl`: cell-major Rush-Larsen gates plus Euler fallback, `batchedSubsteps 1`
- `batched_soa`: state-major SoA Euler, `useSoAEvaluator true`, `batchedSubsteps 100`

The default run is intentionally short:

```bash
./run_bueno_orovio_batched_comparison.sh
```

Useful overrides:

```bash
N_RUNS=10 CONTROL_END_TIME=5 ./run_bueno_orovio_batched_comparison.sh
MODES="cpu batched_soa" N_RUNS=20 ./run_bueno_orovio_batched_comparison.sh
```

Results are written to:

```text
comparisonResults/comparison_metrics.json
comparisonResults/plots/
```

The plot directory contains one combined accuracy plot per exported variable,
using CPU as the reference:

```text
combined_Vm_accuracy_vs_cpu.png
combined_s_accuracy_vs_cpu.png
scatter_Vm_wall_time_vs_accuracy.png
scatter_s_wall_time_vs_accuracy.png
```

The scatter plots put mean wall time on the x-axis and maximum absolute
difference from the CPU reference on the y-axis.

This is a host-batched and SoA benchmark for the migration stage. It validates
the architecture and numerical-policy differences before CUDA/HIP/OpenMP-target
backend work is added.
