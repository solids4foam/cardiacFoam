# benchmarkGPU Analysis (2026-06-28)

This report summarizes the slab benchmark campaign for:

- `BuenoOrovio`
- `TNNP`
- `TWorld`

with focus on:

- full runtime
- ionic runtime
- PDE/diffusion runtime
- sampled activation-time accuracy

Machine-readable exports:

- [benchmarkGPU_results_summary.csv](/home/simao/cardiacFoam/tutorials/benchmarkGPU/benchmarkGPU_results_summary.csv:1)
- [benchmarkGPU_results_summary.json](/home/simao/cardiacFoam/tutorials/benchmarkGPU/benchmarkGPU_results_summary.json:1)

## Scope and caveats

- All CPU adaptive MPI timing ladders are now complete for:
  - `BuenoOrovio`
  - `TNNP`
  - `TWorld`
- Clean CPU `np1` Slurm baselines are available for all three models.
- Accuracy is now established for `BuenoOrovio` and `TNNP`, and can be stated
  directly for `TWorld` as well.
- The newly cloned `dt=2e-5` GPU case directories inherited stale
  `scalingResults` from the `dt=1e-5` source directories when they were copied.
  Those inherited GPU+MPI rows are not used here.
  Only the real new `slurm` serial `dt=2e-5` results are used for:
  - `gpu_tnnp_compact_rl10_dt2e-5`
  - `gpu_tworld_compact_rl10_dt2e-5`

## Timing tables

### CPU adaptive scaling

| Case | MPI | Full s | Ionic s | PDE s |
|---|---:|---:|---:|---:|
| `cpu_bueno_adaptive_dt2e-5` | 1 | 3082.03 | 2530.91 | 535.277 |
| `cpu_bueno_adaptive_dt2e-5` | 4 | 764.07 | 688.516 | 69.0882 |
| `cpu_bueno_adaptive_dt2e-5` | 8 | 365.77 | 301.752 | 32.2588 |
| `cpu_bueno_adaptive_dt2e-5` | 16 | 180.38 | 146.517 | 20.0207 |
| `cpu_bueno_adaptive_dt2e-5` | 24 | 120.63 | 100.548 | 14.2877 |
| `cpu_bueno_adaptive_dt2e-5` | 48 | 78.91 | 62.5231 | 14.2146 |
| `cpu_tnnp_adaptive_dt2e-6` | 1 | 18997.8 | 18488.8 | 493.568 |
| `cpu_tnnp_adaptive_dt2e-6` | 4 | 5727.42 | 5636.43 | 72.3533 |
| `cpu_tnnp_adaptive_dt2e-6` | 8 | 3244.55 | 2699.32 | 45.7715 |
| `cpu_tnnp_adaptive_dt2e-6` | 16 | 1760.44 | 1424.68 | 72.6926 |
| `cpu_tnnp_adaptive_dt2e-6` | 24 | 1232.21 | 1030.67 | 104.948 |
| `cpu_tnnp_adaptive_dt2e-6` | 48 | 691.98 | 511.782 | 142.781 |
| `cpu_tworld_adaptive_dt2e-6` | 1 | 38343.0 | 37872.2 | 459.899 |
| `cpu_tworld_adaptive_dt2e-6` | 4 | 11433.5 | 11307.6 | 88.9536 |
| `cpu_tworld_adaptive_dt2e-6` | 8 | 5819.49 | 5252.28 | 93.1461 |
| `cpu_tworld_adaptive_dt2e-6` | 16 | 3201.23 | 2889.36 | 79.8677 |
| `cpu_tworld_adaptive_dt2e-6` | 24 | 2326.70 | 2112.81 | 115.914 |
| `cpu_tworld_adaptive_dt2e-6` | 48 | 1212.29 | 1024.03 | 133.304 |

### GPU serial baselines

| Case | MPI | Full s | Ionic s | PDE s |
|---|---:|---:|---:|---:|
| `gpu_bueno_compact_rl2_dt2e-5` | 1 | 1020.55 | 13.0541 | 992.833 |
| `gpu_bueno_compact_rl5_dt2e-5` | 1 | 955.80 | 24.2272 | 918.099 |
| `gpu_tnnp_compact_rl5_dt2e-6` | 1 | 10405.7 | 1523.59 | 8852.59 |
| `gpu_tnnp_compact_rl10_dt2e-6` | 1 | 11740.4 | 2691.33 | 9017.00 |
| `gpu_tnnp_compact_rl10_dt1e-5` | 1 | 2785.59 | 586.663 | 2179.64 |
| `gpu_tnnp_compact_rl10_dt2e-5` | 1 | 1127.12 | 271.285 | 842.567 |
| `gpu_tworld_compact_rl5_dt2e-6` | 1 | 12211.6 | 3805.32 | 8375.09 |
| `gpu_tworld_compact_rl10_dt1e-5` | 1 | 3260.52 | 1406.01 | 1836.58 |
| `gpu_tworld_compact_rl10_dt2e-5` | 1 | 1553.99 | 704.363 | 834.606 |

### GPU + MPI scaling

| Case | MPI | Full s | Ionic s | PDE s |
|---|---:|---:|---:|---:|
| `gpu_bueno_compact_rl2_dt2e-5` | 8 | 106.61 | 13.8206 | 89.3327 |
| `gpu_bueno_compact_rl2_dt2e-5` | 24 | 58.66 | 13.2347 | 34.4462 |
| `gpu_bueno_compact_rl2_dt2e-5` | 48 | 56.02 | 19.3398 | 21.5633 |
| `gpu_tnnp_compact_rl10_dt1e-5` | 8 | 825.69 | 648.841 | 168.658 |
| `gpu_tnnp_compact_rl10_dt1e-5` | 24 | 867.11 | 767.022 | 66.8303 |
| `gpu_tnnp_compact_rl10_dt1e-5` | 48 | 879.59 | 785.800 | 56.4674 |
| `gpu_tworld_compact_rl10_dt1e-5` | 8 | 1952.07 | 1760.12 | 163.620 |
| `gpu_tworld_compact_rl10_dt1e-5` | 24 | 2175.37 | 2085.45 | 69.0350 |
| `gpu_tworld_compact_rl10_dt1e-5` | 48 | 2138.60 | 1911.90 | 97.0807 |

## Accuracy tables

Accuracy is measured from:

- `postProcessing/Niedererpoints/0/activationTime`
- `postProcessing/Niedererlines/0/activationTime`

Units are seconds.

### BuenoOrovio

CPU reference: `cpu_bueno_adaptive_dt2e-5 np1` Slurm baseline.

| Comparison | Sample | N | Max abs diff | Mean abs diff |
|---|---|---:|---:|---:|
| CPU `np1` vs GPU `rl2` | points | 9 | 0.0012064 | 0.000687238 |
| CPU `np1` vs GPU `rl2` | lines | 101 | 0.0012064 | 0.000587770 |
| CPU `np1` vs GPU `rl5` | points | 9 | 0.0012064 | 0.000687238 |
| CPU `np1` vs GPU `rl5` | lines | 101 | 0.0012064 | 0.000587770 |

Interpretation:

- `rl2` and `rl5` are effectively equivalent for sampled activation times.
- `BuenoOrovio` accuracy is established.

### TNNP

CPU reference: `cpu_tnnp_adaptive_dt2e-6 np1` Slurm baseline.

| Comparison | Sample | N | Max abs diff | Mean abs diff |
|---|---|---:|---:|---:|
| CPU `np1` vs GPU `rl5 dt2e-6` | points | 9 | 0.0003831 | 0.000200407 |
| CPU `np1` vs GPU `rl5 dt2e-6` | lines | 101 | 0.0003831 | 0.000151286 |
| CPU `np1` vs GPU `rl10 dt2e-6` | points | 9 | 0.0003831 | 0.000200407 |
| CPU `np1` vs GPU `rl10 dt2e-6` | lines | 101 | 0.0003831 | 0.000151285 |
| CPU `np1` vs GPU `rl10 dt1e-5` | points | 9 | 0.0001737 | 0.0000847678 |
| CPU `np1` vs GPU `rl10 dt1e-5` | lines | 101 | 0.0001654 | 0.0000486361 |
| CPU `np1` vs GPU `rl10 dt2e-5` | points | 9 | 0.0001967 | 0.0000891822 |
| CPU `np1` vs GPU `rl10 dt2e-5` | lines | 101 | 0.0001217 | 0.0000790247 |

Interpretation:

- Clean Slurm CPU results removed the earlier false zero-activation mismatch.
- `dt1e-5 rl10` is already very accurate.
- `dt2e-5 rl10` remains very accurate on sampled activation times.
- For `TNNP`, the larger PDE timestep looks acceptable on this metric.

### TWorld

CPU reference: `cpu_tworld_adaptive_dt2e-6 np1` Slurm baseline.

| Comparison | Sample | N | Max abs diff | Mean abs diff |
|---|---|---:|---:|---:|
| CPU `np1` vs GPU `rl5 dt2e-6` | points | 9 | 0.0005442 | 0.000275750 |
| CPU `np1` vs GPU `rl5 dt2e-6` | lines | 101 | 0.0005442 | 0.000222296 |
| CPU `np1` vs GPU `rl10 dt1e-5` | points | 9 | 0.0003418 | 0.000167151 |
| CPU `np1` vs GPU `rl10 dt1e-5` | lines | 101 | 0.0003375 | 0.000122277 |
| CPU `np1` vs GPU `rl10 dt2e-5` | points | 9 | 0.0001539 | 0.0000787533 |
| CPU `np1` vs GPU `rl10 dt2e-5` | lines | 101 | 0.0000966 | 0.0000274970 |

Interpretation:

- `TWorld` accuracy is now also established against a clean CPU `np1` baseline.
- `dt1e-5 rl10` was already close.
- `dt2e-5 rl10` is even closer on sampled activation times.

## Speedup summaries

### CPU MPI speedup vs CPU `np1`

| Case | Baseline `np1` s | Best completed s | Best run | Speedup |
|---|---:|---:|---|---:|
| `BuenoOrovio` CPU adaptive | 3082.03 | 78.91 | `np48` | 39.1x |
| `TNNP` CPU adaptive | 18997.8 | 691.98 | `np48` | 27.5x |
| `TWorld` CPU adaptive | 38343.0 | 1212.29 | `np48` | 31.6x |

### GPU serial speedup vs CPU `np1`

| Model | CPU `np1` s | Best GPU serial s | GPU serial run | Speedup |
|---|---:|---:|---|---:|
| `BuenoOrovio` | 3082.03 | 955.80 | `rl5 dt2e-5` | 3.22x |
| `TNNP` | 18997.8 | 1127.12 | `rl10 dt2e-5` | 16.9x |
| `TWorld` | 38343.0 | 1553.99 | `rl10 dt2e-5` | 24.7x |

### GPU + MPI speedup vs GPU `np1`

| Model / run | GPU `np1` s | Best GPU+MPI s | Best run | Speedup |
|---|---:|---:|---|---:|
| `BuenoOrovio rl2 dt2e-5` | 1020.55 | 56.02 | `np48` | 18.2x |
| `TNNP rl10 dt1e-5` | 2785.59 | 825.69 | `np8` | 3.37x |
| `TWorld rl10 dt1e-5` | 3260.52 | 1952.07 | `np8` | 1.67x |

## High-level analysis

### 1. CPU-only MPI scaling

- All three models scale strongly on CPU out to `48` ranks.
- `BuenoOrovio` shows the cleanest MPI scaling curve.
- `TNNP` and `TWorld` remain dominated by ionic cost, but still get large gains
  from strong scaling.

### 2. GPU serial behavior

- `BuenoOrovio` gains a modest end-to-end speedup over CPU because the PDE side
  remains dominant.
- `TNNP` and `TWorld` gain much more once the PDE timestep is increased.
- The new `dt=2e-5 rl10` runs are substantially faster than the older `dt2e-6`
  GPU runs for both `TNNP` and `TWorld`.

### 3. GPU + MPI behavior

- `BuenoOrovio` benefits strongly from hybrid MPI + GPU because PDE cost is the
  dominant part of the serial GPU run.
- `TNNP` improves sharply from GPU `np1` to `np8`, then flattens by `24` and
  `48`.
- `TWorld` shows the same pattern: improvement at `np8`, then no real gain at
  larger rank counts on a single GPU.

This is consistent with a one-GPU hybrid scaling limit:

- MPI reduces PDE time.
- The ionic part remains on one physical GPU.
- Once PDE time is no longer dominant, the GPU/ionic side plus MPI overhead
  becomes the floor.

### 4. Why BuenoOrovio differs from TNNP / TWorld

The macro PDE timestep is different:

- `BuenoOrovio`: `deltaT = 2e-5`
- `TNNP` / `TWorld`: conservative CPU and older GPU runs at `deltaT = 2e-6`

So the heavy models perform about `10x` more macro PDE steps over the same
physical time interval. Same mesh does not imply the same measured PDE wall
time.

### 5. Timestep escalation result

The new `dt=2e-5 rl10` tests matter:

- `TNNP dt2e-5 rl10` is both fast and accurate.
- `TWorld dt2e-5 rl10` is also both fast and accurate.

At the level of sampled activation times, pushing from `dt1e-5 rl10` to
`dt2e-5 rl10` did not degrade these comparisons; it improved them.

## Practical conclusions

- `BuenoOrovio`
  - Accuracy is established.
  - `GPU rl2` is already sufficient.
  - Hybrid MPI + GPU scales very well.

- `TNNP`
  - Accuracy is established.
  - `GPU rl10 dt2e-5` currently looks like the best serial GPU candidate from
    this campaign.
  - Hybrid MPI + GPU helps up to moderate MPI counts, then saturates on one GPU.

- `TWorld`
  - Accuracy is established.
  - `GPU rl10 dt2e-5` is currently the best serial GPU candidate from this
    campaign.
  - Hybrid MPI + GPU helps at `np8`, but larger rank counts do not pay on one
    GPU.

## Remaining gaps

No major gap remains in the current dataset for the models and configurations
that were actually run.

The only caution is interpretive:

- do not use the inherited stale `scalingResults` copied into the new
  `dt=2e-5` case directories for GPU+MPI analysis
- use only the real new `slurm` serial results for those `dt=2e-5` cases
