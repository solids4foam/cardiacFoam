# GPU ionic-model validation: discussion brief

**Snapshot:** 2026-09-29, updated after completion of jobs 9355 and 9450
**Target:** `cardiacFoam` branch `main`, commit `39333ef3`
**GPU:** NVIDIA RTX 4000 Ada, 20 GB; OpenFOAM 2412. The reported CUDA builds
used CUDA 11.5 and GCC 10 for nvcc. Most benchmark jobs used one MPI rank,
one CPU core, and one GPU.

## Status in one minute

- All 11 in-scope tissue models (all single-cell models except Fabbri) have
  completed scalar-reference versus CUDA-batched single-cell runs, a short 2D
  slab, and the 3D **Niederer geometry** matrix. The GPU backend was required
  by the run scripts. These are meaningful implementation checks, but the
  short 15 ms slab runs do not establish full action-potential, recovery, or
  production timestep validity.
- The previous 10 us / 1, 5, 25 substep integrator screen completed for all
  11 models. Two GPU cases failed at 1 substep (Grandi and TWorld); the
  provisional best settings used 5 substeps for both. Several models benefit
  materially from the GPU, but settings and acceptance remain provisional.
- The published Niederer benchmark reproduction is in scope for **TNNP**.
  Its current-main case used the original 10 us step and five substeps, and
  passed the supplied six-probe activation regression. The all-model
  Niederer-shaped runs are geometry comparisons, not benchmark reproductions
  for every model.
- The manufactured-FDA monodomain batched model now has CPU and CUDA
  implementations. Job 9361 completed a fixed-mesh temporal ladder on both
  backends. Both passed the first-order `Vm L2` gate, and CPU-batched versus
  CUDA-batched fields/states matched within `2.0e-14`. Its analytic `Vm`
  error norms agree with the backup scalar Godunov/Euler reference to about
  0.4% at every tested step. This validates the manufactured convergence and
  verification path; it is not pointwise scalar-vs-batched trace equality.
- The GPU backup's temporal reference is
  `tutorials/manufacturedSolutions/monodomainPseudoECG/setupManufacturedFDA/`
  (`TEMPORAL_CONVERGENCE_REPORT.md`), with fixed 2D `N=640`, final `t=0.2`,
  `Vm/u1/u2` L1/L2/Linf errors, and observed order over timestep halvings.
  Its Godunov/Euler method is first order; its SBDF2/backward method is
  second order until the spatial floor. Job 9359 completed the scalar
  SBDF2/backward ladder with `Vm L2` halving orders 1.96983, 2.00954, and
  2.10796 over the declared coarse window (acceptance PASS). It does not
  validate the batched GPU model.
- Job 9358 completed its CPU and CUDA case matrices, but the direct TNNP
  CPU-batched/GPU-batched SBDF2 comparison **failed** the strict parity gate.
  The run script stopped at that failure; a manual comparison of the saved
  TWorld cases failed as well. Both discrepancies remain open.
- Job 9361 completed the manufactured-batched Godunov/Euler ladder. CPU and
  CUDA `Vm L2` halving orders were 0.973278, 0.987154, and 0.995009 over the
  four coarsest points. Across the five levels, CPU-batched/CUDA-batched
  fields and restart states matched within `1.9984e-14` maximum absolute
  difference against a `1e-10` gate. Against the backup scalar Godunov/Euler
  reference `Vm L2` norms, batched errors are higher by only 0.39–0.43% over
  all five levels, with the same first-order trend.

## Implementation changes currently in the target worktree

The GPU port was compared against the detached backup and applied as targeted
fixes while retaining newer `main` model files and local feature work. Current
changes include uploading the current tissue Vm slice when required, honoring
the requested batched Euler versus Rush–Larsen path in CUDA launchers, and
using each model's hot-path support calculation in the shared batched host
path. CUDA build and launcher compatibility fixes include OpenFOAM's `-iquote`
flag handling, checked launch failures, and device-compatible math calls.
These fixes have been exercised by the compiled GPU matrices below; they do
not yet replace the required per-state/parameter audit or longer accuracy
tests.

## Per-model coverage and completed 3D geometry comparison

The completed 3D matrix used 52,500 cells, 15 ms of simulation, and 25 GPU
ionic substeps. Scalar and GPU-batched cases used matched tissue settings.
The scalar reference used RKF45 while the CUDA-batched path used
Rush–Larsen plus Euler; timings are full `cardiacFoam` execution times, not
isolated kernel timings.
`Speedup` is scalar runtime divided by GPU runtime. `Activation p95` compares
solver activation-time fields at 15 ms. The measured geometry is the
Niederer-shaped slab; this is not, by itself, a published benchmark
reproduction.

| Model | Single-cell trace | 2D slab | 3D Niederer geometry | Scalar / GPU (s) | Speedup | Final Vm RMSE (mV) | Activation p95 (ms) |
|---|---|---|---|---:|---:|---:|---:|
| AlievPanfilov | completed | completed | completed | 205.34 / 205.80 | 1.00x | 0.000004 | 0.000000 |
| BuenoOrovio | completed | completed | completed | 274.65 / 79.22 | 3.47x | 0.006426 | 0.000600 |
| Courtemanche | completed | completed | completed | 483.18 / 105.79 | 4.57x | 0.000172 | 0.000010 |
| Gaur | completed | completed | completed | 1914.98 / 230.89 | 8.29x | 0.000333 | 0.000000 |
| Grandi | completed on GPU | completed | completed | 1000.34 / 267.46 | 3.74x | 0.000030 | 0.000100 |
| PerisYague | completed | completed | completed | 460.95 / 99.66 | 4.63x | 0.000149 | 0.000000 |
| Stewart | completed | completed | completed | 934.53 / 120.33 | 7.77x | 0.001410 | 0.000010 |
| TNNP | completed | completed | completed | 1199.76 / 116.47 | 10.30x | 0.000174 | 0.000010 |
| TWorld | completed on GPU | completed | completed | 2042.49 / 310.96 | 6.57x | 0.000520 | 0.000100 |
| ToRORd_dynCl | completed | completed | completed | 2048.76 / 269.44 | 7.60x | 0.011899 | 0.000900 |
| Trovato | completed | completed | completed | 796.89 / 271.19 | 2.94x | 0.000061 | 0.000010 |

The first full single-cell matrix had invalid trace extraction for Grandi and
TWorld. Targeted 0.5 s GPU reruns produced valid traces and APD/current
comparisons for both. A separate CPU-batched run of those two at 10 us
failed with a floating-point exception in Grandi's `pow` calculation; the
scalar references and CUDA runs completed. This needs a controlled CPU-batched
reproduction before attributing it to the equations or to the GPU port.

In the targeted 10 us full-trace comparisons, Gaur, Grandi, TNNP, and TWorld
had APD90 differences of approximately +0.036, −0.009, −0.0034, and
−0.0030 ms, respectively. Their activation shifts were 0.00076–0.00352 ms.
These support the fine-step baseline, while the older full-matrix Gaur result
was materially worse and needs a like-for-like rerun before its discrepancy
is considered closed.

Fabbri is excluded from the tissue matrix because it is an AN-node model. It
has separate single-cell evidence and should be assessed in that context.

## Provisional timestep results

The Rush–Larsen plus Euler screen tested outer tissue steps of 10, 20, and
50 us, with 1, 5, 10, and 25 ionic substeps. It used a 1,500-cell 2D slab,
15 ms horizon, and checked finite traces, activated counts, and p95 activation
shifts against both matched-step and fine scalar references. Pointwise Vm
error was recorded but was not the gate because the upstroke moves in time.

| Model | Fastest passing setting | Runtime (s) | Speedup vs matched scalar RKF45 | Worst activation p95 (ms) |
|---|---|---:|---:|---:|
| AlievPanfilov | 4 us / 1 substep | 5.01 | 1.15x | 0.0341 |
| BuenoOrovio | 20 us / 1 substep | 1.48 | 1.40x | 0.0018 |
| Courtemanche | 10 us / 1 substep | 2.85 | 3.56x | 0.02744 |
| Gaur | 10 us / 1 substep | 3.01 | 13.88x | 0.07212 |
| Grandi | 10 us / 5 substeps | 5.50 | 4.15x | 0.07994 |
| PerisYague | 10 us / 1 substep | 2.67 | 3.54x | 0.0302 |
| Stewart | 10 us / 1 substep | 3.09 | 6.42x | 0.0495 |
| TNNP | 10 us / 1 substep | 3.04 | 8.49x | 0.0454 |
| TWorld | 10 us / 5 substeps | 8.14 | 5.75x | 0.04777 |
| ToRORd_dynCl | 10 us / 1 substep | 3.45 | 12.02x | 0.0442 |
| Trovato | 10 us / 1 substep | 3.48 | 5.15x | 0.04483 |

These are screening choices, not declared production settings. The 1-substep
10 us GPU cases for Grandi and TWorld failed; the tested 5-substep settings
passed the present activation screen. Longer APD, recovery, current, and
conduction-velocity checks are still required. These failures mark an observed
stability boundary in the timestep/substep sweep; by themselves they do not
indicate a Rush--Larsen implementation defect. At fixed tissue `deltaT`, more
ionic substeps primarily refine the ionic update. Changing tissue `deltaT`
also changes tissue integration and splitting error, so the whole sweep is a
practical stable/accurate settings study, not a pure ODE-only error measure.
In the separate extension,
TNNP and TWorld at 20 and 40 us failed the 0.1 ms activation gate against the
fine reference, so increasing tissue `deltaT` is not currently justified.

## SBDF2: what exists and what is being checked

The manufactured check now follows the backup's coarse template: fixed
640×640 2D mesh, final time 0.2, and `dt = 0.025, 0.0125, 0.00625, 0.003125,
0.0015625`. It uses implicit SBDF2 with backward `ddt(Vm)` and the scalar
manufactured model's fixed RKF45 settings. The predeclared criterion is
monotonic `Vm L2` reduction and observed order 1.8–2.3 for the four coarsest
levels. The finest point is reported but not included in the order gate because
it may reach the spatial floor. This is a test of the tissue SBDF2 method, not
a Rush–Larsen test.

Job 9362 completed the same ladder with the manufactured batched model on CPU
and CUDA. It doubled Euler substeps as tissue `deltaT` halved (25, 50, 100,
200, 400), so the ionic Euler error did not mask the second-order tissue
scheme. Both backends passed second-order L2 gates for `Vm`, `u1`, and `u2`.
Observed orders were `Vm`: 1.97109, 2.00964, 2.10618; `u1`: 1.91077,
1.96285, 2.00616; `u2`: 1.96748, 1.98566, 2.00410. CPU/GPU field and saved
state parity was below `2e-14`. Batched `Vm L2` errors were about 1.8–1.9%
above the scalar SBDF2 reference on the coarse levels. The job ran while
9355 shared the node, so its timings are correctness-only.

Job 9450 completed all 30 GPU cases in the fixed-Euler-substep SBDF2 ladder:
1, 5, 25, 100, 400, and 800 substeps at each tissue timestep. With fixed
substep counts, the first three observed `Vm` orders (`deltaT` 0.025 to
0.003125) were 1.601, 1.451, 1.307 for 1 substep; 1.861, 1.800, 1.718 for
5; 1.946, 1.960, 2.001 for 25; 1.964, 1.997, 2.079 for 100; 1.968, 2.006,
2.101 for 400; and 1.969, 2.008, 2.104 for 800. This supports the expected
effect: coarse fixed-substep Euler error lowers apparent order, while 25 or
more substeps recover the SBDF2 second-order trend over these levels. The
finest point was excluded because its apparent orders rise as the fixed-mesh
spatial floor is approached. Comparing final results at each `deltaT`, the
400-to-800-substep `Vm` RMSE was at most 1.45% of the 800-substep analytic
`Vm L2` error (maximum absolute RMSE `2.62e-7`). This suggests a practical
plateau by 400 substeps for this manufactured setup; 800 is a finest-run
reference, not an exact solution. Job 9355 shared the node, so sweep timings
are not performance measurements. Outputs are in
`/tmp/cardiac_manufactured_euler_substeps_20260929/`.

Separately, completed 15 ms, 1,500-cell TNNP/TWorld runs compared Godunov with
SBDF2. Those results show the method changes the solution somewhat, as
expected; they do not show a GPU bug. At 15 ms:

| Model | GPU Godunov vs GPU SBDF2: Vm RMSE / max (mV) | Activation p95 (ms) | Scalar Godunov vs scalar SBDF2: Vm RMSE / max (mV) | Activation p95 (ms) | Scalar SBDF2 / GPU SBDF2 runtime (s) |
|---|---:|---:|---:|---:|---:|
| TNNP | 0.355 / 4.838 | 0.0176 | 0.219 / 3.101 | 0.0108 | 63.22 / 27.64 |
| TWorld | 0.800 / 9.545 | 0.0228 | 0.424 / 5.253 | 0.0116 | 177.83 / 62.24 |

For the last column, the scalar model used RKF45 and the GPU model used
Rush–Larsen plus Euler; the ratios are speed measurements, not backend
equivalence tests. Job 9358 rebuilt once without CUDA and once with CUDA,
then compared CPU-batched and GPU-batched SBDF2 traces with the same model,
integrator, tissue step, and substeps. The TNNP direct parity comparison
failed the strict field gate; the job stopped before running the TWorld
comparison. A manual TWorld comparison against its saved CPU/GPU cases failed
as well. The comparator includes saved state variables, current, Vm, and
activation. This is an unresolved SBDF2 backend discrepancy; it does not
change the separate successful Godunov scalar/GPU comparison.

The manufactured CUDA smoke completed on an 8×8×8 mesh over 50 us. Scalar
and GPU Vm analytic errors matched at printed precision (`L2=5.34119e-7`),
but this case was spatially dominated and did not exercise SBDF2. Its near-zero
error is not evidence of temporal convergence.

Job 9361 now supplies CPU-batched/CUDA-batched convergence and backend-parity
evidence on the fixed 640×640 mesh at final `t=0.2`. It used 25 Euler ionic
substeps and `dt=0.025, 0.0125, 0.00625, 0.003125, 0.0015625`. Both backends
passed monotonic `Vm L2` decrease and the predeclared first-order range on
the four coarsest points. The direct CPU/GPU field and state comparison passed
at every level (worst reported max difference `1.9984e-14`). The scalar
reference report's Godunov/Euler `Vm L2` values are 0.00204545, 0.00104205,
0.000525733, 0.000263789, and 0.000131853. Job 9361's CPU/GPU-batched values
are 0.00205419, 0.00104630, 0.000527827, 0.000264828, and 0.000132370
(0.39–0.43% higher), with essentially identical first-order orders. This
checks analytic error and convergence against the scalar study; it does not
claim pointwise equality of scalar RKF45 and batched Euler trajectories.

## Time controls and benchmark interpretation

- Tissue `deltaT` advances the monodomain solver and is also the total ionic
  interval passed to the model. These controls are coupled.
- `batchedSubsteps` divides that ionic interval; it does not change the tissue
  step. The batched model computes `Im` after the final ionic substep; SBDF2
  also refreshes current after the Vm update. No independent current-update
  interval is present in these paths.
- Thus coarsening tissue `deltaT` changes both diffusion and ionic evolution.
  Holding `deltaT/substeps` constant controls the ionic interval but still
  changes tissue truncation error.
- The earlier matched TNNP control run measured 2 us / 5 substeps at
  64.74 s scalar versus 21.95 s GPU, and 4 us / 10 substeps at 46.21 s versus
  16.59 s GPU. The latter kept the effective ionic step at 0.4 us; relative to
  the 2 us tissue baseline its Vm RMSE was 0.234 mV and activation p95 shift
  was 0.0115 ms. This is evidence that tissue step remains consequential.
- The backup's 16.9x TNNP and 24.7x TWorld results used 420,000 cells for
  55 ms, GPU `deltaT=20 us` with 10 substeps, and scalar RKF45 at `deltaT=2 us`.
  They are not apples-to-apples with the matched-step table here.

## Completed jobs and timing caveat

Jobs 9355, 9356, 9358, 9359, 9361, 9362, and 9450 have completed. Job 9358 completed
all scalar and batched TNNP/TWorld Godunov/SBDF2 cases at `deltaT=2 us`, five
ionic substeps, and 15 ms. Its direct TNNP CPU-batched/GPU-batched SBDF2
comparison failed: at 15 ms the maximum differences were 0.006296 mV in Vm,
62.14 in ionic current, and 6.823 in the saved V state, far above the
`1e-10` gate. The job stopped before TWorld; a manual comparison of the saved
TWorld results also failed, with `Vm_max=1.154981e-2`,
`ionicCurrent_max=153.0751`, and saved `v` state `max=12.24162`. Activation
counts differed by one. The CPU folder's generic filenames contain “gpu” for
historical naming reasons; the CSV correctly marks those runs `completed-CPU`.

Job 9356's tuned 3D Niederer-geometry rerun completed using one MPI rank and
one RTX 4000 Ada. Both used 15 ms, tissue `deltaT=10 us`, and their respective
provisional slab settings (TNNP one substep; TWorld five). TNNP took 1077.94 s
scalar and 62.42 s GPU (17.27x); TWorld took 1830.54 s scalar and 192.21 s GPU
(9.52x). Final activation p95 shifts were 0.00030 ms and 0.00010 ms. These
matched-macro-step timings and short-horizon comparisons support those
per-model candidates on this geometry; the runs are not the published
Niederer benchmark reproduction and do not replace longer APD/state checks.

| Job | Latest state | Purpose |
|---:|---|---|
| 9355 | Completed, exit 0 | 48-rank CPU electromechanics run to 0.8 s with scalar Gaur in myocardium and BuenoOrovio in the embedded potential domain; 38,905 s wall time |
| 9356 | Completed | Tuned 3D TNNP/TWorld scalar/GPU cases with per-model substeps |
| 9358 | Completed; parity failed for both models | CPU-batched versus GPU-batched SBDF2 cases for TNNP and TWorld; TWorld comparison was run separately after the script stopped at TNNP |
| 9359 | Completed; acceptance passed | Coarse N=640 scalar manufactured SBDF2 temporal ladder; orders 1.96983, 2.00954, 2.10796 |
| 9361 | Completed; convergence and CPU/GPU parity passed | Matched CPU-batched/CUDA-batched manufactured Godunov/Euler ladder; analytic errors align within 0.43% with the backup scalar reference |
| 9362 | Completed; both order gates and CPU/GPU parity passed | Matched manufactured SBDF2/backward ladder; 25–400 Euler substeps, second-order `Vm/u1/u2` convergence |
| 9450 | Completed, 30/30 cases | CUDA SBDF2 temporal ladders at fixed Euler substeps 1, 5, 25, 100, 400, 800; Euler error lowers observed order at low substep counts and 400-to-800 changes are small in `Vm` relative to analytic error |

Job 9355 was a separate 48-rank CPU electromechanics run on xenosim, completed
to 0.8 s with exit status 0. Its configuration uses scalar Gaur in myocardium
and BuenoOrovio in the embedded potential domain; it does not validate the
batched or CUDA Gaur implementations. It overlapped CPU/GPU timing work, so
those timings may have node-level CPU contention. The completed 3D all-model
matrix above ran earlier and remains the cleaner speed comparison.
A single `nvidia-smi` snapshot showing 0% utilization is not evidence that
earlier CUDA cases fell back to CPU; their run logs explicitly record CUDA
device selection.

## Remaining work, in discussion order

1. **Finish SBDF2 evidence.** Investigate the TNNP and TWorld CPU/GPU-batched
   SBDF2 discrepancies from job 9358 and rerun both parity checks after the fix.
   The scalar manufactured SBDF2 ladder passed. Report backend state/current/Vm deltas separately
   from Godunov/SBDF2 differences. The manufactured convergence order answers
   whether the tissue scheme is second order; it does not qualify a cell
   model's Euler/Rush–Larsen ODE accuracy.
2. **Recheck tuned 3D timings without contention.** Repeat TNNP/TWorld one at a
   time when no other job is using the GPU. Keep one rank, mesh, end time,
   output settings, and each model's chosen `deltaT`/substeps matched between
   scalar and GPU. Use the fine scalar case for accuracy, not just the
   matched-step scalar case.
3. **Extend single-cell acceptance.** For all 11 models, compare full traces,
   peak/upstroke, APD90, relevant gates/ions, and integrated ionic current at
   a baseline and at selected coarser steps. Investigate the high-step Grandi
   and TWorld CPU-batched floating-point exceptions and the full-trace
   parameter/tissue assumptions.
4. **Extend slab acceptance.** Run longer propagation and recovery cases,
   compare activation maps, conduction velocity, APD, and current balance, and
   refine mesh separately from timestep. The current slab is an early
   propagation screen only.
5. **Close TNNP benchmark documentation.** The supplied current-main Niederer
   TNNP case has completed and passed its six-probe activation regression.
   Record any remaining paper-to-case differences if a full published
   reproduction claim is needed. Keep the all-model geometry comparisons
   separate from that TNNP benchmark result.
6. **Complete state and parameter audit.** Compare scalar and batch mappings
   for every state, parameter override, initialization option, tissue/sex
   selection, and current sign/unit. Current smoke coverage does not exercise
   every configurable branch.
7. **Profile the GPU.** Separate host/device transfers, ionic kernels,
   diffusion solve, synchronization, and output over increasing cell counts.
   Current timings are whole-solver elapsed time and do not yet identify the
   one-GPU bottleneck.

## Evidence locations

- 11-model Niederer-geometry matrix:
  `/tmp/cardiac_niederer_all_models_20260928/summary.csv`
- 12-model 2D slab matrix:
  `/tmp/cardiac_gpu_slab_all_20260927/summary.csv`
- Full single-cell GPU matrix:
  `/tmp/cardiac_gpu_singlecell_full_all_dt10us_20260927/summary.csv`
- Targeted Grandi/TWorld full traces:
  `/tmp/cardiac_gpu_singlecell_grandi_tworld_baseline_20260928/summary.csv`
- 11-model substep/timestep sweep:
  `/tmp/cardiac_gpu_integrator_sweep_20260928_nonfabbri/best_settings.md`
- SBDF2/Godunov comparison:
  `/tmp/cardiac_gpu_sbdf2_baseline_20260928/summary.csv`
- Manufactured smoke and failed-fineness convergence attempt:
  `/tmp/cardiac_manufactured_cuda_smoke_20260928/summary.csv` and
  `/tmp/cardiac_manufactured_convergence_20260928/summary.csv`
- Batched manufactured Godunov and SBDF2 convergence/parity:
  `/tmp/cardiac_manufactured_batched_temporal_20260928/` and
  `/tmp/cardiac_manufactured_batched_sbdf2_20260928/`
