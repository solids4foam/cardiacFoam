# GPU ionic-model validation: discussion brief

**Snapshot:** 2026-09-30, updated after completion of jobs 9510, 9512, 9513, 9516, and 9519
**Target:** `cardiacFoam` branch `main`, validation base commit `6900b946`;
the parity implementation was committed as `222c8bfc`.
**GPU:** NVIDIA RTX 4000 Ada, 20 GB; OpenFOAM 2412. The reported CUDA builds
used CUDA 11.5 and GCC 10 for nvcc. Most benchmark jobs used one MPI rank,
one CPU core, and one GPU.

**Current CPU-batched/CUDA-batched parity:** all 11 non-Fabbri models pass
the 15 ms, 1,500-cell 2D slab backend comparison at `deltaT=2 us` and five
Rush--Larsen substeps under both Godunov and SBDF2. The fresh four-model
rerun (BuenoOrovio, Courtemanche, PerisYague, ToRORd_dynCl) closes the four
strict field mismatches recorded by older jobs 9473/9475. See the final
section for the criteria and reports. This confirms backend parity at this
tested setting; it does not establish full-duration accuracy or safe larger
time steps.

Fabbri is outside that tissue-model matrix because it is an AV-node pacemaker,
but its CPU-batched/CUDA-batched single-cell parity has now been checked
separately (job 9519). With zero stimulus amplitude over 0–2 s, the final
second matched exactly at saved precision for 10,001 samples, including 33
nonzero rate columns and all directly compared Vm/state/current columns. The
CUDA run selected device 0. This does not validate AV-tissue coupling or
propagation.

## Equation-alignment changes after job 9481 started

The source now aligns the batched CPU and CUDA equations for reviewed backend
drifts while preserving the scalar CellML-style reference:

- BuenoOrovio's scalar evaluator remains unchanged with hard Heaviside
  thresholds. The CPU-batched hot path now calls the same width-`1e-4`
  smoothed batch evaluator as CUDA for rates, RL support, and `Jion`.
- Courtemanche's scalar CellML evaluator remains unchanged. Its batch
  evaluator now uses the same exact `V == -47.13 mV` special case and
  quotient expression as the scalar source, replacing the batch-only
  `1e-12` near-equality branch. CPU-batched and CUDA use this shared equation.
- PerisYague's batch `Ki` rate no longer includes the optional `IbK` term.
  The scalar rate and published 12-current `Iion` sum omit it, while local
  `gbK` is marked optional and defaults to zero. Any nonzero-`gbK` extension
  would need a consistent voltage-current and concentration-balance model.

The alignment rule is explicit: scalar versus batched differences in model
equations or indexing are defects to investigate and correct; CPU-batched and
CUDA-batched must use the same equations and RL/Euler state choices. The
BuenoOrovio smooth-Heaviside variant is the one documented scalar/batched
equation exception because smoothing is required for the batched solver's
stability. Its scalar CellML source remains the hard-Heaviside reference.

The original job 9481 and the follow-up plans for jobs 9492/9482 are historical
notes from the earlier checkpoint. The targeted equation-alignment changes
were later rebuilt and exercised in jobs 9510, 9512, and 9516 below. The
older 9473-9476 results remain useful as evidence of the original mismatches,
but their four flagged rows are superseded by the fresh 9516 parity run.

Job 9487 has dependency `afterany:9482`. It runs a
50 ms single-cell comparison for BuenoOrovio, Courtemanche, and PerisYague:
first a clean CPU build with scalar versus CPU-batched traces, then a clean
CUDA build with scalar versus CUDA-batched traces. Both use the same case
settings (`deltaT=2 us`, five RL/Euler substeps for batched, identical
stimulus, and all exported variables). The scalar solver remains RKF45, so
these are trajectory comparisons rather than one-step RHS equality tests.
Job 9482's matched CPU-batched/CUDA-batched slab comparison supplies the
direct backend parity check under the same batched integrator. At the latest
queue check, job 9492 was running on xenosim, 9482 was pending on 9492, and
9487 was pending on 9482.

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
- Job 9358's TNNP/TWorld SBDF2 mismatch was subsequently traced to a missing
  CUDA Vm-rate extrapolant path in the model wrappers. Job 9473 then passed
  the strict CPU-batched/CUDA-batched parity gate for seven models; four
  models retained differences at that point. Jobs 9510/9512 corrected the
  Courtemanche, ToRORd, BuenoOrovio, and PerisYague backend paths, and job
  9516 reran those four. All 11 now pass the same strict gate under both
  Godunov and SBDF2 at the documented 2 us / five-substep setting.
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
integrator, tissue step, and substeps. It initially failed for TNNP and
TWorld. The missing Vm-rate path was fixed; see the completed parity follow-up
below. The comparator includes saved state variables, current, Vm, and
activation. This is separate from the scalar RKF45 versus CUDA RL/Euler
comparison and the successful Godunov scalar/GPU comparison.

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

## CPU/GPU-batched SBDF2 follow-up and activation shifts

The missing SBDF2 Vm-rate extrapolant was a real backend integration defect.
All 11 non-Fabbri wrappers now call the same shared CUDA slice-update helper
before each model rate evaluation, including the final current refresh. TNNP
and TWorld previously fused that assignment into their model kernels; that
special path was removed so they now follow the same ordering and time inputs
as the other models. Job 9476 completed the short all-11 SBDF2
CPU-batched/CUDA-batched parity rerun after the final uniform refactor.
TNNP/TWorld still pass the strict gate. Seven models pass strict Vm/state/
current parity overall; the same four models have small residual field
differences. The helper adds a separate kernel launch per ionic substep and
current refresh, so performance has not yet been signed off.

A separate CUDA-kernel audit then found two missing Rush--Larsen state updates:
CPU-batched Courtemanche applies RL to `cajsr_v`, and CPU-batched PerisYague
applies RL to `ryr_v`; both CUDA kernels had left those states in their Euler
remainder. The CUDA lists and Euler exclusions have now been corrected. The
audit also found that CUDA applied RL without the CPU path's invalid-support
fallback. All non-Fabbri CUDA RL updates now fall back to Euler when `tau` is
not finite/positive above `VSMALL` or the steady state is non-finite. This
matches the CPU batched executor's decision rule. These CUDA edits are not yet
compiled or validated; replacement job 9492 will build the current kernels.

Job 9492 is running the shortened 0.08 s
SBDF2 comparison on both the 1,500-cell 2D slab and 52,500-cell Niederer
geometry. It uses `deltaT=10 us`, 25 Rush--Larsen/Euler substeps, one MPI
rank/CPU core, one visible RTX 4000 Ada, and 1 ms field output. It will run
CPU-batched and CUDA-batched versions of every non-Fabbri model, compare each
per-cell trace after activation-time alignment, and retain activation shift as
its own metric. The waveform limits are activation p95 <=0.1 ms, activated
mask mismatch <=max(5 cells, 1%), aligned Vm RMSE outside +/-2 ms of the
upstroke <=0.5 mV, p95 peak shift <=2 mV, and p95 APD90 shift <=max(2 ms,
2% of median CPU APD90). The strict componentwise parity diagnostics remain
reported separately. These jobs are correctness studies, not uncontended
performance measurements; the partition does not advertise a Slurm GPU GRES.
At the latest checkpoint (2026-09-29 17:35 Europe/Dublin), 9492 had completed
all 11 CPU slab cases (4,247 s total case runtime) and the CPU Niederer
AlievPanfilov case (2,345 s). CPU Niederer BuenoOrovio was at
`0.00588/0.08 s`. The CUDA build had not started. The previous 9481 run had
reached `0.38246/0.5 s` in its first CPU/slab AlievPanfilov case when
stopped.

Job 9482 is queued after 9492 for a short all-model
CPU-batched/CUDA-batched Godunov parity run. Job 9487 follows 9482 for the
targeted scalar comparisons. Both will use the current source if it remains
unchanged before they start.

Job 9476 used a 1,500-cell slab, 15 ms, `deltaT=2 us`, and five
Rush--Larsen/Euler substeps. Every model had the same activated-cell mask and
count on CPU and CUDA. The four strict-gate flags were:

| Model | Max Vm difference (mV) | Max current difference | Max state difference | Activation p95 shift (ms) |
|---|---:|---:|---:|---:|
| BuenoOrovio | 0.01579 | 0.07670 | 1.843e-4 | 5.896e-5 |
| Courtemanche | 2.547e-4 | 5.093e-3 | 2.546e-4 | 2.026e-7 |
| PerisYague | 4.213e-4 | 8.400e-5 | 4.213e-4 | 5.384e-9 |
| ToRORd_dynCl | 9.936e-5 | 5.060e-4 | 1.126e-4 | 5.013e-7 |

TNNP and TWorld both passed after switching them to the shared update path;
all other seven strict passes remained passes. These results show that the
refactor did not introduce an observable activation delay in the short run.

Job 9473 compared CPU-batched and CUDA-batched SBDF2 on the same 1,500-cell
2D slab for 15 ms (`deltaT=2 us`, five Rush--Larsen/Euler substeps), using
OpenFOAM 2412 and one MPI rank/CPU core with an RTX 4000 Ada (20 GB, CUDA
12.9). Output was written at 5, 10, and 15 ms. Seven models passed the
strict Vm/state/current gate. Four failed its `1e-10` componentwise gate,
but every model had identical activated-cell masks and counts. At 15 ms:

| Model | Activation p95 shift (ms) | Max activation shift (ms) | Max Vm difference (mV) |
|---|---:|---:|---:|
| AlievPanfilov | 1.91e-14 | 1.01e-13 | 9.02e-13 |
| BuenoOrovio | 5.896e-5 | 7.029e-5 | 0.01579 |
| Courtemanche | 2.026e-7 | 2.446e-7 | 0.0002547 |
| Gaur | 9.02e-14 | 1.70e-13 | 3.14e-12 |
| Grandi | 1.006e-13 | 3.00e-13 | 5.90e-12 |
| PerisYague | 5.384e-9 | 5.963e-9 | 0.0004213 |
| Stewart | 1.91e-14 | 1.01e-13 | 9.20e-12 |
| TNNP | 2.08e-14 | 1.01e-13 | 1.40e-12 |
| TWorld | 1.006e-13 | 1.006e-13 | 1.779e-11 |
| ToRORd_dynCl | 5.013e-7 | 5.232e-7 | 9.936e-5 |
| Trovato | 3.30e-14 | 1.01e-13 | 3.90e-12 |

The remaining strict-gate flags are BuenoOrovio, Courtemanche, PerisYague,
and ToRORd_dynCl. These activation shifts are sub-microsecond, far below a
millisecond-scale waveform lag. The maximum Vm differences are also small in
mV, though they are real and are not waived by the activation result. Job
9475 repeated the four comparisons with Godunov coupling. All four still
missed the same strict componentwise threshold, showing the residual backend
arithmetic differences are not introduced by SBDF2. Godunov activation p95
shifts were 0.000207 ms (BuenoOrovio), 2.00e-7 ms (Courtemanche), 5.39e-9 ms
(PerisYague), and 0.000735 ms (ToRORd_dynCl). SBDF2 shifts were no larger;
for ToRORd_dynCl they were much smaller.

The user's observation about BuenoOrovio's lag applies to the smoothed
batched regularized variant versus the unchanged scalar hard-Heaviside
reference. CPU-batched and CUDA-batched now use the same narrow smoothing;
scalar-versus-batched differences should still be measured and reported.
Earlier runs predate this alignment and cannot quantify its effect. An
upstroke pointwise error alone is not a rejection criterion; compare
activation time, peak, APD, and aligned trace error outside the upstroke. The
present sparse 5 ms outputs do not provide adequate APD or fine-grained
waveform evidence.

These runs establish close SBDF2 backend activation parity for all 11
non-Fabbri models, and strict full-field parity for seven. They do not close
the four residual componentwise differences, validate scalar-vs-batched
waveform accuracy at full action-potential resolution, or establish SBDF2
parity on the 3D Niederer geometry.

## Completed jobs and timing caveat

Jobs 9355, 9356, 9358, 9359, 9361, 9362, 9450, 9473, and 9475 have completed. Job 9358 completed
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
| 9358 | Completed; initial parity failure, superseded | CPU-batched versus GPU-batched SBDF2 cases for TNNP and TWorld; follow-up identified missing Vm-rate extrapolation |
| 9359 | Completed; acceptance passed | Coarse N=640 scalar manufactured SBDF2 temporal ladder; orders 1.96983, 2.00954, 2.10796 |
| 9361 | Completed; convergence and CPU/GPU parity passed | Matched CPU-batched/CUDA-batched manufactured Godunov/Euler ladder; analytic errors align within 0.43% with the backup scalar reference |
| 9362 | Completed; both order gates and CPU/GPU parity passed | Matched manufactured SBDF2/backward ladder; 25–400 Euler substeps, second-order `Vm/u1/u2` convergence |
| 9450 | Completed, 30/30 cases | CUDA SBDF2 temporal ladders at fixed Euler substeps 1, 5, 25, 100, 400, 800; Euler error lowers observed order at low substep counts and 400-to-800 changes are small in `Vm` relative to analytic error |
| 9473 | Completed; historical: 7 strict passes, 4 strict flags | Initial all-model SBDF2 backend matrix before equation alignment |
| 9475 | Completed; historical: 4 strict flags | Initial Godunov controls before equation alignment |
| 9512 | Completed; PASS | Four-model single-cell CPU/CUDA trace and rates parity after equation alignment |
| 9513 | CPU matrix completed; first CUDA attempt stopped at startup I/O | Four-model, two-scheme slab setup; the lazy-allocation rates-copy edge case led to the follow-up |
| 9516 | Completed; 8/8 PASS | Fresh CUDA Godunov/SBDF2 slabs compared with the CPU results from 9513 |

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

1. **Extend physiological accuracy checks.** The strict CPU-batched/CUDA-
   batched field and activation parity differences for the four models flagged
   in jobs 9473/9475 were resolved by the equation-path corrections and passed
   in job 9516 under both coupling schemes. Full APD/recovery, longer slab
   output, and model-specific accuracy versus scalar RKF45 remain open;
   backend equality alone does not qualify a production timestep or the
   Euler/Rush--Larsen ODE error.
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

## Fresh single-cell scalar/CPU/CUDA parity checkpoint (2026-09-29)

Jobs 9504 and 9505 ran the same single-cell check for all 11 non-Fabbri
models. Each model used a 50 ms case, `deltaT=2 us`, five Rush--Larsen
substeps, and 0.1 ms output spacing. Each build compared the scalar RKF45
reference against CPU-batched and CUDA-batched traces. The CUDA cases all log
device 0 selection. The runs used one cell and are correctness checks, not
GPU performance tests. Fabbri remains excluded as the AV-node pacemaker case.

The direct CPU-batched versus CUDA-batched trace comparison has 500 common
samples per model. At the saved seven-decimal precision, Vm, mapped ODE states,
and total ionic current are identical for BuenoOrovio, PerisYague,
AlievPanfilov, Gaur, Grandi, Stewart, TNNP, TWorld, and Trovato. The short
AlievPanfilov case has no -30 mV crossing on either backend.

| Model | CPU/CUDA Vm max difference (mV) | Activation-time difference (ms) | CPU/CUDA total-current max difference | Note |
|---|---:|---:|---:|---|
| BuenoOrovio | 0 | 0 | 0 | Direct backend traces identical at output precision; scalar uses hard Heaviside while batch uses the documented smooth transition |
| Courtemanche | 1.683e-4 | 4.748e-7 | 5.259e-4 | Ki max difference 0.002275; activation-aligned Vm RMSE outside +/-2 ms of upstroke is 9.482e-5 mV |
| PerisYague | 0 | 0 | 0 | Direct backend traces identical at output precision |
| AlievPanfilov | 0 | no crossing | 0 | CPU/CUDA traces identical at output precision |
| Gaur | 0 | 0 | 0 | CPU/CUDA traces identical at output precision |
| Grandi | 0 | 0 | 0 | CPU/CUDA traces identical at output precision |
| Stewart | 0 | 0 | 0 | CPU/CUDA traces identical at output precision |
| TNNP | 0 | 0 | 0 | CPU/CUDA traces identical at output precision |
| TWorld | 0 | 0 | 0 | CPU/CUDA traces identical at output precision |
| ToRORd_dynCl | 0.6302 | 0.003250 | 4.414 | Raw Vm difference peaks during upstroke; after activation alignment Vm RMSE outside +/-2 ms of upstroke is 0.002948 mV; max mapped-state difference is 0.0101 (`INaL_mL`) |
| Trovato | 0 | 0 | 0 | CPU/CUDA traces identical at output precision |

ToRORd's raw upstroke difference is primarily a small timing displacement:
the interpolated CPU/CUDA activation times are 21.997813 and 22.001063 ms,
respectively (3.25 us apart). The case ends at 50 ms, before a reliable APD90
recovery measurement for the long-duration models. These data support close
single-cell backend trajectory parity at this setting; they do not replace
the longer slab, APD, or time-step sensitivity checks.

The scalar-reference comparison is separate from direct CPU/CUDA parity.
For example, Gaur has Vm RMSE 0.06676 mV and max difference 0.8446 mV against
the scalar RKF45 trace on both batched backends. ToRORd's CPU/CUDA scalar
comparison metrics differ modestly in this short case, consistent with its
3.25 us backend activation displacement. APD90 is not available in these
50 ms traces.

**Historical diagnostic output limitation (fixed in job 9512):** this earlier
matrix exported zero CUDA `RATES_*` columns because there was no device-to-host
rates copy. Job 9512 added the requested-rate copy and confirmed nonzero,
CPU-matching rates for BuenoOrovio, Courtemanche, PerisYague, and ToRORd_dynCl.
Job 9516 then exposed and fixed a startup-order edge case: before the first
ionic update, the tissue writer can request rates before lazy CUDA allocation.
The copy is now skipped until CUDA buffers exist. See the final validation
section for the fresh tests.

Results are preserved at
`/tmp/cardiac_equation_alignment_singlecell_9504/` (BuenoOrovio,
Courtemanche, PerisYague) and
`/tmp/cardiac_equation_alignment_singlecell_9505/` (the other eight models).
Job 9504 also required a CUDA macro-continuation correction. This was a
compile-only syntax fix; it did not change model equations.

## Courtemanche and ToRORd backend-path correction (2026-09-30)

The CPU-batched hot paths for Courtemanche and ToRORd were still evaluating
their scalar generated CellML functions, while CUDA evaluated their separate
`*ComputeVariablesBatch` functions. A source audit found no substantive
Courtemanche equation mismatch. The ToRORd batch evaluator had one real typo:
the `a3` numerator used `nao/Knao` where the scalar reference uses `ko/Kko`.
The batch expression was corrected, and the CPU-batched per-cell hot paths for
both models now call the same host/device batch evaluators used by CUDA. The
scalar generated reference paths were left unchanged.

Job 9510 rebuilt CPU and CUDA libraries and ran matched 50 ms single-cell cases
at `deltaT=2 us`, five Rush--Larsen substeps, and 0.1 ms output spacing. Both
backends selected the CUDA device for the device cases. At the saved trace
precision, CPU-batched and CUDA-batched voltage, all mapped state variables,
and total ionic current are identical for both models (500 samples each). The
scalar-reference comparison metrics are also identical across backends:

| Model | Scalar-reference Vm RMSE (mV) | Max (mV) | Activation shift (ms) | Off-upstroke RMSE (mV) | Integrated current error |
|---|---:|---:|---:|---:|---:|
| Courtemanche | 0.021436 | 0.304321 | 0.001137 | 0.005925 | 1.829e-5 |
| ToRORd_dynCl | 0.090211 | 1.465961 | 0.007344 | 0.006952 | 4.781e-4 |

These scalar comparisons use RKF45 as the scalar reference and do not imply
that the scalar and batched integrators should agree exactly. At this stage,
CPU-batched/CUDA-batched single-cell state/current parity passed at the
recorded precision. The rates-export limitation noted in the original run was
fixed and verified in job 9512 below. ToRORd's auxiliary `x1`--`x4`
algebraic intermediates showed small absolute differences; direct state and
current columns matched at the trace precision.

Job 9511 was the original two-model slab follow-up. The expanded four-model
CPU matrix completed in job 9513; its first GPU attempt exposed a startup
rates-copy ordering issue, fixed before the complete GPU retry in job 9516.
The final 8/8 slab comparisons passed; see the completed matrix below.

Single-cell traces and build logs are in
`/tmp/cardiac_equation_alignment_singlecell_9510/`. The focused slab job writes
to `/tmp/cardiac_court_torord_slab_parity_9511/`.

### Fresh four-model parity and CUDA rates export (2026-09-30)

Job 9512 rebuilt the CPU and CUDA libraries and reran 50 ms single-cell cases
for BuenoOrovio, Courtemanche, PerisYague, and ToRORd_dynCl. Each used
`deltaT=2 us`, five Rush--Larsen substeps, and 0.1 ms output spacing. CUDA
single-cell logs confirm device 0 was selected. Direct CPU-batched/CUDA-batched
comparison at the saved seven-decimal precision is exact for Vm, mapped
states, ionic-current columns, and rates in all four models.

The job also confirmed the diagnosis of the zero CUDA rates export: the
`.dat` writer did request `RATES_*` for full-variable output, but no device to
host rates copy existed. `prepareIOAccess` now copies rates when the requested
selection includes rates (including full trace selection). A fresh CUDA build
and output check passed for all four models. CPU and CUDA rates were nonzero
and identical at saved precision:

| Model | Rate columns | Nonzero CPU/CUDA entries | Max CPU/CUDA difference | Direct state/current max difference |
|---|---:|---:|---:|---:|
| BuenoOrovio | 4 | 1341 / 1341 | 0 | 0 |
| Courtemanche | 21 | 9209 / 9209 | 0 | 0 |
| PerisYague | 22 | 10281 / 10281 | 0 | 0 |
| ToRORd_dynCl | 45 | 18011 / 18011 | 0 | 0 |

These are single-cell correctness checks, not speed measurements. Results and
build logs are under `/tmp/cardiac_equation_alignment_singlecell_9512/`;
`compare_singlecell_rates.py` reproduces the exported-rate and direct-column
comparison.

Jobs 9513 and 9516 completed the matching four-model slab matrix. Job 9513
ran the CPU-batched cases; its initial CUDA attempt stopped before ionic
updates because the writer requested rate output before lazy device
allocation. The CUDA rate-copy guard now skips the copy until allocation is
complete, preserving initialized zero rates at startup. Job 9516 rebuilt and
ran the eight GPU cases against the CPU results. All eight reports passed for
BuenoOrovio, Courtemanche, PerisYague, and ToRORd_dynCl under both Godunov
and SBDF2 coupling.

Each case used 1,500 cells, 15 ms, `deltaT=2 us`, and five Rush--Larsen
substeps. At 5, 10, and 15 ms, the acceptance limits were maximum absolute
Vm/state difference `<=1e-10`, current difference `<=1e-10 + 1e-12 *
max(abs(reference current))`, and activation difference `<=1e-10 s`. Across
the eight reports the largest Vm, state, and current differences were
`1.11e-14 mV`, `1.11e-11`, and `5.70e-11`; the largest activation p95 was
`1.01e-16 s`. The active-cell counts matched for each CPU/CUDA comparison.
CUDA logs selected device 0. These were one-rank correctness runs, not
performance measurements.

The retry data and reports are under
`/tmp/cardiac_changed_models_gpu_retry_9516/`; the CPU baselines are under
`/tmp/cardiac_court_torord_slab_parity_9513/cpu/`. The complete reproducible
CPU/CUDA matrix is launched from the repository root with
`sbatch tutorials/gpuValidation/run_courtemanche_torord_slab_backend_parity.slurm`.

## Fabbri AV-node single-cell CPU/CUDA parity (job 9519)

Fabbri has a batched CPU implementation (`FabbriBatched`) and an actual CUDA
implementation (`FabbricompactBatched`). It was excluded from tissue slab
parity because it is an AV-node pacemaker; this focused run validates only the
single-cell backend path.

The run used stimulus amplitude zero, 2 s duration, `deltaT=2 us`, five
Rush--Larsen substeps, and all-variable output every 0.1 ms. CPU-batched and
CUDA-batched traces from 1 to 2 s had 10,001 matching samples. All 33
`RATES_*` columns and directly compared Vm/state/current columns were exactly
equal at saved precision (maximum absolute difference zero); both had 289,628
nonzero rate entries. The CUDA log confirms device 0. Scalar RKF45 and
CPU-batched APD90 were 150.05343 and 150.05470 ms, respectively, a 0.00127 ms
difference. The scalar comparison is contextual because it uses a different
integrator than the batched Rush--Larsen/Euler path.

The simulation logs report about 27 s for scalar CPU, 19 s for CPU-batched,
19 s for scalar in the CUDA build, and 505 s for CUDA-batched. This single-cell
CUDA runtime is dominated by kernel-launch overhead and diagnostic output; it
is not a performance benchmark. Logs, traces, and summaries are under
`/tmp/cardiac_fabbri_singlecell_parity_9519/`. Reproduce from the repository
root on xenosim with:

```bash
sbatch tutorials/gpuValidation/run_fabbri_singlecell_backend_parity.slurm
```
