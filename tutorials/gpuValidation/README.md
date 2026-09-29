# GPU ionic-model implementation and validation

For the concise per-model results, job state, limitations, and open work, see
[`VALIDATION_STATUS_2026-09-28.md`](VALIDATION_STATUS_2026-09-28.md).

This guide is the reproducible validation ladder for the CUDA ionic ODE
backend. The myocardium tissue PDE currently runs on the CPU; whole-solver
timings therefore measure GPU ODE integration coupled to CPU PDE solves and
the transfers between them. They are not GPU-PDE timings.

## Validation ladder and conclusions

Run the checks in this order, preserving each stage's reference and
acceptance criteria:

1. **Build and require a CUDA device.** Use the build commands below on a
   node with a visible NVIDIA GPU. Matrix launchers use `--require-gpu` so a
   CPU fallback cannot be counted as a CUDA result.
2. **Single-cell traces.** Compare the scalar reference and CUDA-batched
   trace, saved states, and currents. The single-cell runner covers all
   registered cases; Fabbri is an AN-node model and is assessed separately
   from the tissue-model matrix.
3. **2D slab and integration controls.** Run finite-trace and propagation
   screens, then vary tissue `deltaT` and ionic substeps. Use Rush–Larsen for
   mapped gates and Euler for other states. Compare against matched-step and
   fine scalar RKF45 references; do not treat a fast, unstable, or inaccurate
   setting as a pass.
4. **3D geometry and benchmark.** Run all in-scope models in the
   Niederer-shaped geometry for cross-model comparison. Treat published
   Niederer reproduction as a separate TNNP check under the specified
   benchmark conditions.
5. **Manufactured convergence and batch parity.** Use analytic errors and
   timestep refinement to verify the manufactured batched/CUDA field path.
   The Godunov/Euler and SBDF2/backward ladders pass their respective
   convergence checks; CPU-batched and CUDA-batched manufactured fields and
   states agree within `2e-14`. The fixed-substep Euler sweep shows the
   expected order reduction at low substep counts and a small 400-to-800
   `Vm` difference relative to analytic error.
6. **Physiological SBDF2 device parity.** Compare CPU-batched and
   CUDA-batched TNNP/TWorld using identical model settings, integration,
   mesh, and timestep. This stage currently fails its strict parity check and
   remains open; manufactured parity does not close this model-specific
   issue.

The evidence supports working CUDA ODE paths and useful speedups in the
tested 3D workloads, but it does not certify every model at production
timesteps. The CPU PDE and ODE/PDE transfers remain part of the measured
runtime. See the status report for model-by-model results and known failures.

## Repository comparison and changes

The target was `cardiacFoam` `main` at `39333ef3`. The historical GPU port was
inspected at `d4f9ed3e` in `cardiacFoAM_GPU_backup`. The backup working tree
also contains later local edits, so individual changes were checked against
the port and the current source. Current `main` already contained scalar,
batched CPU, and CUDA classes for all 12 single-cell model families. The
backup has an untracked `tutorials/benchmarkGPU` collection, including
Niederer-geometry runs for TNNP, TWorld, and BuenoOrovio; that collection
has not been copied as a published benchmark claim. Current `main` retains
its newer TWorld, Trovato, and ToRORd model revisions and file names. Their
batch evaluator headers compare identically to the historical port after
normalizing year labels, the added author line, and C math qualifiers;
current wrapper initialization and parameter handling still need separate
validation.

Targeted fixes in the target repository:

- All 12 CUDA launchers now upload the tissue voltage state slice when the
  full host state is not dirty. The SoA slice comes directly from the current
  host state array, so device-resident non-voltage states are preserved.
- CUDA launchers now choose the requested `euler` or `rushLarsen` kernel.
  Before this fix they always chose Rush–Larsen, including when the default
  `batchedIntegrator euler` was requested.
- The shared host hot path now dispatches to the model's hot-path support
  calculation for models without a cell-specific override. Previously it
  treated ordinary algebraic slots as Rush–Larsen time constants and ionic
  current slots. The eight non-heterogeneous model families used this path.
- The CUDA build converts OpenFOAM's `-iquote` include flag for nvcc. CUDA
  launch checks call `std::abort`, and device math calls use global C math
  functions where needed, notably in Gaur and TWorld.

These changes preserve current model names, state indices, and equations.
The scalar and batched evaluators include the same model-specific `Names.H`
files, but this alone is not a complete state and parameter mapping audit.

## Build and hardware

An initial build against CUDA 13.1 on the xenosim login node linked
`libcudart.so.13`, but that node has no visible NVIDIA device; its runs used
the CPU fallback and are not GPU evidence. The validation runs below used
the allocated xenosim GPU node, an NVIDIA RTX 4000 Ada, OpenFOAM v2412,
CUDA 11.5, and GCC 10 as nvcc's host compiler. The CUDA runtime log
confirmed the GPU backend in each batched run. `NVARCH=75` was the tested
build setting; it is recorded here as used, not as a recommendation for a
different GPU.

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
export CARDIAC_ENABLE_CUDA=1
export CUDA_HOME=/usr
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
export NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
(cd src/ionicModels && wclean libso)
(cd src/ionicModels && wmake libso)
```

Run cases on a node with a visible NVIDIA device. Set `FOAM_SIGFPE=false`
for the CUDA batch while diagnosing; a CPU batched case can otherwise stop
at an invalid intermediate state with a floating-point exception. Use
`--require-gpu` in the supplied matrix scripts so a successful CPU fallback
cannot be counted as a CUDA pass.

## Time controls observed in code

For batched models, each call to `solveODE(t0, deltaT, ...)` converts
`deltaT` to model time and divides it by `batchedSubsteps`. The batched host
path evaluates `Im` after the final ionic substep. The CUDA launchers also
evaluate and download `Im` after their substeps. The myocardium domain calls
the ionic solve once per tissue advance, using the solver's current
`deltaT`. In SBDF2 coupling, it refreshes ionic current after updating tissue
`Vm`; that extra evaluation must be considered when comparing modes. The
explicit tissue algorithm can cap `deltaT` using `maxCo`. No independent
`Im` update interval was found in these paths.

Consequently, changing `deltaT` changes both the tissue step and the ionic
interval. Changing `batchedSubsteps` changes the ionic substep while leaving
the tissue step fixed. The control study below measures those effects; it
does not establish a globally safe larger step.

On xenosim (RTX 4000 Ada, one MPI rank/CPU core, 1,500-cell TNNP slab, 15 ms, scalar RKF45 versus
CUDA Rush--Larsen gates plus Euler for the other states), the measured
single-process wall times and field comparisons were:

| Tissue `deltaT` | Ionic substeps | Effective ionic step | Scalar/GPU runtime (s) | GPU vs scalar Vm RMSE / max (mV) | p95 activation shift (ms) |
|---:|---:|---:|---:|---:|---:|
| 2 us | 5 | 0.4 us | 64.74 / 21.95 | 0.000193 / 0.0020 | 0.000010 |
| 2 us | 1 | 2 us | 69.78 / 14.20 | 0.000959 / 0.0099 | 0.000100 |
| 4 us | 10 | 0.4 us | 46.21 / 16.59 | 0.000192 / 0.0019 | 0.000010 |

All three cases activated the same 622 cells. One ionic step per 2 us tissue
step remained stable for this TNNP case and was about 1.5x faster than its
five-substep run, with a small increase in error. Doubling tissue `deltaT`
while holding the effective ionic step at 0.4 us kept the GPU close to its
matched scalar result, but relative to the 2 us tissue baseline its Vm RMSE
was 0.234376 mV, maximum difference 3.375 mV, and p95 activation shift
0.0115 ms. The changed tissue discretization is therefore a material part of
the result even when ionic substep size is held constant. Timings are full
case wall time, not isolated ionic-kernel timing. See
`/tmp/cardiac_gpu_time_control_tnnp_20260928/{summary.csv,comparisons.txt}`
for the captured records.

To find the fastest acceptable GPU setting per in-scope model, the primary
slab sweep excludes Fabbri and uses Rush--Larsen for mapped gates plus Euler
for the remaining states. It
varies tissue `deltaT` at 1x, 2x, and 5x each model's baseline and tests one,
five, ten, and 25 ionic substeps. Every candidate is compared with scalar
RKF45 at the matching tissue step, the fine scalar reference, and the fine
GPU Rush--Larsen baseline. Since `deltaT` controls both tissue advancement
and the total ionic interval, these comparisons separate ionic integration
error from the combined larger-step effect. Full Euler is available as an
optional control (`--integrators euler rushLarsen`), but it is not the
optimization target. Failed or nonfinite cases do not count as speedups.
The 15 ms sweep screens stability and early propagation; follow it with the
full action-potential and slab acceptance checks before selecting a
production setting. The sweep writes a provisional fastest-setting table to
`best_settings.md`.

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
python3 tutorials/gpuValidation/run_slab_integrator_study.py \
    /tmp/cardiac_gpu_integrator_sweep --require-gpu
```

This study is intended to identify which combinations are both stable and
accurate before using a larger tissue step. It does not treat a faster but
unstable or inaccurate run as a successful setting.
For example, Grandi and TWorld failing at 10 us / one substep while passing
the current screen at five substeps records a measured limit in this sweep;
that observation alone is not evidence of a Rush--Larsen implementation bug.
At fixed tissue `deltaT`, changing ionic substeps mainly probes ionic update
resolution. Changing `deltaT` also changes tissue integration and splitting
error, so this is a practical settings envelope rather than a pure ODE-only
error study. The later SBDF2 comparison tests a different tissue coupling
scheme and must be judged against its own matched and fine references.

## Godunov versus SBDF2 coupling

The separate SBDF2 baseline uses TNNP and TWorld on the 2D slab, with
`deltaT=2e-6 s`, Rush--Larsen plus Euler, and five ionic substeps on the
batched path. Within each backend, Godunov and SBDF2 use the same ionic
method. Across scalar and CUDA backends, however, the ionic methods differ:
the scalar reference uses RKF45 and the CUDA batched case uses
Rush--Larsen/Euler. Therefore the original scalar/GPU SBDF2 comparison does
not isolate device parity. A follow-up rebuilds once without CUDA and once
with CUDA, then compares CPU-batched and GPU-batched SBDF2 at identical
settings:

```bash
sbatch tutorials/gpuValidation/run_sbdf2_batched_parity.slurm
```

The job waits until the current 3D tuned comparison finishes, so CPU and GPU
library rebuilds do not overlap that run. Results are written under
`/tmp/cardiac_gpu_sbdf2_batched_parity_20260928/`. It writes fields at
15-digit precision and compares `Vm`, `ionicCurrent`, activation times, and
every saved state variable. The parity gate is a maximum absolute difference
of `1e-10` for those fields; activation counts must also match. This is a
backend-equivalence gate, separate from the expected difference between
Godunov and SBDF2.

The original run rebuilds `libionicModels.so` with CUDA enabled so the SBDF2
current refresh includes the host synchronization fix in
`batchedIonicModel.H`.

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
export CARDIAC_ENABLE_CUDA=1 CUDA_HOME=/usr NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10 FOAM_SIGFPE=false
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
(cd src/ionicModels && wmake libso)
python3 tutorials/gpuValidation/run_sbdf2_baseline.py \
    /tmp/cardiac_gpu_sbdf2_baseline \
    --models TNNP TWorld --delta-t 2e-6 --substeps 5 \
    --end-time 0.015 --require-gpu
```

It compares scalar Godunov versus SBDF2, CUDA Godunov versus SBDF2, and
scalar versus CUDA for each coupling. Those comparisons characterize the
scheme/backend combinations but do not replace direct CPU-batched versus
CUDA-batched parity. SBDF2 is not considered GPU-validated until that direct
backend comparison passes. The refresh path currently synchronizes device state to host
and computes the refreshed current in a host loop; its correctness and
runtime cost are both part of this baseline.

Both models have completed. Godunov scalar versus CUDA at 15 ms had `Vm`
RMSE/max differences of 0.000193/0.0020 mV for TNNP and 0.000869/0.0106 mV
for TWorld, with matching activation counts. Full-case scalar/GPU runtimes
were 59.28/21.22 s (TNNP) and 164.03/38.84 s (TWorld).

The direct CPU-batched versus CUDA-batched SBDF2 comparison (job 9358) failed
the strict `1e-10` parity gate for both models. At 15 ms, TNNP had
`Vm_max=6.296152e-3`, `ionicCurrent_max=62.13696`, and saved `V` state
`max=6.822647`; TWorld had `Vm_max=1.154981e-2`,
`ionicCurrent_max=153.0751`, and saved `v` state `max=12.24162`. Activation
counts differed by one cell in each case. These are matched CPU/GPU batched
SBDF2 runs and expose an unresolved discrepancy. They are distinct from the
scalar RKF45 versus CUDA RL/Euler comparison below. The original job stopped
after TNNP failed; the saved TWorld cases were compared separately. The Slurm
script now records both model comparisons before returning failure.

At 15 ms, scalar-RKF45 versus CUDA Rush--Larsen/Euler SBDF2 had `Vm`
RMSE/max differences of 0.572/7.785 mV for TNNP and 1.220/14.568 mV for
TWorld, with 1 and 2 activated-cell differences. Since that pairing changes
both the ionic integrator and backend, it does not establish a GPU mismatch.
The Godunov-versus-SBDF2 comparison changes tissue integration and current
refresh timing; different traces there are expected and do not isolate a GPU
error.

The activation metric is computed from the solver-written `activationTime`
field at 5, 10, and 15 ms. It reports activated-cell counts (`activationTime
> 0`) and the 95th percentile of absolute activation-time differences among
cells activated in both cases. The percentile does not include cells that
activated in only one case; their count difference is reported separately.

## Single-cell screen

Run the repository-native screen from the repository root:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
python3 tutorials/gpuValidation/run_single_cell_smoke.py \
    --output-dir /tmp/cardiac-single-cell-screen
```

The script refuses to overwrite an output directory. It writes each
dictionary, log, trace, and `summary.csv`. Use `--require-gpu` on a GPU node;
the summary will flag any batched case that fell back to CPU. The default
screen uses 50 ms, `deltaT=2e-6 s`, five batched ionic substeps, scalar
RKF45, and batched Rush–Larsen. The script records each model's tissue and
stimulus. The scalar and batched cases share those settings. The script
requires Python 3 and NumPy.

Results below are **CPU batched versus scalar**, from this short screen.
RMSE and maximum difference are in mV. Activation shift uses the first
upward crossing of -30 mV; `n/a` means no such crossing in 50 ms. Outside
RMSE excludes 2 ms on either side of the scalar activation crossing.

| Model | Scalar | Batched CPU | GPU | RMSE | Maximum | Outside RMSE | Activation shift (ms) |
|---|---|---|---|---:|---:|---:|---:|
| AlievPanfilov | complete | complete | untested | 0.00001 | 0.00001 | 0.00001 | n/a |
| BuenoOrovio | complete | complete | untested | 0.01503 | 0.25300 | 0.00120 | 0.00056 |
| Courtemanche | complete | complete | untested | 0.02143 | 0.30682 | 0.00608 | 0.00147 |
| Fabbri | complete | complete | untested | 0.00001 | 0.00001 | 0.00001 | n/a |
| Gaur | complete | complete | untested | 0.06757 | 1.06841 | 0.03578 | 0.00281 |
| Grandi | complete | complete | untested | 0.00707 | 0.08083 | 0.00245 | 0.00077 |
| PerisYague | complete | complete | untested | 0.02503 | 0.41525 | 0.00038 | 0.00092 |
| Stewart | complete | complete | untested | 0.01391 | 0.19725 | 0.00750 | 0.00055 |
| TNNP | complete | complete | untested | 0.03945 | 0.71547 | 0.00247 | 0.00205 |
| TWorld | complete | complete | untested | 0.04103 | 0.84504 | 0.00427 | 0.00157 |
| ToRORd_dynCl | complete | complete | untested | 0.09194 | 1.64828 | 0.00680 | 0.00503 |
| Trovato | complete | complete | untested | 0.02351 | 0.47690 | 0.00522 | 0.00112 |

AlievPanfilov's stimulus starts after this 50 ms screen in physical time.
A separate 0.5 s stimulated comparison at `deltaT=1e-5 s` and one ionic
substep gave peak times 0.29134 s (scalar) and 0.29133 s (batched
Rush–Larsen), with 0.01614 mV voltage RMSE. Fabbri was screened without
an imposed stimulus; a 60-unit generic stimulus produced an unsuitable
pacemaker case. The other screen stimuli are provisional and do not yet
establish physiological plausibility or full action-potential duration.

A longer 0.5 s run with the same `deltaT` and substeps completed for all
12 CPU pairs. APD90 below is measured from the first -30 mV upward crossing
to 90% repolarization relative to the pre-stimulus voltage. AlievPanfilov
needed a 1 s run to repolarize. The separate `--all-variables` mode writes
every exported state, rate, and algebraic at 0.1 ms intervals. All 373
declared ODE states appeared under matching scalar/batched field names;
all exported samples were finite. The script records APD90 and normalized
integrated total-current error in `summary.csv`, and per-field errors in
`variable_summary.csv`. The largest state difference below is
scaled by that state's scalar trace range or maximum magnitude, whichever
is larger. Current error is the absolute difference of time-integrated
total ionic current divided by the integral of the absolute scalar current
(`Jion` for BuenoOrovio, `Iion_cm` for the others).

```bash
python3 tutorials/gpuValidation/run_single_cell_smoke.py \
    --output-dir /tmp/cardiac-all-variables --end-time 0.5 --all-variables
python3 tutorials/gpuValidation/run_single_cell_smoke.py \
    --output-dir /tmp/cardiac-aliev-long --models AlievPanfilov --end-time 1.0
```

| Model | ODE states | Scalar APD90 (ms) | APD90 difference (ms) | Largest scaled state difference | Integrated current error |
|---|---:|---:|---:|---:|---:|
| AlievPanfilov | 2 | 304.78 | -0.003 | 0.002% | 0.0007% |
| BuenoOrovio | 4 | 272.47 | 0.003 | 0.045% | 0.0083% |
| Courtemanche | 21 | 243.01 | -0.031 | 0.44% | 0.0004% |
| Fabbri | 33 | 150.05 | 0.001 | 1.37% | 0.0002% |
| Gaur | 29 | 180.98 | 0.037 | 1.00% | 0.0438% |
| Grandi | 41 | 326.55 | -0.009 | 0.16% | 0.0004% |
| PerisYague | 22 | 153.60 | 0.001 | 0.36% | 0.0004% |
| Stewart | 20 | 295.54 | -0.014 | 0.67% | 0.0053% |
| TNNP | 17 | 275.53 | -0.004 | 0.86% | 0.0350% |
| TWorld | 93 | 211.30 | -0.004 | 0.65% | 0.1184% |
| ToRORd_dynCl | 45 | 236.25 | -0.001 | 2.44% | 0.0252% |
| Trovato | 46 | 306.76 | -0.010 | 1.07% | 0.0973% |

The larger instantaneous state and current differences occur around fast
upstrokes. This table does not replace an equation-by-equation Rush–Larsen
audit, a parameter override audit, or time-step convergence. In particular,
AlievPanfilov's recovery coefficient depends on its recovery state, so its
current frozen-coefficient exponential update is an approximation rather
than an exact linear Rush–Larsen update.

The 0.5 s scalar/batched CPU wall times below were measured on xenosim
(AMD EPYC 9684X), OpenFOAM v2412, one process, with per-step trace output.
They include case setup and output and are not a clean kernel benchmark.

| Model | Scalar (s) | Batched CPU (s) |
|---|---:|---:|
| AlievPanfilov | 2.37 | 2.43 |
| BuenoOrovio | 2.46 | 2.58 |
| Courtemanche | 3.19 | 3.45 |
| Fabbri | 3.44 | 3.76 |
| Gaur | 3.55 | 3.79 |
| Grandi | 3.96 | 4.20 |
| PerisYague | 3.14 | 3.45 |
| Stewart | 3.15 | 3.29 |
| TNNP | 3.33 | 3.32 |
| TWorld | 5.38 | 6.17 |
| ToRORd_dynCl | 4.14 | 4.49 |
| Trovato | 3.95 | 4.48 |

## Acceptance gates for the longer validation

The short screen above is diagnostic and has no validated-model verdict.
Before a model receives one, its scalar reference must be stable under at
least two reductions of the integration step. Scalar, batched CPU, and GPU
cases must have identical state initialization, parameters, stimulus, tissue
type, and physical output times. State and parameter indices must match
exactly. Every trace must remain finite, with positive concentrations and
bounded gates where the model equations require them.

For a complete paced action potential, the initial comparison targets are
activation-time difference at most 0.1 ms, peak-voltage difference at most
2 mV, APD90 difference at most 2 ms or 2% (whichever is larger), and voltage
RMSE outside a 2 ms window around the upstroke at most 0.5 mV. These limits
separate a small timing shift from a sustained waveform error and exceed
the trace precision and the short-screen discrepancies above. Relevant
states and currents must also show decreasing error as the ionic step is
refined; current sign and integrated charge must agree with the scalar
reference. A model fails validation if convergence is absent even when a
single run meets the voltage limits. Tissue validation additionally
requires a stable propagating front, activation maps and conduction speed
that converge with time and mesh refinement, and an explicit comparison
with the scalar tissue run. Model-specific tolerances will be recorded
before its full run if scalar reference convergence warrants a change.

## 2D slab and time-control screen

`slab2D/` supplies a 20 mm by 3 mm, 1,500-cell TNNP monodomain slab with
one empty mesh layer. It uses TNNP epicardial cells and a 2 ms external
stimulus. This is a simple coupling test with Niederer-style conductivity,
not the published Niederer benchmark. See its README for commands.

At 15 ms, the scalar and batched CPU baseline both activated 622 cells.
Their `Vm` field RMSE was 0.0002 mV, maximum difference 0.0020 mV, and
95th percentile absolute activation-time difference 0.00001 ms. The
batched run took 59.56 s and the scalar run 64.96 s on xenosim. These are
single process CPU wall times from a short case, not GPU speedups.

| Tissue `deltaT` | Ionic substeps | Effective ionic step | 15 ms `Vm` RMSE versus baseline | Activated cells | Runtime |
|---:|---:|---:|---:|---:|---:|
| 2 us | 5 | 0.4 us | baseline | 622 | 59.56 s |
| 2 us | 1 | 2 us | 0.00077 mV | 622 | 27.36 s |
| 4 us | 10 | 0.4 us | 0.23438 mV | 622 | 50.24 s |

The doubled tissue step changed the solution despite the same effective
ionic step. These three runs are too short to establish long-term stability
or an accurate larger-step setting.

Additional scalar versus batched CPU runs used the same simple slab and
stimulus for the two models highlighted in the historical GPU port. TWorld
used `deltaT=2e-6 s`, five ionic steps; BuenoOrovio used `deltaT=2e-5 s`,
five ionic steps. At 15 ms, TWorld activated 730 cells in each run (voltage
RMSE 0.00087 mV, maximum 0.0106 mV, p95 activation shift 0.00010 ms).
BuenoOrovio initially failed to propagate in the batched case. Inspection
found its batched CPU coupling returned raw `Jion`; the scalar and CUDA
paths use `85.7*Jion`. After fixing the CPU hot-path conversion and
rebuilding, scalar and batched both activated 712 cells (voltage RMSE
0.02099 mV, maximum 0.4556 mV, p95 activation shift 0.00159 ms). These
runs used the CPU backend, and remain short coupling checks rather than
full tissue validation.

All 12 models then completed a paired 5 ms simple-slab screen with five
ionic substeps. The configured tissue step was `2e-6 s` for 11 models
(nominal ionic substep `4e-7 s`) and `2e-5 s` for BuenoOrovio (nominal
ionic substep `4e-6 s`). Each pair used its single-cell tissue type and the
shared external slab stimulus. These are nominal settings: the solver may
shorten the final step to land on `endTime`; each batched call divides the
actual step it receives into five ionic substeps. This implicit case does
not apply the explicit-diffusion `maxCo` cap described above.
Every pair completed on the CPU backend. The table reports activated cell
counts and final-time voltage comparison; it is a smoke screen and the
activation counts are not propagation validation. Fabbri had no threshold
crossing by 5 ms in either case.

| Model | Scalar/batched activated cells | Final Vm RMSE (mV) | Maximum difference (mV) |
|---|---:|---:|---:|
| AlievPanfilov | 141 / 141 | 0.000000 | 0.000010 |
| BuenoOrovio | 177 / 177 | 0.017587 | 0.484080 |
| Courtemanche | 93 / 93 | 0.000030 | 0.000160 |
| Fabbri | 0 / 0 | 0.000005 | 0.000100 |
| Gaur | 156 / 156 | 0.000004 | 0.000100 |
| Grandi | 60 / 60 | 0.000009 | 0.000100 |
| PerisYague | 161 / 161 | 0.000082 | 0.000500 |
| Stewart | 181 / 181 | 0.000822 | 0.005400 |
| TNNP | 159 / 159 | 0.000148 | 0.000700 |
| TWorld | 193 / 193 | 0.000229 | 0.003030 |
| ToRORd_dynCl | 157 / 157 | 0.000023 | 0.000300 |
| Trovato | 148 / 148 | 0.000011 | 0.000100 |

Prepare matched scalar/batched cases for any of the 12 models with
`prepare_slab_case.py`, or run the 5 ms matrix with `run_slab_matrix.py`.
Both generators refuse to overwrite an existing output path. Compare saved
fields with `compare_slab.py`.

## Actual CUDA single-cell validation

The following paired traces were run on the RTX 4000 Ada with the GPU
backend explicitly required. The baseline uses scalar RKF45 versus batched
Rush–Larsen for exponential gate states and Euler for remaining states, with
`deltaT=2e-6 s` and `batchedSubsteps 5` (nominal ionic substep `4e-7 s`). Each
case wrote all exported states, algebraics, and rates every 0.1 ms; all traces
were finite and scalar/batched field names matched. AlievPanfilov was extended
to 1 s because its action potential had not repolarized by 0.5 s. Errors are
computed as described in the acceptance section above; current error is the
percent difference of integrated total current relative to integrated
absolute scalar current.

| Model | GPU status | Activation shift (ms) | Peak shift (ms) | APD90 difference (ms) | Voltage RMSE outside upstroke (mV) | Integrated current error (%) |
|---|---|---:|---:|---:|---:|---:|
| AlievPanfilov | complete, 1 s | 0.00160 | 0.000 | -0.0028 | 0.00087 | 0.00062 |
| BuenoOrovio | complete | 0.00051 | 0.000 | 0.0014 | 0.00077 | 0.0820 |
| Courtemanche | complete | 0.00114 | 0.000 | -0.0367 | 0.00257 | 0.00066 |
| Fabbri | complete, spontaneous | 0.00162 | 0.000 | 0.0013 | 0.00195 | 0.00024 |
| Gaur | complete | 0.00352 | -0.100 | 0.0360 | 0.0267 | 0.0438 |
| Grandi | complete | 0.00076 | 0.000 | -0.0089 | 0.00076 | 0.00035 |
| PerisYague | complete | 0.00083 | 0.000 | 0.0007 | 0.00176 | 0.00038 |
| Stewart | complete | 0.00062 | 0.000 | -0.0144 | 0.00303 | 0.00531 |
| TNNP | complete | 0.00149 | 0.000 | -0.0034 | 0.00083 | 0.0350 |
| TWorld | complete | 0.00098 | 0.000 | -0.0030 | 0.00166 | 0.1184 |
| ToRORd_dynCl | complete | 0.01060 | 0.000 | 0.0855 | 0.0425 | 0.0901 |
| Trovato | complete | 0.00076 | 0.000 | -0.0094 | 0.00166 | 0.0973 |

Every activation and APD difference is inside the declared single-cell
limits. Outside-upstroke voltage errors are below 0.05 mV for every model,
and integrated current errors are below 0.12%. The peak-time shift is shown
separately because it can reflect an upstroke timing difference. Each run's
`summary.csv` also records peak voltage error and the all-state traces.

The full-length 10 µs / one-substep stress run completed with finite traces
for ten models, but Grandi and TWorld became invalid before stimulation.
Matched CPU-batched runs at the same step also failed with a floating-point
exception in `pow`, after calcium/buffer states left their physical domain;
scalar RKF45 remained finite. Both CUDA models completed at the validated
2 µs / five-substep baseline. This identifies an integration stability limit
for Euler-updated concentration/buffer states, rather than an observed
CUDA-only state-index or model-equation discrepancy. Do not use
`deltaT=1e-5 s`, one substep as a general setting for these models.

Fabbri is the AV-node model and fires spontaneously in its single-cell case;
the 5 ms slab screen had no threshold crossing for either backend. Its slab
result is a finite-state/parity check, not a failure to propagate under the
generic stimulus.

## Actual CUDA 15 ms slab propagation

The 11 non-Fabbri models completed matched scalar and actual-GPU-batched runs
on the 20 mm by 3 mm, 1,500-cell simple slab, using each model's tissue type,
the same external pulse, and its listed stable tissue step and five ionic
substeps. The Fabbri slab remains the short stability/parity case described
above. Activated counts match for every included model. These cases check
propagation on the slab geometry; they do not reproduce the published
Niederer benchmark.

| Model | Activated cells, scalar/GPU | Final Vm RMSE (mV) | Maximum difference (mV) | p95 activation shift (ms) | Scalar/GPU runtime (s) | End-to-end speedup |
|---|---:|---:|---:|---:|---:|
| AlievPanfilov | 545 / 545 | 0.000004 | 0.0001 | 0.000000 | 13.92 / 11.46 | 1.21x |
| BuenoOrovio | 712 / 712 | 0.022613 | 0.4273 | 0.001800 | 2.60 / 1.89 | 1.38x |
| Courtemanche | 183 / 183 | 0.000346 | 0.0012 | 0.000000 | 62.41 / 19.39 | 3.22x |
| Gaur | 595 / 595 | 0.000675 | 0.0019 | 0.000000 | 83.67 / 23.15 | 3.61x |
| Grandi | 399 / 399 | 0.000049 | 0.0005 | 0.000100 | 102.02 / 27.20 | 3.75x |
| PerisYague | 666 / 666 | 0.000200 | 0.0005 | 0.000000 | 59.21 / 16.85 | 3.51x |
| Stewart | 805 / 805 | 0.001982 | 0.0063 | 0.000008 | 62.66 / 23.50 | 2.67x |
| TNNP | 622 / 622 | 0.000193 | 0.0020 | 0.000010 | 72.71 / 22.34 | 3.25x |
| TWorld | 730 / 730 | 0.000869 | 0.0106 | 0.000100 | 189.38 / 41.42 | 4.57x |
| ToRORd_dynCl | 569 / 569 | 0.012042 | 0.2023 | 0.000700 | 107.12 / 29.90 | 3.58x |
| Trovato | 557 / 557 | 0.000021 | 0.0002 | 0.000000 | 108.01 / 30.45 | 3.55x |

Runtime is wall time for the full serial case command using one MPI rank and
one CPU core, including setup and output, on xenosim (AMD EPYC 9684X plus RTX
4000 Ada); it is not a kernel-only speedup measurement or a comparison with a
fully parallel CPU run. Reproduce with
`run_slab_matrix.py --models ... --end-time 0.015 --workers 1 --require-gpu`.

## Niederer GPU verification

The current `tutorials/NiedererEtAl2011/NiedererEtAl2011verification` case
was copied to `/tmp` and run with `TNNPcompactBatched`, Rush–Larsen, five
ionic substeps, the case's original 10 µs tissue step, geometry, stimulus,
and 15 ms end time. The mesh has 100 x 15 x 35 = 52,500 cells. The CUDA log
confirmed device 0. The existing Niederer activation regression passed all
6 probe checks with differences exactly 0 at the printed precision and
tolerance `1e-4`. Solver time was 75.25 s on xenosim (RTX 4000 Ada, one host
rank); this includes tissue solves and output and is not a kernel-only time.

Reproduce from the repository root on a GPU allocation:

```bash
cp -a tutorials/NiedererEtAl2011/NiedererEtAl2011verification \
    /tmp/niederer-tnnp-gpu
cd /tmp/niederer-tnnp-gpu
foamDictionary constant/electroProperties \
    -entry monodomainSolverCoeffs.ionicModel -set TNNPcompactBatched
foamDictionary constant/electroProperties \
    -entry monodomainSolverCoeffs.batchedIntegrator -set rushLarsen
foamDictionary constant/electroProperties \
    -entry monodomainSolverCoeffs.batchedSubsteps -set 5
```

Then run the case and check the supplied reference:

```bash
srun -N1 -n1 -c1 -p dev --time=00:30:00 bash -lc '
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
export CARDIAC_ENABLE_CUDA=1 CUDA_HOME=/usr NVARCH=75 FOAM_SIGFPE=false
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
cd /tmp/niederer-tnnp-gpu
./Allrun
bash regression/regressionTest.sh --check-only
'
```

This validates the current TNNP Niederer case against its stored activation
reference. The backup also has historical TWorld and BuenoOrovio variants;
their old reference matrices are not reused as current-main validation.

## Cross-model Niederer-geometry matrix

`run_niederer_model_matrix.py` prepares paired scalar RKF45 and CUDA batched
cases for all 11 in-scope models (Fabbri is omitted because it is an AV-node
automatic pacemaker, not a normal paced ventricular cell). It preserves the
current case's mesh, tissue solver, conductivity, stimulus, 10 us tissue
`deltaT`, and output probes, while setting each model's supported tissue
type. GPU cases use Rush--Larsen gate updates plus Euler for other states and
25 substeps (0.4 us nominal ionic step). This common setting gives a
conservative cross-model comparison. The completed paired matrix is evidence
that all 11 models selected CUDA and completed on this geometry and protocol.
It is not published Niederer benchmark validation; only TNNP has the supplied
benchmark reference.

The scalar Gaur pilot reached only 1.5 ms simulated time in about 3.3 minutes
on xenosim, so a full 15 ms cross-model matrix needs a long allocation. Run it
serially to avoid sharing one GPU across simultaneous jobs:

```bash
srun -N1 -n1 -c1 -p main --time=08:00:00 bash -lc '
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
export CARDIAC_ENABLE_CUDA=1 CUDA_HOME=/usr NVARCH=75 FOAM_SIGFPE=false
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
cd /home/simao/cardiacFoam
python3 tutorials/gpuValidation/run_niederer_model_matrix.py \\
    /tmp/cardiac_niederer_all_models_20260928 --substeps 25
'
```

The runner records per-backend wall times, verifies the CUDA device was
selected, writes `summary.csv`, and compares Vm and activation at 5, 10, and
15 ms. Those times include mesh setup, tissue solves, and output, not only
ionic GPU work. The completed result on xenosim used one MPI rank, one CPU
core, and one RTX 4000 Ada GPU:

This matrix is deliberately a common-setting comparison: every model used the
same 52,500-cell mesh, 15 ms duration, and tissue `deltaT=10 us`; scalar cases
used RKF45, while GPU cases used Rush--Larsen plus Euler with 25 ionic
substeps. It controls the TNNP/TWorld comparison, but neither model was run
at its individually fastest accurate setting. The slab sweep is selecting
candidate settings independently for each model.

| Model | Scalar s | GPU s | Speedup | Vm RMSE at 15 ms (mV) | Activation p95 shift (ms) |
|---|---:|---:|---:|---:|---:|
| AlievPanfilov | 205.34 | 205.80 | 1.00x | 0.000004 | 0.000000 |
| BuenoOrovio | 274.65 | 79.22 | 3.47x | 0.006426 | 0.000600 |
| Courtemanche | 483.18 | 105.79 | 4.57x | 0.000172 | 0.000010 |
| Gaur | 1914.98 | 230.89 | 8.29x | 0.000333 | 0.000000 |
| Grandi | 1000.34 | 267.46 | 3.74x | 0.000030 | 0.000100 |
| PerisYague | 460.95 | 99.66 | 4.63x | 0.000149 | 0.000000 |
| Stewart | 934.53 | 120.33 | 7.77x | 0.001410 | 0.000010 |
| TNNP | 1199.76 | 116.47 | 10.30x | 0.000174 | 0.000010 |
| TWorld | 2042.49 | 310.96 | 6.57x | 0.000520 | 0.000100 |
| ToRORd_dynCl | 2048.76 | 269.44 | 7.60x | 0.011899 | 0.000900 |
| Trovato | 796.89 | 271.19 | 2.94x | 0.000061 | 0.000010 |

All models had matching activated-cell counts at 15 ms except BuenoOrovio,
which differed by one cell (11,170 scalar versus 11,171 GPU). These are short
propagation comparisons; they do not establish full action-potential accuracy
or the largest stable timestep for each model. Raw logs and traces are in
`/tmp/cardiac_niederer_all_models_20260928/summary.csv` and the model subfolders.

## Manufactured CUDA correctness smoke

The scalar and batched model were run on the same 8 x 8 x 8 mesh for 50 us,
with outer `deltaT=10 us`; the CUDA Euler path used five ionic substeps. The
run log selected CUDA device 0. Both cases report identical analytic `Vm`
errors (`L1=3.87484e-7`, `L2=5.34119e-7`, `Linf=2.05871e-6`). The `u1/u2`
errors also remain finite and near `1e-10`; their scalar/GPU differences are
expected because the reference uses RKF45 while the batched model uses fixed
Euler. This is an initial correctness smoke, not temporal-convergence or
performance validation. Reproduce on a GPU node with:

```bash
sbatch tutorials/gpuValidation/run_manufactured_cuda_smoke.slurm
```

The paired case dictionaries, logs, and summary are under
`/tmp/cardiac_manufactured_cuda_smoke_20260928/`.

## Remaining validation

- Audit every model's state and parameter mapping equation by equation,
  including non-default parameter overrides and non-default sex/tissue
  selections. Trace field names, finiteness, default initialization, and
  reference comparisons passed, but every override is not covered.
- Run mesh and longer-time refinement for slab conduction speed, activation
  maps, APD, and recovery. Current 15 ms slabs establish provisional
  activation-time candidates, but do not establish full action-potential
  accuracy or mesh independence.
- The 11-model Rush--Larsen slab sweep has completed at outer `deltaT` values
  10, 20, and 50 us with 1, 5, 10, and 25 GPU substeps. Results are in
  `/tmp/cardiac_gpu_integrator_sweep_20260928_nonfabbri/summary.csv` and
  `best_settings.md`. On the current 2D test, TNNP's provisional choice is
  10 us / 1 substep (3.04 s, 8.49x over matched scalar RKF45); TWorld's is
  10 us / 5 substeps (8.14 s, 5.75x). These pass the current activation
  gates against the fine scalar reference, but still need longer APD/state
  validation.
- The manufactured fixed-mesh ladder has completed under
  `/tmp/cardiac_manufactured_convergence_20260928/`. It used an 8 x 8 x 8
  mesh, 50 us final time, outer `deltaT` values 2.5, 5, and 10 us, and GPU
  Euler substep counts 1, 5, and 25. It compared `Vm`, `u1`, and `u2` error
  norms against the analytic solution. `Vm` L2 error stayed near
  `5.34e-7` across the ladder, indicating a spatial/error floor rather than
  a measurable temporal slope; `u1`/`u2` errors were small, with `u2`
  showing clearer substep sensitivity. CUDA and scalar `Vm` norms matched at printed
  precision. This is a short-horizon parity check, not evidence of second
  order convergence, and it did not exercise SBDF2 coupling.
- The 20 and 40 us TNNP/TWorld slab extension has completed under
  `/tmp/cardiac_gpu_tnnp_tworld_coarse_dt_20260928/`. Matched scalar/GPU
  traces remain close, but neither model passed the current 0.1 ms p95
  activation-shift gate against the fine reference at these coarser steps.
  Thus 10 us is the present candidate for both pending longer validation.
- The per-model 3D Niederer-geometry rerun (job 9356) completed on one
  RTX 4000 Ada and one MPI rank. It used 15 ms, 10 us tissue `deltaT`, and
  1/5 ionic substeps for TNNP/TWorld respectively. TNNP took 1077.94 s scalar
  and 62.42 s GPU (17.27x); TWorld took 1830.54 s scalar and 192.21 s GPU
  (9.52x). At 15 ms their activation p95 shifts were 0.00030 ms and 0.00010
  ms. These are tuned matched-step timings on the Niederer geometry, not a
  reproduction of the published benchmark and not long APD/state validation.
  The prior common-setting 3D matrix used 25 substeps for both and measured
  10.30x for TNNP and 6.57x for TWorld.
- Job 9358 completed all CPU and CUDA scalar/batched TNNP/TWorld
  Godunov/SBDF2 cases at 2 us, five substeps, and 15 ms. Godunov scalar/CUDA
  comparisons are close, but direct CPU-batched/CUDA-batched SBDF2 parity
  failed for both models as quantified above. Job 9359 completed the coarse
  640x640 scalar SBDF2 temporal ladder with acceptance PASS and observed
  `Vm L2` orders 1.96983, 2.00954, and 2.10796.
- Job 9361 completed the manufactured FDA batched model on CPU and CUDA using
  the backup study's 2D N=640 mesh, `t=0.2`, `Vm/u1/u2` L1/L2/Linf metrics,
  Godunov/Euler, 25 substeps, and the coarse timestep ladder. Both backends
  passed the first-order convergence gate. CPU-batched and CUDA-batched
  `Vm`, current, activation, analytic-error fields, and restart states
  matched within `1.9984e-14` maximum absolute difference. Their `Vm L2`
  errors are 0.39–0.43% above the backup scalar Godunov/Euler reference at
  all five steps, with the same first-order trend. This validates manufactured
  convergence/verification and batched CPU/GPU parity. Job timings are
  excluded because a separate 48-rank scalar Gaur job shared the node.
- Job 9362 completed the manufactured SBDF2/backward ladder. CPU-batched and
  CUDA-batched runs used 25, 50, 100, 200, and 400 Euler substeps as `deltaT`
  halved. Both passed second-order L2 convergence gates for `Vm`, `u1`, and
  `u2`; their CPU/GPU fields and saved states matched within `2e-14`. Coarse
  `Vm L2` errors were about 1.8–1.9% above the scalar SBDF2 reference. The
  job shared the node with scalar Gaur electromechanics, so its times are not
  used as performance data. Full ladder details are in
  `SBDF2_MANUFACTURED_TEMPORAL.md`.
- The published Niederer TNNP case is the only benchmark-reproduction
  target. Its current-main case has already run with the original 10 us
  tissue step, five ionic substeps, and stored activation reference; all six
  probe checks passed. The other models' Niederer-shaped runs are geometry
  comparisons, not claims of reproducing that published benchmark.
- Profile transfers and kernels separately from case setup/output, and test
  multiple cell counts and GPU devices before publishing general performance
  claims.

Compilation and CPU fallback runs are not recorded as GPU validation. Each
CUDA matrix command above used `--require-gpu` and checked the runtime device
message.

The backup's June 2026 analysis reported 16.9x for TNNP and 24.7x for TWorld,
but those results used a different workload: 420,000 cells for 55 ms. The GPU
cases used `deltaT=20 us` with 10 ionic substeps; their scalar comparisons used
RKF45 at `deltaT=2 us`. Those measurements support testing coarser settings,
but they are not directly comparable to the current common-setting matrix.
The Niederer case generator now accepts explicit settings. For example, to
compare each backend at a matched 20 us tissue step on a GPU allocation,
using the CUDA environment documented above:

```bash
python3 tutorials/gpuValidation/run_niederer_model_matrix.py \
    /tmp/niederer-tnnp-20us --models TNNP --delta-t 2e-5 \
    --substeps 10 --end-time 0.015
python3 tutorials/gpuValidation/run_niederer_model_matrix.py \
    /tmp/niederer-tworld-20us --models TWorld --delta-t 2e-5 \
    --substeps 10 --end-time 0.015
```

The tuned 3D rerun compares scalar RKF45 and CUDA at the same macro time step;
it does not itself compare against the fine-step scalar reference. The 2D
fine-reference comparison is the basis for the provisional candidates. After
the 3D run, compare activation and propagation against the fine-step 3D
reference before treating its speedup as an accuracy-matched result.

Run the manufactured convergence ladder on the GPU node after the larger
single-GPU sweep, or submit its serialized job:

```bash
sbatch tutorials/gpuValidation/run_manufactured_convergence.slurm
```

By default it writes `/tmp/cardiac_manufactured_convergence_20260928/` and
uses the same 8 x 8 x 8 mesh and 50 us final time as the smoke. The scalar
reference runs once per outer `deltaT`; CUDA Euler runs use 1, 5, and 25
substeps. The `summary.csv` contains exact-solution error norms and runtime
for `Vm`, `u1`, and `u2`. Compare rows at fixed `deltaT` to assess ionic
substep convergence, then fixed substeps across `deltaT` to assess outer
temporal behavior. This small case is a correctness/convergence check, not a
GPU performance benchmark.
