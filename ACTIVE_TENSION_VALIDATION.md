# Active-tension CUDA validation

## Scope

These results cover the three active-tension families that have scalar,
batched-host, and CUDA implementations: `NashPanfilov`, `LandNiederer`, and
`LandNiedererTWorld`.  The tissue PDE and mechanics solver remain host
resident; this is GPU acceleration of the batched active-tension ODEs only.

## Test environment

- Source worktree: `gpu-characterization-20261001` from `a6813e26`.
- OpenFOAM: v2412; CUDA compiler: nvcc with GCC 10 host compiler.
- GPU: NVIDIA RTX 4000 Ada, one rank on the local `dev` Slurm partition.
- Active-tension library:
  `/tmp/cardiac-gpu-characterization-20261001/lib/libactiveTensionModels.so`.

## Controlled single-cell parity

`tutorials/electrophysiologyProtocols/singleCell` was run for 100 ms with a
1 us solver step and the same TWorld ionic drive, stimulus, initial state, and
parameters in each variant.  Scalar models use their existing ODE solver;
batched host and CUDA models use the documented explicit Euler update with one
substep.  The exported `Ta` trace contains 101 samples.

| Model | Scalar vs host-batched max \|Ta\| [kPa] | Host-batched vs CUDA | CUDA device dispatch | Result |
| --- | ---: | ---: | --- | --- |
| NashPanfilov | 0.0067594 | 0 at 101/101 samples | yes | pass |
| LandNiederer | 8.04e-05 | 0 at 101/101 samples | yes | pass |
| LandNiedererTWorld | 5.21e-05 | 0 at 100 written samples | yes | pass |

The scalar/batched values are integration-method differences, not CUDA
differences: every batched host/CUDA pair is exactly equal in the written
trace.  The scalar comparison is acceptable for these traces; it must still
be repeated with timestep-convergence and event metrics before being used as a
physiological accuracy claim.

## Coupled spring-supported slab gate

The `springSupportedSlab` case has 3,360 integration cells and a 10 us coupled
step.  The following are short 0.2 ms (20-step) checks, not the documented
0.25 s mechanics benchmark.

| Model | Host-batched | CUDA | Mechanics | Host/CUDA Ta probes | Result |
| --- | --- | --- | --- | --- | --- |
| NashPanfilovBatched | End | End; device logged | converged every step | identical | pass |
| LandNiedererBatched | End | End; device logged | converged every step | identical | pass |
| LandNiedererTWorldBatched | End | End; device logged | converged every step | identical | pass |

This also verifies the corrected Land-family preconditioning path.  The former
100 x 10 ms Euler preconditioner was unstable at TNNP resting calcium; batched
Land models now default to a maximum 0.1 ms startup increment through
`batchedPreconditioningMaxStep`.

An additional CUDA `LandNiedererTWorldBatched` stress run completed 20 ms
(2,000 coupled steps) in the same 3,360-cell spring slab. It reached `End` in
139.64 s, converged the momentum equation at every step, and wrote both 10 ms
and 20 ms `LandNiedererTWorldBatchedState` restart files. This explicitly
exercises the device-resident state synchronization used by restart/export;
it is a sustained-stability gate, not a CPU/GPU performance comparison.

A matched host-batched 20 ms run also reached `End` (149.43 s). Its complete
2,000-sample `Ta` probe file is byte-identical to CUDA. The serialized restart
states are not byte-identical because host and CUDA expression evaluation use
different floating-point instruction order: across all 23,520 saved values,
the maximum absolute difference is 9.29e-14 at 10 ms and 8.41e-14 at 20 ms
(RMS 4.02e-15 and 3.88e-15). This is accepted roundoff, not a restart-state
mapping or residency error.

The same 20 ms host/CUDA check was completed for `NashPanfilovBatched`.
Both runs reached `End`; the complete `Ta` probe file is byte-identical. The
host and CUDA wall times were 153.98 s and 144.38 s, respectively. Restart
state roundoff at 10 and 20 ms has maximum absolute error 7.11e-15 (RMS
6.22e-16 and 7.79e-16), so it also passes the sustained backend-parity gate.

`LandNiederer` does **not** yet pass this particular long spring-supported
slab. It is not a CUDA or batched regression: the scalar model stopped at
18.12 ms and host-batched at 18.07 ms when the nonlinear solid solver reached
its 1,000-corrector limit. A deliberately relaxed solid relative residual
criterion (`rTol` 0.15 instead of 0.02) still stopped at 18.09 ms, and a
scalar temporal-refinement run (`deltaT=5 us`, original solid tolerances)
stopped at 17.535 ms. These tests rule out a near-threshold convergence setting
and a simple 10-us time-step artifact. This is an unresolved Land-2017/
mechanical-law coupling issue; the 20-step result above remains only its short
parity gate.

At the failed interval the Land probe tension is about 1.0--1.1 kPa, not an
obvious kPa/Pa scale explosion. The model declares `Tref=120 kPa` and its
active output is converted once through the shared `TaScale=1000` kPa-to-Pa
interface. Full-twitch amplitude/timing checks against the Land reference are
still required before it can be accepted in deforming tissue.

## CUDA residency and performance interpretation

CUDA ODE state, rates, and algebraics now remain resident after the first
update.  Per solve, the implementation transfers only current model inputs and
the `Ta` row required by the host coupled solver.  Full state/rate/algebraic
copies are performed only when the existing restart/export I/O interface asks
for them.  A post-change 20-step LandNiedererTWorld CUDA slab reproduced the
host probe file exactly and completed in 2.83 s wall time (including Slurm
launch and complete solver work).

With `batchedCUDAProfile true`, the 3,360-cell / 2,000-step TWorld CUDA run
reported 0.103586 s H2D input, 0.147324 s kernels, 0.0368411 s D2H `Ta`, and
0.542988 s total active-tension wrapper time. This is 51.8, 73.7, 18.4, and
271.5 us per coupled update, respectively. The named device operations take
0.288 s (0.20% of the 141.00 s profiled whole-case time); the wrapper takes
0.39%. Timed CUDA events deliberately synchronize operations, so those
profiling-run wall-time figures must not be substituted for production wall
time. The matched unprofiled TWorld A/B wall time is 149.43 s host versus
139.64 s CUDA (1.070x speedup; 6.6% time reduction), while the PDE and solid
solver remain host-resident.

The recorded 100 ms single-cell CUDA wall times—4.86 s (Nash), 6.91 s (Land),
and 6.88 s (LandTWorld)—include solver startup, host ionic work, I/O, and Slurm
launch. They are not ionic- or active-tension-kernel speedup measurements.

## Reproduction

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
export LD_LIBRARY_PATH=/tmp/cardiac-gpu-characterization-20261001/lib:$LD_LIBRARY_PATH
blockMesh -case tutorials/electrophysiologyProtocols/singleCell
srun -N1 -n1 -p dev cardiacFoam -case tutorials/electrophysiologyProtocols/singleCell
```

Select `NashPanfilovBatched`, `LandNiedererBatched`, or
`LandNiedererTWorldBatched` in the case dictionary.  Leave `batchedUseCUDA`
at its default for CUDA, or set it to `false` for the host-batched reference.
For the coupled gate, follow `springSupportedSlab/Allrun` after selecting the
same model in `constant/electroMechanicalProperties`.

## Remaining work

- Resolve and validate Land-2017 length/rate feedback with the deforming solid;
  do not enable it as a long coupled production case until then.
- Run the complete 0.25 s coupled benchmark and timestep-convergence study
  for models whose spring cases are stable.
- Validate restart continuation across a write/restart boundary on CUDA.
- Complete the requested per-ionic-model single-cell/slab matrix and the
  Niederer benchmark characterization, including its own compute/transfer
  profiling rather than extrapolating this active-tension result.
