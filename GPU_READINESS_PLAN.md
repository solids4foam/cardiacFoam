# GPU ionic-model readiness plan

## Scope and target

The target is the CUDA execution path for the batched ionic and active-tension
models in `cardiacFoam`.  The tissue PDE remains CPU-resident; a successful
result must not be described as a full-PDE GPU port.

This worktree starts from `a6813e26` on branch
`gpu-characterization-20261001`.  The source and untracked outputs in the two
existing repository directories are not modified.

### Active-tension validation route

1. `tutorials/electrophysiologyProtocols/singleCell` is the primary controlled
   scalar/batched/CUDA trace-parity harness.  It already couples
   `LandNiedererTWorld` active tension to the TWorld ionic signal and exports
   `Ta`.
2. `tutorials/electromechanicsProtocols/springSupportedSlab` is the coupled
   validation case after trace parity.  The local solids4foam branch contains
   both the `electroMechanicalLaw` and `solidRobin` support-condition sources
   required by that configuration.
3. CUDA active tension is not validated by compilation: device dispatch,
   state update, and scalar/batched/CUDA trace parity are each required.

## Gates

1. Establish a reproducible CUDA build and verify that the solver actually
   selects a device.  A CUDA-linked library or successful compilation alone is
   insufficient.
2. Inventory all scalar, host-batched, and CUDA ionic model implementations
   and audit state, parameter, tissue, unit, and current mappings.
3. Validate all models in matched single-cell simulations before tissue runs.
3a. Implement and validate matched host-batched/CUDA execution for every
    batched active-tension model before claiming it has GPU support.
4. Document the exact time-control and coupling call paths, then validate them
   in a small 2-D slab.
5. Validate the applicable myocardial models in matched 2-D slabs; use the
   Niederer geometry separately from reproduction of a published benchmark.
6. Profile, optimize only demonstrated bottlenecks, and rerun the affected
   validation matrix.
7. Deliver reproducible commands, results, hardware metadata, and limitations.

## Initial evidence and decisions

- CUDA device nodes were absent at the start of this work.  The documented
  `nvidia-modprobe -u -c=0` recovery path restored visibility of an NVIDIA RTX
  4000 Ada GPU on 2026-10-01.
- Legacy `tutorials/benchmarkGPU` is a source of historical cases and results,
  not an implementation baseline.  Its tree is uncommitted and materially
  diverges from the target branch.
- The current tree provides 12 scalar models, 12 batched wrappers, and 12 CUDA
  launchers.  Fabbri is single-cell-only for this campaign because it is an
  AV-node model.
- Historical local artifacts establish broad scalar/CUDA slab coverage, but
  the tracked record distinguishes strict host/CUDA equality from
  physiologically small arithmetic-order differences. See
  `IONIC_VALIDATION_EVIDENCE.md`; final release evidence must be reproduced
  from this worktree rather than retained only in `/tmp`.

## Acceptance policy (to refine per model before execution)

- Matched host-batched/CUDA runs: same model settings, no invalid values, and
  agreement of high-precision Vm, states, and Iion traces under an
  absolute-plus-relative tolerance recorded with the run.
- Scalar reference/batched runs: compare time-step convergence, activation
  time, resting and peak Vm, APD50/APD90, relevant currents/states, and
  post-upstroke trace error.  A pointwise upstroke difference is interpreted
  with activation-time shift rather than independently deemed a failure.
- Slabs: stable propagation, matching activation masks/maps, CV and
  activation-time agreement, and no backend-specific instability.
- Performance: separately report ionic compute, host/device transfer and
  synchronization, tissue/PDE time, total wall time, GPU model/driver/CUDA,
  rank placement, mesh, duration, and warm-up policy.

## Change log

| Date | Change | Reason |
| --- | --- | --- |
| 2026-10-01 | Replaced the initial “GPU allocation required” assumption. | Legacy runtime notes and a direct test showed that this host uses an unreserved GPU whose missing device nodes can be restored with `nvidia-modprobe`. |
| 2026-10-01 | CUDA smoke path confirmed through `srun -p dev`. | The current TNNP slab logged `using CUDA device 0 of 1` and reached `End`; direct-shell probes remain unreliable on this host. |
| 2026-10-01 | Built the target commit's CUDA ionic library in `/tmp/cardiac-gpu-characterization-20261001/lib`. | This removes ambiguity from the previously installed ionic library timestamp.  The build completed with non-fatal external-header warnings and batched-wrapper member-order warnings. |
| 2026-10-01 | Repaired the active-tension CUDA build wiring and device annotations. | CUDA sources were listed but not compiled; the three kernel translation units also included their batch equations without `__host__ __device__` annotation. |
| 2026-10-01 | Verified the rebuilt target library stack through the physical GPU. | A short TNNP 2-D slab using the isolated target libraries logged `using CUDA device 0 of 1` and reached `End`. |
| 2026-10-01 | Recorded a non-blocking external build limitation. | Rebuilding `cardiacFoam` itself stops at the site PETSc link because `libpetsc.so` leaves `sgemmt_`/`dgemmt_` unresolved.  The target electro libraries nevertheless build and run correctly through the existing driver binary. |
| 2026-10-01 | Expanded scope to active-tension GPU parity. | The requested audit found three scalar/batched/CUDA model families: NashPanfilov, LandNiederer, and LandNiedererTWorld. |
| 2026-10-01 | Reclassified active-tension CUDA as unimplemented, not merely unvalidated. | Its CUDA launch wrappers are compiled but have no call sites. `batchedActiveTensionModel::calculateTension()` always dispatches the host/OpenMP executor, so current batched active-tension runs cannot select or test the GPU. |
| 2026-10-01 | Selected the existing single-cell and spring-supported cases for active-tension validation. | The single-cell case provides an existing `LandNiedererTWorld`/`Ta` trace path; the spring slab supplies the follow-on coupled mechanics test.  Source audit confirmed the local solids4foam branch includes `electroMechanicalLaw` and `solidRobin`. |
| 2026-10-01 | Implemented and tested first active-tension CUDA dispatch. | A 100 ms single-cell `LandNiedererTWorld` run on the RTX 4000 Ada exercised the device path.  An initial 1000x rate-unit defect was isolated against the host-batched path and corrected.  The corrected CUDA and host-batched `Ta` traces are identical at all 100 written samples; scalar/CUDA maximum sampled `Ta` difference is 5.21e-05 kPa. |
| 2026-10-01 | Reproduced and corrected the initial spring-slab batched failure. | Scalar active tension completed 5 ms while host-batched and CUDA failed at the second 10 us step. The batched 100-step (10 ms Euler) preconditioner produced invalid crossbridge states at TNNP resting Cai. It now uses a configurable maximum 0.1 ms Euler step. Fresh 20-step host-batched and CUDA slabs both reached `End`, logged finite mechanics convergence, and wrote identical five-probe `Ta` values at 0.2 ms. |
| 2026-10-01 | Completed controlled 100 ms single-cell parity for NashPanfilov and LandNiederer active tension. | Each used the existing TWorld ionic drive, identical 1 us step, and scalar/host-batched/CUDA variants. In both families host-batched and CUDA `Ta` traces were identical at all 101 samples. Scalar-versus-batched maximum sampled `Ta` differences were 0.0067594 kPa (NashPanfilov) and 8.04e-05 kPa (LandNiederer), attributable to the intentionally different scalar ODE solver and batched Euler integrator; neither is a CUDA discrepancy. |
| 2026-10-01 | Completed the 20-step spring-slab CUDA gate for all three active-tension families. | NashPanfilov, LandNiederer, and LandNiedererTWorld host/CUDA pairs reached `End`, converged mechanics at each step, and wrote identical `Ta` probe files. This is a coupled stability/parity check, not a full mechanics benchmark. |
| 2026-10-01 | Made active-tension ODE data device-resident between coupled updates. | The steady-state path uploads only model inputs and downloads only `Ta`; complete state/rate/algebraic transfer is deferred until restart/export I/O requests it. All three single-cell CUDA traces and a post-change LandNiedererTWorld slab remain exactly equal to their host-batched references. |
| 2026-10-01 | Completed a sustained CUDA active-tension coupled run. | The 3,360-cell LandNiedererTWorld spring slab completed 2,000 coupled steps (20 ms) in 139.64 s with finite mechanics convergence and restart state files at 10 and 20 ms. |
| 2026-10-01 | Completed matched 20 ms host/CUDA active-tension slab parity. | All 2,000 `Ta` probe samples are byte-identical. Restart-state serialization differs only by bounded instruction-order roundoff (maximum 9.29e-14 at 10 ms and 8.41e-14 at 20 ms), so this is not presented as binary restart parity. |
| 2026-10-01 | Completed sustained NashPanfilov host/CUDA coupled parity. | The 3,360-cell 20 ms spring slab reached `End` on both backends; all 2,000 `Ta` probe samples are byte-identical and restart-state disagreement is at most 7.11e-15. |
| 2026-10-01 | Isolated a long LandNiederer spring-slab failure from GPU and batching. | Scalar and host-batched models both hit the nonlinear-solid 1,000-corrector ceiling near 18 ms. Relaxing `rTol` from 0.02 to 0.15 and halving `deltaT` to 5 us did not fix it; the path remains an unresolved model/mechanics-coupling limitation. |
| 2026-10-01 | Added opt-in active-tension CUDA event characterization. | The 3,360-cell, 2,000-step TWorld run measured 0.103586 s H2D, 0.147324 s kernels, 0.0368411 s D2H and 0.542988 s wrapper time. The unprofiled matched host/CUDA whole-case wall times were 149.43/139.64 s (1.070x). |
| 2026-10-01 | Added and tested an opt-in Land-2017 stretch-rate cap. | `maxLambdaRate` defaults to unlimited. At 20 s^-1, scalar, host-batched, and confirmed CUDA runs all complete the 20 ms spring slab that fails uncapped. The CUDA/host-batched coupled tension difference is at most 9.62 Pa (0.90% of sampled peak), so it is a stabilization candidate, not yet a default or physiological acceptance. |
| 2026-10-01 | Completed the full scalar capped-Land spring twitch. | Six ranks completed 250 ms in 430 s with convergence at every step. Five active-tension probes peak at 14.70--15.68 kPa between 83.0 and 123.2 ms; this rules out a gross tension-scale error but does not replace isometric Land-reference validation. |
