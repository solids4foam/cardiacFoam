# Manufactured ionic model on the GPU

## Why this is a useful test

`monodomainFDAManufactured` provides an analytic field reference, so a batched
counterpart can test state import/export, ODE updates, current coupling, and
CUDA execution against convergence norms rather than only against another
numerical implementation. The scalar model remains the reference; a registered
batched counterpart and its CUDA kernels are now implemented in this branch.
The CUDA library compiled, and a small CUDA runtime smoke completed with
matching scalar/GPU `Vm` error norms. That is a device-path smoke test, not a
temporal or spatial convergence validation.

The matching equations, initialization, names, units, exact-field import, and
PDE source must be preserved. With those held fixed, the batched and scalar
implementations should approach the same analytical solution and show the
same spatial convergence order once temporal error is sufficiently small.
That is an expectation to verify, not an automatic consequence of using the
same case. The batched model's Euler substeps can dominate the error and alter
the observed slope if `deltaT` is not refined appropriately.

## Scalar model contract to preserve

The model lives in `src/ionicModels/verificationModels/monodomainFDAManufactured`.
Its state order is `V, u1, u2, u3`; constants are `Cm, Beta, Chi`; and its
algebraic is `Iion`. `Cm=2`, `Chi=3`, and `Beta` is selected by dimension
(`-1.1`, `-5.9`, or `-8.6` for 1D, 2D, or 3D). Its nonzero rates are

```text
du1 = (u1 + u3 - V)^2*u2^2
      + 0.5*(u1 + u3 - V)*u2^2*(V - u3)
du2 = -(u1 + u3 - V)*u2^3
```

`dV=du3=0`, and `Im = Iion/Cm`. These polynomial rates are not gate equations
with Rush--Larsen steady-state/time-constant pairs, so the manufactured GPU
test should use Euler and refine its substeps.

The `monodomainPseudoECG` case retains the
`manufacturedAnisotropicMonodomainVerifier`. Its preprocessing imports the
analytic `Vm,u1,u2,u3` fields; its PDE source balances the analytic reference;
and its postprocessing measures the `Vm,u1,u2` errors. The model must retain
`verificationFamily() == "monodomainFDAManufactured"`, the correct geometric
dimension, and matching field names/order. The optional scalar
`setManufacturedSourceTerm` callback is also used by the separate 1D--3D
coupling workflow; that callback needs its own parity work before claiming
full manufactured-model API parity.

## Validation sequence

1. Confirm all state and current fields in the registered host batch evaluator
   and CUDA Euler kernels match the scalar equations in a small case.
2. Compare scalar RKF45 and actual CUDA Euler at matched ionic substeps.
   Require CUDA device selection in the solver log and finite state/current
   fields. A first paired 3D smoke case has completed; its `Vm`, `u1`, and
   `u2` errors are recorded in
   `/tmp/cardiac_manufactured_cuda_smoke_20260928`.
3. Run the `monodomainPseudoECG` temporal convergence specification
   on a fixed mesh. Compare `Vm`, `u1`, `u2`, and pseudo-ECG norms while
   refining `deltaT` and ionic substeps; verify error decreases and measure
   observed temporal order.
4. Run the spatial convergence study with temporal error kept below spatial
   error, then compare scalar and GPU observed spatial order and norms.
5. Repeat a focused SBDF2 comparison only after the Godunov GPU path passes.

## Temporal-convergence method from the GPU backup

The backup contains a completed scalar manufactured-monodomain study in
`cardiacFoAM_GPU_backup/tutorials/manufacturedSolutions/monodomainPseudoECG/setupManufacturedFDA`.
Its useful method is to hold the mesh fixed, halve the outer `deltaT`, and
measure errors against the analytical solution at the same final time. It also
varies the tissue coupling scheme (`godunov` or `sbdf2`) separately from the
`ddt(Vm)` scheme (`Euler` or `backward`). These are tissue-solver time
discretizations; they do not choose the ionic ODE integrator or its CUDA
substeps.

The backup report's 2D `N=640` data show the intended diagnostic. For
SBDF2/backward, `Vm` L2 error falls from `4.36024e-4` at `deltaT=0.025` to
`6.41272e-6` at `0.003125` with observed order near 2, then reaches a spatial
error floor around `1e-6`. Godunov/Euler stays near first order. This validates
the report's scalar scheme separation and warns against calculating an order
after the fixed-mesh spatial floor dominates. One raw case's metadata and log
are present and agree with the reported error/runtime; the bundled
`run_manifest.json` is a separate dry-run manifest, so it is not provenance
for this study. This is CPU/scalar evidence, not a GPU result.

For a batched/GPU correctness study, split the refinements so they identify
which control is responsible:

1. **ODE-substep convergence:** fix mesh, outer tissue `deltaT`, PDE coupling,
   and output times; refine only the ionic Euler substep count. Compare scalar
   RKF45, batched CPU Euler, and CUDA Euler to the exact solution and to one
   another. This tests state updates and batched/device parity.
2. **Outer temporal convergence:** fix a sufficiently fine ionic substep
   count (or otherwise bound its error well below the measured field error),
   then halve tissue `deltaT` on the same mesh. Compare `Vm`, `u1`, `u2`, and
   pseudo-ECG errors and estimate order only while error decreases
   monotonically and remains above the spatial floor.
3. **Scheme comparison:** after Godunov passes, compare SBDF2 and Godunov with
   the same mesh, outer `deltaT`, ionic substeps, initial fields, and output
   times. Treat `timeCouplingScheme sbdf2` plus the selected `ddt(Vm)` as a
   separate tissue integration path; do not label it as RL or as an ODE
   substep method.

The earlier scalar SBDF2 temporal runner is retained as a separate tissue
scheme check. Job 9361 completed the fixed-mesh manufactured ladder using
`run_manufactured_batched_temporal.py` and
`compare_manufactured_batched_parity.py` for matched CPU-batched and
CUDA-batched Godunov/Euler cases. Both backends passed first-order convergence
gates and matched each other within `1.9984e-14`; their analytic `Vm L2`
errors are 0.39–0.43% above the backup scalar Godunov/Euler reference at the
same five steps. This validates manufactured convergence/verification and
batched CPU/GPU parity, without asserting pointwise equality with scalar
RKF45 trajectories. Job 9362 then completed the SBDF2/backward ladder with
25, 50, 100, 200, and 400 Euler substeps as the tissue step halved. Both
backends passed second-order convergence gates for `Vm`, `u1`, and `u2`; their
fields and saved states matched within `2e-14`. Batched `Vm L2` errors were
about 1.9% above the scalar SBDF2 reference on the coarse levels. See
`SBDF2_MANUFACTURED_TEMPORAL.md` for the detailed results and limits.

The manufactured case is an orthogonal correctness/convergence check; its
small ODE system is not a representative GPU performance workload.

## Current run status (2026-09-29)

The batched host and CUDA implementations are compiled into `libionicModels`.
Job 9361 completed the backup-style 2D `N=640`, final-time `0.2` ladder with
Godunov/Euler and 25 ionic substeps. Both batched backends pass first-order
convergence gates; their `Vm`, current, activation, analytic fields, and
restart states agree to a maximum absolute difference of `1.9984e-14`. Their
analytic `Vm L2` errors are within 0.43% of the backup scalar Godunov/Euler
reference. This establishes manufactured convergence/verification and
batched CPU/GPU parity; it does not establish pointwise equality with scalar
RKF45.
The node also has a separate 48-rank scalar Gaur run, so job 9361 and 9362
timings must not be used as performance measurements.

The source temporal study used to set its error metrics and floor criterion is
`cardiacFoAM_GPU_backup/tutorials/manufacturedSolutions/monodomainPseudoECG/setupManufacturedFDA/TEMPORAL_CONVERGENCE_REPORT.md`.
It is adjacent to, not inside, the backup's `tutorials/benchmarkGPU` runtime
reports. The actual `benchmarkGPU` notes cover slab timing and activation
comparisons; the manufactured report contains the fixed-mesh temporal-order
study.
