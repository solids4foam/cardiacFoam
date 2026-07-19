# Numerical Correctness Audit

Role: Agent 3 — Numerical Correctness Skeptic  
Branch reviewed: `no-frontend-minor-errors` (`f6b79880`)  
Review mode: independent, read-only discovery; this report is the only file written.

## Scope and authority consulted

I read `.agents/skills/cardiacfoam/PROJECT_MEMORY.md`, `README.md`,
`src/electroModels/README.md`, and `src/ionicModels/README.md` before judging
implementation behavior. I then traced the electrophysics advance schemes,
myocardium and conduction-domain advances, monodomain/bidomain kernels, ECG
quadrature, electromechanical sequencing, and scalar/batched active-tension
paths. Generated equation headers were treated as scientific/generated code;
no equation, constant, or production source was edited.

The checkout already contained unrelated modified and untracked user files.
They were not changed. No executable OpenFOAM reproduction was attempted because
this audit context did not establish a sourced, built OpenFOAM/cardiacFoam
runtime; the proposed reproductions below are therefore validation recipes, not
claims of completed dynamic reproduction.

## Confirmed / high-confidence findings

### NC-1 — The PIMPLE strong-coupling loop is consumed recursively, preventing coupling refresh on each corrector

- **Evidence:**
  - `src/electroModels/core/advanceSchemes/pimpleStaggered/pimpleStaggeredElectrophysicsAdvanceScheme.C:67-80` owns `while (pimplePtr->loop())`, refreshes both coupling directions, advances the conduction domains, and calls `myocardium.advance(..., pimplePtr)` inside that loop.
  - `src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.C:512-522` performs the ionic reaction step and dispatches the implicit diffusion solve with the same `pimpleControl` object.
  - `src/electroModels/myocardiumModels/monodomainSolver/monodomainSolver.C:145-155` independently calls `while (pimple.loop())` on that same object. The sibling bidomain implementation has the same nested ownership at `src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C:307-317`.
  - The public configuration description explicitly calls this scheme “strongly coupled with PIMPLE iteration” and suitable for bidirectional coupling: `applications/scripts/driverFoam/openfoam_driver/dict_entries.py:340-343`.
- **Affected data flow/equation:** intended corrector sequence is
  `coupling(n) -> conduction(n+1) -> coupling(n+1) -> reaction/diffusion(n+1)`,
  repeated for every PIMPLE corrector. Observed control flow enters the outer
  `pimple.loop()` once, then the myocardium diffusion kernel consumes subsequent
  `pimple.loop()` states internally while solving
  `chi*Cm*dVm/dt = div(G grad(Vm)) - chi*Cm*Iion + source`. Those inner solves do
  not re-run either coupling preparation or the conduction-domain update.
- **Expected versus observed:** each strong-coupling corrector should recompute
  both domain states and coupling currents. Instead, correctors after the first
  are diffusion-only; depending on `pimpleControl` loop semantics/version, the
  outer loop is exhausted or advances with an already-mutated loop state. Thus
  the advertised strong iteration is not performed.
- **Minimal reproducible configuration:** use a monodomain myocardium plus a
  bidirectional reaction-diffusion PVJ coupler, set
  `electrophysicsAdvanceScheme pimpleStaggeredElectrophysicsAdvanceScheme`,
  `solutionAlgorithm implicit`, and compare `nOuterCorrectors 1` with 2 and 3.
  Instrument/count `prepareConductionCouplings`, `advanceConductionDomains`,
  `prepareMyocardiumCouplings`, and the linear solve. For `N` outer correctors,
  the first three counts should each be `N`; current structure predicts fewer
  coupling refreshes than diffusion solves.
- **Comparison and tolerance:** first use an exact call-count assertion (zero
  tolerance): every requested outer corrector must execute exactly one complete
  coupling/conduction/myocardium sequence. Then compare terminal `Vm1D`, tissue
  `Vm`, PVJ current, and activation times against a reference implementation
  with one owner of the PIMPLE loop. Require `Linf(Vm) <= 1e-8 V`,
  `Linf(activationTime) <= max(1e-8 s, 0.01*dt)`, and relative `L2` PVJ-current
  error `<= 1e-6`; tighten after observing cross-version solver noise.
- **Reference-result impact:** yes. Correcting loop ownership can change coupled
  Purkinje–myocardium trajectories, activation times, and accepted coupled
  reference outputs. Uncoupled implicit monodomain/bidomain results may also
  change in iteration history but should converge to equivalent values.
- **Severity / confidence:** **S1 — High / high**.
- **Minimal remediation:** establish exactly one owner for `pimple.loop()`.
  Prefer the orchestration scheme as owner, because coupling refresh belongs
  there; make the per-kernel overload perform one implicit linear solve per
  corrector. Confirm OpenFOAM v2312–v2512 loop semantics before implementation.
- **Required validation:** a solver-free mock/call-count test for loop ownership;
  coupled 1D–3D regression at 1/2/3 correctors; monodomain and bidomain implicit
  MMS checks; serial and decomposed runs; all supported OpenFOAM versions.

### NC-2 — `LandNiedererBatched` skips the resting-Ca preconditioning applied by the scalar model

- **Evidence:**
  - The electromechanical constructor always calls the virtual hook when Cai is required: `src/electroMechanicalModels/sequentialElectroMechanical/sequentialElectroMechanical.C:87-95`.
  - The scalar model overrides it at `src/activeTensionModels/LandNiederer/LandNiederer.C:186-227`, integrates for 1000 ms at resting Cai, and explicitly explains that the raw states are equilibrium only at Cai=0 and otherwise cause a spurious global tension transient (`:191-195`).
  - The base hook is a no-op at `src/activeTensionModels/activeTensionModel/activeTensionModel.H:465`.
  - `src/activeTensionModels/LandNiedererBatched/LandNiedererBatched.H` declares no `preconditionToRestingState` override, and its constructor initializes every cell directly from the raw generated prototype states at `src/activeTensionModels/LandNiedererBatched/LandNiedererBatched.C:120-127`.
- **Affected data flow/equation:** resting ionic `Cai` is supplied to the active
  tension model; Land–Niederer evolves `Ca_TRPN`, `TmBlocked`, `XW`, and `XS`,
  from which algebraic `Ta` is computed. Scalar starts the physical run from the
  Cai-conditioned state. Batched starts from the Cai=0 generated state, despite
  receiving nonzero resting Cai during the physical run.
- **Expected versus observed:** mathematically equivalent scalar and batched
  models should begin from the same resting equilibrium for identical constants,
  Cai, and lambda. Batched instead exhibits the exact initial transient the
  scalar path was designed to remove.
- **Minimal reproducible configuration:** one integration point, constant
  `Cai=0.0002 mM`, `lambda=1`, no electrical stimulus, default
  `preconditioningTime 1000`, and otherwise identical `LandNiederer` versus
  `LandNiedererBatched` dictionaries. Record the six states and `Ta` at physical
  time zero and over the first 100 ms.
- **Comparison and tolerance:** after applying identical preconditioning,
  require initial state `Linf <= 1e-8` and `|Ta_scalar-Ta_batched| <= 1e-6 kPa`
  at t=0. Across the first 100 ms, compare batched Euler with enough substeps to
  separate integrator error from initialization; target relative `L2(Ta) <= 1e-3`
  and peak absolute difference `<= 1e-3 kPa`.
- **Reference-result impact:** yes, for every electromechanical reference using
  `LandNiedererBatched`; initial tension and early deformation may change.
  Scalar Land–Niederer references should not change.
- **Severity / confidence:** **S1 — High / high** (wrong initial physical state
  and potentially global spurious active stress).
- **Minimal remediation:** implement the same configurable resting-Cai
  preconditioning contract for the batched core, using its declared time-scale
  convention and a numerically adequate substep count. Do not copy final scalar
  states blindly if constants or integrator settings can differ.
- **Required validation:** zero-stimulus resting test, scalar/batched state and
  tension comparison, multi-cell uniformity test, CPU/OpenMP/CUDA parity where
  the backend exists, and rerun electromechanical references with explicit lead
  approval for changed baselines.

### NC-3 — Scalar and batched Land–Niederer apply different stretch-rate inputs

- **Evidence:** scalar computes `(lambda-prevLambda)/dt` and clamps it to
  `[-20,20] s^-1` before conversion to `ms^-1` at
  `src/activeTensionModels/LandNiederer/LandNiederer.C:296-305,322-325`.
  Batched computes the same quotient without the clamp at
  `src/activeTensionModels/LandNiedererBatched/LandNiedererBatched.C:143-155`
  and passes the unbounded value (times `1e-3`) into the identical generated
  equation core at `:173-185`.
- **Affected data flow/equation:** deformation gives `lambda`; finite difference
  gives `lambda_rate`; `AV_lambda_rate` enters `LandNiederer2017computeVariables`
  and therefore the crossbridge/distortion state rates and algebraic tension.
- **Expected versus observed:** the paired implementations should expose the
  same model policy for the same lambda history. A step from `lambda=1.0` to
  `1.1` over `dt=0.001 s` produces `100 s^-1`: scalar supplies `20 s^-1`, batched
  supplies `100 s^-1`. This is a fivefold input difference before integration.
- **Minimal reproducible configuration:** one point with identical preconditioned
  state and constant Cai; call once at lambda 1, then at lambda 1.1 with
  `dt=1e-3 s`. Export `AV_lambda_rate`, state rates, and `Ta` from both variants.
- **Comparison and tolerance:** `AV_lambda_rate` is a contract value and should
  agree to machine precision (`abs <= 1e-14 ms^-1`). With matched explicit
  integration/substeps, require state-rate relative `Linf <= 1e-10` and tension
  relative error `<= 1e-6` for the single step. Also test values just below,
  at, and above both clamp boundaries.
- **Reference-result impact:** yes for batched Land–Niederer simulations that
  experience numerical lambda jumps or rates over 20/s. Smooth low-rate cases
  should be unchanged.
- **Severity / confidence:** **S2 — Medium / high**. The divergence is certain;
  physiological impact depends on whether a run crosses the clamp.
- **Minimal remediation:** centralize or identically implement the stretch-rate
  policy in both variants. Because the clamp itself changes scientific behavior,
  first obtain maintainer approval on whether the scalar clamp is the intended
  contract; do not silently remove it or add it to the batched path.
- **Required validation:** boundary-value unit test for rate policy;
  scalar/batched forced-lambda trace comparison; a realistic electromechanical
  case reporting the maximum unclamped rate; CPU/OpenMP/CUDA equivalence.

## Investigation items, not confirmed defects

- The scalar active-tension models use adaptive OpenFOAM ODE solvers while the
  batched family defaults to one explicit Euler substep. That alone is not a
  defect, but there is no evidence in the reviewed paths of a model-by-model
  error bound establishing that default `batchedSubsteps=1` is adequate. A
  timestep-refinement comparison should be added before claiming scalar/batched
  scientific equivalence.
- Activation time in the voltage solvers is recorded at the end of the step
  when `Vm >= threshold`, rather than interpolating the crossing. This is a
  first-order, `O(dt)` timing convention, not necessarily a bug. Reference tests
  should explicitly permit at most one timestep of timing quantization and avoid
  presenting it as sub-timestep accuracy.

## Prioritization

1. Resolve NC-1 loop ownership before trusting strongly coupled PIMPLE results.
2. Resolve NC-2 before using batched Land–Niederer for reference electromechanics.
3. Decide and test the intended lambda-rate policy (NC-3) before harmonizing it.

