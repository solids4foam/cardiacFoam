# Agent 8 — Test and CI Chaos Audit

## Scope and method

Read-only review of `.github/workflows/`, the top-level build/clean scripts,
`tutorials/Alltest-regression`, every tracked `regressionTest.sh`, the
manufactured-solution entry points, and the driverFoam pytest suite. I also read
the existing discovery reports only to determine whether their important defects
should have been caught, and independently checked each cited test/build path.
No production file was edited. The working tree was already dirty, so no test
that runs or cleans tracked tutorial cases was launched. A solver/OpenFOAM build
was not available in this shell.

## Findings

### TC-01 — the build gate can turn a failed child build into a successful job

- **Location:** `Allwmake:4-5,37-55`; exercised by
  `.github/workflows/buildAndTest.yml:56-72`.
- **Evidence:** the script promises to stop at the first error and enables
  `set -e`, but runs both child builds as `./Allwmake ... 2>&1 | tee
  log.Allwmake` without `pipefail`. Bash therefore returns `tee`'s status. I
  independently reproduced the shell semantics with `set -e; false | tee ...;
  echo MASKED`, which printed `MASKED`. The fallback oracle merely greps logs
  for the literal strings `Error ` and `Stop.`; a nonzero child that emits
  neither reaches the success message. This is contrary to the repository's
  strict regression-script convention, e.g.
  `tutorials/Alltest-regression:1-3` (`set -Eeuo pipefail`).
- **Failure masking:** direct and CI-critical. GitHub Actions only sees the
  top-level exit status. A partial build may also pass later tests by loading
  stale user libraries from the persistent OpenFOAM user installation inside
  the container.
- **Impact:** green CI can certify an incomplete build across every matrix row.
- **Severity / confidence:** **S1 — High / high**.
- **Smallest remediation:** enable `set -o pipefail` (or explicitly check
  `PIPESTATUS[0]`) and retain log scanning only for diagnostics.
- **Smallest regression test:** a shell test that substitutes a child
  `Allwmake` returning 42 with output containing neither grep token, asserts the
  top-level status is nonzero, and asserts the applications child was not run.
  Then perform clean full and forced-lightweight builds.

### TC-02 — the advertised MMS suite is substantially larger than the CI MMS gate

- **Location:** `tutorials/Alltest-regression:23-34` versus
  `tutorials/manufacturedSolutions/Allrun:1-58` and the case directories
  `monodomain1D3D`, `monodomainNonOrthoMMS`, and
  `monodomainTotalLagrangianEM`.
- **Evidence:** CI invokes only `tutorials/Alltest-regression`
  (`buildAndTest.yml:74-87`). Its fixed list includes four MMS cases:
  `bidomain`, `monodomainPseudoECG`, `eikonalECG`, and `bathBidomain`. The
  repository also has executable MMS workflows for coupled 1D-3D,
  non-orthogonal monodomain, and total-Lagrangian electromechanics; none has a
  `regressionTest.sh`, none is in the fixed list, and CI does not invoke the MMS
  aggregate `Allrun`.
- **Local contract:** `README.md:67-68` advertises
  `tutorials/Alltest-regression` as the regression command, while
  `tutorials/manufacturedSolutions/README.md:18` advertises a whole
  manufactured suite. Those surfaces currently imply more verification than
  CI performs.
- **Failure masking:** omission. Compilation succeeds while regressions in
  coupled mapping/forcing, non-orthogonal correction, convergence order, or EM
  MMS remain unobserved.
- **Impact:** scientific solver paths can regress under green CI; notably the
  coupled 1D-3D implicit path is the best existing home for a strong-coupling
  iteration regression.
- **Severity / confidence:** **S2 — Medium / high**.
- **Smallest remediation/test:** add focused, reduced-mesh `regressionTest.sh`
  wrappers for these three workflows. Require quantitative norms/orders and a
  nonzero check count, not artifact existence. Register the first two in both
  modes and total-Lagrangian EM only in full mode.

### TC-03 — arbitrary regression skips are accepted as an all-passed CI result

- **Location:** `tutorials/Alltest-regression:61-84` and
  `tutorials/NiedererEtAl2011/electroMechanicalNiedererEtAl2011/regressionTest.sh:137-146`.
- **Evidence:** exit 77 increments `skips`, but the aggregate exits zero whenever
  `failures == 0`, printing `All regression tests PASSED`; there is no allowed
  skip set, mode-specific expectation, or maximum. The build matrix knows
  whether it is `with-solids4foam` or `lightweight`
  (`buildAndTest.yml:44,79-85`) but does not pass an expectation to the runner.
  The EM case is legitimately unavailable in lightweight mode, but in full mode
  its skip should be a configuration failure, not a pass.
- **Failure masking:** a broken availability probe, renamed library, or future
  test converted to exit 77 can silently remove coverage in all six matrix
  jobs. The summary wording further hides the distinction.
- **Impact:** the full-mode electromechanics claim can be green without executing
  its only quantitative EM regression.
- **Severity / confidence:** **S2 — Medium / high**.
- **Smallest remediation/test:** let the runner accept an explicit expected-skip
  list (EM only for lightweight); fail on any unexpected skip and print
  `PASSED WITH EXPECTED SKIPS`, never `All ... PASSED`. Unit-test one expected
  and one unexpected synthetic exit-77 case.

### TC-04 — cellML2foam has no tests or CI trigger despite generating compiled code

- **Location:** `.github/workflows/driverFoamTests.yml:20-23,41-57` and
  `applications/scripts/cellML2foam/`.
- **Evidence:** the only Python workflow installs and tests driverFoam, and its
  path filter includes only `src/**`, `applications/scripts/driverFoam/**`, and
  the workflow itself. There are no maintained tests under cellML2foam. A PR
  changing only `applications/scripts/cellML2foam/**` neither triggers this
  workflow nor exercises that CLI in the unrestricted build workflow.
- **Failure masking:** complete omission of path placement, subprocess failure,
  partial-write, and generated-header pairing behavior.
- **Impact:** the independently verified PY-01/PY-02 output-path defects could
  ship while all required checks remain green.
- **Severity / confidence:** **S2 — Medium / high**.
- **Smallest remediation/test:** add a small pytest job triggered by
  `applications/scripts/cellML2foam/**`; with a neutral CWD and temporary
  `--outdir`, mock external transformation only, assert both paired headers are
  siblings under the destination and no output leaks to CWD. Add one
  nonzero-subprocess negative test.

### TC-05 — OpenFOAM-backed driver mutations always skip in the Python CI job

- **Location:** `driverFoamTests.yml:36-57`,
  `openfoam_driver/tests/test_apply_overrides.py:145-174`, and
  `openfoam_driver/tests/conftest.py:29`.
- **Evidence:** the Python job does not source/install OpenFOAM and explicitly
  sets `SKIP_ENV_DIAGNOSTICS=1`. Both tests that exercise real
  `foamDictionary` behavior call `pytest.skip` when the executable is absent;
  the remaining suite predominantly tests Python parsing/mocks. This skip is
  expected in that job but no OpenFOAM build-matrix step runs the pytest suite,
  so the real dictionary mutation integration has no CI home.
- **Failure masking:** OpenFOAM-version parser/command differences across the
  advertised v2312/v2412/v2512 matrix are invisible.
- **Impact:** driver overrides can pass unit tests yet fail on actual case
  dictionaries or a supported OpenFOAM version.
- **Severity / confidence:** **S3 — Low / high**.
- **Smallest remediation/test:** run just the two integration tests in one
  sourced OpenFOAM matrix row initially, asserting zero skips (`pytest -rs` plus
  a skip-count guard); expand versions if incompatibilities are found.

## Should the other important findings have been caught?

| Finding | Existing test should catch it? | Independently verified gap / smallest effective regression |
|---|---|---|
| OFCPP-02 build status | **Yes, but does not.** The build workflow is the primary gate. | TC-01: injected child exit status with token-free output. |
| OFCPP-01 / AD-1 stale utility selector | **No.** Regression cases run solvers, not `listCellModelsVariables`; driver catalogue tests do not execute the utility. | Build and invoke the utility on the tracked canonical single-cell and monodomain dictionaries; add a missing-`Coeffs` negative case. |
| OFCPP-03 explicit unbuilt solids path | **No.** CI always points at a prebuilt image path or unsets the variable; it never supplies an initialized unbuilt tree. | Source the resolver in four isolated shell fixtures and assert mode/path for forced lightweight, explicit built, explicit unbuilt, and unset. |
| NC-1 nested PIMPLE ownership | **No.** The coupled Purkinje regression explicitly configures `solutionAlgorithm explicit`; current implicit MMS cases have no 1D-3D coupling. | Solver-level call-count instrumentation for 1/2/3 correctors, then the reduced coupled `monodomain1D3D` MMS regression proposed in TC-02. |
| NC-2/NC-3 LandNiederer batched initialization/rate | **No.** The sole EM regression selects scalar `LandNiederer`; `LandNiedererBatched` is commented out in `electroMechanicalProperties:24-26`. | One-point zero-stimulus scalar/batched resting-state comparison plus forced-lambda values below/at/above the clamp. Require explicit state/Ta tolerances. |
| PY-01/PY-02 cellML2foam paths | **No.** No tests or workflow trigger exist. | TC-04 temporary-directory paired-output/CLI test. |
| AD-2 eikonal heterogeneity drift | **No.** eikonal ECG regression checks its fixed electrode/output references, not equivalence with ionic-region parsing or unsupported modes. | Parameterized contract tests for every supported heterogeneity mode and negative rejection tests for unsupported modes, followed by one cell-by-cell shared-dictionary equivalence test. |
| AD-3 CUDA policy drift | **No.** CI never sets `CARDIAC_ENABLE_CUDA` and uses no GPU runner. | Compile-only CUDA configuration plus mocked/helper policy tests for 0/1/many devices; reserve numerical GPU parity for a GPU-capable periodic job. |

PY-03 (stale dashboard documentation) is not a behavior regression the present
unit suite can infer reliably; a small packaging smoke test enumerating every
README-advertised CLI/import would catch it. OFCPP-04 is an API-design violation
rather than a demonstrated numerical failure and should first receive an
explicit mutable-state API contract test if repaired.

## Additional observations not promoted to findings

- `buildAndTest.yml` does cover both full and lightweight modes on all three
  stated OpenFOAM release families; that matrix is a strong baseline once build
  status and skip semantics are made trustworthy.
- Most maintained regression scripts use strict shell options, copy cases out
  of the tutorial tree, and compare quantitative references. Those are sound
  local conventions.
- The driverFoam unit suite has broad negative/assertion coverage. Its three
  expected environment-dependent skips are visible under `pytest -rs`, but CI's
  `-q` presentation makes them easy to overlook.
- `Allwclean:32-36` removes `tutorialsTest`, whereas the current regression
  runner creates `tutorialsTest-regression` (`Alltest-regression:17`). This is a
  cleanup drift (S3 at most), not a correctness gate, because each regression
  run removes its own destination first.

## Validation performed

- Static line-by-line inspection of workflows, aggregate and per-case regression
  scripts, MMS entry points, skip sites, and workflow path filters.
- Shell reproduction of the pipeline-status masking in TC-01.
- Cross-check of all important discovery findings against actual selected cases,
  algorithms, active-tension model choice, and workflow commands.
- No OpenFOAM build/tutorial run was attempted because the current shell lacks
  that environment and tracked tutorials contain pre-existing user changes.
