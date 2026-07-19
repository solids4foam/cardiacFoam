# Agent 2 — OpenFOAM C++ Maintainer audit

## Scope and method

Read-only review of project-owned C++ and build logic under `src/`,
`applications/solvers/`, `applications/utilities/`, `modules/physicsModel/`, the
top-level `Allwmake`, component `Allwmake` scripts, and
`etc/resolveSolids4Foam.sh`. The architecture contract was established from
`README.md`, `src/electroModels/README.md`,
`src/electroModels/ARCHITECTURE.md`, `src/electroModels/core/ARCHITECTURE.md`,
`src/ionicModels/README.md`, and
`.agents/skills/cardiacfoam/PROJECT_MEMORY.md` before judging implementation.

The worktree was already dirty. This audit did not alter production code or
the user's existing changes. `wmake` was not available in the current shell,
so build/run validations below remain required rather than claimed as run.

## Findings

### OFCPP-01 — `listCellModelsVariables` does not read the canonical selector/coefficients contract

- **Location/symbol:**
  `applications/utilities/listCellModelsVariables/listCellModelsVariables.C:44-75`,
  `electroModelDict`; reporting repeats the stale contract at lines 106-114.
- **Classification:** hand-maintained project utility.
- **Evidence:** `electroModelDict` looks only for the obsolete top-level key
  `electroModel`. When it is absent it returns the entire `electroProperties`
  dictionary (lines 73-75). `main` then passes that result to
  `ionicModel::New` (lines 266-270), whose required `ionicModel` lookup is in
  `src/ionicModels/ionicModel/ionicModel.C`, symbol `ionicModel::New`. In a
  current case, `ionicModel` is nested under `<myocardiumSolver>Coeffs`, so the
  utility fails to find it.
- **Canonical local reference:** `src/electroModels/core/electroModel.C:147-152`
  defines `myocardiumSolver` as the single canonical selector;
  `src/electroModels/core/electrophysiologyModel/electrophysiologyModel.C:64-101`
  reads that key and selects its `<type>Coeffs`; current tracked cases such as
  `tutorials/NiedererEtAl2011/NiedererEtAl2011verification/constant/electroProperties:18-26`
  use exactly that structure. The sibling utility
  `applications/utilities/ionicHeterogeneityProbe/ionicHeterogeneityProbe.C:67-77`
  already implements the current lookup.
- **Impact:** the documented utility cannot inspect ordinary current
  monodomain, bidomain, eikonal, or single-cell cases and terminates on the
  missing flat `ionicModel` entry. Its diagnostic also tells users to add an
  obsolete public key, encouraging configuration drift.
- **Severity/confidence:** **S1 — High**, **high** confidence.
- **Minimal remediation:** make `electroModelDict` require/read
  `myocardiumSolver`, require `<selected>Coeffs`, and return that subdictionary.
  Update `writeHeader` and its diagnostics to use the same terminology. If
  legacy `electroModel` compatibility is intentionally retained, place it
  strictly after the canonical path and label it deprecated; do not treat a
  canonical dictionary as flat.
- **Required validation:** build the utility on the oldest and newest supported
  OpenFOAM versions; run it against one tracked `monodomainSolver` case and the
  tracked `singleCellSolver` heterogeneity-probe case; assert it reports the
  configured ionic model and exits zero. Add a negative case with a missing
  `<type>Coeffs` dictionary and verify a `FatalIOError`/clear dictionary-context
  diagnostic.

### OFCPP-02 — top-level build can report success after a failed child build

- **Location/symbol:** `Allwmake:37-55`; child scripts
  `src/Allwmake:22-30` and `applications/Allwmake:4-5`.
- **Classification:** hand-maintained project build scripts.
- **Evidence:** the top-level script enables `set -e` (line 5), but each child
  build is piped through `tee` (lines 37 and 40) without `set -o pipefail` or a
  `PIPESTATUS` check. Bash therefore uses `tee`'s status for each pipeline, so a
  nonzero child `Allwmake` status is masked. The later success test does not use
  exit status; it searches logs only for the literal strings `Error ` and
  `Stop.` (lines 44-46). A compiler/tool failure with different wording, a
  signal, or an early shell failure can therefore reach the “installation was
  a success” message. The component scripts also do not fail fast, allowing a
  failed library/solver to be followed by later successful commands.
- **Canonical local reference:** repository-owned regression and case scripts,
  for example
  `tutorials/NiedererEtAl2011/NiedererEtAl2011verification/regressionTest.sh:2`
  and
  `applications/scripts/driverFoam/openfoam_driver/scripts/run_case.sh:3`, use
  `set -euo pipefail` for precisely this shell contract. The top-level comment
  `Allwmake:4` explicitly promises “Stop at first error.”
- **Impact:** CI and users may accept a partial or stale installation as a
  successful build. Subsequent runs can load old libraries or fail later in a
  way disconnected from the actual compilation error.
- **Severity/confidence:** **S1 — High**, **high** confidence.
- **Minimal remediation:** enable `pipefail` before the logged pipelines (or
  explicitly capture/check the child status from `PIPESTATUS[0]`) and make
  project-owned component `Allwmake` scripts propagate the first failed
  `wmake`. Keep log scanning only as supplemental diagnostics, not the success
  oracle.
- **Required validation:** in a disposable worktree inject a deterministic
  failing child command that prints neither `Error ` nor `Stop.`; verify
  top-level `./Allwmake` returns nonzero and does not invoke the subsequent
  application build. Then perform clean lightweight and full builds on at least
  OpenFOAM v2312 and v2512.

### OFCPP-03 — an explicitly configured unbuilt solids4foam tree is accepted as a full installation

- **Location/symbol:** `etc/resolveSolids4Foam.sh:20-25,51-62`, explicit
  `SOLIDS4FOAM_INST_DIR` branch.
- **Classification:** hand-maintained project build-mode resolver; the
  `modules/solids4foam` contents themselves are external/submodule content and
  are not proposed for modification.
- **Evidence:** the script defines `_s4fLnHeader` as the representative proof
  that a solids4foam tree is built (lines 23-25), and auto-discovery correctly
  tests it at lines 65-78. In contrast, when `SOLIDS4FOAM_INST_DIR` is already
  set, lines 51-59 test only `_solids4FoamHeader`, the source-tree
  `solidModel.H`. Any initialized but unbuilt solids4foam checkout contains that
  file, so the resolver announces “Using full solids4foam installation” and
  sets full mode even though its own compiled/`lnInclude` prerequisite is not
  met. The `else` diagnostic at lines 60-61 also says “not compiled” despite
  testing a source header that does not establish compilation.
- **Canonical local reference:** the same resolver's auto-discovery contract at
  `etc/resolveSolids4Foam.sh:65-89` requires `_s4fLnHeader` and deliberately
  falls back when the bundled tree is present but unbuilt. `README.md` describes
  full mode as using a solids4foam installation, while project memory states
  that the resolver prefers a **built** install.
- **Impact:** builds with a commonly exported path to a fresh/unbuilt checkout
  select the wrong dependency mode and fail on missing generated include links
  or libraries instead of safely using lightweight mode. This is a full versus
  lightweight portability divergence.
- **Severity/confidence:** **S2 — Medium**, **high** confidence.
- **Minimal remediation:** apply the same `_s4fLnHeader` built-tree predicate to
  explicit `SOLIDS4FOAM_INST_DIR` as to auto-discovery. Preserve the separate
  source-header check only to improve the diagnostic (uninitialized/invalid
  versus present-but-unbuilt).
- **Required validation:** shell-test four isolated environments: unset path
  with no submodule, explicit initialized-but-unbuilt checkout, explicit built
  install, and `FORCE_LIGHTWEIGHT_PHYSICSMODEL=1`. Assert both exported mode
  variables and chosen directory, then run one lightweight and one full build.

### OFCPP-04 — manufactured graph verification mutates ionic state through a const-only I/O API

- **Location/symbol:**
  `src/verificationModels/monodomainVerification/manufacturedGraphVerifier.C:125-142`,
  `manufacturedGraphVerifier::preProcess`; interface at
  `src/ionicModels/ionicModel/ionicModel.H:114-139`.
- **Classification:** verifier is hand-maintained project code. Most concrete
  scalar ionic wrappers expose generated state arrays through hand-maintained
  wrapper headers; equation headers such as `*_20xx.H` are generated outputs
  and should not be edited for this issue.
- **Evidence:** `preProcess` receives a mutable `ionicModel*`, but the only state
  accessor is `ioStatesPtr() const`, returning
  `const PtrList<scalarField>*`. The verifier promises that this is “mutable
  state storage” in its diagnostic (lines 128-132), then removes constness with
  `const_cast` (lines 140-141) and writes three states. By contrast, the same
  base API explicitly distinguishes `ioConstantsPtr() const` from
  `ioMutableConstantsPtr() const` at
  `src/ionicModels/ionicModel/ionicModel.H:131-139`; constant overrides consume
  that mutable hook in `ionicModel.C`, symbol
  `applyIonicConstantOverrides`.
- **Impact:** the verification layer depends on an undocumented assumption that
  every returned state container is actually mutable. A future model may
  legitimately expose const-backed or synchronized storage, making the cast
  invalid or bypassing required upload/dirty-state bookkeeping. Current scalar
  manufactured models likely work because their underlying containers are
  mutable, so this is not classified as a demonstrated numerical defect.
- **Severity/confidence:** **S3 — Low**, **high** confidence that the API
  contract is violated; **medium** confidence of current runtime impact.
- **Minimal remediation:** add an explicit mutable-state hook (mirroring
  `ioMutableConstantsPtr`) or a narrow virtual state-initialization operation,
  implement it only for models that support verifier mutation, and make the
  verifier reject models without that capability. Do not edit generated
  equation files. For batched/GPU implementations, the hook must define host /
  device synchronization rather than exposing a stale mirror.
- **Required validation:** compile all scalar and batched model families; run
  the coupled 1D/3D manufactured verification in serial and parallel; compare
  pre-change/post-change norms bitwise where possible (otherwise relative and
  absolute tolerance `1e-12` for this API-only refactor). Add a capability
  negative test proving a const-only model fails cleanly rather than being
  cast.

## Audited items with no finding

- Runtime registrations for top-level spatial workflows are aligned:
  `electrophysiologyModel` is registered under `monodomainSolver`,
  `bidomainSolver`, and `eikonalSolver`; `singleCellSolver` registers directly
  in the `electroModel` table.
- A source-list sweep found no project-owned, non-template `.C` source omitted
  from its applicable `Make/files` in the reviewed scope.
- `modules/solids4foam` is classified as external/submodule content and was not
  audited for style or proposed for direct repair. `modules/physicsModel` is
  project-owned fallback code and was included in the interface/build review.
- Ionic equation headers (`Model_year.H` and batch equation headers) and CUDA
  kernels derived from those equations are generated/generated-family content;
  no style-only changes are proposed. Wrapper `.H/.C`, factories, metadata,
  and `Make/files` remain project-maintained contracts.

## Priority summary

1. **OFCPP-01 (S1):** current public utility is incompatible with current case
   dictionaries.
2. **OFCPP-02 (S1):** build success reporting is not trustworthy.
3. **OFCPP-03 (S2):** explicit dependency selection diverges from the built-tree
   rule.
4. **OFCPP-04 (S3):** verifier mutation needs an explicit ownership/mutability
   contract before additional model families rely on it.
