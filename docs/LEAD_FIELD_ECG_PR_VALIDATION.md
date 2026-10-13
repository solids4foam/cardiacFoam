# Lead-field ECG PR validation

Tested numerical source: `0c1fe468fd14c478ef776bf189877ab7f28f19b5`.
Feature branch: `codex/lead-field-ecg-integration`.
Dependency/base: PR #59, `feat/mono-bi-eikonal-verification` at
`fd0e548b113c991126580f4262bd4975dd3fe226`.

The PR adds passive monodomain/bidomain recording, uniform-template eikonal
recording, static construction, independent verification and maintained tutorials. It targets the PR #59 branch
so its review diff excludes the parent propagation and wall-discretization work.

## Final executed checks

| Platform/control | Result |
|---|---|
| macOS, native OpenFOAM v2412, DP/Int32, Clang, lightweight | Complete Allwmake passes. |
| Configuration harness on v2412 | 113 checks pass, zero skipped. |
| Bidomain lead/pseudo MMS on v2412 | 8 numerical bounds plus independent dual-cohort/quadrature checks pass. |
| Monodomain lead/pseudo MMS on v2412 | 6 numerical bounds plus independent dual-cohort/quadrature checks pass. |
| Existing eikonal ECG regression on v2412 | 25 checks pass. |
| Existing monodomain pseudo-ECG regression on v2412 | 34 checks pass; references unchanged. |
| Linux ARM64, OpenFOAM v2312, GCC, lightweight | Complete Allwmake passes. |
| Both combined MMS regressions on v2312 | Same 8/6 numerical bounds and independent dual checks pass, no skips. |
| New eikonal dual MMS on v2412 and v2312 | Three numerical bounds, independent surface/reference/source diagnostics, four template controls and three rejection controls pass. |
| Paired native serial eikonal runs with/without reciprocal observation | Written activation field and pseudo-ECG trace are byte-identical. |
| Native serial/six-rank eikonal trace comparison | Maximum differences 1e-10 V reciprocal and 1e-11 raw pseudo. |
| Include hierarchy and tutorial index | Pass; 22 indexed regression cases. |
| Project headers and Markdown | Changed C/H normalizer leaves files unchanged; modified Markdown passes lint. |

The final native sweep runs six maintained case regressions and finishes with
no failures/skips. It repeats the complete configuration harness after SVD,
format-2 topology/filename and restart updates. The mixed-observer checks compare
each recording with its own independent reference, not with the other observer.
The final Linux v2312 sweep runs the monodomain, bidomain and eikonal lead-field
regressions, with no failures/skips. The extended source/gradient diagnostic
checker also passes when run directly against both saved eikonal outputs.

The configuration controls exercise physical field modes, source invariance,
serial/MPI construction and projection, Wilson/reference weights, heart-only
export, static recovery, sampling, restart, and expected input failures.
New failures include duplicate/reserved output filenames, a changed physical
patch type and a myocardial baffle preserving cell geometry. The rounding
control checks a saved folder rounding down from the trace timestamp, preserves
existing rows and rejects a genuinely older restart.

The v2312 compatibility check uses the common three-argument boundary exchange
API. Python 3 was installed in the minimal official OpenFOAM image to execute
the standard-library independent checker; missing Python was an environment
failure, not a numerical pass. The final sweep passes all three lead-field
regressions.

On macOS the rounding hook sources the established RunFunctions setup to
restore library lookup after SIP strips DYLD variables. No tracked tutorial
dictionaries are modified by an ad-hoc runner: maintained scripts generate
local includes or use the harness/private run copies.

## Reproduction

Use an isolated lightweight installation and the normal Allwmake entry point.
Initialize the OpenFOAM environment before enabling strict shell error handling.
On macOS keep the selected installation ahead of both normal and
FOAM_LD_LIBRARY_PATH library lookup, as required by RunFunctions.

```bash
FORCE_LIGHTWEIGHT_PHYSICSMODEL=1 ./Allwmake -j 6
CARDIAC_REGRESSION_BUILD_MODE=lightweight ./tutorials/Alltest-regression \
  --only '^(configurationChecks|manufacturedSolutions/(monodomainLeadFieldECG|bidomainLeadFieldECG|monodomainPseudoECG|eikonalECG|eikonalLeadFieldECG))$'
```

The voltage-based tutorials accept serial/fine/timestep variants through
their own Allclean/Allrun entry points. The eikonal variant forwards the
maintained cube options, including resolution and mesh type. Generated fields/meshes and full logs are not
versioned. Prior two-level scientific evidence remains a local archive;
its monodomain/bidomain source and projection kernels are unchanged by the
uniform-template extension. The existing voltage projection now delegates
to the same validated gradient-integral method used by eikonal recording.

## Review and remaining scope

Earlier independent source review found and resolved topology, filename and
restart issues; the IO/common-version exchange pass reported no P0/P1 blocker.
The PR remains a bounded implementation change. It does not establish a
complete discrete adjoint, general rotated/nonorthogonal conductor accuracy,
clinical lead equivalence, anatomical transfer or regional-template eikonal
recording.
Heart-only mode cannot revalidate an absent torso. Fingerprints do not encode
every coupled transform option; AMI and rotating meshes remain outside the
verified set. Regenerate older format-1 fields.

The repository's full CI matrix, including full solids4foam builds and the
remaining tutorials, is distinct from these executed local checks. The papers
and private torso dataset are not part of the code PR.

## Uniform-template eikonal extension

The [extension guide](EIKONAL_LEAD_FIELD_ECG.md) describes the shared template
source, new observer, independent surface-integral verifier and maintained
cube variant. The eikonal source uses `-U'(t-tau) grad(tau)` and imported
conductivity-weighted sensitivities. It does not use the placeholder Vm.

On the 40³/six-rank cube, 201 common samples give 1.68% normalized RMS error
for the prescribed reciprocal cohort and 2.50% for the selected raw
pseudo-ECG cohort. The surface reference Q64/Q96 difference is about
2e-15 V. These are different recording functionals at one selected resolution,
not a model ranking or a convergence-order estimate. The unit/state fixture
covers all four supported templates independently of the PDE solve.

Native/Linux build logs, complete final sweeps, source-independence snapshots,
serial/MPI comparisons and diagnostic controls are retained locally under
`/tmp/cardiacFoam-lead-field-current-numerics/eikonal`. Initial development
runner failures (missing copied base and container MPI root policy) were
resolved before the final successful sweeps; neither was a numerical failure.
The Linux container needs Python for maintained regression checks and the
usual OpenMPI root-run environment when launched as root.

## Additional coarse eikonal controls

Native serial N=10/20 runs of the maintained eikonal dual variant complete.
Together with the existing N=40 serial case, their reciprocal normalized RMS
errors are 20849.451%, 254.243% and 1.679%. The raw pseudo-ECG E3–E5 cohort
errors are 29178.249%, 253.169% and 2.498%, against its own reference.
Independent exact-source midpoint projection has reciprocal RMS errors
20.4783604, 0.250080444 and 0.00164665117 V. The severe coarse loss persists
without propagation or reconstructed-gradient error. No asymptotic rate is
fitted, N=80 is unexecuted and the N=40 regression bounds are unchanged.

Paper 1 source includes these controls; independent eikonal, earlier lead-field,
dual-observer/shell, headline-results and manuscript-structure validators pass.
Docker HTML/PDF render and inspection of the new tables/equations pass.
Local evidence has source, binary and data hashes; immutable publication
release remains open. See [eikonal details](EIKONAL_LEAD_FIELD_ECG.md).

The [matched torso guide](MATCHED_HEART_TORSO_ECG.md) documents existing chest
sites and proposed torso equivalents for absent limbs. VTK inclusion checks
are geometry evidence only, not production ownership or ECG equivalence tests.
