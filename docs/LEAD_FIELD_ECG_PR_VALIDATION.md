# Lead-field ECG PR validation

Original feature checks below use numerical source
`0c1fe468fd14c478ef776bf189877ab7f28f19b5`.
The idealized application is added in `6b47fb42`; its test sequence and
final-wrapper qualifications are recorded separately below.
Feature branch: `codex/lead-field-ecg-integration`.
Dependency/base: PR #59, `feat/mono-bi-eikonal-verification` at
`04e6bd84cf22f0e0e8510fad0d1c806c2c728d81`, merged after the native application
records. Earlier feature checks used its predecessor `fd0e548b`.

The PR adds passive monodomain/bidomain recording, uniform-template eikonal
recording, static construction, independent verification and maintained tutorials. It targets the PR #59 branch
so its review diff excludes the parent propagation and wall-discretization work.

## Original feature checks

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
fitted and the N=40 regression bounds are unchanged. The subsequently
requested N=80 extension is recorded below.

Paper 1 source includes these controls; independent eikonal, earlier lead-field,
dual-observer/shell, headline-results and manuscript-structure validators pass.
Docker HTML/PDF render and inspection of the new tables/equations pass.
Local evidence has source, binary and data hashes; immutable publication
release remains open. See [eikonal details](EIKONAL_LEAD_FIELD_ECG.md).

The [matched torso guide](MATCHED_HEART_TORSO_ECG.md) documents existing chest
sites and the declared torso equivalents for absent limbs. The subsequent
idealized application controls below exercise native ownership and recording.

## N=80 extension

The unchanged native serial dual case completes on 512,000 cells and meets
its 1e-8 outer criterion in 1128 iterations. The independent checker passes
201 common samples, four template controls, output/verifier agreement,
opposite polarity, separate cohort bounds, source/gradient bounds and surface
reference checks without changing thresholds.

| Diagnostic | N=40 | N=80 |
|---|---|---|
| Reciprocal relative RMS | 1.679% | 0.394% |
| Raw pseudo E3–E5 relative RMS | 2.498% | 0.568% |
| Activation L2 error (s) | 0.000159934 | 0.0000390149 |
| Activation-gradient L2 error (s/m) | 0.00564762 | 0.00203656 |
| Exact-source midpoint reciprocal RMS error (V) | 0.00164665117 | 0.000387544686 |

The production reciprocal RMS error is 0.000386975146 V. Fine-pair functional
reductions of 4.26/4.40 are measured observations, not a settled asymptotic
claim. Effective source/sampling/scheme inputs and deployed binaries match
the earlier run; the mesh parameter is the intended change. No tighter outer-
criterion repeat bounds N=80 iterative error. The local hashed N80 supplement
preserves the original three-level archive; manuscripts/raw data remain outside
this code PR. No production source was changed for this extension.

## Idealized heart/torso application extension

Shared LFS inputs from `feat/bath-transmission-mms` at `f5d8011e` supply the
299482-cell heart/blood/lungs/chest-wall/torso mesh. No bath solver changes
are merged. The opt-in `electroHeartECG` tutorial adds source variants and
simultaneous lead-field/pseudo recording. `replayLeadFieldECG` holds each
stored Vm fixed and invokes production passive recovery/torso sampling with
the same physical inputs and Wilson definitions.

The supplied heart has 122 rank-deficient myocardial face stencils. The
builder's explicit `pointFallback` option extends only those stencils with
point-connected myocardial neighbours. Its default remains unchanged and
the required-rank guard remains active. Serial and six-rank actual-mesh
affine controls both have maximum error 2.88e-11. Both utilities compile on
native v2412 and Linux ARM64 v2312; application execution is native v2412.
Four existing serial/MPI configuration/recovery controls pass with zero skips.
The opt-in tutorial's output/log/source/Wilson contracts pass for both fresh
source variants (two grouped checks each). Final-wrapper positive replay and
deliberate mismatched-body rejection also pass.

Fresh monodomain: eight samples at 5–40 ms, six-lead relative RMS
forward/reciprocal difference 8.686%. Fresh heart-only bidomain: two samples
at 0.1/0.2 ms, difference 7.874%; this is a software-path check only. A separate
older coupled source, with byte-identical mesh but uncertified producer
revision, supplies eight 5–40 ms replay samples, difference 7.741%.
The fresh full-window bidomain attempt was intentionally interrupted before
its first snapshot and is not a completed result.

Online/replay reciprocal recordings agree to output rounding. Independent
checks form Wilson leads from raw torso samples and verify source preservation.
Tightening recovery's absolute/relative fixed-point criteria changes the
monodomain forward trace by at most 3.63e-13 V; it does not explain the 8.686%
difference. This measured finite-mesh discrepancy remains open, without an
invented accuracy bound or an equivalence claim. Review corrected-operator
and gradient consistency and mesh sensitivity before a stronger claim.
See the [application report](IDEALIZED_HEART_ECG_VERIFICATION.md) for scope,
archived-source qualifications and evidence provenance.

After integrating parent update `04e6bd84`, the affected Linux v2312 library
rebuild and independent two-snapshot frozen-source replay pass. Source,
Wilson and online/replay checks pass, with maximum native/Linux forward
difference 1.62e-14 V. This selected integration control does not repeat every
earlier regression or certify Linux source propagation. The parent update
reads sub-mesh dictionaries on versions before 2406; native v2412 operators
are unchanged by that conditional fix.
