# Lead-field ECG PR validation

Tested numerical source: `341695f782d892908d2ca977d600877267af69c8`.
Feature branch: `codex/lead-field-ecg-integration`.
Dependency/base: PR #59, `feat/mono-bi-eikonal-verification` at
`fd0e548b113c991126580f4262bd4975dd3fe226`.

The PR adds passive monodomain/bidomain recording, static construction,
independent verification and maintained tutorials. It targets the PR #59 branch
so its review diff excludes the parent sealed-wall/eikonal work.

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
| Include hierarchy and tutorial index | Pass; 21 indexed regression cases. |
| Project headers and Markdown | Changed C/H normalizer leaves files unchanged; modified Markdown passes lint. |

The final native sweep runs five maintained case regressions and finishes with
no failures/skips. It repeats the complete configuration harness after SVD,
format-2 topology/filename and restart updates. The mixed-observer checks compare
each recording with its own independent reference, not with the other observer.

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
failure, not a numerical pass. The final rerun passes both regressions.

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
  --only '^(configurationChecks|manufacturedSolutions/(monodomainLeadFieldECG|bidomainLeadFieldECG|monodomainPseudoECG|eikonalECG))$'
```

The new tutorials also accept serial/fine/timestep variants through their
own Allclean/Allrun entry points. Generated fields/meshes and full logs are not
versioned. Prior two-level scientific evidence remains a local archive;
its source/projection kernels are unchanged by the final input/IO guards.

## Review and remaining scope

Independent source review found and resolved topology, filename and restart
issues; the final IO/common-version exchange pass reports no P0/P1 blocker.
The PR remains a bounded implementation change. It does not establish a
complete discrete adjoint, general rotated/nonorthogonal conductor accuracy,
clinical lead equivalence, anatomical transfer or eikonal reciprocal recording.
Heart-only mode cannot revalidate an absent torso. Fingerprints do not encode
every coupled transform option; AMI and rotating meshes remain outside the
verified set. Regenerate older format-1 fields.

The repository's full CI matrix, including full solids4foam builds and the
remaining tutorials, is distinct from these executed local checks. The papers
and private torso dataset are not part of the code PR.
