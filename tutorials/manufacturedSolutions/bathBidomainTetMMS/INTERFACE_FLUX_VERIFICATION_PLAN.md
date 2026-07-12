# Bath-Bidomain Interface-Flux Verification Plan

> Final formulation and paper-ready results: see `FINAL_SOLUTION.md`. This file
> is retained as the original gated execution plan.

## Execution status (2026-07-12)

Stages 0--6 are complete. The parallel tetrahedral N=10/20/40 comparison in
Stage 7 completed, but the one-sided interface-flux errors are non-monotone
between N=20 and N=40 for both physics-correct interpolation methods. Stage 7
therefore fails its acceptance gate. The exact-field reconstruction diagnostic
then passed at approximately third order for flux and fourth order for trace
continuity, ruling out the measurement stencil and identifying a localized
numerical-gradient issue. Stage 8 and N=80 remain paused pending inspection of
the assembled interface flux/non-orthogonal correction on N=20/N=40. That
inspection now confirms that the assembled flux, its orthogonal part, its
non-orthogonal correction, and the independently reconstructed intracellular
no-leak condition all worsen at N=40. The next gate is localization of those
outliers against interface-face mesh quality before temporal separation. See
`setup/interfaceStudy/STAGE_REPORT.md`.

The corrected same-mesh serial/parallel audit now passes all 73 current metrics
exactly at both N=20 and N=40. Parallel decomposition and reconstruction are
therefore excluded as explanations for the N=40 reversal, but a serial-and-
parallel production discretization issue remains possible.

The subsequent exact discrete-operator audit localizes the strongest split
error to the artificial myocardium-submesh boundary: the mapped intracellular
source bulk error is stable from N=20 to N=40, while its interface-layer error
nearly doubles together with the global operator interface error. The next
implementation is an experimental matched-stencil intracellular-flux variant;
the default must remain unchanged until structured and tet gates pass.

A controlled gradient screen confirms that least-squares gradients and
non-orthogonal correction are necessary. Limiting correction to 0.5 restores
roughly first-order reconstructed extracellular-flux and intracellular-no-leak
convergence from N=20 to N=40, but makes heart `phiE` non-monotone and leaves
assembled-flux convergence near zero order. Limiting is therefore evidence of
localized correction outliers, not an acceptable final production fix.

An exact-operator matched-stencil experiment reduces the N=40 heart-interface
residual from 1.001 to 0.179 and removes most interface localization relative
to bulk. This is supporting evidence only. The next implementation must add the
heart-submesh `Gi*phiE` matrix implicitly to the global `Ge/bath` system and
pair it with the same-stencil `Gi*Vm` right-hand side behind a non-default
selector.

The first serial tetrahedral run of that implicit prototype is promising:
heart/bath `phiE` converge at about 1.8 order, heart extracellular flux at
about 2.4--2.5 order from N=20 to N=40, and intracellular leakage at about
second order. The old N=40 reversal disappears. Bath-side reconstructed flux
remains near `1e-4`, so structured regression, matched assembled-flux output,
and bath-flux diagnosis are required before parallel implementation or N=80.

Structured regression now passes. The corrected composite assembled-flux
metric is roughly seven times better than the current split at N=40, but still
reverses mildly from N=20 to N=40. A matched-plus-0.5-limiter screen is worse
in assembled flux and is rejected. The remaining defect is specifically in the
global extracellular heart--bath flux; a conservative two-domain formulation
is a defensible alternative experiment only if it shares one interface flux or
is assembled monolithically.

Parallel processor-patch matrix mapping is now implemented. Same-mesh serial
and six-rank parallel runs agree exactly across all 109 metrics at N=10, N=20,
and N=40. N=80 is implementation-ready but requires a distributed-memory host;
the local 24 GiB machine is below the estimated aggregate memory requirement.

N=80 subsequently completed locally on six ranks. Potentials retain fitted
orders near 1.74, flux jump near 1.82, and intracellular leakage near 2.01 over
N=10/20/40/80. The matched assembled flux improves from N=40 to N=80 but has a
low fitted order near 0.19 and about 2.8--2.9% relative RMS error. Independent
heart-side reconstructed flux is non-monotone over the final interval despite
a four-level fitted order near 1.6--1.7. These limitations must be explicit in
the paper claim.

## Goal

Demonstrate the finite-volume property that motivates the bath-bidomain
tetrahedral case:

> `cardiacFoam` transfers extracellular potential and current accurately and
> conservatively across a conformal myocardium--bath material interface on
> non-orthogonal tetrahedral meshes.

The primary study is interface-current verification. Global `Vm`, `phiE`, and
`phiI` norms remain system-level checks but are not sufficient evidence by
themselves.

N=80 is excluded until the selected interface treatment has monotone errors
through N=40 and passes serial/parallel equivalence.

## Current Baseline

- Completed tetrahedral meshes: nominal N=10, 20, 40.
- Standard `checkMesh`: pass at all three resolutions.
- Existing implementation: componentwise, unweighted harmonic face tensor on
  all internal and coupled faces.
- At a heart--bath face, the implementation correctly combines heart
  extracellular conductivity `Ge` with bath conductivity, excluding `Gi`.
- Existing manufactured `phiE` norm combines heart and bath cells.
- No direct one-sided interface-flux, boundary-current, or global-current
  balance metric exists.
- N=80 has not been run.

## Evidence Required

### Primary metrics

1. Volume-weighted L1, L2, and Linf `phiE` error in `myocardium`.
2. Volume-weighted L1, L2, and Linf `phiE` error in `bath`.
3. Potential error on both interface surfaces, x=0 and x=1.
4. Independently reconstructed heart-side extracellular normal flux error.
5. Independently reconstructed bath-side normal flux error.
6. Heart--bath flux jump using outward normals.
7. Integrated interface-current error.
8. Global current-balance residual.

The manufactured interface-current magnitude is `alpha = 0.01`.

### Boundary metrics

- `xMin`: potential error against the exact grounded value zero.
- `xMax`: normal-current error against `alpha = 0.01`.
- `sides`: normal-current leakage against zero.

### Secondary metrics

- Myocardium `Vm` and `phiI` norms.
- Maximum and average non-orthogonality.
- Maximum skewness and minimum volume.
- Low-determinant-cell count from the strict mesh audit.
- `phiE` and `Vm` linear iterations.
- Runtime and processor count.

## Comparison Variants

Every comparison must use identical meshes, time controls, solver tolerances,
and manufactured parameters.

### A. `unweightedHarmonic`

Current implementation:

```text
sigma_f = 2 sigma_P sigma_N / (sigma_P + sigma_N)
```

At a heart--bath face use `Ge` on the heart side and bath conductivity on the
bath side.

### B. `distanceWeightedHarmonic`

Candidate production implementation:

```text
sigma_f = (d_P + d_N)/(d_P/sigma_P + d_N/sigma_N)
```

Use normal owner-to-face and neighbour-to-face distances. Preserve the
heart-side `Ge` versus bath-conductivity physics.

Within a uniform material, retain the normal OpenFOAM interpolation behavior
unless a coefficient-level test justifies another choice.

### C. `naiveLinearSigmaTotal`

Control only: ordinary OpenFOAM linear interpolation of heart `Gi+Ge` against
bath conductivity. This is expected to contaminate the interface with
intracellular conductivity and must not be presented as a valid production
method.

## Stage 0 -- Protect the Existing Evidence

1. Record the Git SHA and OpenFOAM/Gmsh versions.
2. Preserve the completed N=10/20/40 result directories.
3. Hash the three meshes or archive their `polyMesh` directories so every
   interpolation variant uses exactly the same cells and faces.
4. Do not overwrite the structured bath-bidomain reference results.

Acceptance:

- Existing `setup/results/summary.csv` remains reproducible.
- Mesh identity can be checked before every comparison run.

## Stage 1 -- Unequal-Cell Coefficient Test

Create a small automated test for a one-dimensional two-material slab with
unequal owner/neighbor distances.

Analytic reference:

```text
R = d_heart/sigma_e + d_bath/sigma_b
q = deltaPhi/R
```

Required cases:

1. Equal distances: current and distance-weighted harmonic must agree.
2. Unequal distances: distance-weighted harmonic must reproduce the analytic
   series resistance.
3. Naive `sigmaTotal` linear interpolation must demonstrate the effect of
   including `Gi` at the interface.

Suggested files:

```text
src/electroModels/electroDomains/extracellularPotentialDomain/
    extracellularFaceConductivity.H
    tests/test_extracellular_face_conductivity.cpp
```

Keep the coefficient calculation in a small independently testable helper;
do not require a full `cardiacFoam` run for this stage.

Acceptance:

- Analytic unequal-cell flux is reproduced to roundoff by variant B.
- Variant A is shown to agree only in the equal-distance case.
- The test fails if `Gi` enters the physics-correct interface coefficient.

## Stage 2 -- Runtime-Selectable Interface Treatment

Add one dictionary key under `bathPotentialDomain`, for example:

```foam
interfaceConductivityInterpolation distanceWeightedHarmonic;
```

Accepted values:

```text
unweightedHarmonic
distanceWeightedHarmonic
naiveLinearSigmaTotal
```

Implementation location:

```text
src/electroModels/electroDomains/extracellularPotentialDomain/
    extracellularPotentialDomain.H
    extracellularPotentialDomain.C
```

Rules:

- Preserve the current method as an explicit selectable baseline.
- Default behavior must remain unchanged until verification supports changing
  it.
- Print the selected method when `reportSetup yes`.
- Reject unknown selector values with a fatal error listing valid choices.
- Keep processor-coupled face behavior consistent with internal faces.

Acceptance:

- The library builds in the sourced OpenFOAM v2412 environment.
- Existing structured regression passes under the default.
- Each selector is exercised by a focused test.

## Stage 3 -- Region, Interface, and Boundary Diagnostics

Extend the manufactured verification output or add a dedicated reconstructed-
field postprocessor.

Preferred first implementation:

1. Run simulations serially or in parallel.
2. Reconstruct the final parallel field when required.
3. Calculate face diagnostics serially on the reconstructed mesh to avoid
   processor-patch double counting.

Suggested output:

```text
postProcessing/bathBidomainInterfaceMetrics.csv
```

One row per case with columns for:

```text
method, N_nominal, h_heart, h_bath,
L1_phiE_heart, L2_phiE_heart, Linf_phiE_heart,
L1_phiE_bath, L2_phiE_bath, Linf_phiE_bath,
x0_potential_L2, x1_potential_L2,
x0_fluxHeart_L2, x0_fluxBath_L2, x0_fluxJump_L2,
x1_fluxHeart_L2, x1_fluxBath_L2, x1_fluxJump_L2,
xMin_potential_L2, xMax_flux_L2, sides_flux_L2,
interfaceCurrentIntegralError, globalCurrentBalance
```

Norm rules:

- Cell L1/L2 norms: volume weighted.
- Face L1/L2 norms: area weighted.
- Report Linf separately.
- Treat x=0 and x=1 separately before aggregating them.
- Use independently reconstructed one-sided gradients for constitutive-flux
  accuracy; do not use only the single assembled conservative face flux.

Acceptance:

- A structured N=10 case reproduces its known whole-field errors.
- Region cell counts sum to the global count.
- Interface face areas and exterior patch areas are reported.
- Exact manufactured current `alpha=0.01` is included in the output metadata.
- Re-running the diagnostic is deterministic.

## Stage 4 -- Structured Reference Comparison

Run the existing structured bath-bidomain case at N=10, 20, and 40 for all
three interpolation variants.

Purpose:

- isolate material-interface interpolation from non-orthogonality;
- validate the new metrics;
- establish whether variant B preserves the structured second-order result.

Acceptance for candidate production method B:

- Monotone heart and bath `phiE` L2 errors.
- Approximately second-order region L2 convergence on the structured ladder.
- Decreasing one-sided interface-flux errors.
- Flux jump and global balance consistent with solver/discretization tolerance.
- Correct `xMin`, `xMax`, and side-boundary behavior.

Stop if the structured reference fails. Do not proceed to tetrahedral scheme
comparisons until the diagnostic or implementation is corrected.

## Stage 5 -- Tetrahedral N=10 Diagnostic Matrix

Run variants A, B, and C on the existing N=10 tetrahedral mesh.

Example result layout:

```text
setup/interfaceStudy/
    unweightedHarmonic/N10/
    distanceWeightedHarmonic/N10/
    naiveLinearSigmaTotal/N10/
```

For each run archive:

- exact command;
- dictionaries;
- mesh hash;
- solver log;
- final reconstructed fields;
- interface/boundary metrics;
- runtime manifest.

Acceptance:

- All runs use the identical mesh.
- Variant B improves or matches interface-flux accuracy relative to A.
- Variant C demonstrates why physics-aware interface treatment is necessary.
- No method is selected from global `phiE` error alone.

## Stage 6 -- Serial/Parallel Equivalence

Use N=10 and the current structured bath-bidomain decomposition pattern:

```text
decomposePar
mpirun -np 6 cardiacFoam -parallel
reconstructPar
```

Compare serial and parallel:

- heart/bath potential norms;
- interface flux norms;
- boundary metrics;
- global current balance;
- final `Vm` and `phiE` fields.

Acceptance:

- Diagnostic differences are within a declared tight floating-point tolerance.
- Manufactured summary is written once to the shared `postProcessing` path.
- No interface face is omitted or double counted after reconstruction.

## Stage 7 -- Parallel Tetrahedral N=10/20/40 Sweep

Run only the best-supported physics-correct method, normally variant B if the
preceding gates pass.

Default execution:

```text
N=10: 2--4 ranks
N=20: 4--6 ranks
N=40: 6--12 ranks
```

Choose rank counts from available local resources; record them in the result
manifest. Avoid reconstructing every time directory. Reconstruct only the final
fields required by the interface diagnostic.

Acceptance:

- Heart and bath `phiE` L2 errors decrease monotonically.
- Both one-sided interface-current errors decrease monotonically.
- Linf outliers are localized and explained.
- The N=20 to N=40 trend is asymptotic enough to justify another refinement.
- Conservation and accuracy are reported as separate quantities.

## Stage 8 -- Temporal Separation Check

At N=20 or N=40, repeat the selected method with `deltaT/2` while holding the
mesh fixed.

Acceptance:

- Changes in the primary spatial/interface metrics are materially smaller than
  the N-to-2N spatial change.
- If they are not, reduce `deltaT` for the spatial ladder and rerun affected
  cases.

## Stage 9 -- N=80 Authorization Gate

N=80 may run in parallel only when all conditions hold:

1. Structured interface verification passes.
2. Unequal-cell coefficient test passes.
3. Serial/parallel equivalence passes.
4. N=10/20/40 heart and bath L2 errors are monotone.
5. Interface one-sided flux errors are monotone.
6. Global current balance is acceptable and understood.
7. Temporal separation is demonstrated.
8. The expected runtime, memory, rank count, and disk use are recorded.

N=80 is a confirmation point, not a method-selection point.

## Stage 10 -- Paper Products

Generate repository-relative artifacts:

```text
interface_method_comparison.csv
structured_interface_convergence.csv
tet_interface_convergence.csv
boundary_flux_summary.csv
serial_parallel_equivalence.csv
solver_scaling.csv
```

Recommended paper figures:

1. Schematic of the conformal heart--bath face and one-sided fluxes.
2. Heart and bath `phiE` L2 convergence, plotted separately.
3. Interface-current error and flux-jump convergence.
4. Comparison of interpolation variants on identical tetrahedral meshes.
5. Optional spatial map locating N=40 Linf outliers.

The manuscript must distinguish:

- local conservation from constitutive-flux accuracy;
- smooth-domain non-orthogonal accuracy from discontinuous-interface accuracy;
- the general OpenFOAM interpolation description from the custom
  bath-interface treatment;
- verification from physical validation.

## Recommended Execution Order

```text
1. Stage 0: preserve meshes/results
2. Stage 1: unequal-cell coefficient test
3. Stage 2: runtime selector
4. Stage 3: diagnostics
5. Stage 4: structured comparison
6. Stage 5: tet N=10 method comparison
7. Stage 6: serial/parallel equivalence
8. Stage 7: parallel N=10/20/40
9. Stage 8: temporal separation
10. Stage 9: decide on N=80
11. Stage 10: paper tables and figures
```

Do not start a later stage when an earlier acceptance gate fails. Record the
failure and diagnose it using the smallest preceding case.
