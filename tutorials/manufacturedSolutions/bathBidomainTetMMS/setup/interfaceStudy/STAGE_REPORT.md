# Interface-flux verification stage report

> Authoritative final formulation: `../../FINAL_SOLUTION.md`. This report is
> retained as the chronological audit trail and includes rejected controls.

## Gates completed

- Baseline N=10/20/40 evidence and the materialized N=40 mesh were archived
  with SHA-256 hashes.
- The face coefficient was extracted into a testable helper without changing
  the default solver behavior.
- The unequal-cell analytic test passed:
  - equal-distance harmonic formulas agree;
  - distance-weighted harmonic reproduces exact series resistance;
  - unweighted harmonic does not reproduce unequal-distance resistance;
  - physics-correct interface interpolation excludes `Gi`;
  - naive `sigmaTotal` interpolation includes `Gi` and gives a different
    coefficient.
- Runtime selectors build and run:
  - `unweightedHarmonic`;
  - `distanceWeightedHarmonic`;
  - `naiveLinearSigmaTotal`.
- Invalid selector values fail with a diagnostic listing valid values.
- N=10 serial/parallel equivalence passed for all 47 interface metrics.
- Structured N=10/20/40 comparisons completed.
- Parallel tet N=10/20/40 comparisons completed on frozen per-resolution
  meshes for both physics-correct methods.
- The identical one-sided reconstruction was repeated with manufactured exact
  `phiE` samples on the weighted-method N=10/20/40 meshes.

## Structured result

Both harmonic methods are identical on the equal-width structured mesh and
recover second-order heart and bath `phiE` L2 convergence. Reconstructed
interface traces converge at approximately second order. The deliberately
naive `Gi+Ge` control is substantially less accurate and converges only around
first order in the region potential norms.

This validates the region/interface diagnostic and demonstrates why the bath
interface must use heart `Ge`, not heart `Gi+Ge`.

## Tetrahedral result

### Unweighted harmonic

| Metric | N=10 | N=20 | N=40 | rates |
| --- | ---: | ---: | ---: | ---: |
| heart `phiE` L2 | 5.928e-3 | 1.844e-3 | 6.237e-4 | 1.72, 1.59 |
| bath `phiE` L2 | 7.980e-3 | 2.486e-3 | 8.415e-4 | 1.72, 1.59 |
| x=0 trace-jump L2 | 5.264e-3 | 9.364e-4 | 2.068e-4 | 2.55, 2.21 |
| x=1 trace-jump L2 | 5.534e-3 | 9.203e-4 | 2.065e-4 | 2.65, 2.19 |
| x=0 heart-flux L2 | 3.571e-3 | 4.097e-4 | 9.109e-4 | 3.20, -1.17 |
| x=1 heart-flux L2 | 3.818e-3 | 4.039e-4 | 8.193e-4 | 3.32, -1.04 |

### Distance-weighted harmonic

| Metric | N=10 | N=20 | N=40 | rates |
| --- | ---: | ---: | ---: | ---: |
| heart `phiE` L2 | 5.991e-3 | 1.849e-3 | 6.256e-4 | 1.74, 1.59 |
| bath `phiE` L2 | 8.088e-3 | 2.491e-3 | 8.415e-4 | 1.74, 1.59 |
| x=0 trace-jump L2 | 5.369e-3 | 9.645e-4 | 2.134e-4 | 2.54, 2.21 |
| x=1 trace-jump L2 | 5.657e-3 | 9.422e-4 | 2.171e-4 | 2.65, 2.15 |
| x=0 heart-flux L2 | 3.564e-3 | 4.019e-4 | 9.174e-4 | 3.22, -1.21 |
| x=1 heart-flux L2 | 3.815e-3 | 3.952e-4 | 8.821e-4 | 3.35, -1.18 |

Net exterior current remains approximately 2--3e-9 for every tet case.

## Conclusion

Distance weighting is analytically correct for unequal two-point resistance,
but it does not fix the observed tetrahedral convergence behavior. The two
physics-correct methods are nearly indistinguishable on this mesh family, so
the earlier hypothesis that unweighted harmonic interpolation was the dominant
cause is not supported by the N=10/20/40 evidence.

The potential trace mismatch converges cleanly while the numerical one-sided
flux becomes non-monotone at N=40.

## Exact-field reconstruction isolation

| Metric | N=10 | N=20 | N=40 | rates |
| --- | ---: | ---: | ---: | ---: |
| x=0 exact trace-jump L2 | 2.535e-3 | 1.812e-4 | 1.135e-5 | 3.81, 4.00 |
| x=1 exact trace-jump L2 | 2.763e-3 | 1.754e-4 | 1.142e-5 | 3.98, 3.94 |
| x=0 exact heart-flux L2 | 3.532e-3 | 4.950e-4 | 6.228e-5 | 2.84, 2.99 |
| x=1 exact heart-flux L2 | 3.773e-3 | 4.816e-4 | 6.252e-5 | 2.97, 2.95 |

The exact sampled field converges cleanly under the identical mesh-dependent
reconstruction. This rules out the postprocessor stencil as the cause of the
N=40 reversal. The non-monotone localized interface gradient is present in the
numerical `phiE` field even though its region norms and reconstructed traces
continue to converge.

## Bidomain interface identities and assembled-flux isolation

The mathematical audit in `docs/bidomain_foundations.md` requires three
distinct interface checks: extracellular potential continuity, extracellular
current continuity, and intracellular no-leak
`Gi*grad(Vm + phiE).n = 0`. The last condition was missing from the initial
diagnostic and has now been added.

The OpenFOAM Laplacian face flux was also reconstructed with the production
surface tensor and separated into its orthogonal and explicit non-orthogonal
parts. The table reports N=20 to N=40 L2 behavior for the weighted method.

| Metric | N=20 | N=40 | rate |
| --- | ---: | ---: | ---: |
| x=0 independently reconstructed extracellular flux | 4.019e-4 | 9.174e-4 | -1.19 |
| x=0 intracellular leakage | 6.969e-4 | 1.266e-3 | -0.86 |
| x=0 assembled extracellular flux | 5.504e-4 | 2.575e-3 | -2.23 |
| x=0 non-orthogonal correction magnitude | 2.367e-4 | 5.940e-4 | -1.33 |
| x=0 orthogonal-part flux error | 5.163e-4 | 2.450e-3 | -2.25 |
| x=1 independently reconstructed extracellular flux | 3.952e-4 | 8.821e-4 | -1.16 |
| x=1 intracellular leakage | 6.735e-4 | 9.990e-4 | -0.57 |
| x=1 assembled extracellular flux | 5.398e-4 | 2.634e-3 | -2.29 |
| x=1 non-orthogonal correction magnitude | 2.245e-4 | 7.858e-4 | -1.81 |
| x=1 orthogonal-part flux error | 4.993e-4 | 2.527e-3 | -2.34 |

The assembled face flux confirms that the reversal is not an artifact of the
independent reconstruction. Both components worsen, but the orthogonal part
dominates the assembled error. Intracellular no-leak also worsens, showing that
the numerical cancellation between the `Vm` and `phiE` gradients deteriorates
on the N=40 interface layer.

## Same-mesh serial/parallel audit

The original N=10 equivalence harness used an older archived serial CSV and did
not require identical CSV schemas. The harness was corrected to run serial and
parallel solutions consecutively on the same `constant/polyMesh` and to reject
schema differences.

- N=20: all 73 current metrics agree exactly after reconstruction.
- N=40: all 73 current metrics agree exactly after reconstruction.
- The reconstructed interface inventories contain 1888 and 7412 faces,
  respectively, with the expected myocardium/bath zone partition.

This rules out processor decomposition, reconstruction, and diagnostic
owner/neighbour orientation as causes of the N=40 reversal. It does not rule
out a production discretization or boundary-assembly problem shared by serial
and parallel execution.

## Exact discrete operator audit

The exact manufactured fields were inserted into the same global and submesh
operators used by production. Residuals were separated between cells adjacent
to the heart--bath interface and bulk cells.

| Exact-field L2 quantity | N=20 | N=40 | rate |
| --- | ---: | ---: | ---: |
| global/submesh equation residual, heart interface | 4.729e-1 | 1.001 | -1.08 |
| global/submesh equation residual, heart bulk | 8.008e-1 | 1.095 | -0.45 |
| global `phiE` operator error, heart interface | 1.489 | 2.817 | -0.92 |
| global `phiE` operator error, heart bulk | 8.454e-1 | 1.131 | -0.42 |
| mapped `Gi*grad(Vm)` source error, heart interface | 1.603 | 2.958 | -0.88 |
| mapped `Gi*grad(Vm)` source error, heart bulk | 2.792e-1 | 2.836e-1 | -0.02 |

Raw cell-divergence residuals need not converge pointwise on irregular tets,
so these values are diagnostic rather than paper convergence norms. Their
localization is nevertheless clear: the mapped intracellular source remains
stable in the bulk but deteriorates at the artificial myocardium-submesh
boundary. The global operator has a matching interface-localized error. The
current formulation therefore relies on cancellation between boundary-adjacent
corrections evaluated on different meshes and gradient stencils. Regular
structured geometry preserves that cancellation; the tetrahedral sequence does
not preserve it monotonically.

## Gradient-scheme screen

Four schemes were screened at N=20 using the same weighted interface
conductivity. Replacing least-squares gradients with `Gauss linear` was
catastrophically inaccurate, and removing non-orthogonal correction was also
worse. Least-squares plus a correction limited to 0.5 was the only informative
alternative and was carried through N=10/20/40.

| Limited-correction metric | N=10 | N=20 | N=40 | rates |
| --- | ---: | ---: | ---: | ---: |
| heart `phiE` L2 | 4.888e-3 | 8.299e-4 | 1.060e-3 | 2.56, -0.35 |
| bath `phiE` L2 | 7.013e-3 | 1.174e-3 | 8.143e-4 | 2.58, 0.53 |
| x=0 reconstructed extracellular flux L2 | 3.574e-3 | 4.935e-4 | 1.955e-4 | 2.86, 1.34 |
| x=1 reconstructed extracellular flux L2 | 3.817e-3 | 4.845e-4 | 2.073e-4 | 2.98, 1.23 |
| x=0 intracellular leakage L2 | 3.372e-3 | 4.850e-4 | 2.195e-4 | 2.80, 1.14 |
| x=1 intracellular leakage L2 | 3.602e-3 | 4.896e-4 | 2.476e-4 | 2.88, 0.98 |
| x=0 assembled flux L2 | 7.800e-4 | 7.880e-4 | 7.305e-4 | -0.01, 0.11 |
| x=1 assembled flux L2 | 7.738e-4 | 8.087e-4 | 7.179e-4 | -0.06, 0.17 |

Limiting suppresses the localized constitutive-gradient and intracellular-
leakage failures, recovering approximately first-order N=20 to N=40 behavior.
It does not provide an acceptable complete scheme: assembled flux barely
improves and heart potential is non-monotone. Thus the current failure contains
a correction-outlier component, but changing a dictionary scheme alone is not
the final solution.

## Matched-stencil exact-operator experiment

Before changing the solver, the exact-field diagnostic was reformulated as

```text
div(Ge grad(phiE)) = -div(Gi grad(Vm + phiE))
```

using the global conformal mesh for the extracellular operator and one combined
heart-submesh operator for the intracellular current. This directly imposes
zero intracellular flux on the artificial submesh boundary instead of relying
on cancellation between separately corrected `Vm` and `phiE` terms.

| Raw exact residual L2 | N=10 | N=20 | N=40 |
| --- | ---: | ---: | ---: |
| current split, heart interface | 3.372e-1 | 4.733e-1 | 1.001 |
| matched split, heart interface | 1.011e-1 | 1.033e-1 | 1.792e-1 |
| matched split, heart bulk | 1.716e-1 | 2.345e-1 | 3.209e-1 |

Raw pointwise divergence residuals on irregular tets are not convergence norms,
and both matched interface and bulk values grow with refinement. The useful
comparison is localization and magnitude: at N=40 the matched interface
residual is 5.6 times smaller than the current split and is only 0.56 of its
matched bulk value; the current interface residual is 0.91 of current bulk.
This supports, but does not yet validate, matched assembly. The decisive test
requires an implicit solver variant that maps the heart-submesh `Gi*phiE`
matrix into the global extracellular system.

## Experimental implicit matched-submesh convergence

A non-default `matchedSubmesh` prototype now assembles the global
`Ge`/bath matrix, maps the implicit heart-submesh `Gi*phiE` matrix into it, and
uses the same heart-submesh stencil for the `Gi*Vm` right-hand side. The first
tetrahedral convergence ladder gives:

| Metric | N=10 | N=20 | N=40 | rates |
| --- | ---: | ---: | ---: | ---: |
| heart `phiE` L2 | 6.066e-3 | 1.812e-3 | 5.165e-4 | 1.74, 1.81 |
| bath `phiE` L2 | 8.214e-3 | 2.445e-3 | 7.086e-4 | 1.75, 1.79 |
| x=0 heart extracellular-flux L2 | 3.515e-3 | 4.607e-4 | 8.478e-5 | 2.93, 2.44 |
| x=1 heart extracellular-flux L2 | 3.751e-3 | 4.564e-4 | 8.132e-5 | 3.04, 2.49 |
| x=0 intracellular leakage L2 | 3.509e-3 | 5.345e-4 | 1.327e-4 | 2.72, 2.01 |
| x=1 intracellular leakage L2 | 3.752e-3 | 5.272e-4 | 1.355e-4 | 2.83, 1.96 |
| x=0 bath extracellular-flux L2 | 9.450e-5 | 1.020e-4 | 9.458e-5 | -0.11, 0.11 |
| x=1 bath extracellular-flux L2 | 9.015e-5 | 9.935e-5 | 9.720e-5 | -0.14, 0.03 |

The matched assembly removes the previous N=40 reversal in heart/bath
potential, heart-side extracellular flux, and intracellular leakage. This is
the first solver-level evidence that the global/submesh stencil mismatch was a
dominant cause. It is not yet a completed method: bath-side reconstructed flux
remains near a `1e-4` floor, the matched assembled face flux is not yet exposed
by the diagnostic, structured regression is pending, and coupled-patch matrix
mapping was initially pending.

### Structured and composite-flux follow-up

The matched prototype passes the structured N=10/20/40 regression: heart and
bath `phiE` retain second-order L2 convergence and the matched/current results
are equivalent on orthogonal geometry.

The corrected matched assembled interface flux uses the global `Ge`/bath face
matrix; the mapped intracellular submesh flux is zero at the heart boundary.

| Metric | N=10 | N=20 | N=40 | rates |
| --- | ---: | ---: | ---: | ---: |
| x=0 matched assembled flux | 4.680e-4 | 2.942e-4 | 3.444e-4 | 0.67, -0.23 |
| x=1 matched assembled flux | 4.612e-4 | 3.027e-4 | 3.510e-4 | 0.61, -0.21 |

This is about seven times more accurate than the current split at N=40, but it
still has a mild N=20 to N=40 reversal. The matched formulation fixes the
dominant intracellular inconsistency but not all global extracellular
face-gradient error, consistent with the remaining bath-side flux floor.

Combining `matchedSubmesh` with a fixed 0.5 correction limiter was rejected at
N=10/N=20: its assembled-flux error increases from `7.43e-4` to `7.90e-4` and
is worse than matched plus fully corrected. No N=40 run is justified.

### Matched parallel equivalence

The mapped matrix now includes processor-coupled heart-submesh boundary
coefficients in the corresponding global processor patches. Same-mesh serial
and six-rank parallel comparisons pass all 109 current metrics exactly at
N=10, N=20, and N=40. The matched formulation is therefore parallel-safe over
the tested decomposition and mesh ladder.

Before execution, N=80 was conservatively estimated to exceed local memory
based on cell-count scaling. Runtime monitoring superseded that estimate, as
recorded below.

### N=80 confirmation

N=80 completed locally on six ranks in about 91 minutes. The mesh contained
6,831,714 cells and 29,528 internal heart--bath faces. Peak observed aggregate
rank memory remained well below the conservative estimate. Final reconstruction
and numerical/exact diagnostics completed successfully.

| Metric | N=10 | N=20 | N=40 | N=80 | fitted order |
| --- | ---: | ---: | ---: | ---: | ---: |
| heart `phiE` L2 | 6.066e-3 | 1.812e-3 | 5.165e-4 | 1.630e-4 | 1.75 |
| bath `phiE` L2 | 8.214e-3 | 2.445e-3 | 7.086e-4 | 2.241e-4 | 1.74 |
| x=0 reconstructed heart flux L2 | 3.515e-3 | 4.607e-4 | 8.478e-5 | 1.388e-4 | 1.64 |
| x=1 reconstructed heart flux L2 | 3.751e-3 | 4.564e-4 | 8.132e-5 | 1.348e-4 | 1.69 |
| x=0 reconstructed flux jump L2 | 3.558e-3 | 5.493e-4 | 1.191e-4 | 9.040e-5 | 1.81 |
| x=1 reconstructed flux jump L2 | 3.796e-3 | 5.431e-4 | 1.209e-4 | 9.137e-5 | 1.83 |
| x=0 intracellular leakage L2 | 3.509e-3 | 5.345e-4 | 1.327e-4 | 5.368e-5 | 2.01 |
| x=1 intracellular leakage L2 | 3.752e-3 | 5.272e-4 | 1.355e-4 | 5.613e-5 | 2.01 |
| x=0 matched assembled flux L2 | 4.680e-4 | 2.942e-4 | 3.444e-4 | 2.875e-4 | 0.19 |
| x=1 matched assembled flux L2 | 4.612e-4 | 3.027e-4 | 3.510e-4 | 2.828e-4 | 0.19 |

The assembled flux resumes improvement at N=80 and is approximately 2.8--2.9%
of the exact current magnitude. Its four-level fitted order is only about 0.19,
so it should be described as bounded and slowly convergent, not high-order.
The independently reconstructed heart flux is non-monotone from N=40 to N=80,
although its four-level fitted rate remains 1.6--1.7. Exact-field reconstruction
converges at about third order over that interval, proving the N=80 reversal is
in the numerical local gradient rather than the measurement stencil.

Potential, reconstructed flux jump, and intracellular no-leak results are
strong and paper-usable. The manuscript must disclose the local heart-gradient
non-monotonicity and the low assembled-flux order separately from exact local
and global conservation.

## Stop decision

The Stage 7 acceptance gate fails because one-sided interface-flux errors are
not monotone from N=20 to N=40. N=80 remains unauthorized.

The exact-field isolation and assembled-flux diagnostics identify a numerical,
localized interface-gradient problem. Harmonic distance weighting is not its
dominant cause, and the non-orthogonal correction is not its only component.
The next smallest diagnostic is to localize the N=40 outlier faces against
interface non-orthogonality, skewness, interpolation weight, and adjacent-cell
quality. In parallel, a non-default matched-stencil operator variant should
assemble the intracellular `Vm + phiE` contribution with one interface-aware
discrete flux instead of depending on cancellation of separately corrected
global/submesh operators. That variant must pass the structured case before a
tet comparison or production-default change. N=80 remains blocked.
