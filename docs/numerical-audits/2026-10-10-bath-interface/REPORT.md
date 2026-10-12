# Bath-bidomain interface diagnostic audit — 10 October 2026

The tests identify two independent reconstruction inconsistencies: the stock global gradient crossing a conductivity jump contaminates the scalar non-orthogonal correction, and the heart-submesh boundary gradient treats physical face values with a normal-projected displacement. A scratch change to the first path makes the sheared affine solution exact. The second inconsistency exists even with exact cell values and exact face traces. These findings do not establish that every remaining nonlinear error comes from the boundary; an exact-trace oracle alone did not change the coupled solution errors.

## Scope and reproducibility

Source: `feat/bath-transmission-mms`, commit `b9ec31966`, **including the uncommitted transmission MMS**. See `manifest.json` for source hashes. Native OpenFOAM v2412, double precision, 32-bit labels. An isolated lightweight build used copied installed ionic/generic-writer/active-tension libraries. Solver changes were made only in `/private/tmp/bath-mms-audit-20261010`, not in the repository's solver sources or installed solver.

Ran 32 small solver cases: 30 transmission cases and two heart-only cases, plus the generator tests, affine interface utility, leakage utility, and exact-field gradient probes. Raw error summaries and probe logs are in `evidence/`; machine-readable transmission results are in `results.json`. Scratch controls and diagnostic source are archived here.

The elliptic cases use one step of 1e-9 s, with 10 outer correctors and 10 non-orthogonal correctors unless indicated. Vm is effectively frozen at its exact initial value; its tiny one-step error is **not** evidence of spatial accuracy of the Vm advance. Nonlinear convergence orders below are short ladders, not a proof of the asymptotic order.

## 1. Verify the manufactured construction independently

`test_transmission_mms.py` passes with SymPy 1.14.0, mpmath and numpy. Checks cover interface identities, finite-difference agreement of sources, boundedness, periodicity, tangential richness, common-mode invariance, the affine solution, ionic closure, generated-header agreement, and refusal to write a header after a failed check. The logged `forced failure` is the expected negative test, not an MMS failure. The compiled verifier also reports 819 reference samples passing.

This verifies the construction and generated reference at the sampled conditions; it does not verify every FVM operator or quadrature rule.

## 2. Affine production-path test

Potential and current errors are L2, in the case's units.

| Configuration | N | phiE heart | Interface flux | Stock tangential-gradient diagnostic |
|---|---:|---:|---:|---:|
| Orthogonal, leastSquares | 10 | 1.6884e-13 | 1.4225e-12 | 1.6450e-12 |
| Sheared, leastSquares | 10 | 4.1998e-5 | 2.2163e-6 | 4.0022e-3 |
| Sheared, leastSquares | 20 | 2.0583e-5 | 1.1148e-6 | 4.0022e-3 |
| Sheared, leastSquares | 40 | 1.0186e-5 | 5.5908e-7 | 4.0021e-3 |
| Sheared, Gauss linear | 10 | 1.9406e-13 | 6.2694e-13 | 4.4763e-12 |
| Sheared, consistent correction control | 10 | 1.6831e-13 | 7.3682e-13 | 3.9135e-3 |
| Sheared, consistent correction control | 20 | 1.3613e-13 | 1.6171e-12 | 3.9135e-3 |

The current production assembly already uses `interfaceGradient` for the tensor cross-flux. However, `fvm::laplacian(gammaNN, phiE)` with the corrected scheme still obtains its scalar snGrad correction from the stock gradient. The isolated control changes this to an uncorrected implicit normal operator and adds the scalar non-orthogonal flux explicitly using the same jump-aware gradient as the tensor cross-flux:

```text
Fexplicit = SfGammaCorr . interpolate(gJump)
          + gammaNN magSf nonOrthCorrectionVectors . interpolate(gJump)
```

The control also includes both terms in the stored flux diagnostic. It removes the affine failure without changing the mesh, solution, conductivity, time step, or MMS forcing. This is strong causal evidence against the stock scalar correction path.

**Diagnostic limitation:** `gradTangentialError` in the verifier evaluates `fvc::grad(phiE)`, not the gradient actually used in `interfaceGradient`. It remains nonzero in the exact corrected solution. Its value alone cannot establish that the assembled operator has a gradient error. Report stock, reconstructed, and flux-used gradients separately.

## 3. Boundary geometry with exact inputs

`bathGradientProbe` assigns exact cell values and exact fixedValue face traces. It compares the stock gradient against a least-squares fit using actual centre-to-face displacements, and a control fit using normal-projected boundary displacements. These fits use area/distance-squared weights; they are diagnostic reconstructions, not an exact reimplementation of OpenFOAM's leastSquares weights.

For a smooth, deliberately nonperiodic affine heart scalar with a tangential gradient and a constant-tensor conormal-null direction, cells within two topological rings of the cyclic seams are excluded. There is no time step or potential solve in this measurement.

| Sheared heart, exact affine inputs | N=10 wall L2 | N=20 wall L2 | N=40 wall L2 |
|---|---:|---:|---:|
| Stock leastSquares | 0.457530 | 0.457530 | 0.457530 |
| Fit using true Cf-C | 1.66e-15 | 3.12e-15 | 6.63e-15 |
| Fit using projected displacement | 0.366547 | 0.366547 | 0.366547 |

For the nonlinear MMS heart phiE with exact inputs:

| Sheared heart phiE gradient | N=10 wall L2 | N=20 wall L2 | N=40 wall L2 |
|---|---:|---:|---:|
| Stock, exact physical face traces | 1.65840 | 1.55291 | 1.54540 |
| Fit using true Cf-C | 0.678826 | 0.234074 | 0.100161 |
| Fit using projected displacement | 1.45535 | 1.29906 | 1.26372 |

The physical-face-value/normal-displacement mismatch causes an error independent of refinement in the affine test. For an affine field, writing d = d_n n + d_t, treating the face difference as a purely normal sample introduces `grad_t . d_t / d_n`. With fixed mesh shear, this ratio does not vanish with h.

The production heart phiE trace is fixedValue, so this geometry issue is relevant. **Do not blindly replace all fixedGradient geometry with Cf-C:** the conormal Vm patch stores a normal-gradient extrapolation and uses phiE.snGrad as an offset. Its face-value convention, normal derivative, skew correction and gradient reconstruction must be made consistent together.

## 4. Nonlinear orthogonal versus sheared measurements

All are the original solver path with converged correction counts, G1, elliptic, ATBI.

| Mesh | N | phiE heart L2 | phiE bath L2 | Interface flux L2 |
|---|---:|---:|---:|---:|
| Orthogonal | 10 | 0.026581 | 0.022188 | 0.029718 |
| Orthogonal | 20 | 0.0065419 | 0.0052083 | 0.0080417 |
| Orthogonal | 40 | 0.0016440 | 0.0012745 | 0.0021810 |
| Sheared | 10 | 0.044229 | 0.020750 | 0.093555 |
| Sheared | 20 | 0.017898 | 0.0083548 | 0.044796 |
| Sheared | 40 | 0.0064043 | 0.0035712 | 0.018209 |

Observed heart-potential orders: orthogonal 2.023, 1.992; sheared 1.305, 1.483. Interface-flux orders: orthogonal 1.886, 1.882; sheared 1.062, 1.299. More levels are needed to claim a limiting order.

The consistent scalar correction improves the N=20 sheared heart potential to 0.015095, bath potential to 0.0054945, and interface flux to 0.037538. It does not restore general second-order behavior by itself. Switching the whole case to Gauss linear also does not solve G1: at N=10 its heart-potential error is 0.051469, worse than leastSquares.

## 5. Controls against alternative explanations

- **Correction convergence:** 10 versus 20 outer/non-orthogonal counts changes the reported norms by at most 1e-12 at N=10. Counts of one understate the nonlinear error; more iterations do not remove the converged deficit.
- **Parallel mapping:** serial versus four ranks, including ranks without heart cells and interfaces crossing processor boundaries, agrees in all reported norms to 2e-13 at N=10.
- **Short-time temporal control:** T=0.01, N=10, dt=0.001 versus 0.0005. Vm L2 is 9.2753e-4 versus 9.2768e-4; phiE heart L2 is about 0.044545 in both. At this resolution and duration, halving dt does not remove the spatial/interface error. This is not a complete temporal-order study.
- **Exact trace oracle:** replacing the reconstructed heart phiE face values with the analytical trace leaves the elliptic potential errors and the short-time Vm/potential errors unchanged to the precision relevant here. It changes the gradient diagnostics. This rules out inaccurate phiE trace values as a sufficient explanation in these tests; it does not repair the boundary gradient/snGrad geometry.
- **Heart-only bidomain control:** N=10 and N=20, one tiny step, gauge-removed phiE L2 is 0.00417655 and 0.00106007, order about 1.98. Raw potential norms contain a gauge offset and must not be used for this comparison.
- **Older componentwise interpolation:** isolated tangential affine patch test has flux error 0.00227907 at both N=10 and N=20; the strip residual doubles with refinement. The normal-only control passes. Switching this utility to the current conormalHarmonic helper reduces tangential flux error to below 2e-16 and strip residual to below 4e-15. This older defect is already addressed in the current helper.
- **Leakage utility:** on the N=10 short-time fields, zeroGradient reconstruction with raw conductivity gives net summed current -2.2568e-4; applying the sealed face conductivity gives -2.5974e-14. This is a diagnostic comparison of operators on the stored fields, not evidence that the sealed production solver leaks current. Zero net current also does not establish local gradient accuracy.

## Interpretation and next bounded change

The heart-only control has smooth conductivities and orthogonal wall geometry. The bath adds an internal derivative jump and a submesh physical-face trace on skewed cells. Its failure is therefore compatible with the working heart-only bidomain result.

1. Make the scalar non-orthogonal correction, tensor cross-flux, and stored interface flux use one consistent jump-aware reconstruction. The scratch control demonstrates the affine acceptance criterion.
2. Repair the heart boundary reconstruction and normal-derivative conventions together. Use physical Cf-C for actual fixedValue samples; account for tangential skew displacement in the conormal/offset derivative path. Keep the combined intracellular boundary flux sealed.
3. Expose stock and operator-used gradient diagnostics separately. Add the sheared affine exact-field and operator tests as automatic gates; a single-resolution regression can otherwise preserve the inconsistency.
4. Recheck G1 at N=20,40,80 with correction convergence and a short time-step control after each bounded change. Require improvement in interface-flux and interface-band norms as well as volume potentials.

There is no fresh evidence here of a failure of the outer Dirichlet boundary law, the general linear solver, or the parallel cell/face mapping. The identified boundary problem is the **internal heart-bath interface reconstruction**. The nonlinear deficit has not been fully apportioned between all operator paths, and no production fix has been applied.

## Re-running the preserved scratch tests

```bash
cd /private/tmp/bath-mms-audit-20261010
source env.sh
python3 run_cases.py
python3 run_extra.py
BATH_AUDIT_CORRECTION=1 python3 run_correction.py
python3 run_trace.py
python3 run_time_controls.py
bathGradientProbe -case cases/g1_shear_N20
bathGradientProbe -case cases/g0_shear_N20 -affine -linearWall
```

The Allrun scripts refuse to overwrite existing scratch cases with log.blockMesh. Recreate a fresh scratch case or clean only that scratch case before repeating a solve. `BATH_AUDIT_EXACT_TRACE` is an analytical oracle for these default-parameter test cases, not a deployable algorithm. `scratch-diagnostic-controls.patch` must not be applied as a production fix.
