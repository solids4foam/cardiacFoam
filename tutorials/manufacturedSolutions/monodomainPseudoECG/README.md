# manufacturedSolutions/monodomainPseudoECG tutorial

Manufactured-solution verification for the monodomain stack with pseudo-ECG verification.

## Overview

### Stack

- myocardium solver: `monodomainSolver`
- ionic model: `monodomainFDAManufactured`
- field verification: `manufacturedFDAMonodomainVerifier` from `libverificationModels`
- shared analytical oracle: `verificationModels`
- optional manufactured pseudo-ECG verification: `manufacturedPseudoECGVerifier` from `libverificationModels`

### Purpose

This case verifies:

- field convergence against an analytical oracle
- manufactured ionic export variables (`u1`, `u2`, `u3`)
- manufactured pseudo-ECG reference output

### Key Configuration

For this workflow, the ionic model exposes manufactured verification metadata.

- ionic-model-side behavior: `monodomainFDAManufactured`
- analytical oracle: `verificationModels`
- field verification hook: `modelPrePostProcessors`
- field verifier: `verificationModels/monodomainVerification`
- ECG verifier: `verificationModels/ecgVerification`

Typical outputs include:

- manufactured field summaries in `postProcessing/` (`<dim>_<N>_cells.dat`)
- `postProcessing/pseudoECG.dat`
- `postProcessing/manufacturedPseudoECG_ECG.dat`
- `postProcessing/manufacturedPseudoECGSummary_ECG.dat`

Corrected 2026-09-27 (tutorials-are-pointers 5.4b-P): the last two names used
to omit the `_ECG` suffix (the `ecgDomains` sub-dict name, `ECG` here) --
confirmed against a real `blockMesh`/`cardiacFoam` run's actual output
filenames, not assumed.

## Variants & Extensions

### Tetrahedral (unstructured) Mesh Variant

`setup/studies/tetConvergence/` is an activatable overlay of this same case on a genuinely unstructured mesh: identical `constant/` and `system/` dicts (electroProperties, physicsProperties, controlDict, decomposeParDict), except the mesh generator changes and `setup/studies/tetConvergence/fvSchemes` (with `ddtSchemes.default none`/`ddt(Vm) backward` spelled out explicitly) is swapped in for the duration of a tet run and restored on exit.

Corrected 2026-09-26 (plan §5g Q10): this study used to also ship a local `fvSolution` copy, described as "a tighter `nOuterCorrectors`/`nNonOrthogonalCorrectors` pair", but it was byte-identical to the case's own `system/fvSolution` (`cmp`): both already set `nOuterCorrectors 1`/`nNonOrthogonalCorrectors 1`. It carried no override and is deleted; the tet route uses `system/fvSolution` directly. `fvSchemes` remains, since it differs from `system/fvSchemes` (only in formatting/header, not in any scheme value; both already set `gradSchemes.default leastSquares`).

Corrected 2026-09-27 (tutorials-are-pointers 5.4b-P, `manufacturedMonodomainPseudoECG`'s tutorial record): every entry `fvSchemes` sets (`ddtSchemes.default`/`ddt(Vm)`, `gradSchemes.default`, `divSchemes.default`, `laplacianSchemes.default`, `interpolationSchemes.default`, `snGradSchemes.default`) already equals `system/fvSchemes`'s own value -- confirmed by a real tet run using `system/fvSchemes` directly, with no swap. So this file, while not byte-identical to `system/fvSchemes` (Q10 only deletes byte-identical overlays), is never actually installed by the omniD record's tet route: it stays here as a native quirk, unused.

#### Tetrahedral Variant Purpose

Verifies that OpenFOAM's non-orthogonal `Gauss linear corrected` Laplacian scheme sustains its spatial convergence rate on a genuinely unstructured tetrahedral mesh of the unit cube:

- a smoothly perturbed hex mesh folds before reaching high non-orthogonality (~8 deg), whereas tets naturally reach tens of degrees
- cardiac anatomical meshes are predominantly tetrahedral, the topology the framework must actually handle

#### Grid Generation

`setup/studies/tetConvergence/box.geo.template` is a gmsh (OpenCASCADE) unit cube. Resolution is set with `gmsh -3 box.geo.template -setnumber lc <value>` (`lc = 1/N`), meshes with gmsh (legacy msh2 format), and imports via `gmshToFoam`. All six boundary faces lie on axis-aligned planes `x,y,z in {0,1}`, where the manufactured cosine field has zero normal derivative, so the solver's default zeroGradient boundary stays compatible with the exact solution.

Corrected 2026-09-27 (plan §5g Q2/Q3/Q8, tutorials-are-pointers 5.4b-P): this
used to describe a characteristic-length placeholder `__LC__`, substituted by
`setup/studies/tetConvergence/run_mono_tet.sh`. The template now carries a
`DefineConstant[ lc = {0.1, Name "lc"} ]` instead (native `21f7bc82e`), set
directly on the `gmsh` command line; `run_mono_tet.sh` no longer exists.

#### Effective Mesh Spacing and Observed Order

The manufactured verifier back-computes an *effective* spacing `dx = 1/round(cbrt(nCells))` from the total cell count. For an unstructured tet mesh this is the mean cell size and the correct convergence abscissa, `p = log(e_coarse/e_fine) / log(dx_coarse/dx_fine)`, rather than assuming factor-of-two refinement.

Corrected 2026-09-27: this used to also cite
`setup/studies/tetConvergence/summarize_tet.py` computing the observed order
next to `checkMesh`'s max non-orthogonality and max skewness. That script no
longer exists in this study directory; nothing here recomputes it.

#### Tetrahedral Convergence Sweep

Tetrahedral mesh convergence study via driverFOAM:

```bash
driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetConvergence/sweep_tet_generic.json
python3 applications/scripts/paperI_results/aggregate.py tet
```

Runs both `Gauss linear` and `leastSquares` gradient reconstruction across multiple resolutions on unstructured tet mesh.

## Usage

### Manual Execution

```bash
./Allrun
./regressionTest.sh
```

### Driver-Managed Sweeps (Suggested)

Spatial convergence (hex mesh):

```bash
driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/cartesianConvergence/sweep_hex_convergence.json --output-dir .tmp/driverfoam/monodomainPseudoECG-cartesian
python3 applications/scripts/paperI_results/aggregate.py mono_spatial
python3 applications/scripts/paperI_results/aggregate.py pseudo_ecg_spatial
```

Temporal discretization (fixed fine mesh, `dt` refinement). The finest 1D and 2D studies use `N = 640`; 3D sweeps did not reach clean asymptotic regime:

```bash
driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/temporalConvergence/sweep_temporal_convergence.json --output-dir .tmp/driverfoam/monodomainPseudoECG-temporal
```

The temporal spec holds a fixed fine mesh (`N = 640` in 1D/2D and `N = 160`
in 3D) while halving `dt`.  It retains pseudo-ECG outputs and uses explicit
RKF45 baseline tolerances (`absTol=1e-10`, `relTol=1e-8`), so field and
functional temporal effects can be compared from the same runs.

The tetrahedral spatial paths use `dt ~ h^2` and are therefore combined
space--time paths.  Their mesh-fixed `dt/2` controls live in
`setup/studies/tetTemporalControl/`; do not label a tet slope as purely spatial
unless the control change is smaller than the accepted field-error separation.

The checked-in sweep JSON files are the source of truth for these studies.

Tetrahedral mesh variant (example of running a study):

```bash
driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetConvergence/sweep_tet_generic.json
python3 applications/scripts/paperI_results/aggregate.py tet
```

Tetrahedral mesh-fixed timestep controls:

```bash
driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetTemporalControl/sweep_tet_dt_half.json --output-dir .tmp/driverfoam/monodomainPseudoECG-tet-dt-half
```
