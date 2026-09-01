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

- manufactured field summaries in `postProcessing/`
- `postProcessing/pseudoECG.dat`
- `postProcessing/manufacturedPseudoECG.dat`
- `postProcessing/manufacturedPseudoECGSummary.dat`

## Variants & Extensions

### Tetrahedral (unstructured) Mesh Variant

`setup/studies/tetConvergence/` is an activatable overlay of this same case on a genuinely unstructured mesh: identical `constant/` and `system/` dicts (electroProperties, physicsProperties, fvSchemes, controlDict, decomposeParDict), except the mesh generator changes and `setup/studies/tetConvergence/fvSolution` (a tighter `nOuterCorrectors`/`nNonOrthogonalCorrectors` pair) is swapped in for the duration of a tet run and restored on exit. It was formerly the standalone `monodomainTetMMS` tutorial, merged here the same way `eikonalTetMMS` was merged into `eikonalECG`.

#### Tetrahedral Variant Purpose

Verifies that OpenFOAM's non-orthogonal `Gauss linear corrected` Laplacian scheme sustains its spatial convergence rate on a genuinely unstructured tetrahedral mesh of the unit cube:

- a smoothly perturbed hex mesh folds before reaching high non-orthogonality (~8 deg), whereas tets naturally reach tens of degrees
- cardiac anatomical meshes are predominantly tetrahedral, the topology the framework must actually handle

#### Grid Generation

`setup/studies/tetConvergence/box.geo.template` is a gmsh (OpenCASCADE) unit cube with a characteristic length placeholder `__LC__`. `setup/studies/tetConvergence/run_mono_tet.sh` substitutes `lc = 1/N` per resolution, meshes with gmsh (legacy msh2 format), and imports via `gmshToFoam`. All six boundary faces lie on axis-aligned planes `x,y,z in {0,1}`, where the manufactured cosine field has zero normal derivative, so the solver's default zeroGradient boundary stays compatible with the exact solution.

#### Effective Mesh Spacing and Observed Order

The manufactured verifier back-computes an *effective* spacing `dx = 1/round(cbrt(nCells))` from the total cell count. For an unstructured tet mesh this is the mean cell size and the correct convergence abscissa. `setup/studies/tetConvergence/summarize_tet.py` computes the observed order from consecutive `dx` values, `p = log(e_coarse/e_fine) / log(dx_coarse/dx_fine)`, rather than assuming factor-of-two refinement, and reports it next to `checkMesh` max non-orthogonality and max skewness.

#### Tetrahedral Convergence Sweep

Tetrahedral mesh convergence study via driverFOAM:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetConvergence/sweep_tet_generic.json
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
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/cartesianConvergence/sweep_hex_convergence.json --output-dir .tmp/driverfoam/monodomainPseudoECG-cartesian
python3 applications/scripts/paperI_results/aggregate.py mono_spatial
python3 applications/scripts/paperI_results/aggregate.py pseudo_ecg_spatial
```

Temporal discretization (fixed fine mesh, `dt` refinement). The finest 1D and 2D studies use `N = 640`; 3D sweeps did not reach clean asymptotic regime:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/temporalConvergence/sweep_temporal_convergence.json --output-dir .tmp/driverfoam/monodomainPseudoECG-temporal
```

The temporal spec holds a fixed fine mesh (`N = 640` in 1D/2D and `N = 160`
in 3D) while halving `dt`.  It retains pseudo-ECG outputs and uses explicit
RKF45 baseline tolerances (`absTol=1e-10`, `relTol=1e-8`), so field and
functional temporal effects can be compared from the same runs.

The tetrahedral spatial paths use `dt ~ h^2` and are therefore combined
space--time paths.  Their mesh-fixed `dt/2` controls live in
`setup/studies/tetTemporalControl/`; do not label a tet slope as purely spatial
unless the control change is smaller than the accepted field-error separation.

The finest fixed-grid temporal level of every dimensional ladder is repeated
with tighter RKF45 tolerances in `setup/studies/odeToleranceControl/`.  This
is required to attribute the measured complete-advance temporal response to
the configured Lie--Godunov/PDE path rather than an inherited ODE default.

The checked-in sweep JSON files are the source of truth for these studies.

RKF45 tolerance control (one finest fixed-grid temporal case per dimension):

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/odeToleranceControl/sweep_ode_tolerance.json --output-dir .tmp/driverfoam/monodomainPseudoECG-ode-control
```

Tetrahedral mesh variant (example of running a study):

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetConvergence/sweep_tet_generic.json
python3 applications/scripts/paperI_results/aggregate.py tet
```

Tetrahedral mesh-fixed timestep controls:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/tetTemporalControl/sweep_tet_dt_half.json --output-dir .tmp/driverfoam/monodomainPseudoECG-tet-dt-half
```
