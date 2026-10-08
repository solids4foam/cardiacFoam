# manufacturedSolutions/monodomainPseudoECG tutorial

Manufactured-solution verification for the monodomain stack with pseudo-ECG verification.

Origin: the monodomain field problem (`manufacturedFDAMonodomainVerifier`) is the monodomain problem of the FDA regulatory science tool (Pathmanathan & Gray 2014). The rotated-anisotropy extension (`manufacturedAnisotropicMonodomainVerifier`) and the pseudo-ECG problem (`manufacturedPseudoECGVerifier`) are this work.

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

## Variants & Extensions

### Insulated-Wall Case

[`insulatedWall/`](insulatedWall/README.md) is a separate case with walls that are insulated in the conormal sense, `n.G.grad(Vm) = 0` with `n.grad(Vm) != 0`, on a domain periodic in `y` and `z`. The exact solution here satisfies both `n.grad(Vm) = 0` and `n.G.grad(Vm) = 0` on every wall, so it cannot tell the two apart; `insulatedWall/` can. It has its own regression and tetrahedral study.

### Tetrahedral (unstructured) Mesh Variant

`setup/studies/tetConvergence/` is an activatable overlay of this same case on a genuinely unstructured mesh: identical `constant/` and `system/` dicts (electroProperties, physicsProperties, controlDict, decomposeParDict), except the mesh generator changes.

`setup/studies/tetConvergence/fvSchemes` sets every entry to `system/fvSchemes`'s own value and is not installed by the record's tet route; the tet route uses the case's own `system/fvSolution` and `system/fvSchemes`.

#### Tetrahedral Variant Purpose

Verifies that OpenFOAM's non-orthogonal `Gauss linear corrected` Laplacian scheme sustains its spatial convergence rate on a genuinely unstructured tetrahedral mesh of the unit cube:

- a smoothly perturbed hex mesh folds before reaching high non-orthogonality (~8 deg), whereas tets naturally reach tens of degrees
- cardiac anatomical meshes are predominantly tetrahedral, the topology the framework must actually handle

#### Grid Generation

`setup/studies/tetConvergence/box.geo.template` is a gmsh (OpenCASCADE) unit cube. Resolution is set with `gmsh -3 box.geo.template -setnumber lc <value>` (`lc = 1/N`), meshes with gmsh (legacy msh2 format), and imports via `gmshToFoam`. All six boundary faces lie on axis-aligned planes `x,y,z in {0,1}`, where the manufactured cosine field has zero normal derivative, so the solver's default zeroGradient boundary stays compatible with the exact solution.

#### Effective Mesh Spacing and Observed Order

The manufactured verifier back-computes an *effective* spacing `dx = 1/round(cbrt(nCells))` from the total cell count. For an unstructured tet mesh this is the mean cell size and the correct convergence abscissa, `p = log(e_coarse/e_fine) / log(dx_coarse/dx_fine)`, rather than assuming factor-of-two refinement.

#### Tetrahedral Convergence Sweep

Tetrahedral mesh convergence study via omnidriver. Spec: `setup/studies/tetConvergence/sweep_tet_generic.json`.

```bash
[omnidriver command to run]
python3 applications/scripts/paperI_results/aggregate.py tet
```

Runs both `Gauss linear` and `leastSquares` gradient reconstruction across multiple resolutions on unstructured tet mesh.

## Usage

### Manual Execution

```bash
./Allrun
./regression/regressionTest.sh
```

### Omnidriver-Managed Sweeps (Suggested)

Spatial convergence (hex mesh). Spec: `setup/studies/cartesianConvergence/sweep_hex_convergence.json`.

```bash
[omnidriver command to run]
python3 applications/scripts/paperI_results/aggregate.py mono_spatial
python3 applications/scripts/paperI_results/aggregate.py pseudo_ecg_spatial
```

Temporal discretization (fixed fine mesh, `dt` refinement). The finest 1D and 2D studies use `N = 640`; 3D sweeps did not reach clean asymptotic regime. Spec: `setup/studies/temporalConvergence/sweep_temporal_convergence.json`.

```bash
[omnidriver command to run]
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

Tetrahedral mesh variant (example of running a study). Spec: `setup/studies/tetConvergence/sweep_tet_generic.json`.

```bash
[omnidriver command to run]
python3 applications/scripts/paperI_results/aggregate.py tet
```

Tetrahedral mesh-fixed timestep controls. Spec: `setup/studies/tetTemporalControl/sweep_tet_dt_half.json`.

```bash
[omnidriver command to run]
```
