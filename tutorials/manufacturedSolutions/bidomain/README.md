# manufacturedSolutions/bidomain tutorial

Manufactured-solution verification for the one-mesh bidomain stack.

Origin: the bidomain problem of the FDA regulatory science tool (Pathmanathan & Gray 2014).

## Overview

### Stack

- myocardium solver: `bidomainSolver`
- ionic model: `bidomainFDAManufactured`
- field verification: `manufacturedFDABidomainVerifier` from `libverificationModels`
- shared analytical oracle: `verificationModels`

### Purpose

This case verifies:

- field convergence against an analytical oracle
- manufactured ionic export variables (`u1`, `u2`, `u3`)
- bidomain potentials (`phiE`, `phiI`) on a single mesh

### Key Configuration

For this workflow, the ionic model exposes manufactured verification metadata.

- ionic-model-side behavior: `bidomainFDAManufactured`
- analytical oracle: `verificationModels`
- field verification hook: `modelPrePostProcessors`
- field verifier: `verificationModels/bidomainVerification`

Typical outputs include:

- manufactured field summaries in `postProcessing/`
- `Vm`, `phiE`, `phiI`

## Variants & Extensions

### Insulated-Wall Case

[`insulatedWall/`](insulatedWall/README.md) is a separate case with unequal anisotropy and walls that are insulated in the conormal sense, `n.G_i.grad(phi_i) = 0` and `n.G_e.grad(phi_e) = 0`, on a domain periodic in `y` and `z`. The diagonal conductivities here make `n.grad(u) = 0` and `n.G.grad(u) = 0` coincide on the walls; in `insulatedWall/` they differ. It has its own regression and tetrahedral study.

### Tetrahedral (unstructured) Mesh Variant

`setup/studies/tetConvergence/` holds this case's own tetrahedral-mesh overlay, co-located with the study that drives it: a unit-cube Delaunay mesh (`box.geo.template`, gmsh OpenCASCADE, characteristic length set by a `DefineConstant`, overridden with `gmsh -setnumber lc <value>`).

The tet route uses the case's own `system/fvSchemes` directly.

#### Gradient-Scheme Convergence Sweep

Runs both `Gauss linear` and `leastSquares` gradient reconstruction across `N = 10, 20, 40, 80` for the coupled `Vm`/`phiE` (gauge-shifted) bidomain system. The primary ladder uses the case's two outer correctors with `1e-15` linear tolerance and explicit RKF45 controls. Its rates remain combined space--time evidence until the mesh-fixed `dt/2` controls show that the temporal perturbation is smaller than the resolved field-error separation.

#### Corrector Study

`setup/studies/corrector/sweep_corrector_study.json` produces a same-mesh sensitivity screen that separates two solver-loop controls the segregated bidomain equations expose on this tetrahedral family:

- an outer sweep, which repeats the coupled `phiE -> Vm` block
- an equation-level non-orthogonal reassembly, which resolves each corrected equation before advancing to the next block

The four reported variants (`baseline`, `outer2`, `nonorth1`, `combined`) cross `nOuterCorrectors = 1,2` with `nNonOrthogonalCorrectors = 0,1` on the `N=10,20,40` Delaunay meshes. `setup/studies/correctorN80/` repeats those variants at the finest level. Together these are same-mesh, same-time-step iteration-sensitivity controls that test solver-loop convergence rather than spatial resolution.

## Usage

### Manual Execution

```bash
./Allrun
./regression/regressionTest.sh
```

### Omnidriver-Managed Sweeps (Suggested)

Cartesian spatial convergence. Spec: `setup/studies/cartesianConvergence/sweep_hex_convergence.json`.

```bash
[omnidriver command to run]
python3 applications/scripts/paperI_results/aggregate.py bidomain_cartesian
```

Temporal convergence. Spec: `setup/studies/temporalConvergence/sweep_temporal_convergence.json`.

```bash
[omnidriver command to run]
python3 applications/scripts/paperI_results/aggregate.py bidomain_temporal
```

Tetrahedral gradient-scheme and corrector studies. Specs: `setup/studies/tetConvergence/sweep_tet_generic.json`, `setup/studies/corrector/sweep_corrector_study.json`, `setup/studies/linearToleranceControl/sweep_tet_phi_tolerance.json`, `setup/studies/tetTemporalControl/sweep_tet_dt_half.json`, `setup/studies/correctorN80/sweep_corrector_n80.json`.

```bash
[omnidriver command to run]
python3 applications/scripts/paperI_results/aggregate.py bidomain_tet_generic
```

The complete rerun matrix contains 54 cases: Cartesian spatial (12), primary
tetrahedral reconstruction (8), fixed-grid temporal (8), ODE (2), loose
`phiE|phiI` tolerance (4), mesh-fixed tetrahedral `dt/2` (4), and corrector
controls (12 plus 4 at `N=80`). Run the listed omnidriver-managed specifications;
do not rely on an aggregate wrapper unless it has been versioned with the
release.
