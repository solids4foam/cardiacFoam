# manufacturedSolutions/eikonalECG tutorial

Manufactured-solution verification for the eikonal activation-time solve and its template ECG calculation.

## Overview

### Stack

- myocardium solver: `eikonalSolver`
- field verification: `manufacturedEikonalVerifier` from `libverificationModels`
- ECG verification: `manufacturedEikonalECGVerifier` from `libverificationModels`
- shared analytical oracle: `verificationModels`

### Purpose

This case verifies:

- activation-time (`psi`) convergence against an analytical oracle
  ($\tau(x) = \exp(k \cdot x)$), compared point-by-point at each cell centre
  (no quadrature needed -- the manufactured solution is exact there)
- the eikonal template ECG, `Vm(x,t) = U(t - psi(x))`, sampled internally and
  written to `postProcessing/eikonalECG.dat`, checked against a
  Gauss-Legendre quadrature reference (`referenceQuadratureOrder = 96`) of
  the manufactured volume integral

### Key Configuration

- solver-side behavior: `eikonalSolverCoeffs` (`constant/electroProperties`)
- analytical oracle: `verificationModels`
- field verifier: `verificationModels/eikonalVerification`
- ECG verifier: `verificationModels/ecgVerification`

Typical outputs include:

- `postProcessing/manufacturedEikonalActivationTime.dat`
- `postProcessing/eikonalECG.dat`
- `postProcessing/manufacturedEikonalECG.dat`
- `postProcessing/manufacturedEikonalECGSummary.dat`

## Variants & Extensions

### Tetrahedral (unstructured) Mesh Variant

`setup/studies/tetConvergence/` is an activatable overlay of this same case
on a genuinely unstructured mesh: identical `constant/` and `system/` dicts,
except the mesh generator changes and `setup/studies/tetConvergence/fvSolution`
is swapped in for the duration of a tet run and restored on exit. It was
formerly the standalone `eikonalTetMMS` tutorial, merged in here (the 5
shared dicts were byte-identical, so this case stayed the canonical hex case
unchanged; only the tet-specific `box.geo.template` and `fvSolution` became
the overlay). `box.geo.template` and its `fvSolution` overlay are co-located
with the study that drives them (`setup/studies/tetConvergence/`) rather than
under `setup/studies/tetConvergence/`, matching `bidomain/setup/studies/tetConvergence/`
and `monodomainPseudoECG`'s own tet variant.

#### Tet Convergence

`setup/studies/tetConvergence/box.geo.template` is a unit-cube gmsh
(OpenCASCADE, Delaunay) template with a characteristic-length placeholder
`__LC__`, instantiated per resolution by the driverFOAM tutorial
(`manufactured_eikonal_ecg.py`'s `render_tet_geo`). Every tet study in this
tutorial (`tetConvergence/sweep_tet_generic.json`, `errorLocalisation/`,
`gradientVerification/`, `gradient_reconstruction/`) renders from this one
template -- there is no separate mesh geometry variant. All six boundary
faces lie on the axis-aligned planes `x,y,z in {0,1}`, where the manufactured
cosine field has zero normal derivative, so the solver's default
zeroGradient (no-flux) boundary stays compatible with the exact solution.

#### Solved-Field Bulk/Boundary Error Decomposition

`setup/studies/errorLocalisation/` sweeps the tet `N` ladder with
`writeErrorField` enabled to separate the solved field's boundary-adjacent
error from its interior (bulk) error.

#### Gradient-Operator Studies

`gradientReconstructionOrder` (`applications/test/gradientReconstructionOrder/`)
exercises the gradient reconstruction operator alone against an exact
analytic field. It is driven as a `gradient_reconstruction=True` workflow_dag
step, appended after the case's solve (see `manufactured_eikonal_ecg.py`'s
`_workflow_dag_for`) -- a real driverFOAM sweep, not bash, the same as every
other study here. `setup/studies/gradientVerification/` covers the full
gaussLinear-vs-leastSquares comparison; `setup/studies/gradient_reconstruction/`
restricts that same matrix to the registered `eikonal_gradient_tet`
experiment's leastSquares-only subset. See each folder's `README.md` for the
exact commands.

## Usage

### Driver-Managed Sweeps

Run this verification suite through driverFOAM. The study manifests below are
the supported execution paths; historic direct OpenFOAM commands and the
retired `reproduce_verification.sh` wrapper are not release procedures.

Cartesian spatial convergence (1D/2D/3D):

    applications/scripts/driverFoam/bin/driverFoam sweep-plan \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/cartesianConvergence/sweep_hex_convergence.json
    applications/scripts/driverFoam/bin/driverFoam sweep-run \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/cartesianConvergence/sweep_hex_convergence.json

Tet convergence:

    applications/scripts/driverFoam/bin/driverFoam sweep-plan \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/tetConvergence/sweep_tet_generic.json
    applications/scripts/driverFoam/bin/driverFoam sweep-run \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/tetConvergence/sweep_tet_generic.json

Nonlinear stopping-criterion control (least-squares tet cases; axis and both
rotated configurations, N=10/20/40/80):

    applications/scripts/driverFoam/bin/driverFoam sweep-plan \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/nonlinearControl/sweep_tet_outer_tolerance.json
    applications/scripts/driverFoam/bin/driverFoam sweep-run \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/nonlinearControl/sweep_tet_outer_tolerance.json

Bulk/boundary error decomposition:

    applications/scripts/driverFoam/bin/driverFoam sweep-plan \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/errorLocalisation/sweep_tet_error_localisation.json
    applications/scripts/driverFoam/bin/driverFoam sweep-run \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/errorLocalisation/sweep_tet_error_localisation.json
    python3 tutorials/manufacturedSolutions/eikonalECG/setup/studies/errorLocalisation/aggregate_bulk_boundary.py

Isolated gradient-operator reconstruction (registered `eikonal_gradient_tet` table, leastSquares only):

    applications/scripts/driverFoam/bin/driverFoam sweep-plan \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradient_reconstruction/sweep_gradient_tet.json
    applications/scripts/driverFoam/bin/driverFoam sweep-run \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradient_reconstruction/sweep_gradient_tet.json
    python3 tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradient_reconstruction/aggregate_gradient_reconstruction.py

Full gaussLinear-vs-leastSquares gradient-operator comparison (not a registered table):

    applications/scripts/driverFoam/bin/driverFoam sweep-plan \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradientVerification/sweep_gradient_tet.json
    applications/scripts/driverFoam/bin/driverFoam sweep-run \
        --spec tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradientVerification/sweep_gradient_tet.json
    python3 tutorials/manufacturedSolutions/eikonalECG/setup/studies/gradientVerification/aggregate_gradient_verification.py

## Effective mesh spacing and observed order

For an unstructured unit-cube tet mesh, use the realised cell-count scale
`h_eff = nCells^(-1/3)`. `setup/studies/tetConvergence/summarize_tet.py`
therefore computes the observed order from consecutive realised scales,
`p = log(e_coarse/e_fine) / log(h_coarse/h_fine)`, and reports each actual
refinement ratio next to the requested N, realised cell count, maximum
non-orthogonality, and maximum skewness. Do not assume a factor-of-two ratio
from the nominal Gmsh length parameter alone.

## Error Calculation and Quadrature

It is important to clarify how the errors are evaluated for the different fields in this verification suite:

1. **Activation Times ($\tau$)**: No quadrature is used here. Because the manufactured activation time is a simple analytical function $\tau(x) = \exp(k \cdot x)$, it can be evaluated exactly at any point. To calculate the error, the solver performs a point-by-point comparison between the numerical activation time solved at each OpenFOAM mesh cell's center and the exact analytical mathematical formula evaluated at that identical cell center point.

2. **The ECG Computation (Where Quadrature is Used)**: Unlike the activation time, the ECG signal is defined mathematically as a volume integral over the entire domain. To get the "exact" baseline reference to compare OpenFOAM against, we calculate the integral of the manufactured analytical gradient field. Because this specific multidimensional integral does not have a simple closed-form algebraic solution, the exact "analytical" reference integral is computed using a highly accurate Gauss-Legendre Quadrature (e.g. $q=96$), which integrates the continuous analytical function down to machine precision. The ECG error is the comparison between OpenFOAM's numerical mesh integration (which simply sums up the cell values $\times$ cell volumes) against this near-perfect reference integral computed using the Gauss-Legendre quadrature.
