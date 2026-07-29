# manufacturedSolutions/monodomainPseudoECG tutorial

This is the manufactured-solution verification workflow for the monodomain
stack with pseudo-ECG verification.

## Stack

- myocardium solver: `monodomainSolver`
- ionic model: `monodomainFDAManufactured`
- field verification:
  `manufacturedFDAMonodomainVerifier` from `libverificationModels`
- shared analytical oracle:
  `verificationModels`
- optional manufactured pseudo-ECG verification:
  `pseudoECGManufacturedVerifier` from `libverificationModels`

## Purpose

This case verifies:

- field convergence against an analytical oracle
- manufactured ionic export variables (`u1`, `u2`, `u3`)
- manufactured pseudo-ECG reference output

## Key configuration idea

For this workflow, the ionic model exposes manufactured verification metadata,
and the exact manufactured reference no longer lives inside the ionic-model
folder.

- ionic-model-side behavior: `monodomainFDAManufactured`
- analytical oracle: `verificationModels`
- field verification hook: `modelPrePostProcessors`
- field verifier: `verificationModels/monodomainVerification`
- ECG verifier: `verificationModels/ecgVerification`

## Outputs

Typical outputs include:

- manufactured field summaries in `postProcessing/`
- `postProcessing/pseudoECG.dat`
- `postProcessing/manufacturedPseudoECG.dat`
- `postProcessing/manufacturedPseudoECGSummary.dat`

## Execution

Manual:

```bash
./Allrun
./regressionTest.sh
```

Driver-managed sweeps:

```bash
applications/scripts/driverFoam/bin/driverFoam all --entry manufacturedFDA --config tutorials/manufacturedSolutions/monodomainPseudoECG/setup/driver_config.json
```

After the sweep completes, persist the canonical Paper I convergence tables (reads
`driverPostProcessingArchive_postProcessing/`, writes `setup/results/*.csv`, never
touches the sweep's own output):

```bash
python3 applications/scripts/paperI_results/aggregate.py mono_spatial
python3 applications/scripts/paperI_results/aggregate.py pseudo_ecg_spatial
```

Temporal-discretization sweep (fixed fine mesh, `dt` refinement). The paper reports
the finest 1D and 2D studies (`N = 640`); the 3D temporal sweeps did not reach a clean
asymptotic regime and were dropped, and the coarser 2D `N = 320` study was superseded
by `N = 640`.

```bash
applications/scripts/driverFoam/bin/driverFoam all --entry manufacturedFDA --config tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/temporal1D_N640/config.json
applications/scripts/driverFoam/bin/driverFoam all --entry manufacturedFDA --config tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/temporal2D_N640/config.json
```

Each temporal config holds a fixed fine mesh (`N = 640`) with a `dt` ladder
(`[1.121075e-3, 5.60538e-4, 2.80269e-4, 1.401345e-4]`), disables ECG post-processing,
and isolates the field temporal-order measurement.

## Tetrahedral (unstructured) mesh variant

`setup/mesh/tet/` is an activatable overlay of this same case on a genuinely
unstructured mesh: identical `constant/` and `system/` dicts (electroProperties,
physicsProperties, fvSchemes, controlDict, decomposeParDict), except the mesh
generator changes and `setup/mesh/tet/fvSolution` (a tighter
`nOuterCorrectors`/`nNonOrthogonalCorrectors` pair) is swapped in for the
duration of a tet run and restored on exit. It was formerly the standalone
`monodomainTetMMS` tutorial, merged in here the same way `eikonalTetMMS` was
merged into `eikonalECG`.

### Tetrahedral variant purpose

Verifies that OpenFOAM's non-orthogonal `Gauss linear corrected` Laplacian
scheme sustains its spatial convergence rate on a genuinely unstructured
tetrahedral mesh of the unit cube:

- a smoothly perturbed hex mesh folds before it reaches high
  non-orthogonality (its ceiling is ~8 deg), whereas tets of the cube reach
  tens of degrees naturally; and
- cardiac anatomical meshes are predominantly tetrahedral, so this is the mesh
  topology the framework must actually cope with.

### How the mesh is generated

`setup/mesh/tet/box.geo.template` is a gmsh (OpenCASCADE) unit cube with a
characteristic length placeholder `__LC__`. `setup/mesh/tet/run_tet_sweep.sh`
substitutes `lc = 1/N` per resolution, meshes with gmsh (legacy msh2 format),
and imports via `gmshToFoam`. All six boundary faces lie on the axis-aligned
planes `x,y,z in {0,1}`, where the manufactured cosine field has zero normal
derivative, so the solver's default zeroGradient (no-flux) boundary stays
compatible with the exact solution -- exactly as on the hex mesh. No `0/`
fields or `createPatch` step is needed: cardiacFoam creates `Vm` with a uniform
zeroGradient boundary (`READ_IF_PRESENT`), so gmsh patch naming is irrelevant.

### Effective mesh spacing and observed order

The manufactured verifier assumes a structured mesh and back-computes an
*effective* spacing `dx = 1/round(cbrt(nCells))` from the total cell count.
For an unstructured tet mesh of the unit cube this is the mean cell size, and
it is the correct convergence abscissa. `setup/mesh/tet/summarize_tet.py`
therefore computes the observed order from consecutive `dx` values,
`p = log(e_coarse/e_fine) / log(dx_coarse/dx_fine)`, rather than assuming a
factor-of-two refinement, and reports it next to the `checkMesh` max
non-orthogonality and max skewness so the mesh quality is explicit.

### Gradient-scheme convergence sweep (paper table)

```bash
cd tutorials/manufacturedSolutions/monodomainPseudoECG
bash setup/mesh/tet/run_scheme_study.sh
```

Runs both `Gauss linear` and `leastSquares` gradient reconstruction across
`N = 10, 20, 40, 80` (overridable via `RESOLUTIONS`), writes
`setup/results/scheme_study.csv`, and persists the canonical Paper I table:

```bash
python3 applications/scripts/paperI_results/aggregate.py tet
```

### Mesh-quality-only sweep

```bash
RESOLUTIONS="10 20 40" ENDTIME=0.02 bash setup/mesh/tet/run_tet_sweep.sh
```

- `RESOLUTIONS` -- space-separated nominal cells-per-side (default `10 20 40`).
- `ENDTIME` -- integration window (default `0.02`; the backward `ddt` scheme is
  O(dt^2), so a short window keeps the temporal error subdominant and tet runs
  tractable). Set `ENDTIME=0.2` to match the hex baseline exactly for the paper
  table.

Requires OpenFOAM sourced and `gmsh` on `PATH`. Results land in
`setup/results/<N>/`, and a combined `setup/results/summary.csv` is written at
the end.
