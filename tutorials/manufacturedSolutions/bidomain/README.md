# manufacturedSolutions/bidomain tutorial

This tutorial is the one-mesh manufactured-solution bidomain check.

This is the manufactured-solution verification workflow for the one-mesh
bidomain stack.

## Stack

- myocardium solver: `bidomainSolver`
- ionic model: `bidomainFDAManufactured`
- field verification:
  `manufacturedFDABidomainVerifier` from `libverificationModels`
- shared analytical oracle:
  `verificationModels`

## Purpose

This case verifies:

- field convergence against an analytical oracle
- manufactured ionic export variables (`u1`, `u2`, `u3`)
- bidomain potentials (`phiE`, `phiI`) on a single mesh

## Key configuration idea

For this workflow, the ionic model exposes manufactured verification metadata,
and the exact manufactured reference no longer lives inside the ionic-model
folder.

- ionic-model-side behavior: `bidomainFDAManufactured`
- analytical oracle: `verificationModels`
- field verification hook: `modelPrePostProcessors`
- field verifier: `verificationModels/bidomainVerification`

## Outputs

Typical outputs include:

- manufactured field summaries in `postProcessing/`
- `Vm`
- `phiE`
- `phiI`

## Execution

Manual:

```bash
blockMesh -dict system/blockMeshDict.1D
./Allrun
./regressionTest.sh
```

Driver-managed sweeps (Cartesian spatial convergence):

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/bidomain/setup/studies/cartesianConvergence/sweep_hex_convergence.json
python3 applications/scripts/paperI_results/aggregate.py bidomain_cartesian
```

Temporal convergence:

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/bidomain/setup/studies/temporalConvergence/sweep_temporal_convergence.json
python3 applications/scripts/paperI_results/aggregate.py bidomain_temporal
```

Or run every registered experiment for this case through the normalized registry:

```bash
./reproduce_verification.sh bidomain_cartesian bidomain_temporal bidomain_tet_generic
```

## Tetrahedral (unstructured) mesh variant

`setup/studies/tetConvergence/` holds this case's own tetrahedral-mesh
overlay, co-located with the study that drives it: a unit-cube Delaunay mesh
(`box.geo.template`, gmsh OpenCASCADE, characteristic length placeholder
`__LC__`) and an `fvSchemes` copy with `gradSchemes.default` forced to
`leastSquares`. The geometry and gradient-scheme override are byte-identical
to `monodomainPseudoECG`'s and `eikonalECG`'s own tet overlays -- all three
manufactured-solution families refine on the same unit cube -- but are kept
as a local copy here rather than referenced across cases, matching how every
merged tet overlay in this repo is scoped to its own case.

### Gradient-scheme convergence sweep (paper table)

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/bidomain/setup/studies/tetConvergence/sweep_tet_generic.json
python3 applications/scripts/paperI_results/aggregate.py bidomain_tet_generic
```

Runs both `Gauss linear` and `leastSquares` gradient reconstruction across
`N = 10, 20, 40, 80` for the coupled `Vm`/`phiE` (gauge-shifted) bidomain
system and persists the canonical Paper I table.

`nOuterCorrectors` is left at this case's own default (2): the corrector
study below already established that two outer sweeps are within 1% of the
fully converged block, so the reported spatial order isn't iteration-error
limited.

### Corrector study purpose

`setup/studies/corrector/sweep_corrector_study.json` produces
`@tbl-bidomain-corrector-sensitivity`: a same-mesh sensitivity screen that
separates two solver-loop controls the segregated bidomain equations expose
on this tetrahedral family --

- an outer sweep, which repeats the coupled `phiE -> Vm` block; and
- an equation-level non-orthogonal reassembly, which resolves each corrected
  equation before advancing to the next block.

The four reported variants (`baseline`, `outer2`, `nonorth1`, `combined`)
cross `nOuterCorrectors = 1,2` with `nNonOrthogonalCorrectors = 0,1` on the
`N=10,20,40` Delaunay meshes -- a same-mesh, same-time-step
iteration-sensitivity screen, not an additional spatial convergence study.

### Running the corrector study

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec tutorials/manufacturedSolutions/bidomain/setup/studies/corrector/sweep_corrector_study.json
```

See `setup/studies/corrector/README.md` for how the variants map to
`n_outer_correctors`/`n_nonorthogonal_correctors` overrides.
