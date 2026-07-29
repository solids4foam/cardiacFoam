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

Driver-managed sweeps:

```bash
applications/scripts/driverFoam/bin/driverFoam all --entry manufacturedFDABidomain --config tutorials/manufacturedSolutions/bidomain/setup/driver_config.json
```

After the sweep completes, persist the canonical Paper I convergence table:

```bash
python3 applications/scripts/paperI_results/aggregate.py bidomain
```

## Tetrahedral (unstructured) mesh variant

`setup/mesh/tet/` holds this case's own tetrahedral-mesh overlay: a unit-cube
Delaunay mesh (`box.geo.template`, gmsh OpenCASCADE, characteristic length
placeholder `__LC__`) and an `fvSchemes` copy with `gradSchemes.default`
forced to `leastSquares`. The geometry and gradient-scheme override are
byte-identical to `monodomainPseudoECG`'s and `eikonalECG`'s own tet
overlays -- all three manufactured-solution families refine on the same unit
cube -- but are kept as a local copy here rather than referenced across
cases, matching how every merged tet overlay in this repo is scoped to its
own case.

### Gradient-scheme convergence sweep (paper table)

```bash
cd tutorials/manufacturedSolutions/bidomain
bash setup/mesh/tet/run_scheme_study.sh
```

Runs both `Gauss linear` and `leastSquares` gradient reconstruction across
`N = 10, 20, 40, 80` (overridable via `RESOLUTIONS`) for the coupled `Vm`/
`phiE` (gauge-shifted) bidomain system, writes
`setup/results/scheme_study.csv`, and persists the canonical Paper I table:

```bash
python3 applications/scripts/paperI_results/aggregate.py bidomain_tet
```

`nOuterCorrectors` is left at this case's own default (2): the corrector
study below already established that two outer sweeps are within 1% of the
fully converged block, so the reported spatial order isn't iteration-error
limited.

### Corrector study purpose

`setup/mesh/tet/studies/corrector/run_corrector_study.sh` produces
`@tbl-bidomain-corrector-sensitivity`: a same-mesh sensitivity screen that
separates two solver-loop controls the segregated bidomain equations expose
on this tetrahedral family --

- an outer sweep, which repeats the coupled `phiE -> Vm` block; and
- an equation-level non-orthogonal reassembly, which resolves each corrected
  equation before advancing to the next block.

The four reported variants (`baseline`, `outer2`, `nonorth1`, `combined`)
cross `nOuterCorrectors = 1,2` with `nNonOrthogonalCorrectors = 0,1` on the
`N=10,20,40` Delaunay meshes; `run_corrector_study.sh` also defines
`outer3`/`outer4`/`outer8`/`outer16` variants used for exploratory screening
but not reported in the paper table. These are same-mesh, same-time-step
iteration-sensitivity controls, not an additional spatial convergence study.

### Running the corrector study

```bash
cd tutorials/manufacturedSolutions/bidomain
bash setup/mesh/tet/studies/corrector/run_corrector_study.sh
```

Override `RESOLUTIONS`, `VARIANTS`, `RESULTS_DIR`, or `KEEP_WORK` through the
environment for focused reruns. `MESH_MODE=ortho` instead exercises the
case's own `blockMeshDict.3D` ladder (the configuration behind the reported
Cartesian `bidomain` convergence table) rather than the tet overlay.
