# Idealized biventricular electromechanics

## What it does

This case uses cardiacFoam's sequential electromechanical model on the idealized biventricular heart: monodomain (or eikonal) electrophysiology with the TNNP ionic model in an `electro` region, the Land--Niederer (`LandNiedererTWorld`) active-tension model, and a total-Lagrangian neo-Hookean solid in a `solid` region. It paces the apex and runs to 20 ms. Mechanics uses the rule-based fibre field `f0`, copied from `../mesh/0/fiber`.

- **Exercises:** the electromechanics coupling on an anatomical mesh, with the fibre field interpolated onto the solid's faces (`interpolateFibreField -region solid` makes `f0f`). It is the only regression that runs electromechanics on an anatomical mesh and the only one that runs `interpolateFibreField`.
- **Use it for:** a starting point for whole-heart electromechanics (add cavity pressures and a pericardial or basal spring, see Mechanical boundary conditions), and for checking a solids4foam build against an anatomical case.

## Regression

Tier: `full`. Markers: `needs-solids4foam` (exits 77 in `lightweight` build mode).

`regression/regressionTest.sh` runs `./Allrun parallel` (the `monodomain` variant; 6 subdomains per region, `simple` strips along x) with `endTime` set to 0.01 instead of the tutorial's 0.02. It compares `regression/electroMechHeart.reference` (10 rows, values frozen from an earlier run of the code): `activationTime` at the apex probe (3.2 ms, tolerance 1e-5 s) and at the mid-wall probe (still -1); the displacement `D` at both probes at 5 and 10 ms (apex x and z, mid-wall y; tolerance 1e-6 m, 2e-6 m for the mid-wall); and the apical active tension `Ta` at 5 and 10 ms (tolerance 10 Pa). The `eikonal` variant is not in the regression. This is the slowest regression of the suite.

```bash
./regression/regressionTest.sh
```

## Mesh and activation

The shared mesh in `../mesh/` is already scaled uniformly by `0.02` from the raw geometry, giving an approximately `80 x 100 x 80 mm`, `226 mL` idealized adult-heart-scale mesh in SI metres. `Allrun` copies it in as-is — no scaling happens here. The apicobasal coordinate has its minimum at the anatomical apex, even though Cartesian x increases toward it. The external stimulus starts at `2 ms` and occupies a verified 26-cell apical box:

```text
x = [0.0756, 0.0800] m
y = [-0.0040, 0.0020] m
z = [0.0030, 0.0060] m
```

`monodomainSolverCoeffs.ionicHeterogeneity` (`mode namedRegions;`) classifies cells into `endocardialCells`/`mCells`/`epicardialCells` regions (`range (0 0.3)`/`range (0.3 0.7)`/`range (0.7 1)`) from a field that must be 0 at the endocardium and 1 at the epicardium — the opposite of this mesh's `tm` convention. `system/electro/setExprFieldsDict` derives `t = 1 - tm` for it to read; `Allrun` runs `setExprFields -region electro` for this.

## Variants

`./Allrun` runs the `monodomain` variant; `./Allrun eikonal` runs the `eikonal` variant. `Allrun` links `constant/electro/electroProperties` and `constant/electroMechanicalProperties` to the chosen variant's `.monodomain` or `.eikonal` files.

- `monodomain`: `monodomainSolver` tissue; the active tension reads `Cai` from the TNNP ionic model in each cell.
- `eikonal`: `eikonalSolver` tissue with the same mesh, conductivity, `c0` as `../electroHeart` and the same apical box, where the activation time is 0. The activation is solved once in the first step. With `activeTensionSignal templates;` the active tension reads `Cai` from TNNP single-cell templates, one per `ionicHeterogeneity` region, at `t - activationTime`. The templates' time zero is their upstroke. `stim_start 0` makes the first template sample TNNP's initial state, the resting state the monodomain variant starts from. The eikonal variant activates the apex at 0 ms, the monodomain variant at 3.2 ms (stimulus at 2 ms).

## Mechanical boundary conditions

`BASE` has zero displacement. `EPI`, `ENDO_LV`, and `ENDO_RV` are traction-free in this first, unloaded activation/contraction case. This is an intentional benchmark simplification: it provides a stable reference configuration and removes rigid-body modes, but omits cavity pressure and pericardial restraint.

For a more physiological pumping simulation, replace the zero endocardial tractions with LV/RV cavity pressures and replace the fixed base with a basal spring constraint; add normal, sliding pericardial restraint at `EPI`.

## Parameter scope

The isotropic passive law (`E=100 kPa`, `nu=0.3`) is a first-integration-test approximation. `monodomainSolverCoeffs` sets `conductivitySource field;`, reading the anisotropic conductivity tensor `Conductivity` (`0/electro/Conductivity`, a `volSymmTensorField`) rather than a uniform value. Active mechanics use the rule-based fibre field `f0`; the passive mechanics remain an isotropic neo-Hookean approximation. Review and calibrate passive, active, cavity-pressure, and pericardial parameters before using this case for physiological predictions.

## Run prerequisites

This needs cardiacFoam compiled with solids4foam, including `libelectroMechanicalModels`. From this directory, run `./Allrun` (or `./Allrun eikonal`, and `parallel` for either). The mesh and the `fiber`/`sheet`/`tm`/`Conductivity` fields are not committed in this case directory — `../prepareHeartMesh` writes them in `../mesh/` from the heart/torso mesh `../mesh/heartTorso` and `Allrun` copies them in each run (the `fiber` field becomes this case's `f0`). `f0f` is made from `f0` by `interpolateFibreField -region solid`, with the solid region's `skewCorrected linear` interpolation; `t` is derived from `tm`. `Allclean` removes those copies and the derived `t` along with the generated regional meshes and run output.

## Notes on dictionary entries

**Stimulus.** The apical stimulus box in `constant/electro/electroProperties`
contains 26 cell centres at the mesh apex, where the supplied `apicobasal`
coordinate is smallest. All values are SI: the shared mesh in `../mesh/` is
already in adult-heart-scale metres and `Allrun` copies it as-is. The ionic
and mechanical states remain at rest for 2 ms before the apical pacing pulse
begins.

**Run length.** `endTime` is 0.02 s with `deltaT` 1e-5 s; the write interval
is 5 ms.

**Active tension.** `0/solid/Ta` is the initial active tension. The
electromechanical model registers and solves `Ta`, but the solid model builds
its stiffness while it is constructed, before that field exists, so
`electroMechanicalLaw` reads `Ta` from this file. Both
`constant/electroMechanicalProperties.*` select the scalar
`LandNiedererTWorld` model, integrated with the `ODESolver` settings (`RKF45`,
`absTol 1e-8`, `relTol 1e-6`); the `eikonal` variant adds
`activeTensionSignal templates` and its `activeTensionTemplates` (see
Variants). The batched versions of the active-tension models are listed in
[`src/activeTensionModels/README.md`](../../../src/activeTensionModels/README.md).

**Solid solver.** In `constant/solid/solidProperties`, `rTol` is the relative
residual (converged when residual < rTol x initial residual; the solver
plateaus at ~12 % of the initial residual because of under-relaxation, which
0.2 accepts), `aTol` the absolute residual, and `sTol` the step norm
(converged when deltaX < sTol x xNorm). `0/solid/D`: the basal cut plane
anchors this simplified, unloaded ventricle.

**Mesh and decomposition.** `system/decomposeParDict` uses six contiguous
base-to-apex strips along global x, matching the regional electro and solid
decomposition dictionaries. The top-level `system/fvSchemes` and
`system/fvSolution` are dummies for the parent region and are not used
directly. `system/electro/setExprFieldsDict` preloads `tm` so the expression
reads an owned, mutable field.

**Probes.** `Taprobes` and `CaiProbes` sample the active tension and `Cai` at
the apex (inside the stimulus region, `(0.077286 -0.000616 0.004096)`) and the
mid-wall (later activation, `(0.024482 0.017758 -0.001906)`).
