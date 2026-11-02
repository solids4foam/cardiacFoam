# Idealized biventricular electromechanics

This case uses cardiacFoam's sequential electromechanical model, the TNNP ionic model, and the Land--Niederer active-tension model. Mechanics uses the rule-based fibre field `f0`, copied from `../mesh/0/fiber`.

## Mesh and activation

The shared mesh in `../mesh/` is already scaled uniformly by `0.02` from the raw geometry, giving an approximately `80 x 100 x 80 mm`, `226 mL` idealized adult-heart-scale mesh in SI metres. `Allrun` copies it in as-is — no scaling happens here. The apicobasal coordinate has its minimum at the anatomical apex, even though Cartesian x increases toward it. The external stimulus starts at `2 ms` and occupies a verified 26-cell apical box:

```text
x = [0.0756, 0.0800] m
y = [-0.0040, 0.0020] m
z = [0.0030, 0.0060] m
```

`monodomainSolverCoeffs.ionicHeterogeneity` (`mode namedRegions;`) classifies cells into `endocardialCells`/`mCells`/`epicardialCells` regions (`range (0 0.3)`/`range (0.3 0.7)`/`range (0.7 1)`) from a field that must be 0 at the endocardium and 1 at the epicardium — the opposite of this mesh's `tm` convention. `system/electro/setExprFieldsDict` derives `t = 1 - tm` for it to read; `Allrun` runs `setExprFields -region electro` for this.

## Mechanical boundary conditions

`BASE` has zero displacement. `EPI`, `ENDO_LV`, and `ENDO_RV` are traction-free in this first, unloaded activation/contraction case. This is an intentional benchmark simplification: it provides a stable reference configuration and removes rigid-body modes, but omits cavity pressure and pericardial restraint.

For a more physiological pumping simulation, replace the zero endocardial tractions with LV/RV cavity pressures and replace the fixed base with a basal spring constraint; add normal, sliding pericardial restraint at `EPI`.

## Parameter scope

The isotropic passive law (`E=100 kPa`, `nu=0.3`) is a first-integration-test approximation. `monodomainSolverCoeffs` sets `conductivitySource field;`, reading the anisotropic conductivity tensor `Conductivity` (`0/electro/Conductivity`, a `volSymmTensorField`) rather than a uniform value. Active mechanics use the rule-based fibre field `f0`; the passive mechanics remain an isotropic neo-Hookean approximation. Review and calibrate passive, active, cavity-pressure, and pericardial parameters before using this case for physiological predictions.

## Run prerequisites

This needs cardiacFoam compiled with solids4foam, including `libelectroMechanicalModels`. From this directory, run `./Allrun` (or `./Allrun parallel`). The mesh and the `fiber`/`sheet`/`tm`/`Conductivity` fields are not committed in this case directory — they live in `../mesh/` and `Allrun` copies them in each run (the `fiber` field becomes this case's `f0`); `t` is derived from `tm`, not copied. `Allclean` removes those copies and the derived `t` along with the generated regional meshes and run output.

## Notes on dictionary entries

**Stimulus.** The apical stimulus box in `constant/electro/electroProperties`
contains 26 cell centres at the mesh apex, where the supplied `apicobasal`
coordinate is smallest. All values are SI: `Allrun` scales the mesh to
adult-heart dimensions. The ionic and mechanical states remain at rest for 2 ms
before the apical pacing pulse begins.

**Run length.** `endTime` covers complete ventricular activation and the
electromechanical response; the write interval is 5 ms.

**Active tension.** `0/solid/Ta` is the initial active tension. The
electromechanical model registers and solves `Ta`, but the solid model builds
its stiffness while it is constructed, before that field exists, so
`electroMechanicalLaw` reads `Ta` from this file. In
`constant/electroMechanicalProperties` the `Batched` variants of the
active-tension models (`LandNiedererTWorldBatched`, `NashPanfilovBatched`) use
GPU/OpenMP vectorisation with `batchedSubsteps` explicit substeps per
mechanics step and `batchedParallelCells` for CPU multi-threading. Batched
active tension always integrates with explicit Euler (its states are not stiff,
the fastest time constant being ~0.7 ms), so it has no `batchedIntegrator`
key; that key belongs to the ionic model only. `preconditioningTime` (default
1000 ms, 0 disables) overrides the pre-conditioning time; the resting `Cai` is
queried from the ionic model's initial state. The `ODESolver` settings apply
to the scalar `LandNiedererTWorld` only.

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
