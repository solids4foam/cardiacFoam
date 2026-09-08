# Idealized biventricular electromechanics

This case adapts `NiedererEtAl2011/electroMechanicalNiedererEtAl2011` to the supplied idealized biventricular mesh. It uses cardiacFoam's sequential electromechanical model, the TNNP ionic model, and the Land--Niederer active-tension model. Mechanics uses the rule-based fibre field `f0`, shared with `electrophysiologyHeart` via `../mesh/0/fiber`.

## Mesh and activation

The shared mesh in `../mesh/` is already scaled uniformly by `0.02` from the raw geometry, giving an approximately `80 x 100 x 80 mm`, `226 mL` idealized adult-heart-scale mesh in SI metres. `Allrun` copies it in as-is — no scaling happens here. The apicobasal coordinate has its minimum at the anatomical apex, even though Cartesian x increases toward it. The external stimulus starts at `2 ms` and occupies a verified 26-cell apical box:

```text
x = [0.0756, 0.0800] m
y = [-0.0040, 0.0020] m
z = [0.0030, 0.0060] m
```

`monodomainSolverCoeffs.ionicHeterogeneity` classifies cells into `endocardialCells`/`mCells`/`epicardialCells` bands (`endoMInterface 0.3`, `mEpiInterface 0.7`) from a field that must be 0 at the endocardium and 1 at the epicardium — the opposite of this mesh's `tm` convention. `system/electro/setExprFieldsDict` derives `t = 1 - tm` for it to read; `Allrun` runs `setExprFields -region electro` for this.

## Mechanical boundary conditions

`BASE` has zero displacement. `EPI`, `ENDO_LV`, and `ENDO_RV` are traction-free in this first, unloaded activation/contraction case. This is an intentional benchmark simplification, not an in-vivo representation: it provides a stable reference configuration and removes rigid-body modes, but omits cavity pressure and pericardial restraint.

For a more physiological pumping simulation, replace the zero endocardial tractions with LV/RV cavity pressures and replace the fixed base with a basal spring constraint; add normal, sliding pericardial restraint at `EPI`.

## Parameter scope

The inherited isotropic passive law (`E=100 kPa`, `nu=0.3`) is useful for a first integration test, but is not a calibrated human biventricular model. `monodomainSolverCoeffs` sets `conductivitySource field;`, reading the anisotropic conductivity tensor `Conductivity` (`0/electro/Conductivity`, a `volSymmTensorField`) rather than a uniform value. Active mechanics use the rule-based fibre field `f0`; the passive mechanics remain an isotropic neo-Hookean approximation. Review and calibrate passive, active, cavity-pressure, and pericardial parameters before using this case for physiological predictions.

## Run prerequisites

This needs cardiacFoam compiled with solids4foam, including `libelectroMechanicalModels`. From this directory, run `./Allrun` (or `./Allrun parallel`). The mesh and the shared `fiber`/`sheet`/`tm`/`Conductivity` fields are not committed in this case directory — they live once in the sibling `../mesh/` directory (shared with `electrophysiologyHeart`, to avoid tracking the same ~10MB mesh twice in git) and `Allrun` copies them in each run (the shared `fiber` field becomes this case's `f0`); `t` is derived from `tm`, not copied. `Allclean` removes those copies and the derived `t` along with the generated regional meshes and run output.
