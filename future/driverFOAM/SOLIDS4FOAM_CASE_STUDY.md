# solids4foam External-Plugin Case Study

## Why solids4foam is the right next project

solids4foam is a strong external-validity case because it is close enough to
exercise driverFOAM's intended OpenFOAM abstractions and different enough to
expose cardiac coupling:

- it is an OpenFOAM toolbox for solid mechanics and fluid–solid interaction;
- it uses the standard OpenFOAM case layout and `Allrun`/`Allclean` conventions;
- it uses a single `solids4Foam` executable with runtime-selected physics and
  model dictionaries;
- its project-specific dictionaries include `physicsProperties`,
  `solidProperties`, `fluidProperties`, and `fsiProperties` rather than
  cardiacFoam's `electroProperties` vocabulary;
- it spans solids, fluids, FSI, thermo-mechanical cases, nonlinear materials,
  contact, and multiple OpenFOAM forks;
- it already has smoke tests and reference-based regression tests;
- it is mature research software with a 2018 methods paper and a 2025 JOSS
  paper [R12, R13].

The local sibling checkout inspected on 2026-08-15 was
`../solids4foam`, branch `feature-electromechanical-fix`, commit `b956d9bc`,
tracking the official `solids4foam/solids4foam` repository. It contained:

- 65 tutorial `Allrun` scripts;
- 59 `physicsProperties` files;
- 54 `solidProperties` files;
- 13 `fluidProperties` files;
- 10 `fsiProperties` files;
- 11 cases in `tutorials/Alltest-regression`.

Official solids4foam documentation distinguishes smoke tests—which only check
that a case advances—from regression tests that compare selected predictions
against stored references [R14]. This distinction should be retained in the
driverFOAM study.

## Claim tested

> A separately packaged solids4foam plugin can plan, execute, sweep, resume,
> audit, and collect artifacts from representative cases using a released,
> unmodified driverFOAM core.

This supports **OpenFOAM project generality**. It does not support simulator
generality beyond OpenFOAM.

## Independence requirements

The study is invalid as an external-plugin demonstration if the plugin is added
inside `openfoam_driver/plugins/` or if driverFOAM core is edited until the cases
pass.

Required structure:

```text
driverfoam-solids4foam-plugin/       # separate repository/distribution
├── pyproject.toml
├── src/driverfoam_solids4foam/
│   ├── plugin.py
│   ├── plugin.yaml
│   ├── dictionaries.py
│   ├── case_introspection.py
│   ├── validation.py
│   ├── artifacts.py
│   ├── sweep.py
│   └── tutorials.py
├── tests/
└── README.md
```

The distribution registers an entry point in the
`driverfoam.plugins` group. Its CI must install a tagged driverFOAM release, not
an editable core checkout containing study-specific patches.

## Proposed plugin responsibilities

| Capability | solids4foam implementation |
|---|---|
| Identity/profile | Stable plugin ID and profile describing supported case files and C++ source roots |
| Case recognition | `constant/physicsProperties` plus `system/controlDict` application or solids4foam-specific dictionaries |
| Solver commands | `solids4Foam`; explicitly list helper utilities only when needed |
| Dictionary catalog | `physicsProperties`, `solidProperties`, `fluidProperties`, `fsiProperties`, plus relevant control/system entries |
| Case introspection | Resolve analysis type (`solid`, `fluid`, `fluidSolidInteraction`) and selected runtime models |
| Semantic validation | Require the matching model dictionary; validate known model/algorithm combinations without pretending to prove mechanics correctness |
| Workflow | Ingest or model `Allrun` steps such as mesh generation, format conversion, solver run, and regression/post-processing |
| Artifacts | Displacement/point-displacement fields, stress/strain outputs, forces, interface histories, logs, and regression summaries as appropriate |
| Provenance | Case dictionaries, custom case source, `solids4Foam` executable, solids4FoamModels library, optional PETSc/preCICE dependencies |
| Sweep materialization | Material parameters, load values, mesh resolution, time-step/tolerance, and selected algorithm/model |

## RunDocument issue that this case should expose

RunDocument v2 currently requires the cardiac-shaped sections `anatomy`,
`physics`, `stimulus`, and `solver`. A solids4foam plugin could temporarily place
values into that shape, but doing so would demonstrate compatibility—not a clean
project-neutral contract.

Use the study in two stages:

1. **Boundary probe:** implement the plugin without core changes and document
   every place where RunDocument v2 or a compatibility fallback is unnatural.
2. **Generic-contract acceptance:** after a versioned generic RunDocument
   envelope/plugin-owned config schema exists, rerun the same plugin unchanged
   except for its declared schema migration.

No non-cardiac explicit v2 plugin should call a cardiac compatibility fallback.
Instrument this as a test.

## Selected case ladder

Start small and increase scientific/architectural difficulty.

| Stage | Candidate | Purpose | Required comparison |
|---|---|---|---|
| S1 | `solids/linearElasticity/patchTest` | Minimal solid case and exact/near-exact patch behavior | Native `Allrun`, solver completion, displacement/stress regression |
| S2 | `solids/linearElasticity/plateHole` | Established linear-elastic regression case | Existing `regressionTest.sh` criteria plus driver artifact manifest |
| S3 | `solids/hyperelasticity/rigidRotation/rotatingSphere` | Nonlinear constitutive/large-motion path | Existing regression values and resume after injected interruption |
| S4 | `solids/thermoelasticity/slabCooling` | Coupled thermal-solid dictionaries and fields | Existing regression plus a small parameter sweep |
| S5 | `fluidSolidInteraction/3dTube` or `HronTurekFsi3` | Multi-region/coupled workflow | Native regression/benchmark observables and multi-step provenance |

The first release need not support all 65 tutorials. It must clearly declare its
supported subset and reject or classify unsupported cases honestly.

## Test matrix

| Test | Native condition | driverFOAM condition | Success measure |
|---|---|---|---|
| Discovery | Locate/categorize case manually | `foamctl describe --plugin solids4foam` | Correct case kind, dictionaries, workflow, and model metadata |
| Plan | Inspect case and run scripts | strict plan | No cardiac fields; correct required files/commands/artifacts |
| Execute | Canonical case `Allrun` | normalized driver workflow | Same executable inputs and regression-pass result |
| Sweep | Hand shell loop | driver sweep | Same case matrix and values; complete manifest |
| Resume | Restart manually after process kill | driver state/resume | Correct restart point; no stale outputs accepted |
| Provenance | Manual notes/git SHA | driver snapshot | Records case, plugin, solver/library, workflow, and environment identity |
| Portability | Native supported OpenFOAM variants | plugin CI matrix | Declared versions pass without conditional core edits |

## Metrics

- driverFOAM core files changed: target **0**;
- plugin source/test lines and person-hours;
- supported case fraction and explicit unsupported classifications;
- plan/run/sweep/resume success rate;
- native-versus-driver numerical/regression agreement;
- planning and execution overhead;
- defect-detection precision/recall on solids4foam-specific mutations;
- provenance fields captured/missed;
- compatibility-fallback calls: target **0** for an explicit v2 plugin;
- OpenFOAM versions/distributions passing CI.

## solids4foam-specific fault injections

- `physicsProperties` type disagrees with available model dictionary;
- missing `solidProperties`, `fluidProperties`, or `fsiProperties`;
- unknown `solidModel`, `fluidModel`, interface/coupling scheme, or material law;
- missing region-specific `fvSolution`/`fvSchemes`;
- missing case-compiled boundary-condition library;
- unavailable PETSc, preCICE, or solids4foam runtime library;
- `Allrun` calls an undeclared helper command;
- solver finishes but required displacement/force/interface artifact is absent;
- solver or `libsolids4FoamModels` changes before resume;
- unsupported OpenFOAM fork/version is selected.

## Scientific restraint

The plugin can validate that a selected model exists, required dictionaries are
present, and an existing regression criterion passes. It cannot infer that a
mesh is adequate, a constitutive law is physically appropriate, contact is
well-posed, or an FSI solution is converged unless explicit domain validators
and benchmark criteria encode those properties.

## Expected publication value

A successful solids4foam study would be more persuasive than adding another
cardiac tutorial because it tests:

- a different scientific domain;
- different project dictionaries and artifacts;
- runtime-selected solid/fluid/FSI models;
- nonlinear and coupled workflows;
- external packaging and entry-point discovery;
- the difference between generic OpenFOAM mechanisms and residual cardiac
  assumptions.

References are listed in [`REFERENCES.md`](REFERENCES.md).

