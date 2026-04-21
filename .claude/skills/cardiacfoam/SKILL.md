---
name: cardiacFoam
description: >
  Use this skill for any task involving the cardiacFoam repository: adding new
  ionic models, editing electro solvers, modifying build files, writing patches,
  following OpenFOAM runtime selection table patterns, or maintaining coding style.
  Trigger whenever the user mentions cardiacFoam, electroModel, ionicModel,
  monoDomain, eikonal, ODE ionic models, electromechanical coupling, wmake, or
  requests C++ changes in this codebase.
---

# Codex Guidelines for cardiacFoam

This document defines how automated coding changes should be made in this repository.

## 1) Core Principles

- Make the smallest correct change.
- Preserve existing architecture and naming patterns.
- Prefer consistency with nearby code over introducing new style variants.
- Avoid broad refactors unless explicitly requested.
- Keep compatibility across supported OpenFOAM versions in mind (v2312–v2512).

## 2) Coding Style Rules

### C++ style

- Follow existing OpenFOAM/cardiacFoam style in surrounding files.
- Use the same indentation, brace style, comment style, and naming conventions as local code.
- Keep lines and expressions readable; avoid clever/condensed code.
- Prefer explicit, local, maintainable changes over abstraction-heavy rewrites.
- Add comments only when behavior is non-obvious; do not add redundant comments.

### File/header conventions

- Preserve existing license/header block format in C++ files.
- Keep include ordering consistent with nearby files.
- Do not change copyright headers unless explicitly asked.

### Ionic model file conventions

Each ionic model follows a three-file pattern:
- `{Model}.H/.C` — OpenFOAM wrapper class (inherits `ionicModel`)
- `{Model}_{year}.H` — Self-contained ODE equations (no OpenFOAM headers)
- `{Model}_{year}Names.H` — Variable name strings for I/O

Preserve this separation when adding or editing ionic models.

### Scripts and docs

- Match existing shell script style in `Allwmake`, `Allrun`, `Allclean`, etc.
- Keep Markdown concise, practical, and repository-specific.
- All Markdown files must pass the lint checker, e.g., `markdownlint README.md`.

## 3) OpenFOAM Conventions to Follow

- Respect OpenFOAM runtime type conventions:
  - `TypeName("...")` in class headers.
  - Registration via `addToRunTimeSelectionTable(...)` in source files.
- Preserve dictionary-driven behavior and runtime configurability.
- Keep OpenFOAM version compatibility intact (v2312–v2512).
- Avoid introducing dependencies that break existing wmake workflows.

## 4) Runtime Selection Tables (How They Work)

cardiacFoam relies on OpenFOAM runtime selection tables to instantiate models
from dictionaries at runtime. Two main selection hierarchies exist:

### electroModel hierarchy

1. `electroModel` (base) declares the selection table.
2. Derived class provides `TypeName("...")`.
3. Derived class registers with `addToRunTimeSelectionTable(electroModel, Derived, dictionary)`.
4. User selects via `electroModel <typeName>;` in `constant/cardiacProperties`.

Registered types:
- `monoDomainElectro` — full monodomain reaction-diffusion PDE
- `eikonalDiffusionElectro` — fast activation-time eikonal solver
- `singleCellElectro` — 0D single-cell ODE solver

### ionicModel hierarchy

1. `ionicModel` (base, also inherits `ODESystem`) declares the selection table.
2. Derived class provides `TypeName("...")`.
3. Derived class registers with `addToRunTimeSelectionTable(ionicModel, Derived, dictionary)`.
4. User selects via `ionicModel <typeName>;` in the case dictionary.

Registered types (12 models):
`BuenoOrovio`, `TNNP`, `Courtemanche`, `Grandi`, `ORd`, `Stewart`,
`Trovato`, `ToRORd_dynCl`, `Fabbri`, `Gaur`, `AlievPanfilov`, `tmanufacturedFDA`

Rules when adding new runtime-selectable classes:
- Add `TypeName("...")` in header.
- Add `addToRunTimeSelectionTable(...)` in source.
- Add source file to `Make/files` in the relevant library.
- Ensure dictionary `type` string exactly matches `TypeName` value.

## 5) Dependency and Build Architecture

cardiacFoam has two operating modes, resolved automatically at build time by
`etc/resolveSolids4Foam.sh`:

- **Full mode** (with solids4foam): electrophysiology + solid mechanics + FSI.
  Uses `modules/solids4foam/` submodule or an external `SOLIDS4FOAM_INST_DIR`.
- **Lightweight mode** (EP-only): uses `modules/physicsModel/` fallback.
  Activated when solids4foam is not found.

Build sequence:
1. `Allwmake` calls `etc/resolveSolids4Foam.sh` to detect mode.
2. `src/Allwmake`: compiles `genericWriter` → `ionicModels` → `electroModels`.
3. `applications/Allwmake`: compiles `cardiacFoam` solver, then utilities.

When touching build files:
- Library source lists live in `src/{library}/Make/files`.
- Solver/utility source lists live in `applications/{type}/{name}/Make/files`.
- Do not add compile-time `#ifdef` blocks unless the existing code already uses them.

## 6) Key Interfaces to Implement

### Adding a new ionicModel

Derived classes must implement:
- `solveODE(...)` — advance ionic state variables with the ODE solver
- `nEqns()` — number of ODEs in the system
- `derivatives(...)` — RHS of ODE system (ODESystem interface)
- `calculateCurrent(...)` — return net ionic current
- `signal(label, CouplingSignal)` — return `Vm`, `Act`, or `Cai` for electromechanical coupling

### Adding a new electroModel

Derived classes must implement:
- `evolve()` — main PDE solver loop per timestep
- `read()` — read from `electroProperties`
- `writeFields(const Time&)` — output fields
- `end()` — cleanup

## 7) Minimise Changes

- Only modify files necessary for the requested task.
- Do not reformat unrelated code.
- Do not rename symbols/files unless required.
- Do not alter behavior outside requested scope.
- Prefer targeted edits over cleanup passes.

Before finalizing, verify:
- Build impact is localized.
- No unrelated files changed.
- No compatibility regressions introduced by style-only edits.

## 8) Change Delivery Format (Mandatory)

All code changes must be returned as **git patches**.

- Provide changes in patch form (`git diff`/unified diff) suitable for application.
- Keep patches focused and reviewable.
- If multiple concerns are required, split into logical commits/patches.
- Do not provide only prose summaries when code edits were requested.

## 9) Practical Workflow for Codex

1. Read nearby code and follow local patterns.
2. Implement minimal patch.
3. Update runtime registration and build lists if adding a new class.
4. Run or describe relevant checks/tests.
5. Return changes as patch-oriented output.

## 10) What to Avoid

- Large-scale refactors without explicit request.
- API redesigns when a local fix is sufficient.
- Introducing new style conventions inconsistent with repository norms.
- Silent behavioral changes not documented in the patch summary.
- Adding OpenFOAM-version `#ifdef` guards unless existing code already uses them.

## 11) Reference Files (Commonly Relevant)

- Main solver: `applications/solvers/cardiacFoam/cardiacFoam.C`
- Core model base classes:
  - `src/electroModels/electroModel/electroModel.H`
  - `src/ionicModels/ionicModel/ionicModel.H`
  - `src/couplingModels/electromechanicsFeedbackProvider.H`
- Example ionic model implementation:
  - `src/ionicModels/BuenoOrovio/BuenoOrovio.H`
  - `src/ionicModels/BuenoOrovio/BuenoOrovio.C`
  - `src/ionicModels/BuenoOrovio/BuenoOrovio_2008.H`
  - `src/ionicModels/BuenoOrovio/BuenoOrovio_2008Names.H`
- Example electro model implementation:
  - `src/electroModels/monoDomainElectro/monoDomainElectro.H`
  - `src/electroModels/monoDomainElectro/monoDomainElectro.C`
- Build lists:
  - `src/electroModels/Make/files`
  - `src/ionicModels/Make/files`
  - `src/genericWriter/Make/files`
- Build scripts: `Allwmake`, `src/Allwmake`, `applications/Allwmake`
- Dependency resolver: `etc/resolveSolids4Foam.sh`
- Benchmark tutorial: `tutorials/NiedererEtAl2011/`
