# Repository Archaeologist: convention and ownership map

Audit basis: branch `no-frontend-minor-errors`, commit `f6b798807e8b8c7b685766d6928c16dc5997e150`. This is a read-only discovery report; no production file was edited. The worktree was already substantially dirty, including a modified `modules/solids4foam` submodule and untracked research/tutorial material, so untracked files and local modifications were not used as sole authorities. Architecture was sampled across the solver, every `src/` library family, utilities, both Python systems, build/resolver scripts, CI, dictionaries, tutorials, and component documentation.

## Canonical architecture and ownership

| Area | Owner / boundary | Canonical evidence |
|---|---|---|
| Top-level execution | `cardiacFoam` owns the time loop; a runtime-selected `physicsModel` owns physics lifecycle. | `applications/solvers/cardiacFoam/cardiacFoam.C:39-80` |
| Build-mode boundary | `etc/resolveSolids4Foam.sh` selects full external solids4foam or the project-owned lightweight fallback; electromechanics is full-mode only. | `etc/resolveSolids4Foam.sh:12-25,27-49,65-95`; `src/Allwmake:22-31` |
| Electro orchestration | `core/` orchestrates; `electroDomains/` owns long-lived domain state; myocardium/conduction/ECG directories own kernels; `electroCouplers/` owns staged coupling. | `src/electroModels/README.md:32-101` |
| Public electro selection | `constant/electroProperties` key `myocardiumSolver`; spatial names map to `electrophysiologyModel`, while `singleCellSolver` is a direct `electroModel`. | `src/electroModels/core/electroModel.C:147-180`; `src/electroModels/core/electrophysiologyModel/electrophysiologyModel.C:29-57` |
| Ionic layer | `ionicModel/` owns selection, ODE/metadata/I/O and heterogeneity contracts; model directories own wrappers and equation/metadata families. | `src/ionicModels/README.md:51-113`; `src/ionicModels/Make/files:1-40` |
| Active tension | `activeTensionModel` owns selection and common point loop; concrete and batched directories own kernels. | `src/activeTensionModels/README.md:26-46`; `src/activeTensionModels/NashPanfilov/NashPanfilov.H:46-145` |
| Electromechanics | Full-mode `electroMechanicalModel` composes electro and solid models; coupling policy belongs here, not in numerical kernels. | `src/electroMechanicalModels/README.md:16-47`; `src/electroMechanicalModels/electroMechanicalModel/electroMechanicalModel.H:55-78` |
| Lightweight compatibility | `modules/physicsModel` is project-owned fallback API compatibility code, deliberately minimal; it is not solids4foam itself. | `modules/physicsModel/README.md:1-40`; `modules/physicsModel/src/solids4FoamModels/physicsModel/physicsModel.H:51-81` |
| Python automation | Packaged `openfoam_driver` owns typed planning/runtime/spec/catalog/schema/test contracts; tutorial setup JSON is an input/output contract. | `applications/scripts/driverFoam/pyproject.toml:5-58`; `applications/scripts/driverFoam/openfoam_driver/core/runtime/models.py:36-150` |
| Model generation | `cellML2foam` owns transformation/templates; generated equations still require documented physiological review. | `applications/scripts/cellML2foam/cellML2Foam.py:59-90`; `applications/scripts/cellML2foam/src/pipeline.py:43-70,73-126` |
| Tutorials | Cases own `Allrun`/`Allclean`, local numerical regressions and reference inputs; driver specs/defaults encode reusable orchestration. | `tutorials/README.md:32-56`; `.agents/skills/cardiacfoam/PROJECT_MEMORY.md:121-133` |

## Canonical conventions

### C++ and OpenFOAM

- Hand-maintained C++ uses the cardiacFoam GPL banner, `Class`/`Description`/`SourceFiles` metadata in headers, include guards, namespace `Foam`, OpenFOAM section rulers, and paired `.H`/`.C` files. Representative references: `applications/solvers/cardiacFoam/cardiacFoam.C:1-39`, `src/ionicModels/BuenoOrovio/BuenoOrovio.H:1-50`.
- Names are OpenFOAM idiomatic: lower-camel class names for framework concepts (`electroModel`, `myocardiumSolver`), model/public runtime names preserve scientific proper names (`BuenoOrovio`), data members end in `_`, and accessors are lower camel case. Runtime names, selector strings, and `<type>Coeffs` dictionaries are compatibility surfaces.
- Runtime-selectable bases declare/define a table and derived types register in the correct table. Ordinary derived models use local `OverrideTypeName("...")` plus `defineTypeNameAndDebug` and `addToRunTimeSelectionTable`; aliases use `addNamedToRunTimeSelectionTable`. Canonical examples: `src/ionicModels/BuenoOrovio/BuenoOrovio.H:140-143`, `.C:35-42`; `src/electroModels/core/electrophysiologyModel/electrophysiologyModel.C:29-57`.
- Ownership follows OpenFOAM types (`autoPtr`, `tmp`, `PtrList`) and factories return `autoPtr`. Virtual base destructors are explicit. Copy prevention varies historically between private declarations and `= delete`; follow the nearest maintained sibling rather than modernizing a whole family.
- Formatting authority is local family style, not a global formatter. Newer code commonly uses four-space continuation/body indentation (`BuenoOrovio.C:84-121`), while older OpenFOAM-derived headers use deeply nested section indentation (`electroMechanicalModel.H:60-150`). Preserve neighboring layout and avoid cross-family reformatting.
- Public configuration is dictionary-driven, with OpenFOAM `FoamFile` headers, semicolon-terminated entries, aligned values, dimensional quantities, and comments that state SI units. Canonical template: `tutorials/template/constant/electroProperties:8-16,41-73`; time units: `tutorials/template/system/controlDict:18-47`.
- Build ownership is explicit: every compiled source belongs in its component `Make/files`, dependencies in `Make/options`; optional CUDA source names come from `Make/files-gpu` through `src/ionicModels/Make/files:37-40`.

### Python

- `openfoam_driver` is the canonical maintained Python package: Python `>=3.10`, absolute standard-library imports before relative package imports, `from __future__ import annotations`, `pathlib.Path`, frozen dataclasses for durable value contracts, built-in generic typing, explicit return annotations, JSON-compatible shapes, and exceptions/diagnostics at library boundaries. Evidence: `pyproject.toml:5-16`; `openfoam_driver/cli.py:28-67`; `core/runtime/models.py:28-49,94-150`.
- CLI entry points return integer status and emit structured JSON where machine consumption is intended; workflow success derives from explicit state, not merely process exit code (`openfoam_driver/cli.py:70-126`). Tests live beside the package under `openfoam_driver/tests/`, and CI installs the package and runs that suite (`.github/workflows/driverFoamTests.yml:41-57`).
- `cellML2foam` is a separate, older script-style subsystem: executable entry script, local `src` path injection, loose typing, direct printing/exits, staged subprocess pipeline (`cellML2Foam.py:1,30-49,95-146`; `src/pipeline.py:43-70,86-125`). Its style is not canonical for new `openfoam_driver` modules.

### Shell, workflows, config, and docs

- Build scripts are Bash, relocate to their own directory, fail explicitly, quote environment-derived paths in maintained resolver code, and source OpenFOAM argument helpers when present (`Allwmake:1-40`; `etc/resolveSolids4Foam.sh:12-18,44-78`). Case scripts conventionally use `Allrun`/`Allclean`; local numerical checks use `regressionTest.sh` or `runRegressionTest.sh` (`tutorials/README.md:32-39`).
- CI uses pinned major action versions, matrixed supported OpenFOAM containers, and exercises both full and lightweight builds (`.github/workflows/buildAndTest.yml:19-44,50-87`). Python contract tests are a separate source-triggered gate (`.github/workflows/driverFoamTests.yml:15-57`).
- JSON schemas/manifests and dictionary keys are public contracts, not formatting targets. The package data declaration establishes the packaged schema and allowlist sources (`applications/scripts/driverFoam/pyproject.toml:48-54`).
- Documentation is layered: root README explains capability/layout/build; component README explains ownership/runtime names; architecture docs explain detailed interfaces; tutorial README explains runnable workflows. Code and build files outrank prose if they disagree (`.agents/skills/cardiacfoam/PROJECT_MEMORY.md:5-18`). Markdown lint is deliberately permissive (`.markdownlint.yaml:1-5`), so concise local consistency is the relevant subjective standard.

## Generated, external, and hand-maintained map

| Classification | Paths / rule | Treatment |
|---|---|---|
| Hand-maintained project code | Most `.C/.H` wrappers and bases in `src/`; solver/utilities; `openfoam_driver`; build/CI scripts; READMEs; tutorial dictionaries/scripts. | Edit only at owning component; preserve nearby convention and public keys. |
| Generator/template code | `applications/scripts/cellML2foam/src/derivedClass_templates/`; transformation/mapping pipeline under `applications/scripts/cellML2foam/src/`. | Prefer changing templates/generator contracts and regenerating; validate physiological semantics. |
| Generated or generator-derived ionic code | Model equation headers such as `src/ionicModels/BuenoOrovio/BuenoOrovio_2008.H` (explicitly “Converted from CellML” at lines 18-20), `*_<year>Names.H`, and batched `*_<year>Batch.H`; some have acquired hand-maintained dependency/metadata additions. | Do not normalize en masse. Establish provenance per model before edits; equation changes require scientific equivalence tests. Wrappers `{Model}.H/.C` are hand-maintained. |
| Manufactured/reference code | `src/ionicModels/verificationModels/`, `src/verificationModels/`, tutorial case inputs/reference results. | Project-owned but correctness-sensitive; do not infer formula equality from style/family similarity. |
| External/submodule | `modules/solids4foam` is a gitlink (recorded at `0bd882172db292c29bf41c4233d61cfa5f116168`, locally modified). | Never normalize or patch during this audit; only project-owned interfaces may change with approval. |
| Project-owned fallback | `modules/physicsModel/` (tracked ordinary files, despite `modules/` location). | Treat as compatibility-sensitive hand-maintained code and compare against both supported modes. |
| Tutorial/reference input | `tutorials/**/constant`, `system`, `0`, graphs/meshes, reference outputs; setup JSON/manifests. | Public examples and regression inputs; numerical/config changes are behavioral, not cleanup. |

## Likely style-drift clusters (investigation targets, not confirmed defects)

1. **C++ layout generations:** `electroModel.H:45-95` uses compact two-space class indentation while `electroMechanicalModel.H:55-100` and fallback `physicsModel.H:51-92` retain older deeply nested OpenFOAM layout. This is visible drift but is subjective unless mixed within one edited family; do not globally reformat.
2. **Copy-control generations:** `electroModel.H:77-81` uses undeﬁned private declarations, while `BuenoOrovio.H:74-78` and `NashPanfilov.H:51-52` use `= delete`. Both are established. Maintainer judgment is needed on minimum supported compiler/OpenFOAM constraints before any convergence.
3. **Generated-wrapper families:** scalar `{Model}.H/.C`, generated equation/Names headers, batched wrappers, Batch headers, and optional CUDA sources are repeated across 12 model pairs (`src/ionicModels/Make/files:11-38`). Header guards already reflect historical names in at least `BuenoOrovioBatched.H:28-29` (`BuenoOrovioGPU_H`), making this a high-yield drift scan, but direct mass repair is unjustified without generator/provenance analysis.
4. **Python maturity split:** packaged/typed `openfoam_driver` versus script-oriented `cellML2foam`. Examples include typed `Path`/dataclass contracts (`core/runtime/models.py:49-99`) versus runtime `sys.path` mutation and untyped CLI functions (`cellML2Foam.py:34-49`). Treat as subsystem boundaries; only defects or shared-contract drift justify convergence.
5. **Shell strictness/quoting:** top-level `Allwmake:2,37-40` mixes unquoted script-directory expansion and pipelines without an explicit `pipefail`; the newer resolver quotes path construction (`resolveSolids4Foam.sh:16-18`). Equivalent build/case scripts should be scanned for failure-propagation divergence, not cosmetically rewritten.
6. **Documentation list duplication:** compiled model/runtime/tutorial lists recur in `README.md:37-60`, `src/ionicModels/README.md:114-152`, `IONIC_MODEL_ARCHITECTURE.md:9-22`, `Make/files`, Python catalogues, and tutorial docs. Existing driver drift guards demonstrate this has caused real drift (`.github/workflows/driverFoamTests.yml:1-13`). Prefer generated/contract-checked lists over manual synchronized edits.
7. **Tutorial naming/history:** root README describes group `singleCellprotocols/` (`README.md:51-58`), whereas the tracked layout and tutorial architecture use `electrophysiologyProtocols/` (`tutorials/README.md:7-15`). This is a concrete documentation-drift candidate for Agent 9, not a style convention.
8. **OpenFOAM banner versions in configuration:** template `electroProperties` says v2412 (`tutorials/template/constant/electroProperties:1-6`) while template `controlDict` says v1912 (`tutorials/template/system/controlDict:1-6`). Banner version may be provenance rather than supported-version contract; maintainer judgment is required before normalization.

## Ambiguous conventions requiring maintainer judgment

- Whether `= delete`, `override`, `final`, `nullptr`, and initializer-list returns are mandatory for new C++ or merely accepted in newer families; local code supports multiple generations and compatibility spans OpenFOAM v2312-v2512 (`Allwmake:13-19`).
- Whether generated equation and Names headers are reproducible outputs or curated generated-derived sources. The CellML CLI explicitly requires manual semantic cleanup (`cellML2Foam.py:59-90`), so a blanket “never edit generated files” rule would be inaccurate.
- Whether Python license banners are required for every new module/test. They recur in production modules (`openfoam_driver/cli.py:1-26`) but tests and package metadata follow lighter conventions.
- Whether shell strict mode should be `set -e`, `set -eu`, or `set -euo pipefail`. Existing authoritative scripts establish only partial consistency.
- Whether tutorial OpenFOAM banner `Version` values should track the oldest supported, authoring, or current release; executable compatibility is controlled elsewhere.
- Whether JSON files under tutorial `setup/` are authored inputs or generated snapshots varies by filename (`driver_config`, `run_manifest`, `artifacts_manifest`); ownership should be documented per artifact before automated rewriting.

## Proposed machine-checkable rules

These rules encode observed contracts and should be scoped to hand-maintained files unless stated otherwise.

1. Parse all runtime registration macros and assert selector names exist in the correct `Make/files`, Python catalogues, dictionary allowlists, template comments/examples, and component model lists. Extend the existing drift-guard approach cited by `.github/workflows/driverFoamTests.yml:1-13`.
2. For every compiled scalar ionic model, verify wrapper `.H/.C`, runtime registration, generated/derived Names metadata, expected batched partner, and optional CUDA entry consistency; allow explicit exceptions for verification models.
3. Validate `FoamFile.object`, directory location, required selector key, matching `<runtimeName>Coeffs`, semicolons, and dimensional entry syntax in tracked canonical/tutorial dictionaries.
4. Require every project-owned compiled `.C`/`.cu` to appear exactly once in its owning `Make/files`/`Make/files-gpu`, and reject missing paths. Exclude the external submodule.
5. Compile/import all project-owned Python, run JSON Schema validation for tracked manifests/configs of declared schema type, and enforce deterministic catalogue-generation tests.
6. Shell-check project-owned build/case scripts for syntax, unsafe unquoted path expansion, and pipelines whose upstream failure is not propagated; use an allowlist for established OpenFOAM idioms.
7. Check that Markdown links and referenced paths exist and that root/component model/tutorial lists match executable catalogues. Apply the committed markdownlint configuration rather than a new style profile.
8. Check C++ include guards match a documented transformation of the filename for newly added hand-maintained headers; report legacy exceptions (such as `BuenoOrovioBatched.H:28-29`) without forcing bulk changes.

## Subjective/local-review rules

- Match indentation, brace placement, section comments, include ordering, and line wrapping in the nearest equivalent maintained sibling; never use one directory to reformat another family.
- Preserve OpenFOAM ownership/error idioms when they clearly communicate lifetime and dictionary context; do not modernize to standard-library constructs for appearance.
- Comments should explain units, contracts, numerical ordering, provenance, or non-obvious behavior; avoid comments that merely restate code.
- Keep component READMEs concise and architecture-specific, and update the nearest owning document when architecture or a public contract materially changes.
- Prefer the smallest patch at the source of truth. Formatting-only drift is S4 and should remain untouched unless adjacent to an approved functional repair.
- Treat any equation, unit, initialization, solver order, tolerance, field layout, or reference-data change as scientific behavior requiring an explicit hypothesis and numerical regression, regardless of how stylistic the diff appears.

## Archaeologist conclusion

The repository's stable convention is **family-local OpenFOAM style plus explicit runtime/build/config contracts**, not uniform formatting. The most valuable preventive controls are contract extraction across registrations/build/catalogues/docs, generated-family completeness checks, and shell failure-propagation checks. Global C++ or Python formatting would erase useful provenance and create high-risk noise without addressing the demonstrated drift surfaces.
