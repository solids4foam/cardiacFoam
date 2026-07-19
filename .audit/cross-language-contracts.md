# Cross-Language Contract Audit (Agent 5)

Date: 2026-07-11  
Scope: read-only discovery across C++, `Make/files`, Python, shell, OpenFOAM dictionaries, JSON Schema, tutorials, generated metadata, documentation, and full/lightweight build modes. No production files were edited.

## Executive result

Four confirmed cross-language inconsistencies were found: three S2 and one S3. The strongest existing contract checks are healthy: the targeted driver suite completed with **38 tests and 122 subtests passing**, including ionic C++-header/catalogue synchronization, catalogue audit, capability-manifest tests, Run-document schema/model parity, and runtime-selection enum guards. The passing tests also expose the remaining gaps: they check variable tuples and selected examples, but not the complete runtime-name/documentation set, the complete active-tension runtime table, semantic field categories, or build-mode availability.

Status vocabulary below follows the brief: **consistent**, **inconsistent**, **absent-required**, **intentionally absent**, and **uncertain**.

## Contract matrix

| Public/cross-language entity | C++ runtime / producer | Build | Python catalogue / registry | Dictionaries / schema | Tutorials / generated metadata | Documentation | Mode contract | Result |
|---|---|---|---|---|---|---|---|---|
| Top-level EP physics type `electroModel` | Registered in `src/electroModels/core/electroModel.C:34` | `src/Allwmake:26` | Dictionary catalogue describes it | `physicsProperties` examples use `type electroModel` | Present | Root/component READMEs present | Both modes | **consistent** |
| Myocardium selectors `monodomainSolver`, `bidomainSolver`, `eikonalSolver` | Named `electroModel` registrations in `electrophysiologyModel.C:35-57`; kernel registrations in sibling solver `.C` files | `libelectroModels` | `dict_entries.py:263-274`; RTST drift guard | `myocardiumSolver` key | Widely used | Architecture docs describe layered dispatch | Both modes | **consistent** |
| `singleCellSolver` | Direct `electroModel` registration at `singleCellSolver.C:41` | `libelectroModels` | Catalogue/defaults use exact name | `myocardiumSolver singleCellSolver` | Registered tutorial | Documented | Both modes | **consistent** |
| Scalar ionic models (12) | Exact runtime registrations and `OverrideTypeName` values | `src/ionicModels/Make/files:11-22` | Exact keys in `IONIC_MODEL_CATALOG` | Used as `ionicModel` enum values | Names metadata parsed from `*_Names.H` | `src/ionicModels/README.md:120-131` | Both modes | **consistent** |
| Manufactured ionic models (3) | Registered under exact names | `src/ionicModels/Make/files:6-8` | Present, deliberately semantic variable metadata | Used by MMS tutorials | Raw enum sync intentionally excluded | Listed at README lines 148-152 | Both modes | **consistent**, with intentional generated-header exclusion |
| Batched ionic runtime models (12) | Only `*compactBatched` subclasses are registered; e.g. `TNNPBatched.C:106-113` | Wrapper sources at `Make/files:24-35`; CUDA sources conditional | `BATCHED_MODELS` uses `*compactBatched` at `ionic_model_catalog.py:485-498` | Tutorials use `TWorldcompactBatched` | Catalogue derives pair metadata | README advertises unregistered `*Batched` names at lines 133-146 and 196-202 | CPU fallback exists; CUDA kernels conditional | **inconsistent** (F1) |
| Active-tension production models (3 scalar + 3 batched) | All six registered | `src/activeTensionModels/Make/files:7-15` | Six exact catalogue keys | `dict_entries.py:1413` enum | Used by single-cell/EM cases | Listed in component README | Library built in both modes; EM consumer full-only | **consistent** |
| Active-tension MMS model `ManufacturedElectromechanics` | Registered at `ManufacturedElectromechanics.C:32-39` | `Make/files:9` | Missing from `ACTIVE_TENSION_MODEL_CATALOG` | Present in enum at `dict_entries.py:1413` | Used at tutorial `constant/electroMechanicalProperties:22` | Documented by MMS tutorial and component tree | Full EM workflow only | **absent-required** (F2) |
| Ionic state/algebraic/constant names | Generated `*_Names.H` is source | Compiled with model | Catalogue exact tuples | Export selections consume names | Regeneration parser + audit tests | Catalogue claims exactness | Both modes | **consistent** for non-manufactured models; manufactured absence intentional |
| Ionic recommended exports | C++ names define admissible tokens | N/A | Hand-curated; legacy exceptions explicitly allowed | Used as planning fallback | Tests enforce minimum content | Documented indirectly | Both modes | **uncertain** for explicitly allowed legacy aliases; no new defect asserted |
| Active-tension state/algebraic/constant names | Model headers/producers | Compiled | Hand-copied catalogue | Capability manifest consumes them | No complete runtime/header parity test found | Documentation is model-level | Both modes for library | **uncertain**; F2 proves registry completeness is not guarded |
| Capability-manifest electro fields | Solvers/model I/O produce fields | N/A | Fixed names plus ionic metadata | Strict planner consumes as allowlist | Tests cover examples only | Described as field names | Both modes | **inconsistent**: species values are inserted (F3a) |
| Capability-manifest solid fields | Solid/EM solver produces `Ta`, `lambda` | EM wrapper omitted in lightweight at `src/Allwmake:29-31` | Spatial solver alone causes solid fields | Strict planner consumes as allowlist | Tests cover only active-tension-positive case | Manifest claims actual exposed fields | Not correctly represented | **inconsistent** (F3b) |
| Run document v2 | `RunDocument` Python model | Packaged data | Model/serializer present | Two identical `run-document.json` copies | Round-trip and equality test | README references v2 | Mode-neutral | **consistent**; duplicate unavoidable packaging copy is tested |
| Dictionary-key catalogue | C++ lookup scanner is approximate | N/A | `dict_entries.py` is public planning catalogue | JSON allowlist carries known scanner exceptions | Scanner/RTST tests exist | Driver README documents allowlist | Both modes | **uncertain**: `dict_key_allowlist.json:91-104` explicitly exempts 13 `stale_paths`; do not treat as defects without parser refinement |
| Registered tutorial/display metadata | Runtime registry and display catalogue cross-checked | N/A | `REGISTERED_TUTORIALS` + `TUTORIALS` | Run document accepts entries | Ten registered tutorials | Tutorials README | No availability flag | **inconsistent** for full-only EM entry (F4) |
| Full/lightweight build choice | Resolver exports `USE_LIGHTWEIGHT_PHYSICSMODEL`; `src/Allwmake` gates EM library | Authoritative shell behavior | Driver registries do not model selected mode | No Run-document capability field | Full-only tutorial still advertised | Root README states limitation | Diverges in Python metadata | **inconsistent** (F4) |

## Confirmed findings

### F1 — Documentation publishes 12 ionic selector names that are not registered (S2, high confidence)

**Evidence.** The public “Batched (SoA) models” list names `AlievPanfilovBatched` through `TWorldBatched` at `src/ionicModels/README.md:133-146`, and repeats those names as heterogeneity-capable selectors at lines 196-202. The current runtime pattern instead defines both wrapper and compact classes but registers only the compact class; for example `src/ionicModels/TNNPBatched/TNNPBatched.C:106-113` defines both type names and registers `TNNPcompactBatched`, while `TNNPBatched.H:146` and `:205` show the distinct strings. The executable Python catalogue correctly lists all 12 `*compactBatched` names at `applications/scripts/driverFoam/openfoam_driver/ionic_model_catalog.py:485-498`. The build file compiles wrapper source directories (`src/ionicModels/Make/files:24-35`), which does not make wrapper type names runtime-selectable.

**Violated contract/canonical example.** README model lists are presented as available dictionary selections; scalar entries at README lines 120-131 match their runtime and catalogue keys exactly. `TWorldcompactBatched` tutorial dictionaries also demonstrate the actual selector spelling.

**Impact.** A user copying any documented batched name into `ionicModel` reaches the OpenFOAM unknown-runtime-type error even though the relevant source was compiled. The heterogeneity section additionally assigns capabilities to names that cannot be selected.

**Minimal remediation.** Update the runtime-model list and heterogeneity examples to the 12 exact `*compactBatched` identifiers. Preserve directory/class implementation terminology separately if the wrapper classes need documentation. Do not add alias registrations without an explicit compatibility design.

**Required validation / automated contract test.** Add a static test that parses the `ionicModel` registration class from each compiled `.C`, resolves its `OverrideTypeName`, and compares the complete set with `IONIC_MODEL_CATALOG` and a machine-readable model list used to generate the README. At minimum, assert that every backticked item in the README’s “Available models” subsections is a registered runtime key.

### F2 — `ManufacturedElectromechanics` is executable and accepted by dictionaries but absent from active-tension introspection (S2, high confidence)

**Evidence.** `src/activeTensionModels/verificationModels/ManufacturedElectromechanics/ManufacturedElectromechanics.C:32-39` defines and registers `ManufacturedElectromechanics` in the `activeTensionModel` table, and `src/activeTensionModels/Make/files:9` compiles it. The driver’s dictionary enum includes it at `applications/scripts/driverFoam/openfoam_driver/dict_entries.py:1408-1414`. The MMS tutorial selects it at `tutorials/manufacturedSolutions/monodomainTotalLagrangianEM/constant/electroMechanicalProperties:22`. However, `ACTIVE_TENSION_MODEL_CATALOG` ends after six production entries at `active_tension_catalog.py:83-178`; introspection serializes only that mapping at `introspection.py:162-168`.

**Violated contract/canonical example.** All six other compiled/registered active-tension types have exact catalogue entries. The ionic catalogue likewise includes its three manufactured runtime models instead of hiding verification types.

**Impact.** `describe` and exported active-tension metadata falsely report an incomplete runtime surface. `build_capability_manifest()` cannot add the manufactured model’s state/algebraic outputs because lookup at `capability_manifest.py:130-137` returns no entry, so strict planning can warn on legitimate sampled MMS fields.

**Minimal remediation.** Add a catalogue entry sourced from `ManufacturedElectromechanics_2026Names.H`, classifying it as manufactured/verification metadata if a new flag is useful. Keep equations untouched.

**Required validation / automated contract test.** Parse compiled `activeTensionModel` registrations and assert exact key equality with `ACTIVE_TENSION_MODEL_CATALOG`, with a documented exclusion mechanism only if a runtime type is intentionally non-public. Add an introspection assertion for `ManufacturedElectromechanics` and its exported variables.

### F3 — Capability metadata mixes biological labels with fields and invents mechanics fields (S2, high confidence)

This is one public contract with two independent defects in the same producer.

**F3a evidence.** `applications/scripts/driverFoam/openfoam_driver/capability_manifest.py:101-107` promises `samplable_fields` are field names exposed by the solver. After adding states, algebraics, and exports at lines 116-119, line **120** also adds `ionic_entry.species`. Catalogue species values are labels such as `human`, `pig`, and `generic`, not C++ output fields.

**F3b evidence.** Lines **122-129** state that a solid region exists for every solver other than `singleCellSolver` and add `Ta` and `lambda`. Ordinary monodomain, bidomain, and eikonal EP cases do not thereby own a mechanics region. The build contract explicitly compiles `electroMechanicalModels` only when `USE_LIGHTWEIGHT_PHYSICSMODEL != 1` at `src/Allwmake:29-31`; the root README lines 7-8 and 70 state that mechanics workflows require full mode. No mode or actual electromechanical top-level selection enters `build_capability_manifest()`.

**Violated contract/canonical example.** `_ELECTRO_SOLVER_FIELDS` and `_SOLID_SOLVER_FIELDS` at lines 45-46 are field-name sets, while catalogue `species` is explicitly biological metadata. The existing single-cell test at `tests/test_capability_manifest.py:26-35` correctly expects no solid region; that same producer contract must be based on actual region/model selection, not merely “not single cell.”

**Impact.** Strict planning passes this manifest to sampled-field validation at `strict_planning.py:333-339`. It can therefore treat `human`/`pig`/`generic`, or `Ta`/`lambda` in a plain EP case, as legitimate sample targets instead of diagnosing them. The solver/function object may silently drop them, exactly the failure the manifest says it prevents.

**Minimal remediation.** Remove `electro.update(ionic_entry.species)`. Populate solid fields only when a solid/electromechanical region is positively resolved (or pass an explicit mode/physics capability into the function); an active-tension model in a single-cell trace should not automatically imply an OpenFOAM solid region without checking the producer contract.

**Required validation / automated contract test.** Assert that no catalogue `species` or `cardiac_region` label appears in `samplable_fields`; test plain monodomain/bidomain/eikonal cases with no EM selection produce an empty solid set; test a positively resolved full electromechanics case produces `Ta`/`lambda`; and test lightweight mode cannot advertise a solid region. Include a strict-planning negative test for sampling `human` and `lambda` from a plain EP case.

### F4 — Driver advertises a full-only tutorial without representing build-mode availability (S3, high confidence)

**Evidence.** The shell source of truth gates `electroMechanicalModels` at `src/Allwmake:29-31`, using the mode exported by `etc/resolveSolids4Foam.sh:27-48` and `:51-95`. Documentation explicitly says electromechanical tutorials cannot run electro-only (`README.md:70`; `tutorials/README.md:22`). Nevertheless `manufacturedMonodomainTotalLagrangianEM` is unconditionally present in `REGISTERED_TUTORIALS` at `core/runtime/registry.py:87-98` and in display metadata at `tutorials_display.py:158-170`. Neither entry contains `requires_full_mode`/required-library metadata, and the Run-document schema has no environment-capability representation for it.

**Violated contract/canonical example.** The registry/display layer claims to describe runnable entries, while the build resolver is the authoritative availability contract. The display entry already has semantic tags and is the natural place to expose requirements rather than removing the tutorial.

**Impact.** In a successful lightweight installation, discovery and planning present this workflow alongside runnable EP workflows; failure is deferred until missing full-mode libraries/types are encountered. This is misleading public behavior and weakens autonomous planning.

**Minimal remediation.** Add explicit required capabilities (for example `requires_full_mode` or required libraries) to registered-entry metadata, surface them in `describe`/Run documents, and have preflight reject or mark unavailable when the resolver/environment indicates lightweight mode. Do not hide the tutorial; retain it with an actionable reason.

**Required validation / automated contract test.** Run registry discovery under simulated `USE_LIGHTWEIGHT_PHYSICSMODEL=1` and full mode. Assert the EM tutorial is marked unavailable with the documented requirement in lightweight mode and available in full mode; assert all other tutorials declare or inherit a mode contract.

## Existing guards and residual gaps

- **Healthy:** `test_ionic_catalog_contract.py` and `test_ionic_catalog_audit.py` compare non-manufactured C++ Names enums with catalogue tuples and deliberately exclude manufactured semantic mappings.
- **Healthy:** `test_run_model.py` asserts the source and packaged Run-document schemas are equal and exercises model round trips. The two schema files were byte-identical (`cmp` exit 0).
- **Healthy:** `drift_guards/test_rtst_enum_contract.py` covers selected dictionary runtime enums.
- **Gap:** no complete set-equality test joins `Make/files` -> registration -> runtime name -> Python catalogue -> README for ionic types.
- **Gap:** no equivalent complete active-tension registration/catalogue test; F2 is the direct consequence.
- **Gap:** capability tests assert presence of selected valid names but not absence of metadata-category contamination or false regions.
- **Gap:** tutorial registry/display equality does not encode run-mode availability.
- **Caution:** `dict_key_allowlist.json:91-104` labels 13 catalogue paths as `stale_paths`. Because the scanner is documented as approximate, these are investigation debt rather than 13 confirmed defects. A future parser-aware guard should replace exemptions one by one.

## Priority order

1. F3: incorrect planning accept-surface can conceal invalid output requests.
2. F1: public copy/paste runtime names fail deterministically.
3. F2: executable MMS model is missing from machine-readable introspection.
4. F4: model build-mode requirements before launch.

All proposed remediations are metadata, documentation, or validation changes; none requires changing scientific equations, solver ordering, units, or generated numerical kernels.
