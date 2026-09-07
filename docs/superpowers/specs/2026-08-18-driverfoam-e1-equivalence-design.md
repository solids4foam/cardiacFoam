# E1 — Numerical Equivalence: Design

Date: 2026-08-18
Branch: `ep-work-onto-main`
Roadmap items: P3.1, P3.2 (`future/driverFOAM/ROADMAP.md`)
Status: design, not yet planned or implemented

## 1. Context

Phases 1 and 2 of the driverFOAM roadmap are complete. Phase 3 produces the
deterministic evidence the rest of the paper depends on. E1 is first because
every other Phase 3 experiment (E2 fault corpus, E4 overhead, the legacy
16-experiment migration) needs the same underlying capability: stage a case,
run it more than one way, compare what came out.

This spec covers E1 only. It also carries two enabling items that E1 forces
into the open (§7 W1, §7 W6) and explicitly defers a third (§4).

**Framing that emerged during design and should survive into the plan:** most of
what E1 needs already exists inside driverFOAM and is either unwired or
under-used. The artifact reconciler is written and tested but has no production
callers (§7 W1); staging, reference parsing, and the strict run path are all in
`regression_equivalence` (§5); the error norms are computed in C++ by the solver
(§3). E1 is therefore predominantly a *find-and-wire* exercise, with two genuinely
new pieces: the minimal Python runner (§7 W2) and the comparator (§7 W4). Plans
derived from this spec should look for existing capability before adding any.

## 2. Research question, and the claim we are actually making

E1 asks: **does driverFOAM alter numerical results?**

**Configuration is held constant by construction.** Dictionary state is an
*input the agent chooses*, not a property under test. The driver mutates
dictionaries when asked to and runs the case as authored when not; both are
correct, and which one happens is determined by the request. E1 therefore
stages one tree, applies no overrides, and runs it three ways. All three
conditions execute identical inputs.

That makes E1's claim precise: **routing an identical case through driverFOAM's
orchestration does not change the numbers.** The candidate mechanisms by which
it *could* are real and none of them are dictionary questions — parallel
decomposition, command ordering, working directory, environment, and restart or
resume handling can all move results.

Dictionary *synthesis* fidelity is a separate question and is not E1's job. It
is already partly covered in the existing harness by
`regression_equivalence/round_trip.py::electro_build_parse_fixpoint`, the
build/parse idempotence check reported as the `idempotent` column of
`python -m openfoam_driver.tests.regression_equivalence`.

The evidence is stratified:

- **Primary — execution-context equivalence.** Same staged tree, same command
  sequence, same environment, plus a guard asserting the driver changed nothing
  it was not asked to change. Solver-free, therefore cheap enough to cover
  *every* case in the registry. This is what makes the claim general rather
  than anecdotal.
- **Confirmation — solver-backed artifact comparison** on a small
  representative set. A reviewer asking "does driverFOAM alter numerical
  results?" wants to see numbers agree, not only an argument that they must.

This stratification shrinks the expensive part of E1 from "six families in
triplicate" to a confirmation set, and widens the cheap part to full registry
coverage.

## 3. Key finding that shapes the design: norms are solver-emitted

`src/verificationModels/verificationUtils.H` computes the error norms in C++,
not Python. Three `computeNorms()` overloads:

| Overload | Returns | Line |
|---|---|---|
| Mesh-free (graph/node data) | parallel L1, L2 (cell-count-averaged), L-inf | :125 |
| `fvMesh` cell-centred | volume-weighted L1/L2 + unweighted L-inf | :177 |
| Bulk/boundary split | volume-weighted L2, partitioned | :249 |

Seven verifier models consume it — `monodomainVerification`,
`bidomainVerification`, `bathBidomainVerification`, `eikonalVerification`,
`ecgVerification`, `coupledVerification`, `electromechanicsVerification` —
covering all six of E1's required case families. They write to
`mesh_.time().globalPath()/"postProcessing"/*.dat`
(`manufacturedEikonalVerifier.C:259`).

The tutorials' Python (`post_processing_manufactured.py` and siblings) only
*parses and plots* those `.dat` files. A raw OpenFOAM field reader already
exists too — `monodomainPseudoECG/setup/mesh/tet/volume_weighted_norm.py`,
which self-validates against the verifier's own `.dat`.

**Consequence:** E1 must not write a field reader, and must not re-derive
norms. Norms are a solver output, read like any other artifact.

## 4. Non-goals

- **No raw-field reader.** Superseded by §3.
- **No invocation of tutorial analysis scripts, in either direction.** E1
  compares artifacts; it does not run analysis to produce them.
- **No internalization of post-processing into driverFOAM.** The existing seam
  is correct and stays: `openfoam_driver/postprocessing/driver.py` defines
  `PostprocessTask(module_relpath)` and loads a case-owned module via
  `importlib`. The driver owns *when analysis runs, where output lands, and
  house style*; the case and the agent own *what is computed*. This mirrors the
  call already made for OpenFOAM function objects.
- **Not E4's sweep-scaling axis** (1/10/100/1000). E1 captures per-condition
  timing as a by-product (§7 W3); only the scaling campaign is deferred.
- **Not E2's fault corpus, and not E3's interruption harness.** E2 reuses this
  spec's staging and comparison code. E3 is an independent harness that
  neither blocks nor is blocked by E1.
- **Not E7.** The analysis-script provenance hole (an external post-processing
  script is currently untracked, so a result CSV's producer is unverifiable) is
  real and is closed by hashing the script as an input. That belongs to E7 and
  is recorded here only so it is not lost.

## 5. Conditions

All three run from one staged input tree, reusing
`regression_equivalence/dual_run.py::_stage_tutorials_root`, which copies a case
into a throwaway tutorials root and strips `postProcessing`,
`workflow_state.json`, and `workflow_logs` so outputs are provably fresh.

| ID | Condition | Status |
|---|---|---|
| A | Canonical direct `Allrun` | exists in each case |
| B | Minimal Python runner | **new, must be written** |
| C | `foamctl run --strict` | exists (`dual_run.py::_drive_agent`) |

**Condition B is the one component with no existing code to lean on.** It is a
deliberately naive subprocess runner (~80 lines): read the case's command
sequence, run each step, capture rc/stdout/stderr, stop on failure. No dict
mutation, no validation, no planning. Its sole purpose is to separate
*driverFOAM's* effect from *any Python orchestration's* effect. Without it, a
difference between A and C cannot be attributed. It must not import
`openfoam_driver`.

### Relationship to the existing harness

`regression_equivalence/dual_run.py` compares **driver output against the
committed `.reference` file** — its own docstring states "the hand-authored path
is not re-run". That answers a different question: it measures whether the
driver still hits historical numbers, which conflates driver effects with any
solver drift since the reference was committed. E1 needs conditions run
*against each other* in one environment at one commit. Reusable from that
module: `_stage_tutorials_root`, `parse_columnar_reference`,
`read_series_value`, `values_agree`, `solver_available`. Not reusable:
`verify_reproduction`'s comparison axis.

**Doc bug to fix in passing:** `dual_run.py`'s module docstring says the strict
run "applies its dict overrides". It does not — `run --strict` plans
(non-mutating) and executes; `apply_case` has no call site in the run path.
Harmless in practice because `dual_run` passes no overrides, but the sentence is
misleading and should be corrected.

## 6. Comparison layers

### L1 — execution-context equivalence (primary, solver-free)

All three conditions run the same staged tree with no overrides applied, so
there is nothing for the driver to materialize and nothing to reconcile between
divergent input states. `strict_plan()` is non-mutating by declaration
(`strict_planning.py:348`), and the run path never calls `apply_case` or
`materialize_case` — those have call sites only in the sweep path
(`sweep_runner.py:147`). Condition C therefore executes the committed
dictionaries unchanged, exactly as A and B do.

Captured per condition, before the solver runs:

- **Input-tree hash.** Per-file sha256 plus a tree digest over the staged case,
  excluding known-volatile paths. Establishes that all three conditions really
  did start from the same inputs.
- **Unchanged-dictionary guard.** Assert the driver modified no dictionary it
  was not asked to modify. This is a guard, not the headline evidence: it fails
  loudly if the no-override invariant ever stops holding.
- **Command sequence.** The ordered list of commands each condition issues,
  with argv and cwd. This is the primary object of comparison — orchestration
  is what differs between conditions.
- **Environment manifest.** Solver binary path and hash, `WM_PROJECT_DIR`,
  OpenFOAM version, relevant environment variables, decomposition settings,
  machine identity.

A divergence here localizes the defect at the point it occurs, rather than
inferring it from a downstream number.

### L2 — artifact reconciliation and content (confirmation, solver-backed)

For the confirmation set only. The comparison object already exists: E1
compares `ReconciliationReport`s (§7 W1) across conditions rather than defining
a new format. Compared:

- **Reconciliation** across conditions: same artifact ids matched, same
  resolved paths, same status, same undeclared set.
- **Content** of artifacts declared tabular: the verifier `.dat` norm files and
  any function-object samples, compared numerically against the frozen
  tolerances (§8).

L2 never opens a tutorial analysis script and never re-derives an observable.

#### 6.1 The common ruler, and why it is sound

Conditions A and B have no reconciler — reconciliation is a driver feature. But
`reconcile_artifacts(case_root, predicted)` is a pure function of a directory
plus a prediction list, so E1 applies **the driver's own reconciler to all three
output trees as a common ruler**.

This is not circular, because the artifact set is **statically determined**. The
solver writes nothing dynamically: time directories follow `controlDict`'s
write controls, verifier `.dat` names are compiled in, and function-object
output lands in `postProcessing/<name>/` where `<name>` is a declaration the
driver itself holds. `niederer_2012.py` shows both ends of the mapping in one
place — the spec holds the function-object name, derives the `postProcessing`
subdirectory from it, and constructs the output filename itself
(`f"{solver}_{ionic_model}_{tissue}_points_DT{dt_tag}_DX{dx_tag}.csv"`).

Prediction is therefore derived from the same declarations that *cause* the
output, not guessed against it. Prediction defects are ordinary bugs in a
deterministic mapping — the empty `export ()` list falling through to catalog
defaults (fixed in `f5f935c2`) is the canonical example — and unit tests are the
right instrument for them, not output sweeping.

**Known limitation.** `#includeFunc` declarations are skipped by the
function-object scanner (`specs/function_object_fields.py:83`, with a test
pinning the behaviour), so their outputs would not be predicted. No tutorial in
the repository currently uses `#includeFunc`, so this is unexercised. It is
recorded here as a stated limitation rather than a work item.

## 7. Work packages

### W1 — Wire the existing artifact reconciler

**The capability already exists and is unwired.** This is the single most
important scoping fact in this spec: E1 is largely an exercise in finding and
connecting driver capability that is already written, not in building new
machinery.

`core/runtime/reconciler.py` provides
`reconcile_artifacts(case_root, predicted, case_id=None) -> ReconciliationReport`:

- `predicted_count`, `matched_count`, `missing_count`
- per predicted artifact: `artifact_id`, `predicted_path`, `status`
  (`matched|missing`), `matched_files` as
  `[{path, kind, size_bytes | entries}]`, `optional`
- `case_id`, with `_substitute_case_id` expanding `{case_id}` so per-case
  outputs are attributable across a sweep
- `_glob_under` resolving a `path_pattern` to the concrete files that landed

Its own inline comment states the intent plainly: the per-match entry count
exists "so an agent inspecting the realized manifest knows whether the matched
dir is non-empty". It is fully unit-tested in `tests/core/test_reconciler.py`
and has **no production callers**.

Two other artifact mechanisms are wired but are not substitutes:

| Mechanism | Location | What it gives | Why it is not enough |
|---|---|---|---|
| Prediction | `core/runtime/artifacts.py:161` | `path_pattern`, `format`, `variables`, `produced_by`, `time_indexed`, `optional` | Pre-run. Patterns, not resolved files. |
| Per-step gate | `core/runtime/workflow_runner.py:297-320` | Fails a step with `missing_artifacts` | Boolean, per-step, keeps no record. |

Work, in order of size:

1. **Wire `reconcile_artifacts` into the strict run path and persist the
   report** into the run document or run manifest. This is the bulk of W1.
2. **Add `sha256` beside `size_bytes`** for file matches. Size alone is weak
   provenance and cannot support E1's content comparison.
3. **Optional: undeclared-match detection.** The report is one entry per
   *predicted* artifact, so a file appearing without having been predicted is
   invisible. Per §6.1 this is **not** a check on solver output, whose artifact
   set is deterministic; at most it is a cheap safety net on the non-solver
   steps (`decomposePar` writing `processor*/`, mesh utilities, the case-owned
   postprocessing module), most of which is ignore-list material. Cut this
   first if the plan needs trimming.

Description fields stay in `expectedArtifacts`. The reconciliation report
references `artifact_id` and adds only resolution, existence, and integrity;
`format`, `variables`, and `produced_by` are not duplicated. Single source of
truth for description, separate record for observation.

**Invariant: the driver describes, the agent interprets.** The manifest states
"columnar `.dat`, columns `[Field, L1-error, L2-error, Linf-error]`, produced by
step `solve`, sha256 ...". It never states whether a value is good, whether the
run converged, or which artifact matters. The moment the manifest judges
values, domain semantics re-enter driverFOAM through a new door.

This is also what makes E1 generalize: because the manifest is domain-neutral,
the same comparison works on a real anatomical case where no verifier model and
no manufactured-solution post-processor exist.

### W2 — Minimal Python runner (condition B)

Per §5. Standalone, no `openfoam_driver` import.

### W3 — Three-condition orchestration

Stage once, run A/B/C, capture L1 evidence per condition, capture the observed
manifest per condition. Captures wall-time, peak RSS, and total output bytes
per condition from day one — roughly twenty lines that hand E4 its
realistic-solver-case rows without a second campaign.

### W4 — Comparator and report

Smaller than first scoped, because W1 supplies the comparison object. Consumes
W3's L1 records plus one `ReconciliationReport` per condition, and emits
per-case comparison CSVs and a summary. Failure output must localize: which
layer, which artifact id, which file, which key or column.

### W5 — `equivalence_protocol.yaml`

Per §8.

### W6 — `plotting_common.py` genericity tidy-up

Independent of the rest; listed here because E1's argument depends on
`openfoam_driver/postprocessing/` being genuinely generic. Audit result: the
directory and path handling *is* clean — `load_csv_folder(folder, glob_pattern,
drop_columns=...)` parameterizes everything, and its only domain assumption is
dropping ParaView export columns, which is OpenFOAM-generic. `driver.py`, the
actual orchestration seam, has zero cardiac vocabulary.

The residue is three helpers in `plotting_common.py`:

| Item | Non-generic aspect | Callers | Action |
|---|---|---|---|
| `rename_cardiacfoam_trace()` :96 | Hardcodes the literal `" cardiacFoam"` and a `", ΔT="` separator — solver-brand-specific, not merely cardiac | 1 (Niederer `line_postProcessing.py`) | Move to the case |
| `parse_model_and_cell()` :74 | Body is generic; the signature encodes electrophysiology ("cell" = ionic cell model) | 1 (`PATHOS/.../singleCellinteractivePlots.py`, untracked) | Rename to `parse_two_part_stem(first_map=, second_map=)`; no behaviour change |
| `extract_dx_dt()` :39 | Hardcodes `DX`/`DT` tokens and mm/ms scaling (`/10`, `/10**(len-1)`) | 5 | Parameterize pattern and scale, or keep and document as a stated convention |

Docstring-only occurrences — `plot_builder.py:98-100`, `:407`, `:515`,
`table_writer.py:48` — are examples inside module docstrings and are a cosmetic
swap.

## 8. Tolerances and preregistration

`EXPERIMENTAL_EVIDENCE_TABLE.md` requires that thresholds "must not be selected
after seeing the final comparison". Tolerances are therefore frozen in
`equivalence_protocol.yaml` **before** the confirmation set runs.

Most rows need no pilot. The committed `.reference` files already carry
per-point tolerances, authored long before E1 existed, so adopting them is
preregistration-safe by construction. `equivalence_protocol.yaml` transcribes
them with per-row provenance: which reference file, which commit, why that
value.

Only observables with no existing tolerance — chiefly the norm rows and any
convergence-order comparison — require a short pilot. The pilot covers *those
rows only*, is run on the cheapest cases, and its output is committed before
the confirmation set runs.

## 9. Component boundaries

| Component | Purpose | Interface | Depends on |
|---|---|---|---|
| Stager | One tree, three conditions | `stage(case) -> (root, case_path)` | existing `_stage_tutorials_root` |
| Condition runners | Execute A / B / C | `run(case_path) -> RunRecord` | subprocess; C uses `foamctl` |
| Evidence capture | L1 records | `capture(case_path, env) -> InputEvidence` | dict parsing, hashing |
| Reconciler (existing) | Observed artifacts | `reconcile_artifacts(case_root, predicted, case_id) -> ReconciliationReport` | already written; needs wiring + sha256 + undeclared set (W1) |
| Comparator | L1 + L2 across conditions | `compare([RunRecord], protocol) -> Report` | `equivalence_protocol.yaml` |

Each is independently testable: the stager, evidence capture, and comparator are
solver-free and unit-testable anywhere. Only the condition runners need a built
cardiacFoam, and they follow the existing `solver_available()` self-skip
convention so the module imports and its pure helpers test without a solver.

## 10. Testing

- Unit: comparator against synthetic divergent/identical record pairs; evidence
  capture against fixture trees. The reconciler already has unit coverage in
  `tests/core/test_reconciler.py`; extend it for sha256 and the undeclared set
  rather than starting a new test module.
- Integration (solver-gated): full three-condition run on the cheapest case.
- Negative: inject a known command-sequence divergence and assert L1 localizes
  it to the correct step; trip the unchanged-dictionary guard and assert it
  names the file and key; inject a known artifact divergence and assert L2
  localizes it to the correct file and column.
- CI: the solver-free layers run on PR; the confirmation set runs on the
  scheduled or release matrix, matching the roadmap's CI-stratification item.

## 11. Open questions

**Closed during design.**

- *Sweep-coordinate binding.* Answered: `ReconciliationReport.case_id` plus
  `_substitute_case_id`, which expands `{case_id}` in a `path_pattern` to the
  literal case id, or to `*` for whole-run reconciliations. No new work.
- *Materialization mechanism.* Dissolved. E1 applies no overrides, so there is
  nothing to materialize. `run --strict` never calls `apply_case` or
  `materialize_case` — those exist only in the sweep path — so condition C
  executes the committed dictionaries unchanged. The earlier question of
  whether a size-one `sweep-plan` writes the same dicts as `run --strict` is
  moot.

- *Confirmation-set membership.* Answered: the nine cases in
  `tutorials/Alltest-regression`'s `REGRESSION_TESTS`, which are exactly the
  nine rows of `regression_equivalence/registry.py`. Verified passing in
  lightweight mode, where `electroMechanicalNiedererEtAl2011` is an expected
  skip (exit 77), leaving eight live cases.
- *Condition A command extraction.* Answered, but **not** as a constant. Each
  case's `regression/regressionTest.sh` runs `./Allclean` then `./Allrun`, and
  all nine carry both as executables — but **six of the nine pass `parallel`**
  (`decomposePar` + `runParallel` + `reconstructPar`); only singleCell,
  bidomain, and bathBidomain are serial. Condition A must therefore be read
  per-case from that script. Hardcoding the serial form runs a different
  simulation from the one the regression certifies, and parallel decomposition
  genuinely moves numbers — this is the precise false-positive mechanism
  Open Question 4 warns about, and it was hit in practice during execution.

Open:

1. **Volatile-path exclusion list** for the input-tree hash. Needs to be
   enumerated empirically and committed, not guessed.
2. **Where the reconciliation report is persisted** — run document versus run
   manifest. Affects whether W1 needs a `run-document.json` schema change or
   only a manifest addition.
3. **Decomposition handling.** Whether the three conditions decompose
   identically, since parallel decomposition is one of the few mechanisms that
   can genuinely move numbers. May belong in the environment manifest rather
   than as a separate check.
4. **Condition B's depth.** As specced, B invokes the same `./Allclean` /
   `./Allrun` scripts as A, differing only in the invoking process. That
   isolates Python's subprocess environment and cwd handling but not command
   decomposition. A stronger variant would have B invoke the individual solver
   and utility commands *inside* `Allrun`; that needs a per-case command list
   and is deferred.

## 12. Sequencing

E1 (this spec) -> E2 fault corpus (reuses stager and comparator) -> legacy
16-experiment migration -> E4 scaling campaign. E3 (interruption and recovery)
is an independent harness and can slot in at any point. W6 is independent of
all of them and can be done at any time.
