# dict_key_allowlist.json — What's Actually In It, And Why

**Status: audited and closed out.** What looked like a large backlog turned
out to be almost entirely already-correct catalog coverage hidden behind
scanner blind spots, plus a handful of dead C++ keys (now deleted) and two
genuine gaps (now resolved — one catalogued, one deliberately left out of
scope with the reasoning documented). See §4d for the final tally. This
document exists so a future session doesn't have to re-derive any of this
from scratch.

## 1. What this file is

`applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/dict_key_allowlist.json`
is a waiver list for the strict dict-key scanner
(`openfoam_driver/scripts/_dict_keys_scanner.py`). The scanner regex-scans
every `.C`/`.H` file under `src/` for OpenFOAM dictionary-read call sites
(`dict.lookup(...)`, `dict.get<T>(...)`, `dict.found(...)`, etc.), and
compares what it finds against `dict_entries_catalog.py`'s `driver_path`
entries. Anything the scanner finds that isn't resolvable to a catalog entry
is "drift"; `test_strict_dict_key_scanner_allowlist_is_current` (in
`openfoam_driver/tests/core/test_strict_planning.py`) fails strict planning
if there's drift that isn't listed in this file.

The file has three arrays:
- `absent_keys` (90 entries as of this writing) — bare key names found by
  the scanner with no matching catalog `driver_path`.
- `stale_paths` (10) — scanner blind spots: constructs the regex can't
  follow (`dimensioned<T>(name, dict)`, `spec.dictionaryEntry` indirection,
  dicts passed into upstream OpenFOAM library code).
- `unmatched_subdicts` (6) — sub-dictionary groupings the scanner reports
  as flat keys, mostly the `<scope>` placeholder problem.

**The mechanism for detecting a stale waiver already exists and works**:
`strict_dict_key_report()`'s `unused_allowlist` field lists any waiver entry
the scanner no longer emits (key was removed/renamed from `src/`, or — as
happened here — someone properly cataloged it). Ran it this session: as of
the commits below, `unused_allowlist == []`. Every one of the 90
`absent_keys` was individually re-verified by hand (grep for the literal
`"key"` string in `src/`) to have a real, current read site — see section 3.

## 2. What triggered this audit

While fixing an unrelated `purkinjeNiedererEtAl2011` regression bug (see
`future/PVJ_KERNEL_MAGNITUDE_FOLLOWUP.md`), a pre-existing uncommitted diff
was noticed adding `useGraphPrePopulation` to `absent_keys` — i.e. someone
had waived a brand new dict key instead of cataloging it. Checked what the
key actually does
(`src/electroModels/electroDomains/myocardiumDomain/eikonalMyocardiumDomain.C`,
`.H`): a `Switch`, default `true`, controlling whether the eikonal
activation-time field is warm-started via a parallel-safe Bellman-Ford relay
from the stimulus/graph seed cells before the nonlinear solve, instead of
being left at the `GREAT` placeholder. Nothing about it resisted proper
cataloging — wrote a full `DictEntry` for it in `dict_entries_catalog.py`
(the `eikonal_diffusion` group), which made the waiver redundant, so it was
removed. Both changes verified via
`test_strict_dict_key_scanner_allowlist_is_current` and committed
(`37478cd8`, catalog entry; the allowlist removal landed as part of the same
change since removing the redundant line reverted that file to its
committed `HEAD` state exactly).

That raised the obvious question for the other 89: **are any of them
actually there for a good reason, or is this just backlog wearing a
waiver's clothes?**

## 3. The full audit — what each key does, and whether it belongs here

Every one of the 90 has a real, current dict-read call site in `src/` (hand
grepped, not just trusted from the scanner's aggregate `status: ok`). Below,
grouped by area, with a one-line description of what each does.

### Purkinje graph / conduction-system topology
(`constant/<graphFile>`, read by `conductionSystemDomain.C`)
- `graphFile` — name of the graph file in `constant/` to load
- `points` — 1D Purkinje node coordinates
- `rootNode` — index of the graph's root/origin node
- `pvjNodes` / `pvjLocations` — indices and 3D coordinates of the
  Purkinje-ventricular junction terminals
- `pvjResistances` — per-junction resistance overrides (falls back to
  `rPvj` if absent). **Resolved permanently out of scope — see §4d.**
- `purkinjeConductivity` — 1D Purkinje fibre conductivity
- `node` — root-stimulus target node index
- `vm1DRest` — resting Vm for the 1D Purkinje solve
- `rootStimulus` (subdict) — root-node pacing config, containing
  `startTime`, `duration`, `intensity`
- `escapeInterval` — escape-rhythm pacing interval (restitution-eikonal 1D
  solver)
- `apdNominal` — nominal APD used by the restitution model
- `referenceConductance` — normalization conductance for the
  restitution-eikonal solver
- `useEdgeConductance` — toggle whether edge-wise conductance is used
- `stimulus` (subdict) / `sites` — spatial stimulus sites for the
  restitution-eikonal 1D solver

### PVJ coupling
(`domainCouplings.<name>`, `pvjCoupler`/`reactionDiffusionPvjCoupler`)
- `electroDomainCoupler` — coupler-type selector
- `conductionNetworkDomain` — which `conductionNetworkDomains` entry this
  coupling targets
- `couplingMode` — unidirectional vs bidirectional
- `debugCoupling` — verbose PVJ coupling logging
- `rPvj` — junction resistance
- `pvjKernel`, `pvjCouplingScheme`, `pvjRadius` — **already properly
  cataloged**; still flagged by the scanner's regex (dynamic-path/subdict
  construct it can't resolve), matching the documented blind-spot pattern
- `coupling` (subdict) — generic per-domain coupling block in the system
  builder

### ECG / torso
- `ecgDomains` (subdict), `ecgSolver` — ECG sub-model selection
  (`pseudoECG`/`torsoECG`/`eikonalECG`)
- `electrodePositions` — named electrode coordinates
- `torsoSurface` — STL path for the torso geometry, used by `torsoECG`
  (`src/genericWriter/ecgModelIO.C`)
- `groundPatches` — patches used as the electrical ground reference
- `surfaceCurrentPatches` — patches where surface current is applied/read
- `sigmaExtracellular` — extracellular conductivity used by
  `eikonalECG`/`pseudoECGSolver`
- `start` / `end` / `deltaT` (inside a `sampling` subdict) — ECG output
  sampling window
- `xMin` / `xMax` — bounding coordinates picking out the single
  ground/current patch in `bathECGManufacturedVerifier`
- `manufacturedBidomain` (subdict) — manufactured-solution override block
  for bidomain ECG
- `manufacturedEikonalECG` (subdict) — manufactured-solution override for
  `eikonalECG`
- `checkQuadratureOrders`, `referenceQuadratureOrder` — Gauss-Legendre
  quadrature settings for manufactured ECG reference integrals

### Ionic models / tissue heterogeneity
- `Vm`, `Iion` — **legitimate exclusion**: field-name string comparisons in
  I/O code, confirmed false positives, never dictionary keys
- `ionicConstantOverrides` (subdict), `ionicHeterogeneity` (subdict) —
  per-region ionic parameter overrides
- `global` — the `<scope>` placeholder value for whole-mesh overrides
- `regions` (subdict), `range`, `apexBaseBands` — heterogeneity region
  definitions (apex-to-base banding, value ranges)
- `endocardialCells`, `mCells`, `epicardialCells` — the three transmural
  tissue-layer names
- `scale`, `set` — per-parameter override operators inside a
  heterogeneity/override scope
- `baseline` — reference/default value inside heterogeneity blocks
- `outputSuffix`, `outputVariables` (subdict), `ionic` (subdict) —
  ionic-model output configuration
- `ODESolver` (subdict) — gating-variable ODE integrator config; the
  `absTol`/`relTol`/`maxSteps`/`solver` leaves live in `stale_paths`, not
  here, because they're read by upstream OpenFOAM's own `ODESolver::New`
  factory, not by any code in this repo

### Active tension / electromechanics + their manufactured verifiers
- `activeTension` (subdict), `constants` (subdict), `initialStates`
  (subdict) — active-tension model config blocks (e.g. `LandNiederer`)
- `preconditioningTime` — ms of ODE preconditioning run at constant
  Ca_i/lambda before the real solve starts, to equilibrate `LandNiederer`'s
  active-tension state variables and avoid a spurious global Ta transient
  (see `LandNiederer::preconditionToRestingState`). **Correction (see
  section 4a): NOT out-of-scope for this catalog** — `activeTensionModel::New`
  is called both from `sequentialElectroMechanical.C` (passes
  `electroMechanicalProperties()`, out of scope) and `singleCellSolver.C`
  (passes `electroProperties()` directly, in scope). Catalog-eligible under
  the `electroProperties.activeTensionModel`/`singleCellSolver` path.
- `electroMechanicalModel` — top-level EM model type selector
- `electromechanicalVerificationModel` (subdict) — manufactured-EM verifier
  selector
- `initializeFields` — **legitimate exclusion**: real key, but belongs to
  `electroMechanicalProperties`, a different dict, deliberately out of
  scope for this `electroProperties` catalog
- `TaScale`, `gamma` — **legitimate exclusion**, same reason as
  `initializeFields`
- `Tmax`, `V0`, `amplitude` — manufactured-electromechanics reference-
  solution parameters (same dict-scope caveat likely applies; not
  double-checked against `electroMechanicalProperties` vs `electroProperties`
  scoping the way `TaScale`/`gamma`/`initializeFields` were)

### Bidomain/bath verification output
- `outputFile` — output CSV/file path for
  `manufacturedFDABidomainVerifier`/`manufacturedFDABathBidomainVerifier`
- `bathPotentialDomain` (subdict) — bath/torso potential-domain config for
  bidomain
- `conductionNetworkDomains` (subdict), `conductionSystemDomain`,
  `conductionSystemSolver` — top-level names/type-selectors tying a
  `conductionNetworkDomains` entry to its solver
- `verificationModel` (subdict) — generic manufactured-verifier selector
  shared by several domain types
- `ecgVerificationModel` — manufactured verifier selector specific to ECG
- `singleCellStimulus`, `externalStimulus` (subdicts) — S1/S2 protocol and
  spatial stimulus config for single-cell/tissue runs
- `nNonOrthogonalCorrectors` — **correction (see section 4a): out of scope,
  not a legacy fallback.** Read via
  `baseMesh.solutionDict().subOrEmptyDict("PIMPLE").lookupOrDefault<label>("nNonOrthogonalCorrectors", 0)`
  in `extracellularPotentialDomain.C` — this comes from `system/fvSolution`'s
  `PIMPLE` block, the exact same key any stock OpenFOAM PIMPLE/SIMPLE solver
  reads. It was never an `electroProperties` key to begin with, so this
  catalog is the wrong place for it regardless of how often it's set —
  same bucket as the `absTol`/`relTol`/`maxSteps`/`solver` ODESolver
  blind spots two sections up.
- `manufactured` (subdict) — generic manufactured-solution override block,
  reused by several verifiers

## 4. Verdict: is there a real reason these have to stay waived?

Mostly no — but the real split isn't the A/B/C/D-by-topic grouping above.
See `future/DICT_CATALOG_NORMALIZED_PATTERN.md` for the full methodology;
the short version is below.

### 4a. Two corrections from an initial pass (don't repeat these mistakes)

An earlier pass through this list called `preconditioningTime` and
`nNonOrthogonalCorrectors` "legacy fallbacks" (optional keys with defaults,
rarely overridden, therefore assumed vestigial). Both were wrong, for two
different reasons worth internalizing:

- **`preconditioningTime`** is not dead — it's a real, purposeful
  electromechanics-initialization step (see the corrected description
  above). It looked out-of-scope/rare only because its containing class
  (`LandNiederer`, via `activeTensionModel::New`) is called from *two*
  different properties files, and only one of its two call sites was
  checked at first. **Lesson: grep every call site of the reading
  function/constructor, not just the first one you find, before deciding a
  key is out of scope.**
- **`nNonOrthogonalCorrectors`** isn't a fallback of *this* catalog at
  all — it's read from `system/fvSolution`'s `PIMPLE` block, standard
  OpenFOAM numerics control, not `electroProperties`. It looked like a
  rarely-set optional key for the same superficial reason (`lookupOrDefault`
  with a default), but the real tell was which **dict object** it's read
  from, not how often it's overridden. **Lesson: "has a default, rarely
  overridden" is not evidence of staleness — check which properties file
  actually owns the key before concluding anything about how commonly
  it's used.**

### 4b. The bigger split: structural containers vs. leaf values

Cross-checked against driverFOAM's actual schema
(`openfoam_driver/core/contracts/dictionary.py`) and dict builder
(`openfoam_driver/specs/dict_builder.py`): `DictEntry` is **leaf-only** —
there is no container/block type in the schema at all. Roughly a third of
the 90 `absent_keys` are pure structural/subdict names the scanner flags
because it sees `dict.subDict("X")`/`dict.found("X")` in the C++ and has no
way to know `X` is just a path segment, not a settable value:

`ODESolver`, `activeTension`, `apexBaseBands`, `bathPotentialDomain`,
`conductionNetworkDomains`, `constants`, `coupling`, `ecgDomains`,
`electrodePositions`, `externalStimulus`, `global`, `groundPatches`,
`initialStates`, `ionic`, `ionicConstantOverrides`, `ionicHeterogeneity`,
`manufactured`, `manufacturedBidomain`, `manufacturedEikonalECG`,
`outputVariables`, `regions`, `rootStimulus`, `scale`, `set`,
`singleCellStimulus`, `stimulus`, `surfaceCurrentPatches`,
`verificationModel` (28 keys, plus `electromechanicalVerificationModel`,
already out-of-scope anyway).

**These should never get their own `DictEntry`.** It's not that they're
hard to catalog — `dict_builder.py`'s `_set_nested`/`_serialize_block`
literally cannot read a value back out of a `driver_path` that terminates
at a container; the builder only ever writes a value at the final segment
of a leaf entry's path. Cataloging the container itself would be
documenting something the code has no way to act on. What actually needs
cataloging is the *leaves already nested inside them* — several of which
(`rootStimulus`'s `startTime`/`duration`/`intensity`/`node`, for example)
are separately already on the genuine-backlog list below.

The remaining ~52 keys not in the container list above — `torsoSurface`,
`pvjResistances`, `apdNominal`, `escapeInterval`, `deltaT`, `outputFile`,
`rPvj`, `xMin`/`xMax`, and most of the rest — are genuine leaf values: real,
single-purpose, one or two clear read sites, no ambiguity about
`driver_path`. Structurally identical to `useGraphPrePopulation` before it
was cataloged. Nobody has written the `DictEntry` yet; nothing prevents it.

A handful need more design thought before a `DictEntry` is written, not
because they're containers, but because the *leaf itself* is a generic
name reused across multiple different parent scopes (`range` inside
different heterogeneity region types, `node`/`start`/`end` used generically
in more than one subdict) — seeing full details worked out in
`future/DICT_CATALOG_NORMALIZED_PATTERN.md` section 4's decision table.

**Bottom line**: this file is a todo list wearing a waiver's clothes for
most of its entries, a set of names that structurally can never become
catalog entries (the containers) for about a third of them, and a handful
of genuinely out-of-scope/false-positive/upstream-owned keys for the rest.
Worth working through the genuine-leaf group incrementally (easy wins, same
pattern as `useGraphPrePopulation`), leaving the container names
permanently un-catalogued (their children are the real target), and the
false-positive/out-of-scope/upstream groups alone.

## 4c. Four keys resolved a different way: deleted, not catalogued

Digging into what each of `report` (eikonalECG), `reportElectrodeLookup`
(torsoECG), `profileBatchedModel` (batchedIonicModel.H), and `probeNodes`
(conductionSystemDomain) actually did on closer inspection turned up a
third resolution besides "catalog it" and "leave it waived": all four
were either fully unwired (`profileBatchedModel` gated a summary-printing
method with zero callers anywhere in the repo) or genuinely pointless as
optional switches (`report`/`reportElectrodeLookup` gated cheap
diagnostic `Info` output with no reason to ever disable it;
`probeNodes` let a case restrict expensive per-node ionic export, but
every real case already got the "restrict to nothing" default of "export
everything" since nothing ever set it).

Rather than write `DictEntry`s for these, they were removed from the C++
source entirely — see commit `5c465ac4`
("refactor: remove unused/unwired debug toggles and dead profiling
machinery"), which also traced `profileBatchedModel`'s removal down
through `BatchedKernelExecutor` in `batchedKernelExecution.H` (its only
construction site, so its whole per-thread timing subsystem was
provably dead too). Regression-verified for three of the four
(`probeNodes`: `purkinjeNiedererEtAl2011`; `profileBatchedModel`'s hot
path: `singleCell`; `report`: `eikonalECG`) with no reference changes
needed. `reportElectrodeLookup` has no tutorial exercising `torsoECG` in
the regression suite, so it's compile-verified only.

Their `dict_key_allowlist.json` waivers are removed accordingly (the
scanner no longer finds them at all, confirmed via `unused_allowlist`
before deletion) — not because they got catalogued, but because they no
longer exist to catalog. This is a fourth outcome worth watching for when
working through the rest of this backlog: some "uncatalogued key" findings
turn out to be a code-cleanup opportunity, not a documentation gap.

## 4d. The "easy" and "hard" buckets: mostly already done, not backlog

Went to implement the ~21 "easy" leaf keys identified above and found there
was almost nothing to do: 20 of the 21 were **already** fully and correctly
catalogued in `dict_entries_catalog.py`, verified by grepping each
`driver_path` directly rather than trusting the scanner's `absent_keys`
flag. `xMin`/`xMax`, for example, are already covered by a generic
`bathPotentialDomain.groundPatches.<patch>` /
`.surfaceCurrentPatches.<patch>` dynamic entry, complete with the real
mutual-exclusivity gotcha (a patch can't appear in both maps) already
documented as a constraint. What's flagging these as "absent" is purely the
scanner's regex being unable to resolve `<name>`/`<patch>`/`<scope>` dynamic
placeholders against literal strings in the C++ — the same blind spot as
`pvjKernel`. **Confirmed by direct experiment**: removed `scale`'s waiver
(already catalogued, dynamic path) to test the theory, and
`test_strict_dict_key_scanner_allowlist_is_current` immediately failed —
proving these waivers cannot be deleted just because the key is properly
catalogued, and must stay permanently. Restored before committing anything.

The "hard" bucket (`scale`, `set`, `range`, `node`, `ionicHeterogeneity`,
`verificationModel`, `coupling`, `activeTension`) was *also* mostly already
done, using exactly the generic-placeholder pattern that would have been
recommended for it. `constants` in particular is handled better than a flat
`DictEntry` ever could: there's a dedicated `active_tension_catalog.py`
module (mirroring `ionic_model_catalog.py`) listing each active-tension
model's actual valid constant/state names per model — the right design for
a key whose valid contents vary structurally by model, already built.

**Two genuine gaps were found and closed this round:**

- `constants.<constant_name>`, `initialStates.<state_name>`, and
  `preconditioningTime` (all flat siblings of `activeTensionModel` directly
  inside `<solver>Coeffs` — confirmed via `electroProperties()` resolving to
  `subDict(type + "Coeffs")` in `electroModel.C:109`, so there is no
  `activeTensionModel { }` sub-block, matching the pattern already
  documented on the `activeTensionModel` entry itself) were added to the
  `active_tension` group. `constants`' entry documents a real, non-obvious
  trap: `LandNiederer` re-derives several dependent constants (`AC_fPKA_TnI`,
  `AC_XSSS`, `AC_XWSS`, `AC_A`, `AC_PKAForceMultiplier`, `AC_k_uw`,
  `AC_k_ws`, `AC_k_wu`, `AC_k_su`, `AC_cds`, `AC_cdw`, `AC_ktm_block`)
  immediately after applying overrides, so overriding one of these specific
  derived names directly is a silent no-op — same "listed != overridable"
  caveat already established for the ionic constant catalog in
  `AGENT_GUIDE.md`. Verified these three waivers are now genuinely
  resolvable (not another `scale`-style permanent blind spot) via the same
  `unused_allowlist` check before removing them — the `unmatched_subdicts`
  matching path apparently *can* resolve container-level dynamic entries
  even though the `absent_keys` leaf-level path can't, which is why
  `preconditioningTime` (flat, no placeholder) came off cleanly while
  `constants`/`initialStates`' *other* remaining `absent_keys` waivers
  stayed (see next point).

- `pvjResistances` (and its siblings `pvjNodes`, `pvjLocations`, `points`,
  `rootNode`) — **resolved as permanently out of scope, not backlog.**
  These are real `dictionary::lookup`/`.get<T>()` reads
  (`conductionGraph.H`, `conductionSystemDomain.C`), structurally identical
  to any `electroProperties` key — just from a different file,
  `constant/<graphFile>` (e.g. `constant/purkinjeGraph`), which this catalog
  module has never addressed at all. That looked like a real gap requiring
  a second catalog module (`electroMechanicalProperties`-style), until
  checking whether driverFOAM ever actually needs to write one: grepped
  every Python file touching `graphFile`/`purkinjeGraph` in
  `openfoam_driver/` and found **no write site anywhere** — driverFOAM only
  ever references a graph file by name via the already-catalogued
  `graphFile` key. The graph file itself is a generated input asset (e.g.
  via `generatePurkinjeTree`), not something driverFOAM constructs or
  mutates. Since the catalog's entire purpose is documenting keys subject
  to *programmatic override*, and nothing here is ever overridden, there is
  no gap to close — a second catalog module would document a capability
  nobody uses. `dict_key_allowlist.json`'s own description was updated
  (finding 5) to document this reasoning permanently, alongside the
  existing false-positive/out-of-scope/blind-spot findings.

Net effect of this round: two real leaf-key gaps closed
(`constants`/`initialStates`/`preconditioningTime`), one apparent gap
resolved by deciding not to build a second catalog module
(`pvjResistances` and siblings), zero remaining actionable backlog from
either the "easy" or "hard" buckets. What's left in `absent_keys`/
`stale_paths`/`unmatched_subdicts` is, as far as this audit found, either a
permanent scanner blind spot or a deliberate scope boundary — not
undocumented work.

## 5. How to verify progress as this gets worked through

After adding a `DictEntry` for any of these keys, confirm the waiver is now
redundant and remove it in the same change:

```bash
cd applications/scripts/driverFoam
source .venv/bin/activate
python3 -m pytest openfoam_driver/tests/core/test_strict_planning.py \
    -k allowlist_is_current -q
```

If it still fails after removing the waiver line, the `DictEntry`'s
`driver_path` doesn't match what the scanner actually resolves (dynamic
path, subdict nesting, or a placeholder scope) — that's the pattern seen
with `pvjKernel`/`pvjCouplingScheme`/`pvjRadius`, which stay in the
allowlist despite being fully cataloged. Don't force-remove those waivers
without first confirming `unused_allowlist` actually flags them as
redundant.
