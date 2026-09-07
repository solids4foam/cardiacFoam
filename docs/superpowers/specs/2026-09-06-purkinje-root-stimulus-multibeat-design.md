# Multi-beat rootStimulus for the Purkinje conduction graph

## Problem

Electromechanical cases need a steady-state beating state before their
results mean anything. The intended workflow is to pace the heart for
several cycles (~5) through the Purkinje system and check that the
mechanical response has become beat-to-beat repeatable.

The Purkinje root stimulus cannot do this today. `readRootStimulus`
(`src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C:232`)
reads a **scalar** `startTime`, and `assembleAppliedCurrent` (`:589`) applies
it in a single window:

```cpp
if (t0 >= rootStartTime_ && t0 <= (rootStartTime_ + rootDuration_))
{
    appliedCurrent[rootNode_] += rootIntensity_;
}
```

One node, one start time, no cycle length, no repeat. Every `rootStimulus`
block in `tutorials/` uses a single scalar `startTime`. The graph fires once.

Neither existing mechanism substitutes for it:

- **`externalStimulus`** (`src/genericWriter/stimulusIO.C:314`) *is*
  list-capable — `stimulusLocationMinList` + `stimulusStartTimeList` gives N
  beats — but it injects current into 3D tissue boxes, bypassing the
  conduction system entirely. It cannot produce a Purkinje-driven activation
  sequence.
- **`singleCellStimulus`** (`nstim1`, `stim_period_S1`) *is* periodic, but it
  is read by the ionic model itself
  (`src/ionicModels/ionicModel/ionicModel.H:107`) and passed into the cell
  kernel inside `forAll(STATES_, integrationPtI)` (e.g.
  `src/ionicModels/Stewart/Stewart.C:205`). On the graph it would fire
  **every Purkinje node simultaneously**, with no propagation from the root
  — physiologically wrong, and silent. (In a reaction-diffusion myocardium
  domain the same protocol is explicitly fatal, via
  `myocardiumDomain::validateNoIonicStimulusInMonodomain()` at
  `src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.C:493`;
  the graph domain has no such guard.)

**Goal:** let `rootStimulus` fire the same node N times, so a Purkinje-driven
pacing train can drive an electromechanical case toward steady-state beating.

## Non-goals

- **No change to the active-tension preconditioning.** It was considered and
  deliberately left alone; see "Relationship to preconditioning" below.
- **No per-firing variation of `duration`, `intensity`, or `node`.** One
  node, one pulse shape, N firing times. Per-firing amplitude or width lists
  are not needed for a pacing train.
- **No protocol validation in the solver** — no ordering, spacing, overlap
  or sign rules on the firing times. See "Design principle this change
  follows".
- **No guard against `rootStimulus` + `singleCellStimulus`.** Deferred for
  separate discussion; see "Out of scope".
- **No `period`/`nBeats` form.** A 0D protocol shape is not imported into a
  spatial domain; see "Dict API".
- **No multi-root / multi-site pacing.** `rootNode_` stays a single node.
- **No change to `externalStimulus` or to the 3D tissue stimulus path.**

## Architecture

### Data model

`conductionSystemDomain.H:73` currently holds:

```cpp
scalar rootStartTime_;
scalar rootDuration_;
scalar rootIntensity_;
```

`rootStartTime_` becomes `List<scalar> rootStartTimes_`. `rootDuration_`,
`rootIntensity_` and `rootNode_` (`:67`) are unchanged — they apply to every
beat.

The `GREAT` sentinel currently used to mean "no root stimulus configured"
(`:236`, and tested at `:296` and `:589`) is retired: an **empty list** is
the natural representation of "never fires", and removes the sentinel
comparisons.

### Dict API

The 1D graph and the 3D tissue are both **spatial domains that exercise a
stimulus**, so they get the same protocol. The 0D `singleCellStimulus`
protocol is a different thing and deliberately does not extend to either.

`rootStimulus` therefore mirrors `externalStimulus`: a `...List` form for
multiple firings, with the existing scalar form as the single-firing case.

```
rootStimulus
{
    node          0;
    startTimeList (0 0.8 1.6 2.4 3.2);
    duration      0.002;
    intensity     500000.0;
}
```

**One node, N times.** `externalStimulus` sizes its list by *box* count, so
N firings at one site means repeating the box N times. The root stimulus
fires at a single user-identified `node`, so the list is purely temporal.
`duration`, `intensity` and `node` stay scalar and apply to every firing.

**Backward compatibility.** A scalar `startTime` expands to a one-element
list, exactly as `stimulusStartTime` is broadcast across boxes at
`src/genericWriter/stimulusIO.C:329`. Every existing tutorial keeps its
current behaviour.

**No `period`/`nBeats` form.** An earlier draft proposed accepting
`startTime` + `period` + `nBeats` alongside the list, borrowing the shape of
`singleCellStimulus`'s `stim_period_S1`/`nstim1`. That is dropped: it would
import a 0D protocol idea into a spatial domain, and the 3D path has no such
form. A fixed-BCL train is `startTimeList` with evenly spaced entries, which
driverFOAM generates. If a period form is ever wanted it should be added to
`externalStimulus` and `rootStimulus` together, not to one of them.

### Parsing (`readRootStimulus`, `:232`)

1. No `rootStimulus` block → `rootStartTimes_` empty, `duration`/`intensity`
   zero. Unchanged behaviour.
2. `startTimeList` present → read it directly.
3. Otherwise → read scalar `startTime` into a one-element list.

`startTime` remains required when there is no `startTimeList`, via the
existing `rsDict.get<scalar>("startTime")`, so a block with neither key
fails through OpenFOAM's own dictionary machinery rather than through added
checking.

**No protocol validation in the solver.** Start times are not checked for
sign, ordering, spacing, or window overlap. This matches `externalStimulus`
exactly, which validates only `checkSize` and `checkNonNegative` on
durations. Ordering and overlap are configuration-level questions, and
configuration validity is driverFOAM's responsibility — a case that reaches
the solver has already been through the strict pre-flight planner. Adding
solver-side protocol rules would put logic in C++ that has an owner
elsewhere, and would diverge the 1D path from the 3D one it is supposed to
match.

The `GREAT` sentinel meaning "no root stimulus configured" (`:236`, tested at
`:296` and `:589`) is retired; an empty list says the same thing without a
sentinel.

The existing `Info<<` summary (`:253`) reports the firing count and the
first/last times instead of a single `start=`.

### Application (`assembleAppliedCurrent`, `:589`)

```cpp
forAll(rootStartTimes_, beatI)
{
    const scalar tStart = rootStartTimes_[beatI];

    if (t0 < tStart || t0 > (tStart + rootDuration_))
    {
        continue;
    }

    appliedCurrent[rootNode_] += rootIntensity_;
}
```

This is `updateExternalStimulusCurrent`'s loop
(`src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.C:432`)
with the box test dropped, down to the `continue` and the `+=`. There is no
`break`: if two windows contain the same `t0` the current sums, which is the
3D domain's existing semantics for overlapping stimuli, not an oversight.
An empty list matches nothing, replacing the `GREAT` sentinel test.

### Which solvers the list actually reaches

This is the sharpest limit on the feature and it is not obvious from the
dictionary.

`assembleAppliedCurrent` has exactly one caller: `monodomain1DSolver`
(`src/electroModels/conductionSystemModels/monodomain1DSolver/monodomain1DSolver.C:74`,
consumed at `:175`). **Multi-firing pacing therefore only happens on the
monodomain 1D path.**

- **`eikonalSolver1D`** never calls `assembleAppliedCurrent`. It reads
  `domain.activationTime()` (`eikonalSolver1D.C:74`) and runs a single
  Dijkstra sweep, seeding from every node with `Tact >= 0` and propagating
  one wavefront. It is a steady-state activation-map formulation with no
  representation of a second beat at all. A `startTimeList` of five entries
  reaches it only through the one-time seeding below, and the four later
  entries have no effect whatsoever.
- **`restitutionEikonalSolver1D`** also never calls
  `assembleAppliedCurrent`. It does track beat-to-beat state
  (`lastActTime_`, `nextTact_`, refractoriness), but its repeat activations
  come from incoming waves and its own escape interval, not from the root
  stimulus list. Later entries in the list do not pace it either. It has no
  live case in the tree in any event — see "Verification".

So the honest scope is: `startTimeList` is a **monodomain 1D pacing
feature**. For both eikonal solvers it is silently inert past the first
entry. Under the design principle below, catching that mismatch —
`startTimeList` with more than one entry against an eikonal solver — is
driverFOAM's job, not a solver-side guard. It is called out here so the
limitation is recorded rather than discovered.

### Eikonal activation seeding (`initialiseState`, `:296`)

```cpp
if (rootStartTime_ < GREAT && rootStartTime_ <= SMALL)
{
    initialActivationTime[rootNode_] = rootStartTime_;
}
```

This seeds the root's activation time when the stimulus fires at t~0, for
the activation-time solvers. The condition becomes `min(rootStartTimes_)`,
guarded by a non-empty list. Because the parser imposes no ordering, the
minimum is taken explicitly rather than reading element 0.

Behaviour for the scalar case is unchanged. This is also the *only* route by
which `rootStimulus` reaches either eikonal solver, and it is inherently
single-valued — which is why later entries in a `startTimeList` cannot
produce further eikonal activations. Seeding them as initial conditions
would be wrong, not merely unsupported: they are not activated at t=0.

`readRootStimulus` runs at `:542`, before `initialiseState` at `:543`, so
the list is populated when the seeding reads it.

### driverFOAM catalog

`applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/dict_entries_catalog.py:1457`
carries four `rootStimulus` entries (`startTime`, `duration`, `intensity`,
`node`). One is added alongside them, same `phases={'stimulus'}` and
`dynamic_path=True`:

- `rootStimulus.startTimeList` — `value_kind='scalar_list'`, unit `s`

`scalar_list` is already an established kind in the catalog. The `startTime`
entry's description is updated to say it is the single-firing form and is
required only when `startTimeList` is absent.

`test_dict_entries.py:225` asserts each `rootStimulus.<sub>` is documented;
its `sub` list extends to the new key.

Because the solver deliberately does no protocol validation, driverFOAM is
where a malformed pacing train is caught. That work is not in this spec's
scope, but it is the reason the solver can stay simple, and it is why an
invalid protocol is not a silent failure: a case that cannot be planned
cannot be run.

## Relationship to preconditioning

The Land--Niederer active-tension models run a one-shot relaxation at
construction (`preconditionToRestingState`, e.g.
`src/activeTensionModels/LandNiedererTWorld/LandNiedererTWorld.C:311`):
`preconditioningTime` ms of the tension ODE at frozen Ca_i, lambda=1,
lambda_rate=0.

This was examined for conflict with a multi-beat run-in and found to be
**complementary, not competing**, so this spec changes nothing about it.

The shipped initial states
(`src/activeTensionModels/LandNiedererTWorld/LandNiedererTWorld_2025.H:223`)
are `Ca_TRPN=0, TmBlocked=1, XW=XS=ZETAS=ZETAW=0` — a degenerate corner of
state space, not an approximate resting state:

- `XU = 1 - TmBlocked - XW - XS = 0`, so `xb_uw = 0` — zero flux into the
  weak crossbridge state
- `blocking_Ca_factor = Ca_TRPN^(-nperm/2)` with `Ca_TRPN = 0` falls into the
  `1e-12` substitution and the cap at 100 (`:317`)
- `RATES[TmBlocked]` is therefore ~0

The crossbridge cascade is stalled until Ca_TRPN charges up, unblocking
TmBlocked, opening XU, filling XW then XS. With the constants at `:136-149`
and Ca_i = 2e-4 mM the slow mode is TmBlocked unblocking, rate
`ktm_unblock * Ca_TRPN^1.018` ~ 0.0026, i.e. a relaxation of order several
hundred time units. Without preconditioning the first beat lands on
artificially under-activated myocardium.

Preconditioning and a multi-beat run-in perform *the same relaxation*.
Preconditioning does it at frozen Ca_i in ODE-only time — and because
`TNNPinitConsts` (`src/ionicModels/TNNP/TNNP_2004.H:315`) ships **identical
initial STATES for endo, M and epi**, varying only conductances, the resting
Ca_i is uniform and the `uniform` fast path in
`preconditionToRestingState` collapses it to a single cell solve. The run-in
does the identical relaxation inside the full 3D monodomain plus solid
solve. Same physics, orders of magnitude apart in cost.

Preconditioning is therefore kept as a cheap reduction of the initial
artifact, and the pacing train handles what preconditioning cannot reach:
the ionic loading (Ca_SR, Na_i) that only accumulates under repeated beats.

Preconditioning's `lambda=1` was considered as a defect and rejected: at
construction no mechanics has run, `D` is zero, so lambda **is** 1 by
definition and no other value exists to pass. A case restarting from a
loaded state skips preconditioning entirely, because
`sequentialElectroMechanical` gates it on `!activeTensionRestarted`
(`src/electroMechanicalModels/sequentialElectroMechanical/sequentialElectroMechanical.C:155`).

## Coexistence with the 3D tissue stimulus

`externalStimulus` (3D tissue) and `rootStimulus` (1D graph) are independent
and this change does not affect that. They live in separate dictionaries
under separate domains — `externalStimulus` under the myocardium solver
coeffs, feeding `myocardiumDomain`'s source field; `rootStimulus` under
`conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs`, feeding
`conductionSystemDomain::assembleAppliedCurrent`. Neither reads the other,
and no code path couples or guards them, so a case may pace the tissue
directly, pace through the Purkinje root, or do both.

Two cases carry both blocks today, but neither is a running demonstration of
simultaneous stimulation, and neither paces more than once:

- `tutorials/template/constant/electroProperties` — both blocks with real
  values, but it is a documentation template.
- `tutorials/manufacturedSolutions/monodomain1D3D/constant/electroProperties`
  — both blocks present but **zeroed** (`duration 0`, `intensity 0` on
  each), because the manufactured source drives that case.

So the evidence in-tree is that both blocks parse and coexist in one
dictionary, not that simultaneous 1D+3D stimulation has been exercised. The
independence argument above rests on the code paths, not on those cases.

## Verification

There is no C++ unit-test framework in this repository — no `tests/` tree,
no CMake or Catch/GTest harness — and none is introduced here. Verification
uses the two mechanisms the repo already has: tutorial regressions for
solver behaviour, and driverFOAM's Python suite for the catalog. Cases are
driven through driverFOAM per `CLAUDE.md`, not through ad-hoc shell.

### The vehicle: purkinjeNiedererEtAl2011

All solver-side verification uses this one case. It already carries a
regression harness that runs three variants in a single script
(`regression/regressionTest.sh`), each against its own reference, and
`Allrun` selects a variant with `solver=<name>`, copying both
`constant/electroProperties.<name>` and `system/controlDict.<name>`.

Its variants are a 1D-solver x 3D-solver matrix, and **only the 1D choice
determines whether the firing list is reached at all**, because
`assembleAppliedCurrent` has exactly one caller (`monodomain1DSolver`):

| variant | 1D Purkinje | 3D myocardium | reaches the list |
|---|---|---|---|
| `monodomain` (default) | `monodomain1DSolver` | `monodomainSolver` | **yes** |
| `eikonal` | `eikonalSolver1D` | `eikonalSolver` | no — inert |
| `hybrid` | `eikonalSolver1D` | `monodomainSolver` | no — inert |

`hybrid` is the clarifying case: monodomain 3D tissue, but an eikonal 1D
Purkinje, and the list is still inert. The 3D solver has no bearing on it.
None of the three uses `restitutionEikonalSolver1D`, and neither does
anything else that is live: `restitutionEikonalSolver1D` is named in a few
dictionaries in the tree, but there is no working case exercising it. For
3D-to-1D coupling the Purkinje side is eikonal alone. The restitution solver
is therefore out of scope here in the strong sense — there is nothing to
regress against even if it were in scope.

`tutorials/heartSim3D-1D/monodomainHeart` was considered and rejected: it
has no `regression/` directory, so it is not a regression bar, and it is a
whole-heart case where this change would be a small signal in a large run.

### Five claims

**A. A scalar `startTime` fires at the same time, with the same current, as
today.** Fails if the one-element expansion mistimes the window, or if
retiring the `GREAT` sentinel changes when "no stimulus" is detected.

The `monodomain` variant is extended in place for claim B, so its reference
is regenerated and cannot also serve as the frozen bar. Coverage is instead:

- The `eikonal` and `hybrid` variants keep scalar `startTime 0.0` and their
  references stay untouched, exercising the scalar parse path and the
  seeding branch.
- **Before regenerating `purkinjeSlab.reference`, the extended run is
  checked against the *existing* reference at `t = 0.02`.** The first firing
  is unchanged at `0.01`, and a second firing at `0.3` cannot influence
  anything before it, so the state at 20 ms must still match the old
  reference within its existing tolerances. Only then is the reference
  regenerated at the new `endTime`. This preserves the backward-compatible
  bar rigorously without keeping a separate variant for it, and it must
  happen in that order — regenerating first would let the change validate
  itself.

**B. N list entries produce N root firings at the stated times.** The only
genuinely new observable. Fails if the window test is off by one, if only
the first entry fires, or if the list is never read.

The observable is **root-node voltage sampled at two times**. The
`monodomain` variant is extended in place — no new variant, no new
`controlDict`, no function object:

- `constant/electroProperties.monodomain` — `startTime 0.01` becomes
  `startTimeList (0.01 0.3)`.
- The same file's Purkinje-graph `outputVariables` — `export
  ( IcouplingSource )` becomes `export ( Vm IcouplingSource )`. This
  restores a default rather than inventing an output: the graph domain's
  built-in export list is `{"Vm", "IcouplingSource"}`
  (`conductionSystemDomain.C:317`) and the variant currently overrides Vm
  away. With it back, `postProcessing/purkinjeNetwork.dat` carries
  `node<N>_Vm_V` columns over time.
- `system/controlDict.monodomain` — `endTime 0.02` becomes `0.32`, so the
  second firing falls inside the run. At `deltaT 5e-5` that is 400 steps to
  ~6400. The tree already carries far longer list-driven runs on the same
  ionic model — `tutorials/electrophysiologyProtocols/rotorInstability` uses
  `stimulusStartTimeList (0.0 0.45 2.0)` with `endTime 4` at the same
  `deltaT`, i.e. 80,000 steps — so this stays cheap.
- `regression/purkinjeSlab.reference` — regenerated (see claim A for the
  ordering constraint), plus two rows recording root-node Vm shortly after
  each firing.
- `regression/regressionTest.sh` — one extractor reading a node's Vm column
  from `purkinjeNetwork.dat` at a given time. The file is already parsed
  there by `extractFinalPvjValue`; the existing extractors take the last row
  (`END {print value}`) or the first (`exit`), so selecting a row by its
  time column is the one new thing.

Each row records the value actually produced, with a tolerance, exactly like
the existing rows — it asserts no physiological interpretation. A break in
the second injection changes the second value, which is the whole point.

**Why root Vm rather than an activation time.** The two `activationTime`
fields differ, and only one of them latches:

- The **3D myocardium** field updates on every upward threshold crossing
  (`myocardiumDomain.C:461`, no first-activation guard), which is why
  probing it over time shows successive beats.
- The **1D graph** field latches: `monodomain1DSolver.C:208` sets it only
  `if (actTime[i] < 0.0 && Vm[i] >= activationVmThreshold_)`, so a second
  firing never updates it. The repeated-activation mechanism that exists —
  `acceptsRepeatedTerminalActivationTimes()` (`conductionSystemSolver.H:90`,
  default false; consumed at `conductionSystemDomain.C:666`) — covers only
  terminal/PVJ nodes and is overridden true only by
  `restitutionEikonalSolver1D`.

The 3D field would therefore work, but it answers a different question —
whether the second beat propagated through the PVJs and captured the tissue
— and at a 290 ms coupling interval on BuenoOrovio (APD ~300 ms) capture is
marginal. A refractory second beat would fail that check for physiological
reasons while the list mechanism is perfectly correct. Root-node Vm records
the injection itself and is independent of capture, so it tests the thing
this change actually alters.

Nor can the existing `final`-value checks show it: a single snapshot after
both firings does not distinguish two beats from one.

**C. A case with no effective `rootStimulus` still injects nothing.** Fails
if empty-list handling diverges from the old `GREAT` sentinel. Covered by
the `eikonal` and `hybrid` variants, both of which ship `intensity 0.0`.

**D. Eikonal seeding is unchanged for the scalar case.** Fails on `min()`
over an empty list, or a changed seeding condition. The `eikonal` variant
has `startTime 0.0`, so the seeding branch fires and the whole activation
map descends from it — `duration` and `intensity` being zero is irrelevant,
since the eikonal path never reads the applied current. Checked against
`regression/eikonalSlab.reference`, with `hybridSlab.reference` covering the
same seeding under a monodomain 3D tissue.

**E. The catalog documents `startTimeList` and export still works.**
`test_dict_entries.py:225` asserts every `rootStimulus.<sub>` key is
documented; its `sub` list extends to `startTimeList`. The existing catalog
export test must still pass.

C and D are one bar — `eikonalSlab.reference` and `hybridSlab.reference`
must not move. A and B share the extended `monodomain` variant, with A
checked against the old reference at `t = 0.02` *before*
`purkinjeSlab.reference` is regenerated. E is Python. No new variant and no
new `controlDict` are introduced.

## Design principle this change follows

Solver logic in C++ is kept to what is intrinsic to the numerics —
discretisation, matrices, and the failures those produce. Protocol and
configuration validity belong to driverFOAM, which owns dictionary mutation
and the strict pre-flight planner (see `CLAUDE.md`).

This is why the parser above checks nothing about the shape of the pacing
train. An invalid protocol is not a silent failure: it is a case that
driverFOAM will not plan, and therefore a case that never reaches the
solver. Pushing the same rules into C++ would duplicate an existing owner
and add branches that can drift from it.

## Out of scope

- **A `validateNoIonicStimulusInMonodomain()` equivalent for the graph
  domain.** `rootStimulus` together with `singleCellStimulus` would
  double-stimulate the Purkinje tree, the latter firing at every node at
  once. Deferred for a separate discussion, including whether the existing
  3D guard at
  `src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.C:493`
  should be generalised or moved to driverFOAM rather than duplicated. Not
  decided here, and nothing in this change depends on the outcome.
- **Multi-beat activation for the eikonal 1D graph.** See "Which solvers
  the list actually reaches": `eikonalSolver1D` is a single-sweep
  steady-state formulation and cannot represent a second beat, so a
  multi-entry `startTimeList` is inert there. Making the eikonal graph
  accept several beats is a **standalone task**, not a follow-up to this
  one — it is a change to the eikonal formulation itself, not to stimulus
  plumbing, and it does not share code with this change beyond the seeding
  call site. Until it exists, driverFOAM rejecting a multi-entry
  `startTimeList` against an eikonal solver is the interim guard, and that
  rule is also separate work.
- **A steady-state convergence criterion for a paced electromechanical
  run.** "Steady state" over ~5 cycles means beat-to-beat repeatability of a
  named mechanical quantity within a stated tolerance; it is **not** full
  ionic steady state, which for TNNP requires far more beats for Ca_SR and
  Na_i to settle. A run that claims steady state needs that metric defined
  first.
