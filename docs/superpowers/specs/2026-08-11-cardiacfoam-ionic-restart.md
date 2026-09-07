# cardiacFoam: correct ionic-state restart

Date: 2026-08-11
Status: proposed, not scheduled. **Optional** — see "Is this worth building?"
Scope: C++ solver work. Not part of the driverFOAM phase programme.

## The problem

Every ionic-model state and algebraic export is allocated
`IOobject::NO_READ, IOobject::AUTO_WRITE`
(`src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.C:531-553`),
and the string `restart` appears nowhere in `src/` or `applications/solvers/`.

Restarting from a later time directory restores `Vm`, `phiE`, `phiI`,
`activationTime` and `externalStimulusCurrent` — and **silently
re-initialises every gating variable, calcium concentration, and
active-tension state to its `initConsts` default.**

The result is a cell state that is internally inconsistent: `Vm` may sit
mid-repolarisation while its gates are at rest values. The solver does not
complain. The trajectory afterwards is wrong, and wrong in a way that looks
plausible on a plot.

## Why the one-line fix is a trap

The obvious change — flip `NO_READ` to `READ_IF_PRESENT` at
`myocardiumDomain.C:531-553` — would make things **worse**, and this is the
central finding of this spec.

`exportedFieldNames()` returns `variableExport_`, filtered against the model's
state and algebraic names (`src/ionicModels/ionicModel/ionicModel.H:856-866`).
That is a **user-selected subset**, and it returns an empty `wordList` when
`variableExport` is unset — which is the default. So:

- Usually **no** ionic state is on disk at all.
- When some is, it is whichever variables the user chose to export for
  plotting — a presentation choice, not a state vector.

Flipping the read flag would therefore restore the exported subset and default
the rest. A partial restore of a stiff coupled ODE system is worse than a full
reset: a full reset is obviously wrong and shows up immediately, while a
partial one produces a plausible-looking trajectory that is quietly incorrect.
That is the same failure class as a silent stale-resume — the thing this
programme exists to eliminate.

## Design

### 1. Separate checkpoint output from presentation output

`variableExport` stays exactly as it is: a presentation concern, controlling
what a user wants written for plotting and post-processing.

Add an independent checkpoint write that emits the **complete** state vector —
every state variable the model declares, not a subset — on write times. It
must be derived from the model's own declared state (`ioStateNames()` /
`ioNumStates()`), so a newly added model or variable is covered automatically
rather than by a hand-maintained list.

### 2. All-or-nothing restore

Read the complete set back `READ_IF_PRESENT`, and then apply a hard rule:

| On disk | Behaviour |
|---|---|
| none of the state present | initialise from `initConsts` — today's behaviour, unchanged |
| **all** of the state present | restore it; this is a true continuation |
| **some but not all** present | **refuse to start**, naming the missing variables |

The third row is the whole point. Partial state must fail loudly rather than
silently defaulting.

### 3. Backward compatibility

A case with no checkpoint data behaves exactly as today, so nothing regresses
and no existing tutorial changes. The feature is opt-in by the presence of
checkpoint output.

## Acceptance test — write this first

**The equivalence test is the deliverable that matters**, and it should be
written before any flag is touched, so it fails for a reason you can point at:

> Run a case `0 → T` in one go. Run the same case `0 → T/2`, then restart
> `T/2 → T`. Assert the final state agrees to solver tolerance.

Today this fails. When it passes, restart works — and no amount of reading the
code substitutes for it. Its absence is precisely why this defect sat
undetected.

Use a single-cell or small monodomain case so the test is cheap enough to run
in CI, and pick a `T/2` that lands mid-action-potential — restarting during the
plateau or upstroke is where inconsistent gate state shows up most sharply.
Restarting during diastole would likely pass even with the bug, and would be a
vacuous test.

## Related, unresolved

`bodyAndOrgansConductivity` is read `MUST_READ` at the literal start-time
instance (`extracellularPotentialDomain.C:495-506`), **not** through the
backward-searching `findFieldInstance` helper the conductivity tensors use
(`src/electroModels/core/conductivityFieldIO.C:105-138`). Since the utility
writes it once into `0/` and nothing rewrites it, a bath-bidomain restart from
a nonzero time may fail outright with a missing-file fatal.

This was **not established** — OpenFOAM's `fileOperation::filePath` fallback
behaviour could not be settled from source alone. The asymmetry between the two
code paths is the evidence. If this spec is implemented, test it empirically as
part of the same work; if it is not, the question is moot, because restart is
not supported anyway.

## Is this worth building?

**Only if restartable long runs are on the roadmap.** All twelve tutorials use
`startFrom startTime; startTime 0` and run end-to-end. If that is how
production work is done, then driverFOAM's `ionic_state_not_restored` warning
(Phase 2, Task 4) is sufficient on its own, and this spec should stay on the
shelf.

Build it when one of these becomes true:

- runs get long enough that losing progress to a crash is expensive;
- a study needs to branch from a common pre-computed state (e.g. pace to
  steady state once, then fan out over interventions) — this is the strongest
  case, since the alternative is re-paying the steady-state cost per branch;
- restart is wanted for time-parallel or checkpoint-restart scheduling on a
  cluster.

Until then the honest position is: restart is unsupported, and the driver says
so out loud.
