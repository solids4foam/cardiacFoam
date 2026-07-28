# Unify Purkinje `outputVariables` with the ionic-model registry

Deferred plan. Apply against a fresh `main` (the current branch is not
the right base — this is a future-work design note, not an in-flight
patch). Captured 2026-04-27 from a session that established the
present per-consumer split for `outputVariables` in the agent docs.

## Motivation

The Purkinje (1-D conduction-system) `outputVariables` block accepts
a hardcoded six-name set:
`Vm, Iion, activationTime, Icoupling, IcouplingSource, IcouplingCurrent`
([`conductionSystemDomain.C:253-289`](../src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C)).
Anything else is silently dropped — the `else` branch on line 290 has
no warning.

This is artificial. The Purkinje block already owns an `ionicModel`
([`conductionSystemDomain.C:220`](../src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C)),
so every 1-D node integrates a full ionic system. All of the model's
state and algebraic variables exist on every node and could be
exported per-node, exactly as the myocardium-side
`outputVariables.ionic.export` does today via
[`ionicModelIO.C`](../src/genericWriter/ionicModelIO.C). A user who
wants to inspect (say) intracellular Ca along a Purkinje strand has no
way to do it today, even though the data is right there in memory.

The ionic-side path also gives us, for free:

- `Vm` / `Iion` aliasing
  ([`ionicVariableCompatibility.C:112-279`](../src/genericWriter/ionicVariableCompatibility.C))
  so the user-facing names are stable across ionic models that name
  their voltage state `V` (TNNP), `u` (BuenoOrovio), `cell_v` (Gaur).
- A runtime warning on unknown names
  ([`ionicModelIO.C:417`](../src/genericWriter/ionicModelIO.C)),
  instead of today's silent drop.

## The natural split

Two token kinds, two handlers.

| Token kind | Source of truth | Handler |
|---|---|---|
| Network-only — `activationTime`, `Icoupling`, `IcouplingSource`, `IcouplingCurrent` | Computed by conduction-system code; no analog in the ionic registry | Built-in. Stay in `purkinjeModelIO`. |
| Ionic — `Vm`, `Iion`, plus any state/algebraic name registered by the active ionic model | The owned `ionicModel`'s `ioStateNames()` / `ioAlgebraicNames()` | Delegate to `ionicModelIO`'s filter — picks up aliasing + warnings for free. |

`conductionSystemDomain` would then call something like
`purkinjeModelIO::filterTokens(exportVars_, ionicModel)` and get back
two lists (network + ionic-resolved) instead of doing its own if-else
ladder.

## The dispatch home

[`purkinjeModelIO.C`](../src/genericWriter/purkinjeModelIO.C)
already exists (~212 lines) as the sibling to
[`ionicModelIO.C`](../src/genericWriter/ionicModelIO.C) (~767 lines).
It is the right place for the new dispatch:

```cpp
// purkinjeModelIO.H
namespace purkinjeModelIO
{
    struct ResolvedTokens
    {
        wordList networkTokens;   // activationTime, Icoupling*, IcouplingCurrent
        wordList ionicTokens;     // Vm, Iion (aliased), or any registered state/alg name
        wordList unknown;         // names that are neither — for one-shot warning
    };

    ResolvedTokens filterTokens
    (
        const wordList& userExport,
        const ionicModel& cellModel
    );
}
```

`conductionSystemDomain::initialiseOutputControls()` then becomes:

```cpp
const dictionary& ovDict = coeffsDict_.subOrEmptyDict("outputVariables");
const wordList raw = ovDict.getOrDefault<wordList>("export", {"Vm", "Icoupling"});
const auto resolved = purkinjeModelIO::filterTokens(raw, *ionicModelPtr_);
exportVars_   = resolved.networkTokens;     // existing path
ionicExport_  = resolved.ionicTokens;       // new
warnUnknown(resolved.unknown);              // new
```

The existing if-else ladder at lines 253-289 still handles the
network tokens. A new ionic loop appends per-node columns by reusing
the lookup machinery in `ionicModelIO::exportedFieldNamesRef`.

## The main tradeoff: column-count blowup

Per-node × per-variable explodes the column width of
`postProcessing/purkinjeNetwork.dat`. A 5 000-node TNNP Purkinje × 17
states = ~85 000 columns in one CSV. Three options, ranked by
operational impact:

1. **Per-variable file split** *(recommended)*. Write
   `purkinjeNetwork_Vm.dat`, `purkinjeNetwork_Cai.dat`, etc. — one
   file per requested ionic variable, each with `nNodes` columns.
   Network tokens stay in the single `purkinjeNetwork.dat`. Scales
   linearly; users drop unused files; downstream analysis tools
   (pandas, paraview) load only what they need.
2. **Node-subset filter**. Add a sibling key
   `outputNodes (14083 14127 ...)` that gates which nodes get
   exported per ionic variable. Keeps everything in one file. More
   intrusive on the dict surface but cheaper for large graphs where
   only a few probe nodes are interesting.
3. **Single fat file**. Keep one CSV; let the column count grow.
   Operationally painful for clinicians; only viable on small
   (~hundreds-of-nodes) Purkinje networks.

Recommendation: ship (1) as the v1, leave (2) as a follow-on switch
when someone hits a real ergonomic wall. A clinician should weigh in
before this lands — see `docs/agent/cardiacfoam-runbook.md` § Open
decisions, panel content.

## Catalog and documentation consequences

When this refactor lands:

- **`dict_entries.py`**: rewrite the Purkinje
  `outputVariables.export/debug` entries' descriptions to reflect the
  new name universe (network tokens + active-model state/algebraic
  names + `Vm`/`Iion` aliases, with silent-drop replaced by warning).
  Mirror the language already in the ionic-side entries.
- **`run-document-cookbook.md` § Convention 5**: collapse the two
  valid-name universes into one description per layer (network
  built-ins are still Purkinje-specific; the open ionic universe is
  shared with the myocardium side).
- **`validator-gaps.md` § 5**: drop the "no warning at all" line for
  the Purkinje path — that gap closes with this refactor.
- **`heartPurkinje-monodomain.json` example**: optionally extend the
  Purkinje export list to demonstrate an ionic state-variable name,
  e.g. `["Vm", "activationTime", "Icoupling", "Cai"]` for a Gaur
  Purkinje, to show off the new freedom.
- **Per-variable-file convention**: if option (1) wins, document the
  filename pattern (`purkinjeNetwork_<token>.dat`) somewhere
  user-facing; that becomes part of the run output contract.

## Open questions before implementing

1. **Per-variable file vs node-subset** — clinician input wanted on
   what they actually want from a 5 k-node Purkinje export.
2. **Aliasing scope** — do we want a `psi` alias on the Purkinje path
   when the conduction-system solver is `eikonalSolver1D`?
   Currently `psi` is the field name there; making it interchangeable
   with `activationTime` would smooth the eikonal-Purkinje story.
3. **Filter on-or-before construction** — `ionicModelIO` filters
   inside `exportedFieldNamesRef`, which is called per-write. Doing
   the same on the Purkinje path is cheap, but pre-resolving once at
   construction time (the way `conductionSystemDomain` already does)
   lets the warning fire once at startup instead of once per write
   step. Pre-resolution is preferable.
4. **Backwards compatibility** — the silent-drop behaviour is
   established. A hard error on unknown names would be cleaner but
   risks failing existing cases that have stale entries in their
   `outputVariables` blocks. A warning-only first release with a hard-
   error follow-up after one cycle is the safer rollout.

## How to start

1. Read the *current* state of the four files this touches —
   architecture may have moved since this note was written:
   `conductionSystemDomain.C`, `purkinjeModelIO.{H,C}`,
   `ionicModelIO.C`, `ionicVariableCompatibility.C`.
2. Sketch the `purkinjeModelIO::filterTokens` signature and the
   `ResolvedTokens` struct. PR-1 lands the helper with unit tests, no
   call sites.
3. PR-2 swaps the call site in `initialiseOutputControls` and adds
   the ionic loop in `openOutputFile` / the per-write code path.
4. PR-3 (option 1 above) implements per-variable file split.
5. PR-4 catalog/doc updates.

Each PR is independently revertible.
