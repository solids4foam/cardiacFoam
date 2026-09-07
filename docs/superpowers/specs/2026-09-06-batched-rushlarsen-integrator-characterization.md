# Batched Rush-Larsen / integrator selection: ionic vs active-tension — complete characterization

**Date:** 2026-09-06
**Status:** Characterization only. No code changes proposed here, by request.
**Scope:** how "which time-integration scheme applies to a given ODE state" is
decided in the batched (SoA, GPU-facing) ionic and active-tension model
families, and where the two families disagree.

Every claim below is anchored to a specific file and line, read directly.
Where something was not directly verified (e.g. whether plain-Euler is
numerically adequate at the substep counts actually used), that is said
explicitly rather than implied.

---

## 0. Two independent families, no shared base

`ionicModels/` and `activeTensionModels/` each implement their own version of
"batched cell-parallel time-stepping with an optional Rush-Larsen (RL) path
and an optional CUDA kernel." There is no shared header or base class between
them:

```
ionicModels/ionicModel/                    activeTensionModels/activeTensionModel/
  batchedIonicModel.H      (1307 lines)      batchedActiveTensionModel.H  (328 lines)
  batchedIonicCore.H                         batchedActiveTensionCore.H
  batchedIonicCoreCUDA.H                     batchedActiveTensionExecution.H
  batchedKernelExecution.H                   (no CUDA-core split file)
  batchedRushLarsenEntry.H                   (no dispatch-table equivalent)
  ionicModelGPU.H
  ionicHeterogeneity*.H
  gpuMath.H
```

`batchedActiveTensionCore.H` includes only `scalarField.H`; it does not
reference anything under `ionicModels/`. Each family reimplemented the same
pattern independently.

---

## 1. What Rush-Larsen requires, precisely

RL applies to a state whose ODE has (or can be locally approximated by) the
form

```
dy/dt = (y_inf - y) / tau
```

with `y_inf` ("steady state") and `tau` (relaxation time constant) evaluated
at the start of the substep and held fixed across it. For that form, RL
integrates the substep exactly:

```
y(t+dt) = y_inf + (y(t) - y_inf) * exp(-dt/tau)
```

which is unconditionally stable and exact for any `dt`, whereas explicit
Euler on the same state is only conditionally stable/accurate and requires
`dt` small relative to `tau`. This is why RL exists at all for gating-type
variables with fast time constants: it is not an optional refinement, it is
what keeps those specific states well-behaved regardless of substep size.
States whose ODEs are not in this form do not benefit from RL and gain
nothing from having it "available."

Both families' RL-eligibility checks apply exactly this form:

```cpp
// batchedActiveTensionModel.H:224-226 (Foam::batchedActiveTensionModel)
&& tau > VSMALL
&& std::isfinite(steadyState)
&& std::isfinite(tau);
```
```cpp
// batchedRushLarsenEntry.H:163 (Foam::resolveSupportRushLarsenEntry)
if (tau <= minTau || !std::isfinite(tau)) return false;
```

---

## 2. Ionic models: per-state declaration is a static dispatch table

Each ionic batched model builds one `const std::array<batchedRushLarsenEntry,
NUM_STATES>` at namespace scope, once, e.g.
[`TNNPBatched.C:117-134`](../../../src/ionicModels/TNNPBatched/TNNPBatched.C):

```cpp
const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
TNNPRushLarsenDispatch = []()
{
    std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
    t.fill(Foam::rlNone());

    t[m]   = Foam::rlScalarAlgAndSupport(tau_m,   m_inf,   TNNP_BATCH_SUPPORT_tau_m,   TNNP_BATCH_SUPPORT_gInf_m);
    t[h]   = Foam::rlScalarAlgAndSupport(tau_h,   h_inf,   TNNP_BATCH_SUPPORT_tau_h,   TNNP_BATCH_SUPPORT_gInf_h);
    ...
    t[fCa] = Foam::rlNone();   // clamped-rate gate — RL returns false on both paths
    t[g]   = Foam::rlNone();   // clamped-rate gate — RL returns false on both paths

    return t;
}();
```

Every state is accounted for explicitly, including a comment explaining why
the non-eligible ones (`fCa`, `g` — clamped-rate gates, not a clean
exponential form) are excluded.

`batchedRushLarsenEntry` ([`batchedRushLarsenEntry.H:53-70`](../../../src/ionicModels/ionicModel/batchedRushLarsenEntry.H))
is a plain-data struct: for `tau` and for the steady state, it records a
*source* (`algebraic` / `constant` / `literal` / `support`) and an index into
that source. Building an entry (`rlScalarAlgAndSupport(...)`, `rlNone()`,
etc.) is a `constexpr` helper call — the table is fixed at compile time, not
computed per cell or per step.

### 2.1 Two-stage resolution: project once, resolve many times

The table is consumed in two separate passes:

1. **Projection**, once per model-time evaluation, over every entry —
   [`TNNPBatched.C:548-558`](../../../src/ionicModels/TNNPBatched/TNNPBatched.C):
   ```cpp
   for (const auto& entry : TNNPRushLarsenDispatch)
   {
       projectScalarRushLarsenEntryToSupport(entry, CONSTANTS_, algebraics, supportValues);
   }
   ```
   `projectScalarRushLarsenEntryToSupport` ([`batchedRushLarsenEntry.H:178-`](../../../src/ionicModels/ionicModel/batchedRushLarsenEntry.H))
   copies each entry's declared tau/steady-state source value (an algebraic,
   a constant, or a literal) into a flat `supportValues` buffer, once, for
   every state that wants it.

2. **Resolution**, per state, on demand —
   [`resolveSupportRushLarsenEntry`](../../../src/ionicModels/ionicModel/batchedRushLarsenEntry.H:142-176):
   ```cpp
   inline bool resolveSupportRushLarsenEntry(
       const batchedRushLarsenEntry& entry, const Constants&,
       const SupportValues& supportValues, const Scalar minTau,
       Scalar& steadyState, Scalar& tau)
   {
       if (entry.supportTauSource == RLSupportSource::none) return false;
       tau = supportValues[entry.supportTauIndex];
       if (tau <= minTau || !std::isfinite(tau)) return false;
       steadyState = supportValues[entry.supportGInfIndex];
       return true;
   }
   ```
   This single function is shared, unchanged, by all 12 ionic batched models
   (`AlievPanfilov`, `BuenoOrovio`, `Courtemanche`, `Fabbri`, `Gaur`,
   `Grandi`, `PerisYague`, `Stewart`, `TNNP`, `ToRORd_dynCl`, `Trovato`,
   `TWorld`). None of them re-derive "how to get tau/steadyState from an
   entry" — they only build the table.

Per-model overrides (`TNNPBatched::rushLarsenParametersFromHotPathSupport`,
[`TNNPBatched.C:589-`](../../../src/ionicModels/TNNPBatched/TNNPBatched.C:589))
are a thin one-liner: look up `TNNPRushLarsenDispatch[stateI]`, call the
shared resolver.

A per-cell overload exists too (`...FromHotPathSupport(cellI, stateI, ...)`,
[`batchedIonicModel.H:898-919`](../../../src/ionicModels/ionicModel/batchedIonicModel.H:898)),
so heterogeneous-tissue constants (`constants(cellI)` instead of the global
`CONSTANTS_`) can feed the same table/resolver pair.

---

## 3. Ionic models: the global switch runs inside one shared loop

`batchedIntegrator` (`"euler"` default, [`batchedIonicModel.H:187-190`](../../../src/ionicModels/ionicModel/batchedIonicModel.H:187))
becomes `kernel.useExplicitEuler` / `kernel.useRushLarsen`
([`batchedIonicModel.H:1051-1052`](../../../src/ionicModels/ionicModel/batchedIonicModel.H:1051)),
validated with a `FatalIOError` on any other value
([`batchedIonicModel.H:1062-1073`](../../../src/ionicModels/ionicModel/batchedIonicModel.H:1062)).

Critically, both integrator choices run through the **same** per-state loop.
`rushLarsenEnabledForState` ([`batchedIonicModel.H:609-638`](../../../src/ionicModels/ionicModel/batchedIonicModel.H:609)):

```cpp
bool rushLarsenEnabledForState(...) const
{
    if (!kernel.useRushLarsen) return false;
    const bool hasParameters = rushLarsenParametersFromHotPathSupport(...);
    return hasParameters && tau > VSMALL && isfinite(steadyState) && isfinite(tau);
}
```

is called from inside `buildPredictorState` ([line 640](../../../src/ionicModels/ionicModel/batchedIonicModel.H:640))
for **every** state, on **every** integrator setting. Setting
`batchedIntegrator euler` makes this function return `false` unconditionally
(first line), so every state falls back to the Euler branch inside the same
loop — but the loop, and the dispatch-table lookup, still execute. The table
content itself never changes; only whether its answer is honored does.

There is no model-level "does this model support RL at all" flag on the
ionic side — an ionic model whose table is all `rlNone()` and is asked for
`rushLarsen` would simply behave as pure Euler for every state, with no
error. (Not evaluated further here — no ionic model's table was audited for
this case.)

---

## 4. Active-tension models: per-state declaration is a hand-written branch

There is no table and no shared resolver. `rushLarsenParametersForCell` is a
virtual function, overridden per model, whose body is an `if` chain over
state indices, with the tau/steady-state math inlined directly:

[`NashPanfilovBatched.C:146-165`](../../../src/activeTensionModels/NashPanfilovBatched/NashPanfilovBatched.C:146):
```cpp
bool Foam::NashPanfilovBatched::rushLarsenParametersForCell(
    const label cellI, const label stateI,
    const scalarUList& stateValues, const scalarUList& rateValues,
    const scalarUList& algebraicValues, scalar& steadyState, scalar& tau) const
{
    if (stateI == ::Ta)
    {
        const scalar eFactor = algebraicValues[::AV_e];
        if (eFactor > VSMALL)
        {
            tau = 1.0 / eFactor;
            steadyState = CONSTANTS_[AC_kTa] * algebraicValues[::AV_u];
            return true;
        }
    }
    return false;
}
```

[`GoktepeKuhlBatched.C:148-167`](../../../src/activeTensionModels/GoktepeKuhlBatched/GoktepeKuhlBatched.C:148)
is structurally identical (same `Ta`-only branch, same formula shape),
written independently rather than shared.

`rushLarsenParametersForCell` reads `algebraicValues` directly — there is no
separate projection/support-buffer stage analogous to ionic's
`projectScalarRushLarsenEntryToSupport` / `supportValues`.

`LandNiedererBatched`, `LandNiedererTWorldBatched`, `LandNiederer`, and
`LandNiedererTWorld` do not override `rushLarsenParametersForCell` at all —
the base class's implicit "no RL" applies to every one of their states.

---

## 5. Active-tension models: the global switch forks into two separate functions

`batchedIntegrator` (`"euler"` default,
[`batchedActiveTensionModel.C:50`](../../../src/activeTensionModels/activeTensionModel/batchedActiveTensionModel.C:50))
becomes `useRushLarsen_`, validated with `FatalIOError` on any other value
([`batchedActiveTensionModel.C:78-83`](../../../src/activeTensionModels/activeTensionModel/batchedActiveTensionModel.C:78)),
plus a model-level guard not present on the ionic side: `FatalIOError` if
`rushLarsen` is requested but `providesRushLarsenParameters_` is `false`
([`batchedActiveTensionModel.C:86-93`](../../../src/activeTensionModels/activeTensionModel/batchedActiveTensionModel.C:86)).
`providesRushLarsenParameters` is passed by each model's constructor:
`true` for `GoktepeKuhlBatched` ([`GoktepeKuhlBatched.C:65`](../../../src/activeTensionModels/GoktepeKuhlBatched/GoktepeKuhlBatched.C:65))
and `NashPanfilovBatched` ([`NashPanfilovBatched.C:64`](../../../src/activeTensionModels/NashPanfilovBatched/NashPanfilovBatched.C:64)),
omitted (defaults `false`) for both `LandNiedererBatched`
([`LandNiedererBatched.C:99`](../../../src/activeTensionModels/LandNiedererBatched/LandNiedererBatched.C:99))
and `LandNiedererTWorldBatched`
([`LandNiedererTWorldBatched.C:178`](../../../src/activeTensionModels/LandNiedererTWorldBatched/LandNiedererTWorldBatched.C:178)).

Unlike the ionic side, the two integrator choices do **not** share a loop.
`BatchedTensionExecutor::solveCells`
([`batchedActiveTensionExecution.H:140-148`](../../../src/activeTensionModels/activeTensionModel/batchedActiveTensionExecution.H:140))
picks between two entirely different functions:

```cpp
if (useRushLarsen)
{
    backend_.buildPredictorState(cellI, dtSubstep, scratch);
    backend_.applyCorrectorState(scratch);
}
else
{
    backend_.applyExplicitEulerStep(dtSubstep, scratch);
}
```

`applyExplicitEulerStep` ([`batchedActiveTensionModel.H:205-211`](../../../src/activeTensionModels/activeTensionModel/batchedActiveTensionModel.H:205)):

```cpp
void applyExplicitEulerStep(const scalar dtSubstep, CellScratch& scratch) const
{
    forAll(scratch.stateValues, i)
    {
        scratch.stateValues[i] += dtSubstep * scratch.rateValues[i];
    }
}
```

steps every state with plain Euler and **never calls
`rushLarsenParametersForCell`**. Only `buildPredictorState`
([`batchedActiveTensionModel.H:213-238`](../../../src/activeTensionModels/activeTensionModel/batchedActiveTensionModel.H:213))
consults it, per state, inside a `forAll` loop — and `buildPredictorState` is
only ever reached when `useRushLarsen_` is `true`.

Consequence: on the active-tension side, "Euler mode" is not "the RL-aware
loop with RL declined" (as on the ionic side) — it is a different function
that is structurally blind to whatever the per-state dispatch would have
said.

---

## 6. Per-model census: what actually declares RL support today, and why

| Model | `providesRushLarsenParameters` | States with RL declared | Basis |
|---|---|---|---|
| `GoktepeKuhlBatched` | `true` | `Ta` (its only state) | `dTa/dt = e·(kTa·u − Ta)` — exact linear relaxation |
| `NashPanfilovBatched` | `true` | `Ta` (its only state) | same form, same formula |
| `LandNiedererBatched` | `false` (default) | none | not overridden |
| `LandNiedererTWorldBatched` | `false` (default) | none | not overridden |
| `LandNiederer` (non-batched) | n/a | n/a | non-batched model, no RL concept applies |
| `LandNiedererTWorld` (non-batched) | n/a | n/a | same |

### 6.1 The Land 2017 model's 7 states, checked against the RL form

State indices from [`LandNiederer_2017Names.H:32-40`](../../../src/activeTensionModels/LandNiederer/LandNiederer_2017Names.H:32):
`XS, XW, TRPN, TmBlocked, ZETAS, ZETAW, Cd`. Their rate equations, read from
[`LandNiederer_2017.H:321-414`](../../../src/activeTensionModels/LandNiederer/LandNiederer_2017.H:321):

- **`TRPN`** (line 346): `d(trpn)/dt = koff·(Cai/ca50)^n·(1−trpn) − koff·trpn`.
  Expands to `A − B·trpn` with `A, B` depending only on `Cai`, `ca50`, `koff`
  (all frozen across a substep) — linear in the single state `trpn`.
- **`TmBlocked`** (line 372): `= ktm_block·blocking_factor·XU − ktm_unblock·trpn^(nperm/2)·tmBlocked`,
  where `XU = 1 − tmBlocked − xw − xs` (line 282). Substituting `XU`, this is
  linear in `tmBlocked` given `xw`, `xs` frozen.
- **`ZETAS`** (line 386): `= A·dlambda_dt − cds·zetas` — linear in `zetas`
  given the forcing `dlambda_dt` frozen.
- **`ZETAW`** (line 390): same form, linear in `zetaw`.
- **`Cd`** (line 413): `= par_k·(C − cd)/eta` = `(par_k/eta)·(C − cd)` —
  linear in `cd` given `C = lambda − 1` frozen.
- **`XS`** (line 321) and **`XW`** (line 326): each is linear in itself given
  the *other* frozen — e.g. `d(xs)/dt = k_ws·xw − xs·(k_su + gamma_rate)` —
  but the two states appear in each other's equation (via `XU` and via the
  `xb_ws`/`xb_uw` cross terms), so they are a **coupled linear pair**, not
  two independent single-state relaxations.

So 5 of the 7 states (`TRPN`, `TmBlocked`, `ZETAS`, `ZETAW`, `Cd`) are, given
frozen coupling terms, individually in the exact `dy/dt = (y_inf − y)/tau`
form the RL check requires. `XS`/`XW` are linear but mutually coupled. None
of this is implemented — `rushLarsenParametersForCell` is not overridden in
any Land-family file, so all 7 states always take the Euler branch,
regardless of `batchedIntegrator`.

This document does not assess whether Euler is numerically adequate for
these states at the substep counts actually configured (`batchedSubsteps`)
in any given case — that would require a convergence study, not a code read.

---

## 7. Defaults and real-world exposure

- Both families default `batchedIntegrator` to `"euler"`
  ([`batchedIonicModel.H:189`](../../../src/ionicModels/ionicModel/batchedIonicModel.H:189),
  [`batchedActiveTensionModel.C:50`](../../../src/activeTensionModels/activeTensionModel/batchedActiveTensionModel.C:50)).
- Every occurrence of `batchedIntegrator` in `tutorials/` is commented out:
  `tutorials/template/constant/electroProperties:265`,
  `tutorials/NiedererEtAl2011/electroMechanicalNiedererEtAl2011/constant/electroMechanicalProperties:34`,
  `tutorials/electromechanicalIdealizedModel/electroMechanicalIdealizedModel/constant/electroMechanicalProperties:34`.
  No tutorial in the repository sets it to `rushLarsen`.
- driverFOAM's catalog entry
  ([`dict_entries_catalog.py:343-352`](../../../applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/dict_entries_catalog.py:343))
  lists `batchedIntegrator` as an `enum` field with `applicable_when`
  including `LandNiedererBatched` and `LandNiedererTWorldBatched` (models
  with zero RL support — requesting `rushLarsen` on them now raises
  `FatalIOError`) alongside `NashPanfilovBatched`/`GoktepeKuhlBatched`
  (models where `Ta` is only actually RL-integrated if the case dict
  explicitly opts in; the catalog entry's own description does not
  distinguish "optional refinement" from "the only path this state's
  dispatch entry was written for").

Net effect measured from the repository as it stands: no case in this repo
currently runs any active-tension state through Rush-Larsen, including the
two models built with it, because the global default is `euler` and nothing
overrides it.

---

## 8. Side-by-side summary

| | Ionic (`ionicModels/`) | Active tension (`activeTensionModels/`) |
|---|---|---|
| Per-state declaration | static `std::array<batchedRushLarsenEntry, NUM_STATES>`, built once, all states listed explicitly (incl. `rlNone()` with a reason comment) | virtual function per model, `if (stateI == ::X)` branches, no exhaustive listing |
| Resolution logic | one shared `resolveSupportRushLarsenEntry`, reused by all 12 models | inlined per model, duplicated between `GoktepeKuhl`/`NashPanfilov` |
| Data path | two-stage: project sources into a `supportValues` cache once, resolve from that cache per state | single-stage: reads `algebraicValues` directly, no cache layer |
| Per-cell heterogeneity | explicit `...ForCell(cellI, ...)` overload feeding `constants(cellI)` into the same resolver | `rushLarsenParametersForCell` already takes `cellI`, but has no separate global/per-cell split since there's only one function |
| Global switch's effect | gates the return of `rushLarsenEnabledForState`, called from *inside* the one loop that also runs on `euler` | selects between two *different* functions; the `euler` function never calls the per-state dispatch |
| Model-level "no RL at all" guard | none — an all-`rlNone()` table would silently behave as Euler if `rushLarsen` were requested | `providesRushLarsenParameters_` + explicit `FatalIOError` if `rushLarsen` requested without it |
| States with RL declared today | per model (e.g. TNNP: 9 of 12 gates; `fCa`/`g` excluded with a comment) | `Ta` only, in the 2 single-state models; 0 of 7 in the 4-model Land family despite 5 states being in the eligible form |

---

## 9. Explicitly not covered here

- Whether Euler is numerically failing (accuracy/stability) for any specific
  state at the substep counts used in any real case — not measured.
- The ionic side's own gap (no model-level "unsupported RL" guard) was not
  traced further (e.g. whether any ionic model's table is fully `rlNone()`).
- No CUDA (`.cu`) file implements time integration for either family — both
  kernels compute rates/algebraics only; the substep loop and the integrator
  choice described above are host-side. This document does not characterize
  the `.cu` kernels themselves beyond that fact.
