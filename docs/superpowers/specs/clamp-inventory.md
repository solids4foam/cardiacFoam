# SoA-Euler clamp inventory — TNNP, ToRORd_dynCl, Gaur

**Generated:** 2026-05-12
**Purpose:** Per-model list of `(stateName, clampType, bound)` tuples that the next-session SoA-shim factoring will consume. Per `validation_strategy.md §2.2` of the upstream cardiacFoam GPU port: clamps are essential for detailed CellML-generated models because explicit Euler can drift gates outside `[0,1]` or concentrations negative within one substep, and the next call's `log(neg)` or `1/(1+exp(...))` will then NaN.

Classification rubric (recap):

- **gate** — probability/fraction in `[0,1]`, typically Hodgkin–Huxley `(x_inf − x)/tau_x` form with `x_inf = 1/(1+exp(...))`. Clamp `[0.0, 1.0]`.
- **concentration** — ion concentration `> 0`, consumed by Nernst `log(C_o/C_i)` or appears in a denominator. Clamp `>= SMALL`.
- **voltage** — `V` / `Vm`. No clamp.
- **other / structural** — Markov-chain components, CaMK/phosphorylation trapping fractions, dynamic SR-release variables (Jrel-style states); each handled with a per-state rationale below.

The agent that will *use* this list is responsible for picking the exact `SMALL` numeric value (the OpenFOAM `SMALL` constant) and for deciding whether Markov state per-component clamps need a follow-up renormalisation step.

---

## TNNP (NUM_STATES = 17)

Source: `src/ionicModels/TNNP/TNNP_2004Names.H`, `TNNP_2004.H`.

| Index | State | Classification | Clamp |
|---|---|---|---|
| 0 | V | voltage | none |
| 1 | K_i | concentration | `>= SMALL` |
| 2 | Na_i | concentration | `>= SMALL` |
| 3 | Ca_i | concentration | `>= SMALL` |
| 4 | Xr1 | gate | `[0.0, 1.0]` |
| 5 | Xr2 | gate | `[0.0, 1.0]` |
| 6 | Xs | gate | `[0.0, 1.0]` |
| 7 | m | gate | `[0.0, 1.0]` |
| 8 | h | gate | `[0.0, 1.0]` |
| 9 | j | gate | `[0.0, 1.0]` |
| 10 | d | gate | `[0.0, 1.0]` |
| 11 | f | gate | `[0.0, 1.0]` |
| 12 | fCa | gate | `[0.0, 1.0]` |
| 13 | s | gate | `[0.0, 1.0]` |
| 14 | r | gate | `[0.0, 1.0]` |
| 15 | Ca_SR | concentration | `>= SMALL` |
| 16 | g | gate | `[0.0, 1.0]` |

Notes:
- `Ca_i` and `Ca_SR` are consumed in `i_rel`, `g_inf`, and Nernst-style `E_Ca = (RT/2F) log(Ca_o/Ca_i)` — must stay strictly positive.
- `g` and `fCa` are gate-like with conditional update (`RATES = 0` if drifting upward at depolarised V); still bounded in `[0,1]` mathematically.

## ToRORd_dynCl (NUM_STATES = 45)

Source: `src/ionicModels/ToRORd_dynCl/ToRORd_dynCl_2023Names.H`, `ToRORd_dynCl_2023.H`.

| Index | State | Classification | Clamp |
|---|---|---|---|
| 0 | V | voltage | none |
| 1 | CaMKt | other (CaMK trapping fraction, `[0,1]`) | `[0.0, 1.0]` |
| 2 | Nai | concentration | `>= SMALL` |
| 3 | Nass | concentration | `>= SMALL` |
| 4 | Ki | concentration | `>= SMALL` |
| 5 | Kss | concentration | `>= SMALL` |
| 6 | Cass | concentration | `>= SMALL` |
| 7 | Cansr | concentration | `>= SMALL` |
| 8 | Cajsr | concentration | `>= SMALL` |
| 9 | Cai | concentration | `>= SMALL` |
| 10 | Cli | concentration | `>= SMALL` |
| 11 | Clss | concentration | `>= SMALL` |
| 12 | INa_m | gate | `[0.0, 1.0]` |
| 13 | INa_h | gate | `[0.0, 1.0]` |
| 14 | INa_j | gate | `[0.0, 1.0]` |
| 15 | INa_hp | gate (phosphorylated h) | `[0.0, 1.0]` |
| 16 | INa_jp | gate (phosphorylated j) | `[0.0, 1.0]` |
| 17 | INaL_mL | gate | `[0.0, 1.0]` |
| 18 | INaL_hL | gate | `[0.0, 1.0]` |
| 19 | INaL_hLp | gate (phosphorylated) | `[0.0, 1.0]` |
| 20 | Ito_a | gate | `[0.0, 1.0]` |
| 21 | Ito_iF | gate | `[0.0, 1.0]` |
| 22 | Ito_iS | gate | `[0.0, 1.0]` |
| 23 | Ito_ap | gate (phosphorylated) | `[0.0, 1.0]` |
| 24 | Ito_iFp | gate (phosphorylated) | `[0.0, 1.0]` |
| 25 | Ito_iSp | gate (phosphorylated) | `[0.0, 1.0]` |
| 26 | ICaL_d | gate | `[0.0, 1.0]` |
| 27 | ICaL_ff | gate | `[0.0, 1.0]` |
| 28 | ICaL_fs | gate | `[0.0, 1.0]` |
| 29 | ICaL_fcaf | gate | `[0.0, 1.0]` |
| 30 | ICaL_fcas | gate | `[0.0, 1.0]` |
| 31 | ICaL_jca | gate | `[0.0, 1.0]` |
| 32 | ICaL_ffp | gate (phosphorylated) | `[0.0, 1.0]` |
| 33 | ICaL_fcafp | gate (phosphorylated) | `[0.0, 1.0]` |
| 34 | ICaL_nca_ss | gate (Ca-binding occupancy, `[0,1]`) | `[0.0, 1.0]` |
| 35 | ICaL_nca_i | gate (Ca-binding occupancy, `[0,1]`) | `[0.0, 1.0]` |
| 36 | IKr_C1 | other (Markov state, see notes) | `[0.0, 1.0]` |
| 37 | IKr_C2 | other (Markov state, see notes) | `[0.0, 1.0]` |
| 38 | IKr_C3 | other (Markov state, see notes) | `[0.0, 1.0]` |
| 39 | IKr_I | other (Markov state, see notes) | `[0.0, 1.0]` |
| 40 | IKr_O | other (Markov state, see notes) | `[0.0, 1.0]` |
| 41 | IKs_xs1 | gate | `[0.0, 1.0]` |
| 42 | IKs_xs2 | gate | `[0.0, 1.0]` |
| 43 | Jrel_np | other (dynamic SR-release variable) | `[0.0, 1.0]` (see notes) |
| 44 | Jrel_p | other (dynamic SR-release variable, phosphorylated) | `[0.0, 1.0]` (see notes) |

## Gaur (NUM_STATES = 29)

Source: `src/ionicModels/Gaur/Gaur_2021Names.H`, `Gaur_2021.H`.

| Index | State | Classification | Clamp |
|---|---|---|---|
| 0 | cell_v | voltage | none |
| 1 | nai | concentration | `>= SMALL` |
| 2 | nass | concentration | `>= SMALL` |
| 3 | ki | concentration | `>= SMALL` |
| 4 | kss | concentration | `>= SMALL` |
| 5 | cai | concentration | `>= SMALL` |
| 6 | cai2 | concentration | `>= SMALL` |
| 7 | cass | concentration | `>= SMALL` |
| 8 | cansr | concentration | `>= SMALL` |
| 9 | cajsr | concentration | `>= SMALL` |
| 10 | cacsr | concentration | `>= SMALL` |
| 11 | I_Na_m | gate | `[0.0, 1.0]` |
| 12 | I_Na_h | gate | `[0.0, 1.0]` |
| 13 | I_Na_j | gate | `[0.0, 1.0]` |
| 14 | INaL_ml | gate | `[0.0, 1.0]` |
| 15 | INaL_hl | gate | `[0.0, 1.0]` |
| 16 | ICaL_d | gate | `[0.0, 1.0]` |
| 17 | ICaL_fca | gate | `[0.0, 1.0]` |
| 18 | IKr_xr | gate | `[0.0, 1.0]` |
| 19 | IKs_xs1 | gate | `[0.0, 1.0]` |
| 20 | IKs_xs2 | gate | `[0.0, 1.0]` |
| 21 | ITo_aa | gate | `[0.0, 1.0]` |
| 22 | CICR_Jrel2 | other (dynamic SR-release variable) | none (see notes) |
| 23 | CICR_Jrel1 | other (dynamic SR-release variable) | none (see notes) |
| 24 | CaMK_CaMKt | other (CaMK trapping fraction, `[0,1]`) | `[0.0, 1.0]` |
| 25 | CICR_tjsrol | other (countdown timer, see notes) | none |
| 26 | CICR_A | other (availability scaler, `[0, 100]`) | none (see notes) |
| 27 | ICaL_fs | gate | `[0.0, 1.0]` |
| 28 | ICaL_ff | gate | `[0.0, 1.0]` |

---

## Notes / ambiguous cases

- **TNNP `g` and `fCa`** — gate-like Hodgkin–Huxley variables but with a conditional `RATES = 0` rule (no upward drift while `V > -60` mV and `x_inf > x`). Still mathematically bounded in `[0,1]`; the standard `[0,1]` clamp is correct.

- **ToRORd_dynCl `CaMKt`** — CaMK trapping fraction in `[0,1]`. Appears as `(1 - CaMKt)` and bare `CaMKt` in `CaMKb` / `CaMKa`; a negative value would corrupt downstream phosphorylation fractions. Clamp `[0,1]`.

- **ToRORd_dynCl `IKr_C1..C3, I, O`** — five Markov states for the IKr channel. Each individual occupancy is in `[0,1]`, *and* mathematically `C1 + C2 + C3 + I + O = 1`. A per-state `[0,1]` clamp prevents NaN propagation but does not enforce the sum-to-1 constraint; small drift accumulates over many substeps. **Decision (recorded 2026-05-12):** **add the renormalisation pass** after the IKr Markov state update, before the next substep. Cost is ~15 FLOPs per cell per substep (5 reads + 4 adds + 1 div + 5 mults + 5 writes) against an evaluator that's ~10,000+ FLOPs per cell — well under 0.5% overhead. Cannot explode: the per-state `[0,1]` clamp already in this inventory handles negative-drift; renormalisation enforces the joint sum, improving long-time accuracy vs RKF45 without changing short-time correctness. Reference snippet for the SoA-shim implementer:

  ```cpp
  // After the IKr Markov state update, before the next substep:
  const double sum = STATES[IKr_C1] + STATES[IKr_C2] + STATES[IKr_C3]
                   + STATES[IKr_I]  + STATES[IKr_O];
  const double inv = 1.0 / sum;
  STATES[IKr_C1] *= inv;
  STATES[IKr_C2] *= inv;
  STATES[IKr_C3] *= inv;
  STATES[IKr_I]  *= inv;
  STATES[IKr_O]  *= inv;
  ```

  In SoA layout this becomes one fully-coherent pass over five contiguous state slots per cell. On GPU it's a single warp-wide reciprocal — negligible relative to the launch.

- **ToRORd_dynCl `ICaL_nca_ss`, `ICaL_nca_i`** — Ca-binding occupancy fractions, structurally in `[0,1]`. Treating as gates is safe.

- **ToRORd_dynCl `Jrel_np`, `Jrel_p`** — dynamic release-channel state evolving as `(Jrel_inf - Jrel) / tau_rel`. Initial values are tiny (~1e-21), and the steady-state magnitude is small. They are not strictly bounded in `[0,1]` mathematically (sign tracks the sign of `ICaL_ss`, which is negative → `Jrel_inf` is positive but small). **Default classification:** clamp `[0,1]` is a reasonable safety net to prevent runaway; however the implementer should confirm the physical sign convention before committing. Flagged for follow-up.

- **Gaur `CaMK_CaMKt`** — same rationale as ToRORd `CaMKt`; clamp `[0,1]`.

- **Gaur `CICR_Jrel1`, `CICR_Jrel2`** — release flux states; no natural `[0,1]` bound (initial values 9.71e-22 and 8.59e-10 respectively, dynamic range model-dependent). Classification: **other**; **no clamp** — drift tolerance must be assessed via the SoA-Euler validation harness rather than enforced statically. Flagged for follow-up.

- **Gaur `CICR_tjsrol`** — a "time-since" countdown variable initialised to 100 and reset to 0 on a release event (`diff = -tjsrol/0.001` while triggered, else `diff = 1`). Not a gate, not a concentration. **No clamp.**

- **Gaur `CICR_A`** — release-channel availability that decays from 100 to 0 after a release and is never directly consumed by a `log` or division. **No clamp**; explicit Euler drift here is a physical-realism concern, not a NaN concern.

- **Gaur ordering quirk** — `ICaL_fs` and `ICaL_ff` appear *after* the CICR/CaMK block in the enum (indices 27, 28) rather than next to `ICaL_d` and `ICaL_fca`. Pure cosmetic; classification is unaffected.

---

## Summary counts

| Model | total states | gates | concentrations | voltage | other |
|---|---|---|---|---|---|
| TNNP | 17 | 12 | 4 | 1 | 0 |
| ToRORd_dynCl | 45 | 32 | 10 | 1 | 2 (Jrel_np, Jrel_p) + 1 (CaMKt, listed under "other" above but clamped `[0,1]`) + 5 (IKr Markov, clamped `[0,1]`) = 8 |
| Gaur | 29 | 13 | 10 | 1 | 5 (CaMK_CaMKt clamped `[0,1]`; CICR_Jrel1, CICR_Jrel2, CICR_tjsrol, CICR_A unclamped) |

Compact view (gates + "other-clamped-as-gate" merged into the gate column for at-a-glance budget planning):

| Model | total | clamped `[0,1]` | clamped `>= SMALL` | unclamped (V) | unclamped (other) |
|---|---|---|---|---|---|
| TNNP | 17 | 12 | 4 | 1 | 0 |
| ToRORd_dynCl | 45 | 41 | 10 | wait — see below | — |
| Gaur | 29 | 14 | 10 | 1 | 4 |

Recount for ToRORd_dynCl (since the multi-category breakdown matters for the implementer):
- `[0,1]` clamps (gates + CaMKt + IKr Markov + nca + Jrel): 32 gate + 1 CaMKt + 5 IKr + 0 (nca already counted in 32) + 2 Jrel = **34** if Jrel-clamping is accepted, else **32** with Jrel deferred.
- Wait: re-tallying from the per-row table: rows with clamp `[0.0, 1.0]` = rows 1 (CaMKt), 12–35 (24 rows), 36–40 (5 rows), 41–42 (2 rows), 43–44 (2 rows) = 1+24+5+2+2 = 34. Concentration rows = 2–11 (10 rows). Voltage = 1 (row 0). 1+10+34 = 45. ✓

Final corrected compact view:

| Model | total | clamped `[0,1]` | clamped `>= SMALL` | voltage (no clamp) | other (no clamp) |
|---|---|---|---|---|---|
| TNNP | 17 | 12 | 4 | 1 | 0 |
| ToRORd_dynCl | 45 | 34 | 10 | 1 | 0 |
| Gaur | 29 | 14 | 10 | 1 | 4 |

The 4 unclamped "other" rows in Gaur are `CICR_Jrel1`, `CICR_Jrel2`, `CICR_tjsrol`, `CICR_A` — each flagged in the notes section above for the implementer to revisit during SoA-Euler validation.
