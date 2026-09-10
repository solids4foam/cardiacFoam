# eikonalECG

Template-voltage surrogate ECG solver for eikonal activation-time fields.

## How it works

The eikonalECG model computes a pseudo-ECG signal without running a full
reaction-diffusion solve at each time step. Instead it uses:

1. **Activation time field** `ψ(x)` — produced by the eikonal solver.
2. **Precomputed tissue templates** `U(t)` — lookup tables of Vm vs. time
   for each tissue type (endocardial, mid-myocardial, epicardial), obtained
   from reference single-cell simulations.
3. **Analytical gradient**: the surrogate is `Vm(x,t) = U(t − ψ(x))`, i.e. the
   action potential waveform shifted in time by the local activation delay. The
   solver does **not** build a `Vm` field and differentiate it numerically; it
   applies the chain rule directly,

   ```
   grad(Vm) = -dU/ds(t − ψ) * grad(ψ)
   ```

   evaluating the template **derivative** `dU/ds` per cell as a weighted blend
   of the endocardial, mid-myocardial and epicardial templates. Templates are
   stored in mV and scaled to volts.
4. **Lead-field sum** over cells using precomputed anisotropic lead vectors
   `z = (σ · r) V_cell / |r|³`, where `r` is the vector from the electrode to
   the cell centre. The sum is reduced across MPI ranks.

   The same lead-vector construction is used by `pseudoECG`, so the two ECG
   paths are consistent by construction.

After the eikonal solve, ECG computation evaluates templates once per output
time step rather than per solver iteration.

## Verification

Setting the manufactured-template path replaces the tabulated blend with
`manufacturedTemplateDerivative(localTime)`, so the ECG functional can be driven
by a manufactured `dU/ds` for method-of-manufactured-solutions verification.

## Tissue templates (`tissueTemplates.H`)

The file `tissueTemplates.H` is a generated C++ header containing three
static scalar arrays:

| Array | Description |
|-------|-------------|
| `endoTimes` / `endoValues` | Endocardial Vm(t) template |
| `midTimes`  / `midValues`  | Mid-myocardial (M-cell) Vm(t) template |
| `epiTimes`  / `epiValues`  | Epicardial Vm(t) template |

Each template covers approximately **1 second** of simulated time
(10001 points, dt = 0.1 ms), downsampled by a factor of 100 from the
original single-cell ODE output to keep compiler memory usage low.

The template ionic model and stimulus conditions used to generate the file are
recorded in the local generation workflow. Regenerate the header from a
reference `singleCell` tutorial run when template provenance changes.

### When to regenerate

Regenerate `tissueTemplates.H` if:

- A new ionic model is adopted as the reference (current templates are
  TWorldcompactBatched-derived for human ventricular tissue).
- Tissue properties (conductivity, chi, cm) change significantly enough
  to alter AP duration.

For most production runs the existing templates are appropriate. If a run
needs templates that reflect its *own* ionic parameterisation instead —
drug effects, patient-specific overrides, a different ionic model entirely —
see `personalizedTemplates` below rather than regenerating this file.

## electroProperties configuration

Place the `eikonalECG` block inside `eikonalSolverCoeffs.ecgDomains`:

```c++
eikonalSolverCoeffs
{
    // ... conductivity, chi, cm, c0, stimulus ...

    ecgDomains
    {
        ECG
        {
            ecgSolver    eikonalECG;

            // REQUIRED — defines the Vm-response interpolation window.
            // start: beginning of AP window [s] (usually 0)
            // end:   must cover full AP duration:
            //        human ventricle ~ 0.3–0.5 s
            //        must not exceed template duration (~1.0 s)
            // deltaT: interpolation step [s] — 0.001–0.005 s is typical
            sampling
            {
                start  0;
                end    0.5;
                deltaT 0.001;
            }

            electrodePositions
            {
                lead_I   ( 0.1  0.0  0.0 );
                lead_II  ( 0.0  0.1  0.0 );
            }
        }
    }
}
```

### `sampling` block guidance

| Parameter | Typical value | Notes |
|-----------|--------------|-------|
| `start` | `0` | Start of Vm window; usually 0 |
| `end` | `0.3`–`0.5` | Must be ≥ AP duration; do not exceed ~1.0 s |
| `deltaT` | `0.001`–`0.005` | Finer than this has no benefit (template is at 0.1 ms) |

Setting `end` shorter than the AP duration truncates repolarisation in the
reconstructed ECG. Setting it longer than 1.0 s requests template samples
outside the tabulated range; the solver clamps to the last tabulated sample **silently** — no warning is issued.

### Transmural heterogeneity

If `ionicHeterogeneity` is configured in the solver coefficients, the ECG
model uses the same endo/mid/epi transmural weight field to blend between
the three templates. Without heterogeneity all cells use the endocardial
template.

### Personalized templates (`personalizedTemplates`)

By default the three templates above are the compiled, generic
`tissueTemplates.H` arrays — the same for every case, regardless of that
case's own ionic model or parameterisation. Adding a `personalizedTemplates`
block replaces them with templates generated for *this* case: it paces
three single cells (endocardial/mid-myocardial/epicardial anchors) with the
case's own ionic model and captures their Vm(t) response, then uses those
captured traces exactly where the compiled arrays would otherwise be used.

```c++
ecgDomains
{
    ECG
    {
        ecgSolver eikonalECG;

        personalizedTemplates
        {
            // REQUIRED. An eikonal-only case has no ionic model of its own
            // to inherit, so this authoritative single-cell parameterisation
            // must be supplied explicitly. Keep it equal to the ionic model
            // configuration (model, overrides, stimulus) the intended
            // monodomain run would use, or the personalized ECG will not
            // reflect that run's actual physiology.
            ionicModelConfig
            {
                ionicModel TWorldcompactBatched;   // exact registered
                                                     // runtime-selection-table
                                                     // name — e.g. plain
                                                     // "TWorldBatched" has no
                                                     // dictionary constructor
                                                     // and fatals
                tissue     endocardialCells;        // REQUIRED by every
                                                     // ionicModel constructor
                                                     // (ionicSelector::
                                                     // selectTissue() fatals
                                                     // without it) even
                                                     // though this generator's
                                                     // own configureIonicHeterogeneity()
                                                     // call supersedes it per
                                                     // anchor — the value here
                                                     // is a required
                                                     // placeholder, not a
                                                     // meaningful per-anchor
                                                     // selector
                solver     RKF45;
                absTol     1e-6;
                relTol     1e-4;

                // batchedSubsteps subdivides each outer `dt` step (below)
                // into this many internal integration steps. Batched models
                // default to 1 (no subdivision), which can be numerically
                // unstable for stiffer ionic models at typical outer dt
                // values — TWorld needs ~100 substeps at dt=1e-4s to avoid
                // a SIGFPE, matching the ~1e-6s step size known to work
                // elsewhere in this repo; BuenoOrovio tolerates dt=1e-4s
                // with the default of 1. Tune per ionic model.
                batchedSubsteps 100;

                ionicConstantOverrides
                {
                    global { /* model constants, e.g. GKr 0.5; */ }
                }

                singleCellStimulus
                {
                    stim_start      20;    // ms
                    stim_period_S1  1000;  // ms
                    stim_duration   1;     // ms
                    stim_amplitude  60.0;  // sized for the chosen ionic
                                           // model's own scale — too weak
                                           // an amplitude leaves the cell
                                           // sub-threshold and fatals the
                                           // generator's amplitude check
                    nstim1          10;    // silently overwritten to match
                                           // nBeats below
                    nstim2          0;     // MUST be exactly 0 — capture
                                           // needs one unambiguous final
                                           // S1 response, not S2 pacing
                }
            }

            nBeats   10;       // S1 beats (including the captured one) to
                                // pace each anchor to steady state
            duration 0.60;     // s, capture window from the final S1 onset
            dt       1e-4;     // s, ODE integration/capture step
        }

        sampling { /* unchanged, see above */ }
        electrodePositions { /* unchanged, see above */ }
    }
}
```

**Units, once, precisely:** `duration` and `dt` are SI seconds, matching
`sampling`'s own units. `singleCellStimulus`'s timing keys (`stim_start`,
`stim_period_S1`, `stim_duration`) are milliseconds, matching this
codebase's general single-cell stimulus convention. Captured Vm templates
are recorded in millivolts, matching the compiled `tissueTemplates.H`
convention; the mV→V conversion happens exactly once, at the point the
template derivative feeds the eikonal chain rule (mirroring the compiled
path's own conversion) — nowhere else in the personalized path does a unit
conversion occur.

**`transmuralBands`, `namedRegions`, or `cellZoneRegions` are supported.**
The `ionicHeterogeneity` block your solver coefficients already configure
is reused as-is (never re-parsed into a separate scheme). For
`transmuralBands`, three templates are generated (endo/mid/epi), exactly as
before. For `namedRegions`, one template is generated per entry in the
mode's `regions` sub-dictionary, paced at each region's own representative
field value, and blended per cell via the same `namedRegionWeightsAt()`
weighting the monodomain path uses — with the same field name, transition
width, mode, and smoothing keys `transmuralBands` already reads. For
`cellZoneRegions`, one template is generated per entry in the mode's
`regions` sub-dictionary (each naming a mesh `cellZone`), with a crisp
(unblended) per-cell weight of 1.0 for cells in that zone — this mode has
no `field`/`transitionWidth`/`smoothing`/`transitionMode` concept at all,
matching how the monodomain path already treats it. `gradientAxes` is not
supported by any mode; using `personalizedTemplates` with it is rejected
(see timing note below).

**Rejected at `eikonalECG` construction** (before any case/mesh setup, so
these specific misconfigurations fail immediately when the case's dict is
first parsed, not partway through a run): a missing `ionicModelConfig` or
its `ionicModel`/`singleCellStimulus`; `nBeats < 1`; non-positive
`duration`/`dt`; `duration` exceeding one S1 period
(`duration > 1e-3*stim_period_S1`); non-zero `nstim2`; and combining
`personalizedTemplates` with a manufactured-ECG verification configuration.

**Rejected at first `solve()`, not construction** (these two need the
mesh and `constant/electroProperties`, which don't exist yet when the
`eikonalECG` object is constructed): a missing `ionicHeterogeneity` block,
and `mode` other than `transmuralBands`, `namedRegions`, or
`cellZoneRegions`. For `cellZoneRegions` specifically, a mesh cell not
claimed by exactly one region's `cellZone` (unclaimed, or claimed by more
than one) also fatals here, since zone membership is a mesh property, not
a dict property, and so cannot be checked earlier. A case with one of
these problems will construct successfully and only fatal once the solver
actually starts solving.

**Without `personalizedTemplates`, nothing changes:** the compiled
`tissueTemplates.H` arrays, `transmuralBands`-only weighting, and every
other existing behavior of this solver are exactly as documented above —
this is an opt-in, additive path with no effect on cases that don't
configure it.

### Optional keys

```c++
report    true;   // Print sampling progress to stdout (default: true)
```

## Output

The solver writes `postProcessing/eikonalECG.dat` — a time-series file
with one column per electrode. Times correspond to the `sampling.deltaT`
grid, not the solver `deltaT`.
