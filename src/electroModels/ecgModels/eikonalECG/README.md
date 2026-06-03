# eikonalECG

Template-voltage surrogate ECG solver for eikonal activation-time fields.

## How it works

The eikonalECG model computes a pseudo-ECG signal without running a full
reaction-diffusion solve at each time step. Instead it uses:

1. **Activation time field** `ψ(x)` — produced by the eikonal solver.
2. **Precomputed tissue templates** `U(t)` — lookup tables of Vm vs. time
   for each tissue type (endocardial, mid-myocardial, epicardial), obtained
   from reference single-cell simulations.
3. **Surrogate reconstruction**: for every cell, `Vm(x,t) = U(t − ψ(x))`.
   The action potential waveform is simply shifted in time by the local
   activation delay.
4. **Pseudo-ECG integration** over the reconstructed Vm field using
   precomputed lead vectors.

This makes ECG computation essentially free after the eikonal solve — the
templates are evaluated once per output time step, not per solver iteration.

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

The template ionic model and stimulus conditions used to generate the file
are recorded in the generation script (not committed — contact the author
or regenerate from a reference singleCell tutorial run).

### When to regenerate

Regenerate `tissueTemplates.H` if:
- A new ionic model is adopted as the reference (current templates are
  TNNP-based for human ventricular tissue).
- Tissue properties (conductivity, chi, cm) change significantly enough
  to alter AP duration.

For most production runs the existing templates are appropriate.

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

Setting `end` shorter than the AP duration will truncate repolarisation
in the reconstructed ECG. Setting it longer than 1.0 s will cause a
lookup out of range (the solver will warn and clamp).

### Transmural heterogeneity

If `ionicHeterogeneity` is configured in the solver coefficients, the ECG
model uses the same endo/mid/epi transmural weight field to blend between
the three templates. Without heterogeneity all cells use the endocardial
template.

### Optional keys

```c++
report    true;   // Print sampling progress to stdout (default: true)
```

## Output

The solver writes `postProcessing/eikonalECG.dat` — a time-series file
with one column per electrode. Times correspond to the `sampling.deltaT`
grid, not the solver `deltaT`.
