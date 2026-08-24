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

Setting `end` shorter than the AP duration truncates repolarisation in the
reconstructed ECG. Setting it longer than 1.0 s requests template samples
outside the tabulated range; the solver clamps to the last tabulated sample **silently** — no warning is issued.

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
