# ionicHeterogeneityProbe

Validates the `ionicHeterogeneity` parameter map for a `BuenoOrovio` ionic
model by sweeping synthetic transmural-distance samples and checking that the
resulting action potentials are physiologically consistent and vary smoothly
across the transmural gradient.

> [!NOTE]
> This utility supports only the `BuenoOrovio` ionic model.

## What it does

1. Reads `constant/electroProperties` to build the ionic-model dictionary.
2. Reads optional overrides from `constant/ionicHeterogeneityProbe`.
3. Constructs `nSamples` integration points uniformly distributed over the
   transmural coordinate `t ∈ [0, 1]` and applies the `ionicHeterogeneity`
   parameter map.
4. Runs a single-cell ODE simulation for `duration` ms at step size `dt` ms.
5. Writes outputs to `postProcessing/ionicHeterogeneityProbe/`:
   - `Vm_traces.csv` — voltage traces for every sample point
   - `AP_metrics.csv` — per-sample RMP, peak, amplitude, max dV/dt, APD30/50/70/90
   - `smoothness_report.csv` — adjacent-sample APD jumps and waveform RMS
   - `APD_envelope_report.csv` — per-sample APD monotonicity check (optional)

The utility returns exit code `1` if any sample produces an invalid action
potential, if any adjacent-sample transition exceeds the configured thresholds,
or if APD envelope checks fail (when enabled).

## Configuration

Optional `constant/ionicHeterogeneityProbe` dictionary keys:

| Key | Default | Description |
|---|---|---|
| `nSamples` | 101 | Number of transmural sample points |
| `duration` | 1000.0 | Simulation duration in ms |
| `dt` | 0.1 | ODE time step in ms |
| `writeEvery` | 1 | Write every N steps to `Vm_traces.csv` |
| `minAmplitude` | 50.0 | Minimum acceptable AP amplitude (mV) |
| `minPeak` | 0.0 | Minimum acceptable peak voltage (mV) |
| `maxRecoveryRise` | 5.0 | Max post-peak rise allowed (mV) before flagging a secondary rise |
| `checkSecondaryRise` | false | Enable secondary-rise check |
| `maxAPDJump` | 20.0 | Max APD difference (ms) between adjacent samples |
| `maxWaveformRMS` | 10.0 | Max waveform RMS (mV) between adjacent samples |
| `checkAPDEnvelope` | true | Check that each sample's APDs lie within the endo–M–epi envelope |
| `failOnAPDEnvelope` | true | Return exit code 1 on envelope violation |
| `maxAPDBoundTolerance` | 5.0 | Tolerance (ms) for the envelope bound check |
| `endoReferenceT` | 0.0 | Transmural coordinate of the endocardial reference |
| `mCellReferenceT` | midpoint | Transmural coordinate of the M-cell reference |
| `epiReferenceT` | 1.0 | Transmural coordinate of the epicardial reference |

A `singleCellStimulus` sub-dictionary can also be provided to override the
stimulus protocol for the probe run.

## Usage

```bash
ionicHeterogeneityProbe -case <caseDir>
```

## Companion script

`plotIonicHeterogeneityProbe.py` in this directory generates plots from the
CSV outputs written by the utility.
