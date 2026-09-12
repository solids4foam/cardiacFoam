# BuenoOrovio Ionic Heterogeneity Probe

This tutorial runs the meshless `ionicHeterogeneityProbe` utility for the
Bueno-Orovio transmural heterogeneity map.

It samples synthetic transmural distances from `t=0` to `t=1`, applies the
configured endocardial, M-cell, and epicardial bands, and writes one voltage
trace per sample. It is intended to check whether the smoothed transition bands
produce plausible action potentials before using the same heterogeneity settings
in a spatial myocardium case.

The current transmural map uses one-sided transition bands. In other words,
`transitionWidth` starts at each interface and extends toward the next tissue
type; it is not interpreted as `interface +/- transitionWidth`.

For the default dictionary:

- Endo plateau: `0 <= t <= 0.3`
- Endo-M smooth blend: `0.3 < t < 0.4`
- M-cell plateau: `0.4 <= t <= 0.7`
- M-Epi smooth blend: `0.7 < t < 0.8`
- Epi plateau: `0.8 <= t <= 1`

This comes from:

```cpp
ionicHeterogeneity
{
    mode              transmuralBands;
    field             t;
    endoMInterface    0.3;
    mEpiInterface     0.7;
    transitionWidth   0.1;
    smoothing         smoothstep;
}
```

## Run

Build the utility first:

```bash
wmake applications/utilities/ionicHeterogeneityProbe
```

Then run this tutorial:

```bash
./Allrun
```

The utility may return non-zero if any configured smoothness or AP-shape check
fails. CSV files are still written, and `Allrun` still attempts to generate the
plots.

## Configuration

The ionic model and transmural heterogeneity map are configured in:

```text
constant/electroProperties
```

The synthetic sampling, time integration, stimulus, and smoothness tolerances are
configured in:

```text
constant/ionicHeterogeneityProbe
```

The default Bueno-Orovio stimulus amplitude is `0.4`.

## Outputs

The utility writes:

```text
postProcessing/ionicHeterogeneityProbe/Vm_traces.csv
postProcessing/ionicHeterogeneityProbe/AP_metrics.csv
postProcessing/ionicHeterogeneityProbe/smoothness_report.csv
postProcessing/ionicHeterogeneityProbe/APD_envelope_report.csv
```

The plotting script writes:

```text
postProcessing/ionicHeterogeneityProbe/Vm_transition_bands_2D.png
postProcessing/ionicHeterogeneityProbe/Vm_transmural_surface_3D.png
postProcessing/ionicHeterogeneityProbe/APD_transmural_metrics.png
```

To regenerate only the figures:

```bash
python3 ../../../applications/utilities/ionicHeterogeneityProbe/plotIonicHeterogeneityProbe.py \
  postProcessing/ionicHeterogeneityProbe
```

To also write separate overlays of every computed voltage trace in the Endo-M
and M-Epi areas, add `--all-traces`.
