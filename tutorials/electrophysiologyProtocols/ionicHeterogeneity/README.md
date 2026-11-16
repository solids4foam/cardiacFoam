# ionicHeterogeneity

## What it does

Runs the meshless `ionicHeterogeneityProbe` utility on the Bueno-Orovio
transmural heterogeneity map: it samples the transmural coordinate from `t=0`
to `t=1`, applies the endocardial, M-cell and epicardial bands, and writes one
action potential per sample with its metrics. `./Allrun overrides` runs the
same probe with the constant overrides of `constant/ionicConstantOverrides`.

- **Exercises:** `ionicHeterogeneity` in `namedRegions` mode with blended
  transitions, and `ionicConstantOverrides` (a `global` scale, a `set` and a
  `scale` on named regions). No other regression exercises the overrides.
- **Use it for:** checking a heterogeneity map or an override parametrisation
  (for example the ones of `../../idealizedHeart/pathos/ionicPathology`) at
  single-cell cost before running it in tissue.

## Regression

Tier: `standard`.

`regression/regressionTest.sh` runs the probe twice, seconds each:

- as configured: `peak` and `APD90` from `AP_metrics.csv` at five transmural
  samples (the three plateaus and both blends) and `waveformRMS` from
  `smoothness_report.csv` at a plateau and both blends, against
  `regression/ionicHeterogeneity.reference` (13 rows);
- with `./Allrun overrides`: its log must report each of the three overrides,
  and the same metrics are compared against `regression/overrides.reference`
  (13 rows).

The values are frozen from an earlier run. Run it alone with:

```bash
./regression/regressionTest.sh
```

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
    mode              namedRegions;
    field             t;
    transitionMode    blend;
    transitionWidth   0.1;
    smoothing         smoothstep;
    regions
    {
        endocardialCells { range (0 0.3); }
        mCells           { range (0.3 0.7); }
        epicardialCells  { range (0.7 1); }
    }
}
```

## Run

Build the utility first:

```bash
wmake applications/utilities/ionicHeterogeneityProbe
```

Then run this tutorial:

```bash
./Allrun              # the map of constant/electroProperties
./Allrun overrides    # the same, with constant/ionicConstantOverrides
```

`./Allrun overrides` copies `constant/ionicConstantOverrides` to
`local/overrides`, which `constant/electroProperties` includes in
`singleCellSolverCoeffs`; a plain `./Allrun` leaves it empty.

The utility returns non-zero if any configured smoothness or AP-shape check
fails; the CSV files are still written and `Allrun` exits with that status.
The figures are produced by the plotting script described below.

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

To generate the figures:

```bash
python3 ../../../applications/utilities/ionicHeterogeneityProbe/plotIonicHeterogeneityProbe.py \
  postProcessing/ionicHeterogeneityProbe
```

To also write separate overlays of every computed voltage trace in the Endo-M
and M-Epi areas, add `--all-traces`.

