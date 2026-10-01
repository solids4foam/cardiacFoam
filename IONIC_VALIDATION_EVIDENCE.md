# Ionic-model validation evidence carried into this worktree

This records local campaign artifacts produced before the active-tension
changes. The active-tension work does not modify `src/ionicModels`, but these
`/tmp` results must be reproduced from the final commit before a release.

## In-scope inventory

| Model | Scalar | Host-batched | CUDA | Tissue scope |
| --- | --- | --- | --- | --- |
| AlievPanfilov | yes | yes | yes | myocardium |
| BuenoOrovio | yes | yes | yes | myocardium |
| Courtemanche | yes | yes | yes | myocardium |
| Fabbri | yes | yes | yes | single cell (AV node) |
| Gaur | yes | yes | yes | myocardium |
| Grandi | yes | yes | yes | myocardium |
| PerisYague | yes | yes | yes | myocardium |
| Stewart | yes | yes | yes | myocardium |
| TNNP | yes | yes | yes | myocardium |
| TWorld | yes | yes | yes | myocardium |
| ToRORd_dynCl | yes | yes | yes | myocardium |
| Trovato | yes | yes | yes | myocardium |

The batched selection name is `<Model>compactBatched`; CUDA is a runtime
backend, confirmed by the `using CUDA device` log line.

## Scalar versus CUDA slabs

`/tmp/cardiac_niederer_all_models_20260928/summary.csv` records all 11
myocardial models on the 52,500-cell Niederer geometry through 15 ms with 25
batched substeps. Every scalar and CUDA case completed. Activated-cell counts
match at 5 and 10 ms; at 15 ms BuenoOrovio differs by one of 52,500 cells,
with a 0.000600 ms activation p95 shift. The largest recorded 15 ms Vm RMSE
is 0.011899 mV and activation p95 shift is 0.000900 ms (ToRORd_dynCl).

The 1,500-cell 15 ms campaign reports end-to-end scalar/CUDA wall-time ratios
of about 1.21x (AlievPanfilov) to 8.30x (Gaur). These are solver times, not
isolated GPU-kernel speedups.

## Host-batched versus CUDA SBDF2

`/tmp/cardiac_sbdf2_all_model_parity_9476` contains paired 1,500-cell cases
at 5, 10, and 15 ms. All models completed and logged CUDA selection.

| Models | Strict 1e-10 state/current comparison | Interpretation |
| --- | --- | --- |
| AlievPanfilov, Gaur, Grandi, Stewart, TNNP, TWorld, Trovato | pass | roundoff-level backend agreement |
| BuenoOrovio, Courtemanche, PerisYague, ToRORd_dynCl | fail | equal activation counts, but arithmetic-order current/state differences exceed the deliberately strict comparator |

At 15 ms, maximum Vm differences in the latter group are respectively
1.58e-05, 2.55e-07, 4.21e-07, and 9.94e-08 mV. They are not exact backend
parity; the final campaign must use justified quantity-specific
absolute-plus-relative tolerances, and report activation/waveform metrics.

## Fabbri and time controls

Fabbri is AV-node-only. The local 1--2 s single-cell host/CUDA campaign
(`/tmp/cardiac_fabbri_singlecell_parity_9519.log`) has 10,001 samples, 33 rate
columns, and 68 Vm/state/current/rate columns with zero saved-precision
difference.

`/tmp/cardiac_gpu_time_control_tnnp_20260928/summary.csv` confirms the actual
batched relation: ionic step = tissue `deltaT / batchedSubsteps`. The tested
combinations were 2 us/5 (0.4 us ionic), 2 us/1 (2 us), and 4 us/10 (0.4 us).
The tissue solve and `Im` update occur on the tissue step in this path; a
final code audit must confirm this for every supported integrator.

## Required release gates

- Materialize and rerun the matrix from the final commit outside `/tmp`.
- Validate full Niederer duration, timestep sensitivity, and restart paths.
- Add CUDA-event H2D/kernel/D2H timings before reporting acceleration.
