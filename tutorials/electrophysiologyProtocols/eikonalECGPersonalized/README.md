# eikonalECGPersonalized

Small, fast, real (non-manufactured) `eikonalECG` regression fixture
exercising the opt-in `personalizedTemplates` block: a coarse 2cm cube with
a real transmural gradient drives `ionicHeterogeneity`, which the
`eikonalECG` `personalizedTemplates` generator reuses as-is to pace three
single cells (endo/mid/epi) of a real ionic model and build dynamic
action-potential templates, replacing the solver's fixed compiled-in
`tissueTemplates.H` traces.

## Stack

- myocardium solver: `eikonalSolver`
- ionic heterogeneity: `mode transmuralBands` (`field t`, `endoMInterface
  0.3`, `mEpiInterface 0.7`)
- ECG: `ecgDomains.ECG.ecgSolver eikonalECG` with `personalizedTemplates`
  (`ionicModel TWorldcompactBatched`, `batchedSubsteps 100`, 10-beat
  single-cell pacing per tissue band)

`TWorldcompactBatched`'s batched CPU kernel needs `batchedSubsteps 100`
(not the default 1) to avoid a SIGFPE from `TWorld`'s stiffer dynamics
under the generator's outer `dt=1e-4`, and `stim_amplitude 60.0` (the
repo's real-mV/real-current TWorld convention, not a BuenoOrovio-scale
value) to clear `eikonalTemplateGenerator.C`'s `validateTemplate()`
amplitude check — see `constant/electroProperties` for the full note.

## Execution

```bash
./Allrun
```

This is the baseline of a 3-tutorial family exercising the same mesh,
stimulus, and `personalizedTemplates` config against the three
`ionicHeterogeneity` modes — see `../eikonalECGPersonalizedNamedRegions`
and `../eikonalECGPersonalizedCellZoneRegions`.
