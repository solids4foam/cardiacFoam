# replayLeadFieldECG

Compare precomputed reciprocal projection with production passive torso
recovery using the same stored myocardial voltage samples. This is serial
post-processing of a fixed, embedded-heart support mesh. Reconstruct parallel
source outputs first. No ionic step or cardiac time integration is performed.

The independent recovery-owned bidomain provider supplies physical Gi/Ge;
`constant/passiveECGProperties` contains `bidomainSolverCoeffs` and
`bathPotentialDomain` blocks plus `maxCorrectors`, `absoluteTolerance` and
`relativeTolerance`, with a separate `maxLinearResidual` acceptance check.
Source snapshots are full-support `Vm` fields in volts,
restricted through the production myocardial cell map. The observer inputs
and fields come from `constant/leadFieldProperties`.

At each selected time the utility holds Vm fixed, iterates `solvePhiE()` until
the maximum potential update meets the declared criterion, samples production
`torsoECG`, and applies exactly the reciprocal lead's balanced electrode
weights. It records both signals, differences, correction count, update and linear residual.
Vm internal and boundary arrays must remain unchanged. No source field is
written. An unconverged recovery stops rather than supplying a comparison.
The stopping criterion is not an error bound: repeat with tighter controls
when needed to distinguish iterative from recording/discretisation effects.

```bash
replayLeadFieldECG -time 0.005:0.04
# A tighter fixed-point repeat (preserve the original outputs first):
replayLeadFieldECG -absoluteTolerance 1e-12 -relativeTolerance 1e-9
```

`postProcessing/matchedLeadFieldECG.dat` retains every sample;
`matchedTorsoPotentials.dat` retains the raw electrode potentials;
`matchedLeadFieldECGSummary.dat` reports per-channel RMS differences, forward
RMS amplitudes and maximum differences in volts. Each source run's online
pseudo-ECG is separate. A different recording kernel is not the reference
for the bounded torso calculation.

Current matching scope requires an insulated exterior, no independent
extracellular forcing, and `distanceWeightedHarmonic` interface interpolation.
Nonempty ground/current-patch dictionaries and incompatible exterior field
patches are rejected. The reciprocal builder uses this interface construction;
other bath interpolation models need their own reciprocal construction or
separate comparison. Agreement is not assumed: the reciprocal cell-gradient
quadrature is not an exact discrete transpose of the corrected forward solve.

The [idealized-heart ECG tutorial](../../../tutorials/idealizedHeart/electroHeartECG/README.md)
supplies monodomain/bidomain source variants, the shared heart/torso geometry,
six Wilson-referenced leads, and the complete dictionary structure.
