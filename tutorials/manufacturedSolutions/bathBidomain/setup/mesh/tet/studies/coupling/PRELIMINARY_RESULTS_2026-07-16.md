# Preliminary bath coupling-control results

The lightweight build was run serially on the exact meshes and time controls
of the reported tetrahedral bath study. With
`bathPredictorCorrector false`, the stored potential,
reconstructed-current-jump, and intracellular leakage metrics were reproduced
at `N=10,20,40`.

| N | variant | heart $L_2(\phi_e)$ | bath $L_2(\phi_b)$ | $x=0$ assembled-current $L_2$ |
|---:|---|---:|---:|---:|
| 10 | one pass | $6.06607\times10^{-3}$ | $8.21426\times10^{-3}$ | $4.67976\times10^{-4}$ |
| 10 | predictor--corrector | $2.81749\times10^{-3}$ | $4.04936\times10^{-3}$ | $4.64175\times10^{-4}$ |
| 20 | one pass | $1.81197\times10^{-3}$ | $2.44473\times10^{-3}$ | $2.94155\times10^{-4}$ |
| 20 | predictor--corrector | $9.74917\times10^{-4}$ | $1.35486\times10^{-3}$ | $2.90330\times10^{-4}$ |
| 40 | one pass | $5.16474\times10^{-4}$ | $7.08605\times10^{-4}$ | $2.26414\times10^{-4}$ |
| 40 | predictor--corrector | $3.07773\times10^{-4}$ | $4.32837\times10^{-4}$ | $2.26572\times10^{-4}$ |

The predictor--corrector reduces the potential errors by 40--54% across the
three levels. Its effect on the assembled-current error is only -1.30% to
+0.07%. Eight additional global non-orthogonal assemblies change that flux
error by +3.75% and +0.27% at `N=10,20`.

These controls reject the hypothesis that incomplete temporal PDE coupling or
an unconverged deferred non-orthogonal source is the dominant cause of the
nearly stationary assembled-current error. Disabling the predictor--corrector
is nevertheless measurably less accurate for the cell-centred potentials, so
`bathPredictorCorrector` defaults to `true`.

The original `matchedSubmeshStudy/N40` CSV reported the diagnostic method as
`unweightedHarmonic`, whereas the `N=10,20,80` rows and the paper formulation
use `distanceWeightedHarmonic`. Re-evaluating the same `N=40` field with the
paper method gives $2.26414\times10^{-4}$, identical to the archived serial and
parallel distance-weighted diagnostics. The paper table must use that value;
mixing the unweighted `N=40` diagnostic into the distance-weighted ladder
artificially produces a larger non-monotone jump.
